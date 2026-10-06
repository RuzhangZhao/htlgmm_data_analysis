## CML vs htlgmm(no penalty) vs internal-only, on the htlgmm DGP, pZ=10, pW=10.
## Per (family, n): NREP reps; stores per-rep coefficient estimates, SEs, and the
## test metric (AUC logistic / R2 linear) for the 3 methods. Aggregate + plot
## with compare_cml_htlgmm_plot.R. Heavy -> run on the cluster.
## Usage: Rscript compare_cml_htlgmm_sim.R [family] [nid]   (no args = loop all)
suppressMessages({library(MASS); library(glmnet); library(magic); library(mvtnorm); library(pROC); library(htlgmm)})
source("cml_apply.R")

pZ <- 10; pW <- 10
n_list <- c(200, 400, 600, 800, 1000, 1200, 1500, 2000, 2500, 3000)
NREP <- 100
rmult <- 10
ntest <- 1e5
outdir <- "compare_results"; dir.create(outdir, showWarnings = FALSE)

Z_nonnull_index <- 1:pZ
W_nonnull_index <- c(1, 2, 3, 7)
ZWlinklist <- vector("list", pZ)
ZWlinklist[[8]] <- 1; ZWlinklist[[6]] <- 2; ZWlinklist[[5]] <- 3; ZWlinklist[[4]] <- 4
cor_block <- function(s, rho) rho^abs(row(diag(s)) - col(diag(s)))

make_dgp <- function(family) {
    coefZ <- rep(0, pZ); coefW <- rep(0, pW)
    if (family == "logistic") { b0 <- -1.66; coefZ[Z_nonnull_index] <- 0.3/2; coefW[W_nonnull_index] <- 0.25/2 }
    else                      { b0 <- -2.3;  coefZ[Z_nonnull_index] <- 0.3;   coefW[W_nonnull_index] <- 0.25   }
    coefZW <- c(coefZ, coefW)
    
    sim <- function(N) {
        Sigma <- Reduce(adiag, replicate((pZ+pW)/10, cor_block(10, 0.5), simplify = FALSE))
        X <- rmvnorm(N, sigma = Sigma)
        for (i in seq_along(ZWlinklist)) for (j in ZWlinklist[[i]]) {
            jj <- pZ + W_nonnull_index[j]; X[, jj] <- X[, jj] + 0.3 * X[, Z_nonnull_index[i]]
        }
        X <- scale(X)
        if (family == "logistic") y <- rbinom(N, 1, expit_cml(X %*% coefZW + b0))
        else                      y <- X %*% coefZW + b0 + rnorm(N, 0, 3)
        list(X = X, Z = X[, 1:pZ, drop = FALSE], W = X[, (pZ+1):(pZ+pW), drop = FALSE], y = c(y))
    }
    
    list(coefZW = coefZW, b0 = b0, sim = sim)
}

run_cell <- function(family, n) {
    dgp <- make_dgp(family); sim <- dgp$sim
    if (family == "logistic") {
        truth <- c(dgp$b0, dgp$coefZW); q <- 1 + pZ + pW
        z_idx <- 1 + (1:pZ); w_idx <- 1 + pZ + (1:pW); metric_name <- "AUC"
    } else {
        truth <- dgp$coefZW; q <- pZ + pW
        z_idx <- 1:pZ; w_idx <- pZ + (1:pW); metric_name <- "R2"
    }
    
    hfam <- if (family == "logistic") "binomial" else "gaussian"   # htlgmm/glm family name
    set.seed(1); te <- sim(ntest)
    y0 <- if (family == "logistic") te$y else scale(te$y, scale = FALSE)
    
    metric_of <- function(beta)
        if (family == "logistic")
            suppressMessages(as.numeric(pROC::auc(y0, c(expit_cml(te$X %*% beta[-1] + beta[1])), direction = "<")))
    else 1 - sum((y0 - te$X %*% beta)^2) / sum(y0^2)
    
    est <- list(internal=matrix(NA,NREP,q), cml=matrix(NA,NREP,q), htlgmm=matrix(NA,NREP,q))
    se <- list(internal=matrix(NA,NREP,q), cml=matrix(NA,NREP,q), htlgmm=matrix(NA,NREP,q))
    metric <- matrix(NA, NREP, 3, dimnames = list(NULL, c("internal","cml","htlgmm")))
    
    for (r in 1:NREP) {
        set.seed(1000 + r)
        ex <- sim(rmult * n); mn <- sim(n)
        
        if (family == "logistic") {
            eg <- glm(ex$y ~ ex$Z, family = binomial)
            ext_info <- cml_external_info(eg, nExt = rmult * n)
            study_info <- list(list(Coeff = coef(eg)[-1], Covariance = vcov(eg)[-1,-1], Sample_size = rmult*n))
            
            ym <- mn$y; A_arg <- 1
            fI <- glm(mn$y ~ mn$Z + mn$W, family = binomial)
            est$internal[r,] <- coef(fI); se$internal[r,] <- sqrt(diag(vcov(fI)))
            metric[r,"internal"] <- metric_of(coef(fI))
            
            cm <- runCML(mn$y, mn$Z, mn$W, ext_info, y0, te$X, pZ, pW, ExUncertainty = TRUE)
            metric[r,"cml"] <- cm$auc
            
        } else {
            ym <- scale(mn$y, scale = FALSE)
            eg <- lm(scale(ex$y, scale = FALSE) ~ ex$Z)
            ext_info <- cml_external_info_linear(eg, nExt = rmult * n)
            study_info <- list(list(Coeff = coef(eg)[-1], Covariance = vcov(eg)[-1,-1], Sample_size = rmult*n))
            
            A_arg <- NULL
            fI <- lm(c(ym) ~ 0 + mn$Z + mn$W)
            est$internal[r,] <- coef(fI); se$internal[r,] <- sqrt(diag(vcov(fI)))
            metric[r,"internal"] <- metric_of(coef(fI))
            
            cm <- runCML_gaussian(c(ym), mn$Z, mn$W, ext_info, y0, te$X, pZ, pW, ExUncertainty = TRUE)
            metric[r,"cml"] <- cm$R2
        }
        
        est$cml[r,] <- cm$beta
        se$cml[r, cm$selected_vars$position] <- sqrt(pmax(cm$selected_vars$variance, 0))
        
        hg <- tryCatch(htlgmm(y = c(ym), Z = mn$Z, W = mn$W, ext_study_info = study_info,
                              family = hfam, penalty_type = "none", A = A_arg),
                       error = function(e) NULL)
        
        if (!is.null(hg)) {
            est$htlgmm[r,] <- as.numeric(hg$beta)
            metric[r,"htlgmm"] <- metric_of(as.numeric(hg$beta))
            if (!is.null(hg$selected_vars))
                se$htlgmm[r, hg$selected_vars$position] <- sqrt(pmax(hg$selected_vars$variance, 0))
        }
        
        if (r %% 50 == 0) message(family, " n=", n, " rep ", r, "/", NREP)
    }
    
    out <- list(est = est, se = se, metric = metric, truth = truth,
                z_idx = z_idx, w_idx = w_idx, family = family, n = n,
                pZ = pZ, pW = pW, metric_name = metric_name)
    
    saveRDS(out, file.path(outdir, sprintf("compare_%s_n_%d_pZ_%d_pW_%d.rds", family, n, pZ, pW)))
    message("saved ", sprintf("compare_%s_n_%d_pZ_%d_pW_%d.rds", family, n, pZ, pW))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
    run_cell(args[1], n_list[as.numeric(args[2])])
} else {
    for (family in c("linear", "logistic")) for (n in n_list) run_cell(family, n)
}
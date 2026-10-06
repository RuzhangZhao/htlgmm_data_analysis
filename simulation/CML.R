## -------------------------------------------------------------------------
## Helpers
## -------------------------------------------------------------------------
expit <- function(x) 1 / (1 + exp(-x))


## -------------------------------------------------------------------------
## CORE 1: CML for LOGISTIC full model    logit P(Y=1|X,B) = gamma'(1,X,B)
##         external REDUCED model         logit P(Y=1|X)   = beta'(1,X)
## -------------------------------------------------------------------------
## YInt        : internal outcome vector (0/1)
## XInt        : internal X covariates used by the external model (no intercept)
## BInt        : internal-only covariates (not in the external model); NULL ok
## betaHatExt  : external reduced-model coefficients, order (intercept, X)
## gammaHatInt : starting value = internal-only full-model MLE, order (int,X,B)
## tol,maxIter : Newton-Raphson controls
## factor      : step-halving factor in (0,1]; 1 = plain Newton-Raphson
## -------------------------------------------------------------------------
cml_logistic <- function(YInt, XInt, BInt, betaHatExt, gammaHatInt,
                         tol = 1e-8, maxIter = 400, factor = 1) {
    XInt <- as.matrix(XInt)
    Xtilde <- cbind(1, XInt)                         # (1, X)    : n x p
    XBtilde <- if (is.null(BInt)) Xtilde else cbind(Xtilde, as.matrix(BInt))  # (1,X,B): n x q
    
    n <- length(YInt)
    p <- ncol(Xtilde)
    q <- ncol(XBtilde)
    
    gammaHat <- matrix(gammaHatInt, q, 1)
    betaHatExt <- matrix(betaHatExt, p, 1)
    lambda <- matrix(0, p, 1)                        # Lagrange multipliers start 0
    
    estDiff <- 1
    counter <- 0
    
    while (estDiff > tol && counter < maxIter) {
        v1 <- matrix(0, q, 1)
        v2 <- matrix(0, p, 1)
        h11 <- matrix(0, q, q)
        h12 <- matrix(0, p, q)
        h22 <- matrix(0, p, p)
        
        for (i in 1:n) {
            Xi <- matrix(Xtilde[i, ], p, 1)
            XBi <- matrix(XBtilde[i, ], q, 1)
            y <- YInt[i]
            
            pInt <- as.numeric(expit(t(XBi) %*% gammaHat))
            pExt <- as.numeric(expit(t(Xi) %*% betaHatExt))
            diffP <- pInt - pExt
            pq <- pInt * (1 - pInt)
            lamTXi <- as.numeric(t(lambda) %*% Xi)
            denom <- as.numeric(n - lamTXi * diffP)      # profile-likelihood denominator
            
            ## score
            v1 <- v1 + (y - pInt) * XBi + lamTXi * pq * XBi / denom
            v2 <- v2 + diffP * Xi / denom
            
            ## negative Hessian blocks
            xiC <- denom * pq
            h11 <- h11 - pq * (XBi %*% t(XBi)) +
                (xiC * (1 - 2 * pInt) * lamTXi + pq^2 * lamTXi^2) *
                (XBi %*% t(XBi)) / denom^2
            
            h12 <- h12 + (xiC * (Xi %*% t(XBi)) +
                              pq * diffP * lamTXi * (Xi %*% t(XBi))) / denom^2
            
            h22 <- h22 + diffP^2 * (Xi %*% t(Xi)) / denom^2
        }
        
        V <- rbind(v1, v2)
        H <- rbind(cbind(h11, t(h12)), cbind(h12, h22))
        
        old <- matrix(c(gammaHat, lambda), p + q, 1)
        pars <- old - factor * solve(H) %*% V          # Newton-Raphson step
        estDiff <- max(abs(pars - old))
        
        gammaHat <- pars[1:q]
        lambda <- pars[(q + 1):(q + p)]
        counter <- counter + 1
    }
    
    if (counter == maxIter) warning("CML (logistic): max iterations reached without convergence")
    names(gammaHat) <- colnames(XBtilde)
    
    list(gammaHat = gammaHat, lambda = lambda, iter = counter, estDiff = estDiff)
}


## -------------------------------------------------------------------------
## CORE 2: CML for LINEAR full model       E(Y|X,B) = gamma'(1,X,B)
##         external REDUCED model          E(Y|X)   = beta'(1,X)
## -------------------------------------------------------------------------
cml_linear <- function(YInt, XInt, BInt, betaHatExt, gammaHatInt,
                       tol = 1e-8, maxIter = 400, factor = 1) {
    XInt <- as.matrix(XInt)
    Xtilde <- cbind(1, XInt)
    XBtilde <- if (is.null(BInt)) Xtilde else cbind(Xtilde, as.matrix(BInt))
    
    n <- length(YInt)
    p <- ncol(Xtilde)
    q <- ncol(XBtilde)
    
    gammaHat <- matrix(gammaHatInt, q, 1)
    betaHatExt <- matrix(betaHatExt, p, 1)
    lambda <- matrix(0, p, 1)
    
    sigma2 <- mean((YInt - XBtilde %*% gammaHat)^2)      # full-model resid var
    sigma2Ext <- mean((YInt - Xtilde %*% betaHatExt)^2)  # external resid var fixed
    
    estDiff <- 1
    counter <- 0
    
    while (estDiff > tol && counter < maxIter) {
        v1 <- matrix(0, q, 1)
        v2 <- matrix(0, p, 1)
        h11 <- matrix(0, q, q)
        h12 <- matrix(0, p, q)
        h22 <- matrix(0, p, p)
        
        for (i in 1:n) {
            Xi <- matrix(Xtilde[i, ], p, 1)
            XBi <- matrix(XBtilde[i, ], q, 1)
            
            resid <- YInt[i] - as.numeric(t(XBi) %*% gammaHat)
            gap <- as.numeric(t(XBi) %*% gammaHat - t(Xi) %*% betaHatExt)
            denom <- as.numeric(n - sigma2Ext^(-1) * (t(lambda) %*% Xi) * gap)
            
            v1 <- v1 + (1 / sigma2) * XBi * resid +
                sigma2Ext^(-1) * (XBi %*% t(Xi)) %*% lambda / denom
            
            v2 <- v2 + sigma2Ext^(-1) * Xi * gap / denom
            
            h11 <- h11 - (1 / sigma2) * (XBi %*% t(XBi)) +
                sigma2Ext^(-2) * (XBi %*% t(Xi) %*% lambda) %*%
                t(XBi %*% t(Xi) %*% lambda) / denom^2
            
            h12 <- h12 + (denom * sigma2Ext^(-1) * (Xi %*% t(XBi)) +
                              sigma2Ext^(-2) * gap *
                              (Xi %*% t(lambda) %*% Xi %*% t(XBi))) / denom^2
            
            h22 <- h22 + sigma2Ext^(-2) * gap^2 * (Xi %*% t(Xi)) / denom^2
        }
        
        V <- rbind(v1, v2)
        H <- rbind(cbind(h11, t(h12)), cbind(h12, h22))
        
        old <- matrix(c(gammaHat, lambda), p + q, 1)
        pars <- old - factor * solve(H) %*% V
        estDiff <- max(abs(pars - old))
        
        gammaHat <- pars[1:q]
        lambda <- pars[(q + 1):(q + p)]
        counter <- counter + 1
        
        sigma2 <- mean((YInt - XBtilde %*% gammaHat)^2)
    }
    
    if (counter == maxIter) warning("CML (linear): max iterations reached without convergence")
    names(gammaHat) <- colnames(XBtilde)
    
    list(gammaHat = gammaHat, lambda = lambda, sigma2 = sigma2,
         iter = counter, estDiff = estDiff)
}


## -------------------------------------------------------------------------
## Asymptotic variance for the LOGISTIC CML estimator (single external model).
## Var(gamma_CML) = (B + C L^{-1} C')^{-1}, which is <= B^{-1} = Var(gamma_Int).
## If ExUncertainty=TRUE, adds the term from estimating beta on a finite
## external sample (rho = n/m, CovExt = Var(betaHatExt)).
## -------------------------------------------------------------------------
cml_var_logistic <- function(YInt, XInt, BInt, gammaHatInt, betaHatExt,
                             CovExt = NULL, rho = 0, ExUncertainty = FALSE) {
    XInt <- as.matrix(XInt)
    Xtilde <- cbind(1, XInt)
    XBtilde <- if (is.null(BInt)) Xtilde else cbind(Xtilde, as.matrix(BInt))
    
    n <- length(YInt)
    p <- ncol(Xtilde)
    q <- ncol(XBtilde)
    
    gammaHatInt <- matrix(gammaHatInt, q, 1)
    betaHatExt <- matrix(betaHatExt, p, 1)
    
    B <- matrix(0, q, q)
    C <- matrix(0, q, p)
    L <- matrix(0, p, p)
    Q <- matrix(0, p, p)
    
    for (i in 1:n) {
        Xi <- matrix(Xtilde[i, ], p, 1)
        XBi <- matrix(XBtilde[i, ], q, 1)
        
        pInt <- as.numeric(expit(t(XBi) %*% gammaHatInt))
        pExt <- as.numeric(expit(t(Xi) %*% betaHatExt))
        
        B <- B + (XBi %*% t(XBi)) * pInt * (1 - pInt)
        C <- C + (XBi %*% t(Xi)) * pInt * (1 - pInt)
        Q <- Q + (Xi %*% t(Xi)) * pInt * (1 - pInt)
        L <- L + ((pInt - pExt) * Xi) %*% t((pInt - pExt) * Xi)
    }
    
    asyV.I <- solve(B)
    core <- solve(B + C %*% solve(L) %*% t(C))       # CML variance ignoring ext. uncertainty
    asyV.CML <- core
    
    if (ExUncertainty) {
        QVQ <- t(Q) %*% CovExt %*% Q
        asyV.CML <- core + (1 / n) * rho *
            core %*% C %*% solve(L) %*% QVQ %*%
            t(solve(L)) %*% t(C) %*% core
    }
    
    dimnames(asyV.I) <- dimnames(asyV.CML) <- list(colnames(XBtilde), colnames(XBtilde))
    list(asyV.I = asyV.I, asyV.CML = asyV.CML)
}
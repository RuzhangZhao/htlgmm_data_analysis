args = commandArgs(trailingOnly = TRUE)
pZid = as.numeric(args[1])
nid = as.numeric(args[2])
pZ = c(10,40)[pZid]
pW = 150
nlist = c(200,400,600,800,1000,1200,1500,2000,2500,3000)
n = nlist[nid]
print(paste0("pW:",pW,";n:",n))
pacman::p_load(expm,magic,glmnet,cluster,MASS,corpcor,caret,pROC,mvtnorm)
library(htlgmm)
if(0){
    remove.packages("htlgmm")
    devtools::install_github("RuzhangZhao/htlgmm")
}
### parameter setting
coefZ = rep(0,pZ)
coefW = rep(0,pW)
intercept = -2.3
savepath="htlgmm_simulate/linear_pW_150/" #self-defined save path
if(pZ==10){
    Z_nonnull_index=1:10
    coefZ[Z_nonnull_index]=0.3
}else if(pZ == 40){
    Z_nonnull_index=c(1,2,3,14,15,26,27,28,39,40)
    coefZ[Z_nonnull_index]=0.3*1.28
}
W_nonnull_index = c(1,2,3,7,11,16,21,24,26,40,54,57,59,78,80)
coefW[W_nonnull_index]=0.25
coefZW = c(coefZ,coefW)

correct_pos = c(Z_nonnull_index,(pZ+W_nonnull_index))
correct_Zpos = Z_nonnull_index
correct_Wpos = pZ+W_nonnull_index

ZWlinklist=list()
ZWlinklist[[10]]=c(8,9)
ZWlinklist[[8]]=c(1)
ZWlinklist[[6]]=c(2)
ZWlinklist[[5]]=c(3)
ZWlinklist[[4]]=c(4)
ZWlinklist[[3]]=c(10)
ZWlinklist[[2]]=c(14)
ZWlinklist[[1]]=c(11,13)

ZWlinkcolnames=c(paste0("Z",1:pZ),paste0("W",1:pW))

cor_block=function(size,rho){return(rho^abs(row(diag(size))-col(diag(size))))}
simulate_Xy=function(n,pZ=pZ,pW=pW,coefs=coefZW,intercept=intercept,
                     rho=0.5,block_size=10,adj_r=0.3,
                     W_nonnull_index=W_nonnull_index,
                     ZWlinklist=ZWlinklist){
    Sigma_unit=cor_block(block_size,rho)
    Sigma_list=list()
    p=pZ+pW
    for(i in 1:(p/block_size)){Sigma_list[[i]]=Sigma_unit}
    Sigma=Reduce(adiag,Sigma_list)
    X=rmvnorm(n,sigma=Sigma)
    for(i in 1:length(ZWlinklist)){
        if(length(ZWlinklist[[i]])>0){
            for(j in ZWlinklist[[i]]){
                jj=pZ+W_nonnull_index[j]
                X[,jj]=X[,jj]+adj_r*X[,Z_nonnull_index[i]]
            }
        }
    }
    X=scale(X)
    y=X%*%coefs+intercept+rnorm(n,0,3)
    list("X"=X,"y"=y)
}

summary_stat=function(Z,y){
    study_info=list()
    reducedlm=lm(y~.,data = data.frame(y,Z))
    study.m = list(Coeff=reducedlm$coefficients[-1],
                   Covariance=vcov(reducedlm)[-1,-1],Sample_size=nrow(Z))
    study_info[[1]] = study.m
    return(list("lm"=reducedlm,"study_info"=study_info))
}

RC=function(y0,predy){
    1-sum((predy-scale(y0,scale = F))^2)/sum((scale(y0,scale = F))^2)
}

adalasso_inference<-function(res_ada,X_main,y_main){
    beta_ada<-coef(res_ada,s="lambda.min")[-1]
    sigma_hat<-mean((X_main%*%beta_ada-y_main)^2)
    index_nonzero<-which(beta_ada!=0)
    Sigsum_scaled_nonzero<-crossprod(X_main[,index_nonzero])
    inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
    final_v<-diag(inv_Sigsum_scaled_nonzero)*sigma_hat
    
    pval_final<-pchisq(beta_ada[index_nonzero]^2/final_v,1,lower.tail = F)
    pval_final1<-p.adjust(pval_final,method = "BH")
    corrected_pos<-index_nonzero[which(pval_final1<0.05)]
    corrected_pos0<-index_nonzero[which(pval_final<0.05)]
    ada_res<-list()
    ada_res_selected_vars<-list()
    ada_res_selected_vars$position<-index_nonzero
    ada_res_selected_vars$coef<-beta_ada[index_nonzero]
    ada_res_selected_vars$variance<-final_v
    ada_res_selected_vars$pval<-pval_final
    ada_res_selected_vars$FDR_adjust_position<-corrected_pos
    
    ada_res$selected_vars<-ada_res_selected_vars
    ada_res$beta<-beta_ada
    ada_res
}

fdrgmm<-function(ada){
    sum(!ada$selected_vars$FDR_adjust_position%in%correct_pos)/max(1e-8,length(ada$selected_vars$FDR_adjust_position))
}

power<-function(ada){
    sum(ada$selected_vars$FDR_adjust_position%in%correct_pos)/length(correct_pos)
}

powerZ<-function(ada){
    sum(ada$selected_vars$FDR_adjust_position%in%correct_Zpos)/length(correct_Zpos)
}

powerW<-function(ada){
    sum(ada$selected_vars$FDR_adjust_position%in%correct_Wpos)/length(correct_Wpos)
}

cover<-function(ada){
    ans<-rep(F,length(correct_pos))
    for(i in 1:length(correct_pos)){
        pos = correct_pos[i]
        if(ada$beta[pos]!=0){
            idhere<-which(ada$selected_vars$position == pos)
            nonzerovar<-ada$selected_vars$variance[idhere]
            up<-ada$beta[pos]+1.96*sqrt(nonzerovar)
            low<-ada$beta[pos]-1.96*sqrt(nonzerovar)
            if( coefZW[pos]<up &coefZW[pos]>low  ){
                ans[i]<-T
            }
        }
    }
    mean(ans)
}

### Infinity Test Data for Validation
if(!file.exists(paste0(savepath,"linear_test_pZ_",pZ,".rds"))){
    ntest = 10^6
    Xy0=simulate_Xy(ntest,pZ=pZ,pW=pW,coefs=coefZW,intercept=intercept,
                    rho=0.5,block_size=10,adj_r=0.3,
                    W_nonnull_index=W_nonnull_index,
                    ZWlinklist=ZWlinklist)
    X0=Xy0$X
    y0=c(Xy0$y)
    RC(y0,X0%*%coefZW)
    saveRDS(Xy0,paste0(savepath,"linear_test_pZ_",pZ,".rds"))
}else{
    Xy0=readRDS(paste0(savepath,"linear_test_pZ_",pZ,".rds"))
    X0=Xy0$X
    y0=c(Xy0$y)
}

savefilename = paste0(savepath,"linear_n_",n,"_pZ_",pZ,"_r_10.rds")
savefilenamea = paste0(savepath,"linear_n_",n,"_pZ_",pZ,"_r_10ada.rds")
if(!file.exists(savefilename)){
    savefilelist=list()
    savefilelista=list()
}else{
    savefilelist=readRDS(savefilename)
    savefilelista=readRDS(savefilenamea)
}
if(length(savefilelist) == 100){break}

for(i in (length(savefilelist)+1):100){
    set.seed(i)
    seed.use=sample(1:2023,1)
    set.seed(seed.use)
    message(paste0("EPOCH:",i,",seed:",seed.use))
    # external data 10
    Xy_ext=simulate_Xy(n*10,pZ=pZ,pW=pW,coefs=coefZW,intercept=intercept,
                       rho=0.5,block_size=10,adj_r=0.3,
                       W_nonnull_index=W_nonnull_index,
                       ZWlinklist=ZWlinklist)
    X_ext=Xy_ext$X
    y_ext=c(Xy_ext$y)
    y_ext=scale(y_ext,scale = F)
    colnames(X_ext)=ZWlinkcolnames
    sum_res=summary_stat(X_ext[,1:pZ],y_ext)
    glm10=sum_res$lm
    study_info10=sum_res$study_info
    
    # internal data
    Xy_main=simulate_Xy(n,pZ=pZ,pW=pW,coefs=coefZW,intercept=intercept,
                        rho=0.5,block_size=10,adj_r=0.3,
                        W_nonnull_index=W_nonnull_index,
                        ZWlinklist=ZWlinklist)
    X_main=Xy_main$X
    y_main=c(Xy_main$y)
    y_main=scale(y_main,scale = F)
    Z_main=X_main[,1:pZ]
    W_main=X_main[,-c(1:pZ)]
    
    res_lasso=cv.glmnet(x = X_main,y = y_main,family = "gaussian")
    beta_lasso=c(coef(res_lasso,s="lambda.min")[-1])
    res_ridge=cv.glmnet(x = X_main,y = y_main,alpha=0,family = "gaussian")
    www=1/abs(c(coef(res_ridge,s="lambda.min")[-1]))^(1/2)
    res_ada=cv.glmnet(x = X_main,y = y_main,penalty.factor = www,family = "gaussian",)
    
    ## Compute RC for true model
    RC_true_model=RC(y0,X0%*%coefZW)
    RC_true_Z=RC(y0,X0[,1:pZ]%*%coefZ)
    RC_true_W=RC(y0,X0[,-c(1:pZ)]%*%coefW)
    ## Compute RC for main study only model
    RC_main_lasso=RC(y0,X0%*%beta_lasso)
    RC_main_adalasso=RC(y0,X0%*%coef(res_ada,s="lambda.min")[-1])
    res_main_adalasso<-adalasso_inference(res_ada,X_main,y_main)
    runhtlgmm=function(y_main,X_main,study_info,res_lasso,penalty_type){
        RC_ext=RC(y0,X0[,1:pZ]%*%study_info[[1]]$Coeff)
        weight_list<-c(1,1.5,2,4,8,16)-1
        nfolds<-10
        doinference=F
        
        if(penalty_type == "adaptivelasso"){
            doinference=TRUE
        }
        res_htlgmm_lasso_tuneC_mshrink=cv.htlgmm(y = y_main,Z = Z_main,W = W_main,
                                               penalty_type = penalty_type,nfolds = nfolds,
                                               initial_with_type="lasso",
                                               beta_initial = as.vector(coef(res_lasso,s="lambda.min"))[-1],
                                               inference = doinference,
                                               weight_adaptivelasso=www,
                                               ext_study_info = study_info,family = "gaussian",
                                               use_sparseC = F,
                                               tune_weight = T,
                                               weight_list = weight_list,
                                               tune_weight_method="mshrink")
        
        res_htlgmm_lasso_tuneC_ridge=cv.htlgmm(y = y_main,Z = Z_main,W = W_main,
                                               penalty_type = penalty_type,nfolds = nfolds,
                                               initial_with_type="lasso",
                                               weight_adaptivelasso=www,
                                               beta_initial = as.vector(coef(res_lasso,s="lambda.min"))[-1],
                                               inference = doinference,
                                               ext_study_info = study_info,family = "gaussian",
                                               use_sparseC = F,
                                               tune_weight = T,
                                               weight_list = weight_list,
                                               tune_weight_method="ridge")
        
        res_htlgmm_lasso_owGMM=cv.htlgmm(y = y_main,Z = Z_main,W = W_main,
                                          penalty_type = penalty_type,nfolds = nfolds,
                                          beta_initial = as.vector(coef(res_lasso,s="lambda.min"))[-1],
                                          inference = doinference,
                                          weight_adaptivelasso=www,
                                          ext_study_info = study_info,family = "gaussian",use_sparseC = F)
        
        R2_tuneC_mshrink=RC(y0,X0%*%res_htlgmm_lasso_tuneC_mshrink$beta)
        weightC_mshrink = res_htlgmm_lasso_tuneC_mshrink$weight_min
        
        R2_tuneC_ridge=RC(y0,X0%*%res_htlgmm_lasso_tuneC_ridge$beta)
        weightC_ridge = res_htlgmm_lasso_tuneC_ridge$weight_min
        R2_owGMM=RC(y0,X0%*%res_htlgmm_lasso_owGMM$beta)
        
        RC_list=c(mean(y_main),RC_true_model,RC_true_Z,RC_true_W,RC_ext,RC_main_lasso,RC_main_adalasso,
                  R2_tuneC_mshrink,weightC_mshrink,
                  R2_tuneC_ridge,weightC_ridge,
                  R2_owGMM)
        names(RC_list)=c("case_ratio","true_model","true_Z","true_W","external_glm","main_lasso","main_ada",
                         "R2_tuneC_mshrink","weightC_mshrink",
                         "R2_tuneC_ridge","weightC_ridge",
                         "R2_owGMM")
        print(round(RC_list,3))
        
        if(penalty_type == "adaptivelasso"){
            f_inf<-function(res_list,curf){
                sapply(res_list, function(i){
                    curf(i)
                })
            }
           
            res_list<-list(res_main_adalasso,res_htlgmm_lasso_tuneC_mshrink,
                           res_htlgmm_lasso_tuneC_ridge,
                           res_htlgmm_lasso_owGMM)
            
            inf_res=list()
            inf_res$fdr=f_inf(res_list,fdrgmm)
            inf_res$power=f_inf(res_list,power)
            inf_res$cover=f_inf(res_list,cover)
            inf_res$powerZ=f_inf(res_list,powerZ)
            inf_res$powerW=f_inf(res_list,powerW)
            
            res11<-rbind(inf_res$powerW,
                         inf_res$fdr,
                         inf_res$cover)
            colnames(res11)<-c('ada','mshrink',
                               'ridge','owGMM')
            rownames(res11)<-c('pW','fdr',"cv")
            print(res11)
        }
        list("RC"=RC_list)
    }
    savefilelist[[i]]=runhtlgmm(y_main,X_main,study_info10,res_lasso,"lasso")
    savefilelista[[i]]=runhtlgmm(y_main,X_main,study_info10,res_lasso,"adaptivelasso")
    if (i %% 5 == 0){
        saveRDS(savefilelist,savefilename)
        saveRDS(savefilelista,savefilenamea)
    }
}




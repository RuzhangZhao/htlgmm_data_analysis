args <- commandArgs(trailingOnly = TRUE)
disease_id = as.numeric(args[1])
seed.use = as.numeric(args[2])
nfold = 5
namelist<-c("colon","Breast","incAsthma","incCVD","incSTR")
savepath=paste0("UKB_Protein/")
diseasename = namelist[disease_id]
useprs="prs"

riskfactor<-readRDS(paste0("riskfactor/",diseasename,"_risk_zrz",useprs,".rds"))

proteomics0<-readRDS("all_proteomic_imputed.rds")
# find overlapping proteomics and risk factor
proteomics0<-proteomics0[which(rownames(proteomics0)%in%rownames(riskfactor)),]
riskfactor_main<-riskfactor[which(rownames(riskfactor)%in%rownames(proteomics0)),]
# construct main study
riskfactor_main<-riskfactor_main[rownames(proteomics0),]
main_study<-data.frame(riskfactor_main,proteomics0)

print(paste0("Main_study_case:all_case:",sum(riskfactor_main[,1]),":",sum(riskfactor$Y)))

family = "binomial"
summary_stat<-function(X,y,family = "binomial"){
    study_info_scaled<-list()
    Xonlylm<-glm(y~.,data = data.frame(y,X),family = family)
    study.m = list(Coeff=Xonlylm$coefficients[-1],
                   Covariance=vcov(Xonlylm)[-1,-1],Sample_size=nrow(X))
    study_info_scaled[[1]] <- study.m
    return(list("multilm"=Xonlylm,"sum_mul"=study_info_scaled))
}
suppressPackageStartupMessages({
    library(caret)
    library(glmnet)
    library(htlgmm)
    library(pROC)
})

set.seed(seed.use)

# hierarchical sampling to generate different folds
# downsample to keep all cases but less control for more efficient computation
allfoldids<-createFolds(main_study$Y,k = nfold)
allfoldids_down<-allfoldids
allfoldids_down2<-allfoldids
for(i in 1:nfold){
    curfold<-allfoldids_down[[i]]
    id1<-which(main_study$Y[curfold]==1)
    id1half<-sample(id1,round(length(id1)/2))
    id2<-sample(which(main_study$Y[curfold]==0),length(id1)*9)
    id2half<-sample(id2,round(length(id2)/2))
    ids<-c(id1,id2)
    idshalf<-c(id1half,id2half)
    curfold2<-curfold[idshalf]
    curfold<-curfold[ids]
    allfoldids_down[[i]]<-curfold
    allfoldids_down2[[i]]<-curfold2
}
allids_down<-Reduce(c,allfoldids_down)
allids_down2<-Reduce(c,allfoldids_down2)

time_start=Sys.time()

message("Generate summary data")
ext_study<-riskfactor[!rownames(riskfactor)%in%rownames(main_study),]
y_ext=ext_study[,1]

sum_res<-summary_stat(X = scale(ext_study[,-1]),
                      y=ext_study[,1])
multilm<-sum_res$multilm
sum_mul<-sum_res$sum_mul
beta_ext = sum_mul[[1]]$Coeff

auc_func<-function(y0,pred_y){
    suppressMessages(cur_auc<-c(auc(y0,c(htlgmm::expit(pred_y)),direction = "<")))
    cur_auc
}

nrisk=ncol(ext_study)-1
auc_list<-list()
for(fold in 1:nfold){
    print(paste0(diseasename,':Fold',fold,':seed',seed.use))
    test_index<-allfoldids[[fold]]
    # downsample
    test_index_down<-allfoldids_down[[fold]]
    train_index_down<-allids_down[!allids_down%in%test_index_down]
    # downsample half case
    test_index_down2<-allfoldids_down2[[fold]]
    train_index_down2<-allids_down2[!allids_down2%in%test_index_down2]
    
    main_test<-as.matrix(main_study[test_index,])
    
    main_train<-as.matrix(main_study[train_index_down,])
    main_train2<-as.matrix(main_study[train_index_down2,])

    runhtlgmm=function(name,main_train){
        main_y_train<-c(main_train[,1])
        main_X_train<-main_train[,-1]
        prefunc <- preProcess(main_X_train, method = c("center", "scale"))
        main_X_train = predict(prefunc, main_X_train)
        main_y_test = c(main_test[,1])
        main_X_test = as.matrix(predict(prefunc, main_test[,-1]))
        
        main_Z_train<-main_X_train[,1:nrisk]
        main_W_train<-main_X_train[,-c(1:nrisk)]
        
        pro_lasso<-cv.glmnet(x=main_X_train,
                             y=main_y_train,
                             family = family,
                             type.measure="auc")
        pred1<-predict(pro_lasso,newx=main_X_test,s="lambda.min")
        auc_lasso<-auc_func(main_y_test,pred1)
        
        pred1<-main_X_test[,1:nrisk]%*%beta_ext
        auc_ext<-auc_func(main_y_test,pred1)
        res<-c(auc_lasso,auc_ext)
        names(res)<-paste0(name,c("lasso","ext"))
        message(paste0(diseasename,"+",name))
        print(res)
        
        htf<-function(htname,res_lasso,initial_with_type){
            
            res_htlgmm_lasso_tuneC=cv.htlgmm(y = y_main,Z = main_Z_train,W = main_W_train,
                                              penalty_type = "lasso",
                                              initial_with_type="lasso",
                                              beta_initial = as.vector(coef(res_lasso,s="lambda.min")),
                                              inference = F,
                                              study_info = sum_mul,family = "binomial",
                                              use_sparseC = F,
                                              tune_weight = T,
                                              tune_weight_method="mshrink")
            
            auc_tuneC<-auc_func(main_y_test,main_X_test%*%res_htlgmm_lasso_tuneC$beta[-1])
            weightC = res_htlgmm_lasso_tuneC$weight_min
            res<-c(auc_tuneC,weightC)
            names(res)<-paste0(htname,c("tuneC","weight"))
            res
        }
        c(res,htf(paste0(name,"htl_"),pro_lasso))
    }
    res2<-runhtlgmm("half_",main_train2)
    print(round(res2,4))
    
    res<-runhtlgmm("",main_train)
    print(round(res,4))
    
    auc_list[[fold]]<-c(res,res2)
    finalfilename<-paste0(savepath,diseasename,"/ukb_protein_analysis_",useprs,"_seed",seed.use,".rds")
    saveRDS(auc_list,finalfilename)
}

print("---------------------------------")
print(round(Reduce("+",auc_list)/nfold,4))

time_diff = Sys.time() - time_start
hours_diff = as.numeric(time_diff, units = "hours")
print("hours")
print(hours_diff)



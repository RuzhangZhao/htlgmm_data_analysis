
allukb=readRDS("/dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/ukbCombined_032022_mortality_covid.rds")
# Exclude Dropout
exclude_sample<-read.table("/dcs05/legacy-dcl01-arking/data/UK_Biobank/active/pheno/withdraw17731_285_20230821.txt")$V1
allukb1<-allukb[which(!allukb$f.eid%in%exclude_sample),] # 502357
# Exclude Related Individuals 
# Data-Field 22020: Used in genetic principal components
allukb1<-allukb1[which(allukb1$f.22020=="Yes"),] # 407001

# Data-Field 31: Sex
# Data-Field 22001: Genetic sex
allukb1<-allukb1[which(allukb1$f.31.0.0==allukb1$f.22001.0.0),]

# Data-Field 21000: Ethnic background
allukb1<-allukb1[which(allukb1$f.21000.0.0 %in% c("White","British","Irish","Any other white background")),] # 381358

allukb2<-allukb1[allukb1$f.31.0.0 == "Female",]

# Data-Field 21022: age
# Data-Field 21001: BMI
# Data-Field 20116 Description: Smoking status
# Data-Field 20161 Description: Pack years of smoking
# Data-Field 1558 Description: Alcohol intake frequency.
index_riskf<-c("f.eid","f.21022.0.0","f.21001.0.0","f.50.0.0","f.20116.0.0","f.20161.0.0","f.1558.0.0")
names_riskf<-c("Y","age","BMI","Height","SmokingStatus","PackYrs","AlcoholFreq")
target_ukb<-allukb2[,index_riskf]
colnames(target_ukb)<-names_riskf
rownames(target_ukb)<-allukb2$f.eid

# Data-Field 20110 Description: Illnesses of mother

illness_mother<-allukb2[,which(grepl("f.20110.0.",colnames(allukb2)))]
target_ukb$mother_breast_cancer<-ifelse(rowSums(illness_mother=="Breast cancer",na.rm=T)>0,1,0)
illness_silbings<-allukb2[,which(grepl("f.20111.0.",colnames(allukb2)))]
target_ukb$silbings_breast_cancer<-ifelse(rowSums(illness_silbings=="Breast cancer",na.rm=T)>0,1,0)


supp_ukb<-allukb2[,c("f.2724.0.0","f.2734.0.0","f.2714.0.0","f.2784.0.0","f.2814.0.0","f.2834.0.0","f.2774.0.0")]
names(supp_ukb)<-c("menopause","Numbirth","age_menarche","oral_contraceptive","HRT","Bilateral_oophorectomy","stillbirth")
target_ukb<-cbind(target_ukb,supp_ukb)

target_ukb$prv_bbd<-ifelse(rownames(target_ukb)%in%bbd_prv_eid,1,0)

# case definition 
target_ukb$Y<-ifelse(rownames(target_ukb)%in%bc_inc_eid,1,0)
target_ukb<-target_ukb[which(!rownames(target_ukb)%in%bc_prv_eid),]

## filtering or imputation 
target_ukb<-target_ukb[rowSums(is.na(target_ukb[,c("age","BMI","Height","SmokingStatus","AlcoholFreq")]))==0,]
target_ukb<-target_ukb[target_ukb$SmokingStatus!="Prefer not to answer",]
target_ukb<-target_ukb[target_ukb$AlcoholFreq!="Prefer not to answer",]

target_ukb<-target_ukb[(target_ukb$menopause == "Yes")|(target_ukb$age>50),]

# Smoking imputation 
smoke<-rep(0,length(target_ukb$SmokingStatus))
smoke[target_ukb$SmokingStatus == "Previous"] = 1
smoke[target_ukb$SmokingStatus == "Current"] = 2
target_ukb$SmokingStatus<-ifelse(smoke==1,1,0)
target_ukb$SmokingStatus_Current<-ifelse(smoke==2,1,0)

target_ukb$PackYrs[which(is.na(target_ukb$PackYrs)&smoke == 0)]=0
target_ukb$PackYrs[which((is.na(target_ukb$PackYrs)&smoke == 1))]<-mean(target_ukb$PackYrs[which((!is.na(target_ukb$PackYrs)&smoke == 1))])
target_ukb$PackYrs[which((is.na(target_ukb$PackYrs)&smoke == 2))]<-mean(target_ukb$PackYrs[which((!is.na(target_ukb$PackYrs)&smoke == 2))])


alchol<-rep(1,length(target_ukb$AlcoholFreq))
alchol[target_ukb$AlcoholFreq == "Three or four times a week"] = 2 
alchol[target_ukb$AlcoholFreq == "Once or twice a week"] = 3
alchol[target_ukb$AlcoholFreq == "One to three times a month"] = 4
alchol[target_ukb$AlcoholFreq == "Special occasions only"] = 5
alchol[target_ukb$AlcoholFreq == "Never"] = 6
target_ukb$AlcoholFreq<-ifelse(alchol==1,1,0)
target_ukb$AlcoholFreq2<-ifelse(alchol==2,1,0)
target_ukb$AlcoholFreq3<-ifelse(alchol==3,1,0)
target_ukb$AlcoholFreq4<-ifelse(alchol==4,1,0)
target_ukb$AlcoholFreq5<-ifelse(alchol==5,1,0)

target_ukb<-target_ukb[!is.na(target_ukb$Numbirth),]
target_ukb$Numbirth[target_ukb$Numbirth>5]<-5

target_ukb<-target_ukb[target_ukb$HRT%in%c("Yes","No"),]
target_ukb$HRT<-ifelse(target_ukb$HRT=='Yes',1,0)
target_ukb<-target_ukb[target_ukb$Bilateral_oophorectomy%in%c("Yes","No"),]
target_ukb$Bilateral_oophorectomy<-ifelse(target_ukb$Bilateral_oophorectomy=='Yes',1,0)
target_ukb<-target_ukb[target_ukb$stillbirth%in%c("Yes","No"),]
target_ukb$stillbirth<-ifelse(target_ukb$stillbirth=='Yes',1,0)
target_ukb<-target_ukb[target_ukb$oral_contraceptive%in%c("Yes","No"),]
target_ukb$oral_contraceptive<-ifelse(target_ukb$oral_contraceptive=='Yes',1,0)
target_ukb<-target_ukb[target_ukb$age_menarche>0,]

target_ukb<-target_ukb[,-which(colnames(target_ukb)%in%c("hysterectomy","menopause"))]




colMeans(is.na(target_ukb))


for(i in 1:ncol(target_ukb)){
    if(sum(is.na(target_ukb[,i]))>0){
        target_ukb[is.na(target_ukb[,i]),i]<-mean(target_ukb[!is.na(target_ukb[,i]),i])
    }
}

target_ukb<-target_ukb[rowSums(is.na(target_ukb))==0,]
table(target_ukb$Y)

diseasename="Breast"
saveRDS(target_ukb,paste0("/dcs04/nilanjan/data/rzhao/UKB/",diseasename,"_risk_zrz.rds"))

## prs
library(data.table)
prsscore<-fread("/dcs05/legacy-dcl01-arking/data/UK_Biobank/active/chatterjeelab_UKB/ruzhang_122123/standard_prs_040124_participant.tsv")
target_ukb<-target_ukb[rownames(target_ukb)%in%prsscore$eid,]
prsscore<-prsscore[match(rownames(target_ukb),prsscore$eid),]
target_ukb<-cbind(target_ukb,prsscore$p26220)
colnames(target_ukb)[ncol(target_ukb)]<-"prs"
target_ukb<-target_ukb[!is.na(target_ukb$prs),]
# 26210  Standard PRS for asthma (AST)
# 26218  Standard PRS for bowel cancer (CRC)
# 26220  Standard PRS for breast cancer (BC)
# 26223  Standard PRS for cardiovascular disease (CVD)
# 26248  Standard PRS for ischaemic stroke (ISS)

diseasename="Breast"
saveRDS(target_ukb,paste0("/dcs04/nilanjan/data/rzhao/UKB/",diseasename,"_risk_zrzprs.rds"))

library(data.table)
library(tidyverse)
# For the cancer, the column is changed to 40005/40006

allukb=readRDS("/dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Phenotype/Downloads/ukbCombined_032022_mortality_covid.rds")
cancer_diag_time=allukb[,c("f.eid","f.53.0.0",colnames(allukb)[grepl("f.40005.",colnames(allukb))])]
cancer_diag=allukb[,c("f.eid",colnames(allukb)[grepl("f.40006.",colnames(allukb))])]
colnames(cancer_diag)=c("eid",0:(ncol(cancer_diag)-2))
colnames(cancer_diag_time)=c("eid","p53_i0", 0:(ncol(cancer_diag)-2))

# only difference between diag and diag_time columns:
# the 2nd column of diag_time is entry time
# 2:end of diag are corresponding to 3:end of diag_time
cancer_diag_long = as.data.frame(pivot_longer(cancer_diag,cols = 2:ncol(cancer_diag),names_to = "diagnosis_instance",values_to = "diagnose"))
cancer_diag_time_long = as.data.frame(pivot_longer(cancer_diag_time,cols = 3:ncol(cancer_diag_time),names_to = "diagnosis_instance",values_to = "diag_time"))
rm(cancer_diag,cancer_diag_time)
ids=which(!is.na(cancer_diag_time_long$diag_time)) # merge the two
cancer_diag_long = cancer_diag_long[ids,]
cancer_diag_time_long = cancer_diag_time_long[ids,]
cancer_diag_merge = cbind(cancer_diag_long,cancer_diag_time_long[,c(2,4)])
cancer_diag_merge$instance_type = ifelse(cancer_diag_merge$p53_i0 > cancer_diag_merge$diag_time,'prevalent','incident')
cancer_diag_merge = cancer_diag_merge[!is.na(cancer_diag_merge$instance_type),]


## Breast Cancer
diag_target = cancer_diag_merge[grepl('C50|D05',cancer_diag_merge$diagnose),]
diag_target$diag = 'breast_cancer'

## add in self-report
self_entry = allukb[,c('f.eid',paste0('f.20001.0.',0:5))]
self_entry_eid=self_entry$f.eid[rowSums(self_entry==1002,na.rm=T)>0]
diag_target$self_entry = ifelse(diag_target$eid %in% self_entry_eid,'Yes','No')

# 1. not in self report at entry time AND not in prevalent cases
# 2. in self report at follow time OR in incident cases
bc_prv_eid=unique(diag_target$eid[diag_target$self_entry=="Yes"|
                                       diag_target$instance_type=="prevalent"])

bc_inc_eid=unique(diag_target$eid[diag_target$instance_type=="incident"&
                                       grepl('C50',diag_target$diagnose)&
                                       !(diag_target$eid%in%bc_prv_eid)])


## Colorectal Cancer
diag_target = cancer_diag_merge[grepl('C18|C19|C20',cancer_diag_merge$diagnose),]
diag_target$diag = 'colorectal_cancer'

## add in self-report
self_entry = allukb[,c('f.eid',paste0('f.20001.0.',0:5))]
self_entry_eid=self_entry$f.eid[rowSums(self_entry%in%c(1020,1022,1023),na.rm=T)>0]
diag_target$self_entry = ifelse(diag_target$eid %in% self_entry_eid,'Yes','No')

# 1. not in self report at entry time AND not in prevalent cases
# 2. in self report at follow time OR in incident cases
colon_prv_eid=unique(diag_target$eid[diag_target$self_entry=="Yes"|
                                      diag_target$instance_type=="prevalent"])

colon_inc_eid=unique(diag_target$eid[diag_target$instance_type=="incident"&
                                      grepl('C18|C19|C20',diag_target$diagnose)&
                                      !(diag_target$eid%in%colon_prv_eid)])





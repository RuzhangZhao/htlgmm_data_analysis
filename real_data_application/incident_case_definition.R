library(data.table)
library(tidyverse)

allukb=readRDS("/dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Phenotype/Downloads/ukbCombined_032022_mortality_covid.rds")
diag_time=allukb[,c("f.eid","f.53.0.0",colnames(allukb)[grepl("f.41280.",colnames(allukb))])]
diag=allukb[,c("f.eid",colnames(allukb)[grepl("f.41270.",colnames(allukb))])]
colnames(diag)=c("eid",0:(ncol(diag)-2))
colnames(diag_time)=c("eid","p53_i0", 0:(ncol(diag)-2))

# only difference between diag and diag_time columns:
# the 2nd column of diag_time is entry time
# 2:end of diag are corresponding to 3:end of diag_time
diag_long = as.data.frame(pivot_longer(diag,cols = 2:ncol(diag),names_to = "diagnosis_instance",values_to = "diagnose"))
diag_time_long = as.data.frame(pivot_longer(diag_time,cols = 3:ncol(diag_time),names_to = "diagnosis_instance",values_to = "diag_time"))
rm(diag,diag_time)
ids=which(!is.na(diag_time_long$diag_time)) # merge the two
diag_long = diag_long[ids,]
diag_time_long = diag_time_long[ids,]
diag_merge = cbind(diag_long,diag_time_long[,c(2,4)])
diag_merge$instance_type = ifelse(diag_merge$p53_i0 > diag_merge$diag_time,'prevalent','incident')
diag_merge = diag_merge[!is.na(diag_merge$instance_type),]

## Stroke
diag_target = diag_merge[grepl('I63|I64|I61|I60',diag_merge$diagnose),]
diag_target$diag = 'stroke'

## add in self-report
self = allukb[,c('f.eid','f.4056.0.0')]
self_entry = self[,c(1,2)]
self_entry_eid=self_entry$f.eid[!is.na(self[,2])]

#Data-Field 6150 Description: Vascular/heart problems diagnosed by doctor
self = allukb[,which(grepl("f.6150.0.",colnames(allukb)))]
self_entry_eid=c(self_entry_eid,
                 allukb$f.eid[rowSums(self=='Stroke',na.rm=T)>0])

self = allukb[,which(grepl("f.20002.0.",colnames(allukb)))]
# 1081 stroke 
# 1082 transient ischaemic attack (tia)
# 1083 subdural haemorrhage/haematoma
# 1086 subarachnoid haemorrhage
self_entry_eid=c(self_entry_eid,
                 allukb$f.eid[rowSums(sapply(c(1081,1082,1083,1086,1583), function(i){
                     rowSums(self==i,na.rm=T)}))>0])
diag_target$self_entry = ifelse(diag_target$eid %in% self_entry_eid,'Yes','No')

stroke_prv_eid=unique(diag_target$eid[diag_target$self_entry=="Yes"|
                                          diag_target$instance_type=="prevalent"])
stroke_inc_eid=unique(diag_target$eid[diag_target$instance_type=="incident"&
                                          grepl('I63',diag_target$diagnose)&
                                          !(diag_target$eid%in%stroke_prv_eid)])

## Asthma
diag_target = diag_merge[grepl('J45|J46|J30|L20',diag_merge$diagnose),]
diag_target$diag = 'asthma'

## add in self-report

self = allukb[,c('f.eid','f.6152.0.0','f.6152.0.1','f.6152.0.2','f.6152.0.3','f.6152.0.4')]
self_entry_eid=self$f.eid[rowSums(self=='Asthma',na.rm=T)>0]
self_entry_eid=c(self_entry_eid,self$f.eid[rowSums(self=='Hayfever, allergic rhinitis or eczema',na.rm=T)>0])

self = allukb[,c('f.eid','f.22127.0.0')]
self_entry_eid=c(self_entry_eid,self$f.eid[rowSums(self=='Yes',na.rm=T)>0])

self = allukb[,which(grepl("f.20002.0.",colnames(allukb)))]
self_entry_eid=c(self_entry_eid,allukb$f.eid[rowSums(self=="1111",na.rm = T)>0])

diag_target$self_entry = ifelse(diag_target$eid %in% self_entry_eid,'Yes','No')

asthma_prv_eid=unique(diag_target$eid[diag_target$self_entry=="Yes"|
                                          diag_target$instance_type=="prevalent"])
asthma_inc_eid=unique(diag_target$eid[diag_target$instance_type=="incident"&
                                          grepl('J45|J46',diag_target$diagnose)&
                                          !(diag_target$eid%in%asthma_prv_eid)])
length(asthma_inc_eid)


# Atrial fibrillation and flutter

diag_target = diag_merge[grepl('I48',diag_merge$diagnose),]
diag_target$diag = 'atrial_fibrillation'
self = allukb[,which(grepl("f.20002.0.",colnames(allukb)))]
self_entry_eid=allukb$f.eid[rowSums(self=="1471",na.rm = T)>0]
diag_target$self_entry = ifelse(diag_target$eid %in% self_entry_eid,'Yes','No')

af_prv_eid=unique(diag_target$eid[diag_target$self_entry=="Yes"|
                                      diag_target$instance_type=="prevalent"])

# Hypertension
diag_target = diag_merge[grepl('I10|I15',diag_merge$diagnose),]
diag_target$diag = 'hypertension'

#Data-Field 6150 Description: Vascular/heart problems diagnosed by doctor
self = allukb[,which(grepl("f.6150.0.",colnames(allukb)))]
self_entry_eid=c(self_entry_eid,
                 allukb$f.eid[rowSums(self=='High blood pressure',na.rm=T)>0])

self = allukb[,which(grepl("f.20002.0.",colnames(allukb)))]
self_entry_eid=c(self_entry_eid,
                 allukb$f.eid[rowSums(sapply(c(1065,1072), function(i){
                     rowSums(self==i,na.rm=T)}))>0])
diag_target$self_entry = ifelse(diag_target$eid %in% self_entry_eid,'Yes','No')

hypertension_prv_eid=unique(diag_target$eid[diag_target$self_entry=="Yes"|
                                                diag_target$instance_type=="prevalent"])


# CVD
diag_target = diag_merge[grepl('I20|I21|I22|I23|I24|I25',diag_merge$diagnose),]
diag_target$diag = 'cvd'

#Data-Field: 3894 Description: Age heart attack diagnosed
#Data-Field: 3627 Description: Age angina diagnosed
self = allukb[,c('f.eid','f.3894.0.0','f.3627.0.0')]
self_entry = self[,c(1,2)]
self_entry_eid=self_entry$f.eid[which(self_entry[,2]>0)]

#Data-Field 6150 Description: Vascular/heart problems diagnosed by doctor
self = allukb[,which(grepl("f.6150.0.",colnames(allukb)))]
self_entry_eid=c(self_entry_eid,
                 allukb$f.eid[rowSums(self=='Heart attack',na.rm=T)>0])
self_entry_eid=c(self_entry_eid,
                 allukb$f.eid[rowSums(self=='Angina',na.rm=T)>0])

self = allukb[,which(grepl("f.20002.0.",colnames(allukb)))]
self_entry_eid=c(self_entry_eid,
                 allukb$f.eid[rowSums(sapply(c(1066,1075,1076), function(i){
                     rowSums(self==i,na.rm=T)}))>0])
diag_target$self_entry = ifelse(diag_target$eid %in% self_entry_eid,'Yes','No')

cvd_prv_eid=unique(diag_target$eid[diag_target$self_entry=="Yes"|
                                       diag_target$instance_type=="prevalent"])

cvd_inc_eid=unique(diag_target$eid[diag_target$instance_type=="incident"&
                                       grepl('I21|I22|I23|I24|I252|I251|I255|I256|I257|I258|I259',diag_target$diagnose)&
                                       !(diag_target$eid%in%cvd_prv_eid)])

## crohn

diag_target = diag_merge[grepl('K50|K51',diag_merge$diagnose),]
diag_target$diag = 'crohn'
self = allukb[,which(grepl("f.20002.0.",colnames(allukb)))]
self_entry_eid=allukb$f.eid[rowSums(self==1462,na.rm=T)>0]
self_entry_eid=c(self_entry_eid,allukb$f.eid[rowSums(self==1463,na.rm=T)>0])
diag_target$self_entry = ifelse(diag_target$eid %in% self_entry_eid,'Yes','No')
crohn_prv_eid=unique(diag_target$eid[diag_target$self_entry=="Yes"|
                                         diag_target$instance_type=="prevalent"])
length(crohn_prv_eid)

## colon polyps
diag_target = diag_merge[grepl('K635',diag_merge$diagnose),]
diag_target$diag = 'polyps'
self = allukb[,which(grepl("f.20002.0.",colnames(allukb)))]
self_entry_eid=allukb$f.eid[rowSums(self==1460,na.rm=T)>0]
diag_target$self_entry = ifelse(diag_target$eid %in% self_entry_eid,'Yes','No')
polyps_prv_eid=unique(diag_target$eid[diag_target$self_entry=="Yes"|
                                          diag_target$instance_type=="prevalent"])
length(polyps_prv_eid)

## benign breast disease 
diag_target = diag_merge[grepl('D24|N60',diag_merge$diagnose),]
diag_target$diag = 'benignbreast'
self = allukb[,which(grepl("f.20002.0.",colnames(allukb)))]
self_entry_eid=allukb$f.eid[rowSums(sapply(c(1364,1366,1367,1560,1666), function(i){
    rowSums(self==i,na.rm=T)}))>0]
diag_target$self_entry = ifelse(diag_target$eid %in% self_entry_eid,'Yes','No')
bbd_prv_eid=unique(diag_target$eid[diag_target$self_entry=="Yes"|
                                       diag_target$instance_type=="prevalent"])
length(bbd_prv_eid)

##数据来源
##https://ascopubs.org/doi/suppl/10.1200/CCI.17.00096

OV_drug <- read.table(file = './data/OV/OV_drug.txt',header = T,sep = '\t')

OV_drug <- subset(OV_drug,subset = Chemotherapy.number.of.lines.of.therapy>0)
OV_drug$OS.time <- as.numeric(OV_drug$total_days_overall_survival)/29
OV_drug$FFI.time <- as.numeric(OV_drug$Days.off.platinum.prior.to.recurrence.1st.line)/29
OV_drug$FFI.time[OV_drug$FFI.time<0]=0

OV_drug$drug.name <- OV_drug$X1st_regimen
OV_drug$primary_therapy_outcome <- OV_drug$primary_therapy_outcome_success
OV_drug$FFI.event <- OV_drug$X1st_chemo_regimen_clinical_benefit_days_outcome
OV_drug$OS.event <- as.numeric(OV_drug$outcome_overall_survival_censoring)

DFS <- (OV_drug$days_to_tumor_progression !='null') | (OV_drug$days_to_tumor_recurrence !='null')
OV_drug$DFS.event <- ifelse(DFS,1,0)

progres_days <- ifelse(OV_drug$days_to_tumor_progression=='null','0',OV_drug$days_to_tumor_progression)
progres_days <- as.numeric(progres_days)
recurence_dys <- ifelse(OV_drug$days_to_tumor_recurrence=='null','0',OV_drug$days_to_tumor_recurrence)
recurence_dys <- as.numeric(recurence_dys)
last_FollowUP <- as.numeric(OV_drug$days_to_last_followup)
re_pro <- cbind(progres_days,recurence_dys,last_FollowUP)
OV_drug$DFS.time <- apply(re_pro,1,function(x){
  if(max(x[1:2])==0){
    return(x[3])
  }else{
    return(max(x[1:2]))
  }
})
OV_drug$DFS.time <- OV_drug$DFS.time/29

OV_drug$FFI.event <- as.character(OV_drug$FFI.event)
ggdensity(OV_drug,x = 'FFI.time',fill = 'FFI.event')
OV_drug$FFI.event <- as.numeric(OV_drug$FFI.event)



OV_drug$platinum_response <- ifelse(OV_drug$FFI.event==0 & OV_drug$primary_therapy_outcome=='COMPLETE RESPONSE','Sensitive',
                                    ifelse(OV_drug$FFI.event==0 & OV_drug$primary_therapy_outcome !='COMPLETE RESPONSE','unknow',
                                    ifelse(OV_drug$FFI.event==1 & OV_drug$FFI.time > 7.67,'Sensitive','Resistant')      
                                           )
                                    )
OV_drug$platinum_response <- ifelse(OV_drug$primary_therapy_outcome=='COMPLETE RESPONSE','Sensitive',
                                    ifelse(OV_drug$primary_therapy_outcome !='COMPLETE RESPONSE' & OV_drug$FFI.time > 7.67,'Sensitive','Resistant')
)

OV_drug$lastFollowup <- as.numeric(OV_drug$days_to_last_followup)/29

max(as.numeric(OV_drug[which(OV_drug$FFI.event==0 & OV_drug$FFI.time<8),]$days_to_last_followup))

cliData_drug <- OV_drug[,c('bcr_patient_barcode','outcome_overall_survival_censoring','age_at_initial_pathologic_diagnosis','anatomic_organ_subdivision',
                           'histological_type','person_neoplasm_cancer_status','primary_therapy_outcome','race',
                           'site_of_tumor_first_recurrence','tumor_grade','tumor_stage','tumor_residual_disease',
                           'OS.time','FFI.time','drug.name','FFI.event','DFS.event','DFS.time','platinum_response')]
names(cliData_drug)[c(1,2,3)] <- c('patient_id','OS.event','age')
cliData_drug$OS.event <- as.numeric(cliData_drug$OS.event)
cliData_drug$age <- as.numeric(cliData_drug$age)

save(cliData_drug,file = './data/bulk/cliData_drug.RData')

####
load('./data/Pancancer/CN_score_merge.RData')
#load('./data/Pancancer/clinicData_survival.RData')
load(file = './data/OV/cliData_drug.RData')
load('./data/Pancancer/subtypes.RData')

subtypes_OV <- as.data.frame(subset(subtypes,cancer.type=='OVCA'))
subtypes_OV$pan.samplesID <- substr(x = subtypes_OV$pan.samplesID,start = 1,stop = 12)
subtypes_OV <- subtypes_OV[,c('pan.samplesID','Subtype_mRNA')]





library(survival)
clinical_merge <- merge.data.frame(cliData_drug,CN_score_merge,by.x = 'patient_id',by.y = 'patient_barcode',all.x = T)
dataMerge_subset <-  clinical_merge
dataMerge_subset$HRD_status <- ifelse(dataMerge_subset$HRD>=41,'HR_deficiency','HR_proficiency')

load('./data/OV/Cosmic_signature.RData')
#Cosmic_signature <- t(Cosmic_signature)
Cosmic_signature <- data.frame(sampleID=rownames(Cosmic_signature),Cosmic_signature)
Cosmic_signature$sampleID <- gsub(pattern = '\\.',replacement = '-',x = Cosmic_signature$sampleID)
Cosmic_signature$sampleID <- substr(x = Cosmic_signature$sampleID,start = 1,stop = 15)
dataMerge_subset <- merge.data.frame(dataMerge_subset,Cosmic_signature,by.x = 'aliquot_barcode',by.y = 'sampleID',all.x = T)

data_explor.continuous2continuous(dataMerge_subset,x_dimension = 'HRD',y_dimension = 'COSMIC_3',plot = T,linear_fitting = T)





#dataMerge_subset$HRD_status <- ifelse(dataMerge_subset$HRD>=63,'HR_deficiency',ifelse(dataMerge_subset$HRD<31,'HR_proficiency','middle'))


dataMerge_subset <- merge.data.frame(dataMerge_subset,subtypes_OV,by.x = 'patient_id',by.y = 'pan.samplesID',all.x = T)

ESTIMATE_OV <- read.table(file = './data/RawData/ESTIMATE_OV.txt',header = T,sep = '\t')
ESTIMATE_OV$class <- substr(ESTIMATE_OV$ID,start = 14,stop = 15)
ESTIMATE_OV$ID <- substr(ESTIMATE_OV$ID,start = 1,stop = 15)

dataMerge_subset <- merge.data.frame(dataMerge_subset,ESTIMATE_OV,by.x = 'aliquot_barcode',by.y = 'ID',all.x = T)
table(dataMerge_subset$tumor_grade)
dataMerge_subset$age <- as.numeric(dataMerge_subset$age)

dataMerge_subset$Non.silent.per.Mb <- ifelse(dataMerge_subset$Non.silent.per.Mb>10,10,dataMerge_subset$Non.silent.per.Mb)

dataMerge_subset$TMB <- dataMerge_subset$Non.silent.per.Mb


save(dataMerge_subset,file='./data/bulk/dataMerge_subset.RData')
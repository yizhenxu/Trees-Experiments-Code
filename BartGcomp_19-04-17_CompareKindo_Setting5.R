### Diabetic Readmission Data
library(data.table)
#datpath = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\MPBART Computing\\UCI Diabetes Readmission"
#d = fread(paste0(datpath,"\\diabetic_data.csv"))
d = fread("diabetic_data.csv")

#website: https://www.hindawi.com/journals/bmri/2014/781670/
#variable summary: https://www.hindawi.com/journals/bmri/2014/781670/tab1/
head(d)

table(d$Y)
d[d=="?"] = NA

### Age as numeric 0~9
#return those that matches (content); 
#".?" means the single character represented by "." appears at most 1 time; 
#"\\-{1}" hyphen appears exactly 1 time;
#"\\d+" numbers, appear at least 1 time;
d$age = gsub("\\[(\\d).?\\-{1}\\d+\\)","\\1",d$age); d$age = as.numeric(d$age)
unique(d$age)

### Weight as numeric 
i1 = which(!is.na(d$weight) & d$weight!=">200")
d$weight[i1] = gsub("\\[(\\d*)\\-{1}\\d+\\)","\\1",d$weight[i1])
d$weight[d$weight == ">200"] = "200"
d$weight = as.numeric(d$weight)
unique(d$weight)

### Admission Type - unchanged 
tmp = table(d$readmitted,d$admission_type_id,useNA = "always")
for(j in 1:ncol(tmp)) tmp[,j] = tmp[,j]/sum(tmp[,j])
tmp
table(d$admission_type_id,useNA = "always")

### Discharge Disposition

table(d$discharge_disposition_id,useNA = "always")
table(d$readmitted,d$discharge_disposition_id,useNA = "always")
tmp = table(d$readmitted,d$discharge_disposition_id)
for(j in 1:ncol(tmp)) tmp[,j] = tmp[,j]/sum(tmp[,j])
tmp

a = d$discharge_disposition_id
a[a %in% 1] = 1 # home wo care
a[a %in% c(2,3,4,24,27)] = 2 # care facility
a[a==5] = 3# other inpatient care
a[a %in% c(6,7,8)] = 4# home w care
a[a==9] = 5# here impatient care 
a[a==10] = 6# newborn?
a[a==11] = 7# expired (all readmission >30 days)
a[a %in% c(13,14,19,20,21)] = 8 # hospice
a[a %in% c(12,16,17)] = 9 # outpatient
a[a==15] = 10 # swing bed
a[a==22] = 11 # drug rehab
a[a==28] = 12 # psychiatric
a[a==18] = 13 # NULL
a[a==23] = 14 # long term care hospital
a[a==25] = 15 # Not mapped
table(a,useNA = "always")
d$discharge_disposition_id = a

### admission_source_id

table(d$admission_source_id,useNA = "always")
tmp = table(d$readmitted,d$admission_source_id,useNA = "always")
for(j in 1:ncol(tmp)) tmp[,j] = tmp[,j]/sum(tmp[,j])
tmp

a = d$admission_source_id
a[a %in% c(1:3)] = 1 # referral
a[a %in% c(4:6,10,22,25)] = 2 # transferred
a[a %in% c(11:14,23,24)] = 3 # baby
a[a %in% c(8,9,15,17,21)] = 4 # NA
a[a==7] = 5 # emergency
a[a == 20] =6 # not mapped

table(a,useNA = "always")
d$admission_source_id = a

### Payer code
a = d$payer_code
table(a,useNA = "always")
a[a=="FR"] = NA
d$payer_code = a

### Medical specialty -> vector of nl
a = d$medical_specialty
tmp = table(a,useNA = "always")
sort(tmp)
nl = c("Neuro","Surg","Ped","Cardio","Onc","Ortho","Emerg","Inter","General","Psy","Uro","Neph","Radio","Pulm","Endo","Gas")
for(i in 1:length(nl)){
  vn = nl[i]
  alln = grep(vn,names(tmp), perl=TRUE, value=TRUE)
  d[,(vn) := 1*(a %in% alln)]
}

d[, medical_specialty := NULL]

### Diagnosis (Primary, Secondary, Additional Secondary)
#The following abbreviations are used for particular icd9 codes:
#circulatory for icd9: 390-459, 785, 
#digestive for icd9: 520-579, 787, 
#genitourinary for icd9: 580-629, 788, 
#diabetes for icd9: 250.xx, 
#injury for icd9: 800-999, 
#musculoskeletal for icd9: 710-739, 
#neoplasms for icd9: 140-239, 
#respiratory for icd9: 460-519, 786, and 
#other for otherwise.
dl = c()
for(vn in paste0("diag_",1:3)){
  a = d[, get(vn)]
  d[,(paste0(vn,"_1")):= 1*(a %in% c(390:459, 785))] #circulatory
  d[,(paste0(vn,"_2")):= 1*(a %in% c(520:579, 787))] #digestive
  d[,(paste0(vn,"_3")):= 1*(a %in% c(580:629, 788))] #genitourinary
  ll = grep("250",unique(a), perl=TRUE, value=TRUE)
  d[,(paste0(vn,"_4")):= 1*(a %in% ll)] #diabetes
  d[,(paste0(vn,"_5")):= 1*(a %in% 800:999)] #injury
  d[,(paste0(vn,"_6")):= 1*(a %in% 710:739)] #musculoskeletal
  d[,(paste0(vn,"_7")):= 1*(a %in% 140:239)] #neoplasms
  d[,(paste0(vn,"_8")):= 1*(a %in% c(460:519, 786))] #respiratory
  dl = c(dl, paste0(vn,"_",1:8))
}

d[,(paste0("diag_",1:3)) := NULL]

if(0){ # check the created diag vars, diag_i_j
  colnames(d)
  for(j in 51:74){
    a = d[,j,with=F]
    print(table(a))
  }
}

### Glucose serum test result
a = d$max_glu_serum
table(a,useNA="always")
a[a=="None"] = 0
a[a=="Norm"] = 1
a[a==">200"] = 2
a[a==">300"] = 3
d[, max_glu_serum := as.numeric(a)]
d[, max_glu_serum_NA := 1*(max_glu_serum==0) ]

### A1c test result
a = d$A1Cresult
table(a,useNA="always")
a[a=="None"] = 0
a[a=="Norm"] = 1
a[a==">7"] = 2
a[a==">8"] = 3
d[, A1Cresult := as.numeric(a)]
d[, A1Cresult_NA := 1*(A1Cresult==0) ]

colnames(d) = gsub("-", "", colnames(d), fixed = TRUE)
### Create indicators from all categorical covariates
factorize =  c("gender","race","payer_code","admission_type_id","discharge_disposition_id",
               "admission_source_id","payer_code","change",
               "diabetesMed", "metformin", "repaglinide","nateglinide",
               "chlorpropamide","glimepiride", "acetohexamide",
               "glipizide","glyburide","tolbutamide","pioglitazone","rosiglitazone",
               "acarbose","miglitol", "troglitazone","tolazamide" ,"examide",
               "citoglipton", "insulin","glyburide-metformin",
               "glipizide-metformin","glimepiride-pioglitazone", "metformin-rosiglitazone",
               "metformin-pioglitazone" )
factorize = gsub("-", "", factorize , fixed = TRUE)
#factorize[which(!(factorize %in% colnames(d)))]
dd = factorizedata(d, factorize)
dd = createind(dd,c=100) #abandon binary levels of less than c one's

numericX = c("age","weight","time_in_hospital","num_lab_procedures","num_procedures", 
             "num_medications", "number_outpatient","number_emergency",
             "number_inpatient", "number_diagnoses")
dd = createmis(dd,numericX)

#check = apply(dd,2,function(x) length(unique(x)))
#sum(check==1)

yind = which(colnames(dd)=="readmitted")
colnames(dd)[yind] = "y"

a = dd$y; a[a=="NO"] = 1; a[a=="<30"] = 2; a[a==">30"] = 3; dd$y = as.numeric(a)

nd = colnames(dd); xvec= nd[-which(nd %in% c("encounter_id","patient_nbr","readmitted"))]
fml = paste(c("y ~ ",paste0(xvec[-length(xvec)]," + "),last(xvec)),collapse = "")
show(fml)
dim(dd)

set.seed(8)
rnd = sample(1:nrow(dd)); ntrd = ceiling(nrow(dd)*0.7)
trd = dd[rnd[1:ntrd],]; ted = dd[rnd[(ntrd+1):nrow(dd)],]





### Multinomial Probit BART Settings

set.seed(99)
nd = 2000 # posterior draws
nb = 500 # burn in
nt = 100 # number of trees
p=length(unique(trd$y)) # number of outcome levels

# seed to start with
CatSeed = 888888



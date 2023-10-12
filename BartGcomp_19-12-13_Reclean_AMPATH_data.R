path = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART"
path = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART"
path = "/Volumes/GoogleDrive/My Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART"

setwd(paste0(path,"/makeGridData/my working datasets"))

library(readstata13)
library(lubridate)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(data.table)
library(zoo)
library(nnet)

source(paste0(path,"/BartGcomp_19-12-13_Reclean_AMPATH_data_Functions.R"))

### dt1 cleaning steps in 17-11-09_selection cleaning lay grid R script
### 5a. carry forward measurement before adult initial to the adult initial record for baseline variables
### 5b. carry backward nearest measurement within a month after adult initial to the adult initial record
### 14. define initial height, weight, whostage, married, age
### 15. lay grid of 200 days; *do this in ATHENA due to memory limit*
### Delete those records after dod grid in grid laying

# merged data of demog and visits, done in STATA, with abnormal death eliminated
dat = read.dta13("./denom_all.dta")

colnames(dat)
var_enrol = c("state_province", "initclinic", 
              "male","yrenrol","ageatenroll")
var_time  = c("arvstart","apptdate","returnvisitdate","death","transfer", 
              "numadmeas","arvperfectad", "clinic","weight","height",  "cd4", 
              "whostage", "onarv","married", "age", "viralload", 
              "visittype","lost2fup" )

dat = dat[,c("patient_id","dod","enroldate","hivstatus",var_enrol,var_time)] # rearrange the order of the variables

dt = setDT(dat)

# load 2 csv files
invisit = fread("hivvisit2_ET.csv")
indemog = fread("hivdemog2_travel&tribe.csv")

# merge travel time and tribe
setkey(indemog, patient_id)
setkey(dt, patient_id)
dt = merge(dt, indemog, all.x = TRUE)

# merge encounter_type
invisit$apptdate = as.Date(invisit$apptdate, "%d%b%Y")
dt[,tmpid:=paste0(patient_id,"_",apptdate)]
invisit[,tmpid:=paste0(patient_id,"_",apptdate)]
invisit[, patient_id:= NULL]
invisit[, apptdate:=NULL]
setkey(dt,tmpid)
setkey(invisit, tmpid)
dt = merge(dt, invisit, all.x = TRUE)

# define record number 
dt[, num := 1:.N, by=patient_id]

### 0. order by patient id and record index TO REDUCE POSSIBLE MISTAKES IN THE FOLLOWING CLEANING PROCESS
dt = dt[order(patient_id,num)]

save(dt, file = "AMPATH_12-13-2019_dt.RData")
#load("AMPATH_12-13-2019_dt.RData")
#dc = copy(dt); remove(invisit, indemog, dat)
length(unique(dt$patient_id))#96045 the very raw
dim(dt)#1687415 

### 1. add "AMPATH transfer form" to the transfer variable 2: 1 for transfer within ampath, 2 for transfer out of ampath
dt[, lastnum := num[.N], by = patient_id]
see = dt[encounter_type =="AMPATH transfer form",] # either the last or second to the last record, hence imply xfer OUT
dt$transfer[dt$encounter_type =="AMPATH transfer form"] = 2

### 2. set extreme values of weight, height, and CD4 to NA; use within-patient median for processing weight and height

dt$cd4[(dt$cd4 > 2000) & !is.na(dt$cd4)] = NA

dt[, wmed := median(weight,na.rm=T), by = patient_id]
dt[, hmed := median(height,na.rm=T), by = patient_id]

dt$weight[xor((dt$weight < dt$wmed-25) , (dt$weight > dt$wmed+25)) & !is.na(dt$weight)] = NA
dt$height[xor((dt$height < dt$hmed-20) , (dt$height > dt$hmed+20)) & !is.na(dt$height)] = NA

dt[, wmed := NULL]; dt[,hmed:=NULL]

### 3.1 clean state province record with text parsing
a = dt$state_province#table(a)
dt$state_province = a

### 3.2 Relevel hivstatus
a1 = dt$hivstatus
a = rep(NA,length(a1))
a[a1 %in% c("POSITIVE","1 Pos PCR")] = 1
a[a1 == "NEGATIVE"] = 0 #table(a) 0: 7678 1: 873722 
dt$hivstatus = a

### 4. selection criterior: 

### 4.1 patients with only one adult initial
tmp = dt[, sum(encounter_type == "adult initial"), by=patient_id]
keepid_1init = tmp$patient_id[tmp$V1 == 1] # those IDs with only 1 adult init
dt1 = dt[patient_id %in% keepid_1init,]

dim(dt1)#1653119  
length(unique(dt1$patient_id))# 92880

remove(keepid_1init,dt)

### 4.2 delete patients whose returnvisitdate < apptdate (RVD were scheduled backward)
# jth record: RVD(j) - APPT(j), clinician's expectation at APPT(j) in days for APPT(j+1), j=1,...,Nj 
dt1[, RVDgap := as.numeric(returnvisitdate - apptdate), by=patient_id]

## Flag records with RVD scheduled backward than the appt date
summary(dt1$RVDgap)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#-3564.00    28.00    42.00    48.34    56.00   730.00   121477 

tmp = dt1[RVDgap <= 0,]
dt1 = dt1[!(patient_id %in% tmp$patient_id),] # 91992 ppl left here

dim(dt1)#1651236 
length(unique(dt1$patient_id))# 92807
 
### 4.4 delete patients with [apptdate - dod > 30]
tmp = dt1[!is.na(dod) & (apptdate - dod > 30) ] # 215 records from 91 patients problematic
dt1 = dt1[!(patient_id %in% tmp$patient_id),] # 91957 ppl left here

dim(dt1)#1650031 
length(unique(dt1$patient_id))# 92716

### Delete those records after dod grid in grid laying
# tmp = dt1[!is.na(dod) & (apptdate - dod > 0)] ; dim(tmp) # 35 records

### 5.1 set transfer from 2 to NA if date of date of death is not missing
dt1[, .N, transfer]
dt1[!is.na(transfer) & transfer==2 & !is.na(dod), transfer:=NA, ] # transfer==2: 3974 -> 3931

### 5.2 Delete records after the FIRST tranfer out (2) 
## Mostly for those who xferred out and then re-entered
## xfer out that's NOT ENDPOINT is definitely informative -- patient returns either very late or very early

# Define new_num being 0 at adult initial
dt1[, ai := which(encounter_type=="adult initial"), by = patient_id]
dt1$new_num = dt1$num - dt1$ai; dt1[, ai:= NULL]
# Redefine lastnum based on new_num
dt1[, lastnum:= new_num[.N],by=patient_id]
#dt1[,.N,by=transfer] # 3931 xfer out records in total
sum(!is.na(dt1$transfer) & dt1$transfer==2 & dt1$new_num != dt1$lastnum) # 957 non-endpoint xfer out records
# Delete records after the FIRST tranfer out (2)  
dt1[, xfer2 := 1*(!is.na(transfer) & transfer==2)];dt1[xfer2==0, xfer2 := NA]
dt1[,.N,by=xfer2] 
dt1[, cfxfer2 := na.locf(xfer2,na.rm = F), by = patient_id]; dt1[is.na(cfxfer2), cfxfer2:=0]
dt1[, cfxfer2_prev := c(NA, cfxfer2[-.N]), by = patient_id ]
dt1[, first_xfer2 := cfxfer2 - cfxfer2_prev ]
dt1[,xfer2:=NULL]; dt1[, cfxfer2_prev :=NULL]
#dt1[,.N, c("first_xfer2","cfxfer2")]
#   first_xfer2 cfxfer2       N
#1:          NA       0   92712 # no xfer out and a first obs -> KEEP
#2:           0       0 1544091 # no xfer out yet up to this point -> KEEP 
#3:           1       1    3825 # first xfer out in follow up -> KEEP 
#4:           0       1    9399 # records after a first xfer out -> delete
#5:          NA       1       4 # first obs being xfer out -> delete; ppl xferred out at baseline - possibly input error
dt1[is.na(first_xfer2) , first_xfer2 := 0]
dt1[,xfer_del:= 1*(first_xfer2==0 & cfxfer2 == 1) ]
dt1 = dt1[xfer_del == 0,] 
dt1[, c("cfxfer2" ,"first_xfer2","xfer_del"):=NULL]

### 5.3 Delete patients with transfer out before adult initial
tmp = dt1[!is.na(transfer) & transfer == 2 & new_num < 0, ] # 1 record problematic, from a patient with 3 records left
dt1 = dt1[!(patient_id %in% tmp$patient_id),] # 91957 ppl left here

dim(dt1)#1640625 
length(unique(dt1$patient_id))#  92711

## Be careful that the last transfer record might not have RVD
dt1[transfer==2,.N,by=is.na(returnvisitdate)]
#   is.na    N
#1: FALSE 3139
#2:  TRUE  685

## Are all xfer out records the endpoint? - Yes.
## Update lastnum
dt1[, lastnum:= new_num[.N],by=patient_id] 

## DEFINE TIME
dt1[, time := as.numeric(apptdate - apptdate[encounter_type=="adult initial"]), by=patient_id]

dt2 = copy(dt1)
###########################################################
# clean the global environment and keep only list / dataframe 
rm(list = ls.str(mode = 'character'))
rm(list = ls.str(mode = 'numeric'))
remove(see,tmp)

save(dt1, file = "AMPATH_12-13-2019_dt1.RData")

# treatment switch; dt1 onarv has no missingness
dt1[, trt_switch := onarv - c(0,onarv[-.N]), by = patient_id] 
#table(dt1$trt_switch,useNA = "always") # 1: start trt, 0: trt status unchanged, -1: switch off trt
#-1       0       1       <NA> 
# 9354 1553542   77729       0
# cd4_obs: indicator of cd4 measurement in (0, 200]; cd4_first: indicator of a first cd4 measurement in (0, 200]
dt1[, cd4_obs := 1*(!is.na(cd4)), by = patient_id] 
dt1[, cd4_first := cumsum(cd4_obs), by=patient_id]
dt1[, cd4_first := 1.0*which(cd4_first==1)[1], by=patient_id]

#table(dt1$cd4_obs)
#0       1 
#1450776  189849 
#table(dt1$cd4_first)

if(0){
  ### DECRIPTION 1. in (0,200], CDF of time for the first CD4 and treatment initiation
  dt1[, ind := 1:(.N), by=patient_id]
  
  a = dt1$time[dt1$trt_switch == 1]
  b = dt1$time[dt1$ind == dt1$cd4_first]
  
  a1 = a[a>0 & a<=200 & !is.na(a)]; b1 = b[b>0 & b<=200 & !is.na(b)]
  p1=p2 = rep(NA,200)
  for(i in 1:200){
    p1[i] = mean(a1<=i,na.rm=T); p2[i] = mean(b1<=i,na.rm=T)
  }
  #black: Init Trt; red: First CD4
  png(filename = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft6\\DCDFTime.png",  width = 5, height = 4, units = 'in', res = 300)
  plot(p1,type="l",xlab="Day", ylab="Conditional CDF",ylim = range(c(p1,p2)),lwd = 2); lines(p2,col="red",lwd = 2)
  abline(v=12,col="blue",lwd = 2)
  legend(55, 0.3, legend=c("Measurement Time of the First CD4", "Time of Treatment Initiation","Day 12"),
         col=c("red", "black","blue"), lty=1, cex=0.8,lwd = 2)
  dev.off()
  
  p2[12]; 1-p1[12]
  #[1] 0.514606
  #[1] 0.9687763
}

##################################################################################
##################################################################################

### 5.4 Delete patients with trt initation before adult initial
tmp = dt1[,sum(time[trt_switch==1] < 0), by = patient_id] #number of pre-baseline trt initiations
table(tmp$V1) #0: 88841  1: 3867     2: 3 ; three people had 2 trt initations before baseline
dt1 = dt1[patient_id %in% tmp$patient_id[tmp$V1==0],]
dim(dt1)#1572255 
length(unique(dt1$patient_id))# 88841

if(0){
  ### DECRIPTION 2. Action proportion at baseline Day 0
  # treatment initation at Day 0: proportion 19.6%
  a = dt1$time[dt1$trt_switch == 1]; mean(a==0)
  # first cd4 measurement at Day 0: proportion 72.6% 
  b = dt1$time[dt1$ind == dt1$cd4_first]; mean(b==0,na.rm=T)
}

### 6.1 Carry forward/backward baseline variables V between (-200, 200] to time 0

V = c("male","age","married","yrenrol","TravelTime","whostage","weight","height","viralload")

gap = c(-199,200); anchor = 0

for(i in 1:length(V)){
  var = V[i]
  a = dt1[time == 0, get(var)]
  varinit = paste0(var,"init")
  if(sum(is.na(a)) != 0){
    w1 = movefun(dt1, var, gap, anchor)
    dt1[,tmp := get(var)]; dt1[time == 0 , tmp :=  w1]
    dt1[, (varinit) := rep(tmp[time==0],.N), by=patient_id]
    dt1[,tmp:=NULL]
  } else { # no missingness
    check = dt1[, length(unique(get(var))), by = patient_id]
    if(length(unique(check$V1))>1) # changes over time
      dt1[, (varinit) := rep(get(var)[time==0],.N), by=patient_id]
  }
}

### 6.2 Carry forward (-200,0] CD4 to time 0 and put into new covariate CD4_preB
w1 = movefun(dt1, "cd4", c(-199,0), 0)
dt1[,tmp := cd4]; dt1[time == 0 , tmp :=  w1]
dt1[, cd4_preB := rep(tmp[time==0],.N), by=patient_id]
dt1[,tmp:=NULL]

### 6.3 The last CD4 in (0,12] carried to time 0 and put into cd40 (what is the frequency of cd4 measurement during (0,30]?)
tmp = copy(dt1)
tmp[, tag1 := 1*(time >= 0 & time <= 12)]
tmp = tmp[tag1 == 1, c("patient_id", "time" , "cd4"), with=F]; tmp[time == 0, cd4 := NA] # NA to time==0 because day 0 is not included
tmp[, cd40 := tail(cd4[!is.na(cd4)],1), by= patient_id]
w1 = tmp[time == 0, cd40]

dt1[,tmp := cd4]; dt1[time == 0 , tmp :=  w1]
dt1[, cd4_0 := rep(tmp[time==0],.N), by=patient_id]
dt1[,tmp:=NULL]

### 6.4 The last onarv in [0,12] carried to time 0 and put into onarv1
tmp = copy(dt1)
tmp[, tag1 := 1*(time >= 0 & time <= 12)]
tmp = tmp[tag1 == 1, c("patient_id", "time","onarv"), with=F]
tmp[, ao := tail(onarv,1), by= patient_id]
w1 = tmp[time == 0, ao]

dt1[,tmp := onarv]; dt1[time == 0 , tmp :=  w1]
dt1[, onarv_0 := rep(tmp[time==0],.N), by=patient_id]
dt1[,tmp:=NULL]

#save(dt1, file = "AMPATH_12-13-2019_dt1_1.RData")
#load("AMPATH_12-13-2019_dt1_1.RData")

### 6.5 Delete records before time 0 and delete patients with only the baseline records
summary(dt1$time)
dt2 = dt1[time >= 0,]
dt2 = dt2[order(patient_id,time)]

dt2[, num := 0:((.N)-1), by=patient_id]
dt2[, lastnum:= num[.N],by=patient_id]

dt2 = dt2[lastnum > 0,]
dim(dt2)# 1521029
length(unique(dt2$patient_id))#  79504

### 7. Tag encounter/measurement records, last records, and last encounters
## DEFINE MEASUREMENT POINTS AND ENCOUNTER PROCESS
# RECORD_TYPE =  1 AS ENCOUNTER, 0 AS NON-ENCOUNTER MEASUREMENT
dt2[, record_type := 1*(encounter_type %in% c("adult initial","adult return"))]
dt2[, .N, by=record_type]# 0  210514; 1 1310515

# define LAST ENCOUNTER & LAST RECORD
# LENC SHOULD BE THE LAST ENCOUNTER (RECORD_TYPE==1) BUT NOT THE LAST RECORD!!!!
dt2[, whichlast := last(which(record_type==1)),by=patient_id] # index of the last encounter
dt2[, lenc := ifelse(!is.na(whichlast), c(rep(0,whichlast[.N]-1),1,rep(0,.N - whichlast[.N])), 0 ) ,by=patient_id]
dt2[, whichlast:= NULL]

dt2[, lrec := c(rep(0,.N-1),1),by = patient_id]

sum(dt2$lrec)#79504
sum(dt2$lenc)#79504

if(0){
  ### DECRIPTION 3. Summary table of covariates 
  colnames(dt2)
  #Vinit = c("male","ageinit","marriedinit","yrenrol","TravelTime","whostageinit","weightinit","heightinit","viralloadinit")
  varlist =c("gender","ageinit","marriedinit","yrenrol","TravelTime","whostageinit","weightinit","heightinit","viralloadinit",
             "cd4_preB","cd4_0","onarv_0")
  tmp = dt2[time == 0, ]
  tab = TabCovSumm1(varlist, tmp)
  
  ### DESCRIPTION 4. Total number of records
  visitnum =  dt2[,num[.N],by=patient_id]
  table(visitnum$V1[visitnum$V1<21])
  pdf("visitnum.pdf", width =5.5, height = 4)
  hist(visitnum$V1, main = "Total Number of Records Within Patient",xlab = "Number of Records")
  dev.off()
  
  ### DESCRIPTION 5. Total number of care encounters
  encnum =  dt2[, sum(record_type),by=patient_id]
  table(encnum$V1[encnum$V1<21])
  pdf("encnum.pdf", width =5.5, height = 4)
  hist(encnum$V1, main = "Total Number of Care Encounters Within Patient",xlab = "Number of Records")
  dev.off()

  ### DESCRIPTION 6: gap between RVD and Record dates
  dt2[,Rj_Aj := returnvisitdate - apptdate]
  
  pdf("Rj_Aj.pdf", width =5.5, height = 4)
  hist(as.numeric(dt2$Rj_Aj),xlim = c(0,200),xlab = "Days",breaks = 1000,main = "Clinician-scheduled Expected Gap Time")
  dev.off()
  
  mean(as.numeric(dt2$Rj_Aj) < 200, na.rm=TRUE) # 99.8% smaller than 200
  mean(!is.na(dt2$Rj_Aj)) # 95.3% non-missing
  max(dt2$Rj_Aj,na.rm=T) # max 730  
  
}
##################################################################################
##################################################################################

### 15. lay grids of 200 days
### 16. define death_outcome; mark all records in a grid where dod falls in, indicator lgrid_death
### 17. define xfer_outcome; mark the grid record after the grid of xfer==2, as having lgrid_xfer == 1 

### define engagement: outcome not death or xfer, outcome := 1*(totrec>1), by = tmpgridid
source(paste0(path,"/BartGcomp_19-12-13_Reclean_AMPATH_data_Athena.R"))

rm(list = ls.str(mode = 'character'))
rm(list = ls.str(mode = 'numeric'))

###########################################################

load("AMPATH_12-13-2019_dt4.RData")

#dat = dt6[,c("patient_id","time","apptdate","onarv_cf","trt_switch","ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit","new_cupd","new_wupd","new_vupd","outcome","log_weight","log_cd4plus1","vl_0","log_vlplus1"),with=F]
#dat = dt6[,c("patient_id","time","onarv_cf","trt_switch","ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit","new_cupd","new_wupd","new_vupd","outcome","log_weight","log_cd4plus1","vl_0","log_vlplus1"),with=F]

#theV = c("ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit")
#theX = c("new_wupd","log_weight",
#         "new_cupd","log_cd4plus1",
#         "new_vupd","vl_0","log_vlplus1")


### SETTING UP VARIABLES
# create new ID
num = dt4[, .N,by=patient_id]
dt4$myID = rep(1:length(num[[2]]), times = num[[2]])

theA = "onarv_cf"

theY = "outcome"

theX = c("log_cd4plus1", "new_cupd","log_cd4preBplus1", "new_cd4_cf")

theV = c("male","ageinit","marriedinit","yrenrol","TravelTime","whostageinit","weightinit","heightinit",
            "log_vlplus1_init","vl_0")

newdat = dt4[, c("patient_id","myID","time",theY, theA, theX, theV), with=F]

#source("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\Functions_new.R")
source(paste0(path,"/Functions_new.R"))

#0 1: onarv_cf male new_cupd new_wupd new_vupd outcome vl_0
#0 1 NA: marriedinit
#-1 0 1: trt_switch
#cont: yrenrol 
#cont NA: TravelTime, whostageinit
### FACTORIZE CATEGORICAL VARIABLES WITH <= c 
newdat = factorizedata(newdat,c = 4) # the 4 is not arbitrary...


### create missing category (for factors)/ (RX,R) (for continuous variables)
newdat = createmis(newdat)

### create indicator variables for factors
newdat = createind(newdat, keep = c(theA,theY,theX)) # variables NOT in 'keep' discards the original factor variable when its indicators are generated


theV = c("ageinit",
         "male_1",
         "yrenrol",
         "TravelTime_RX","TravelTime_R",
         "whostageinit_RX","whostageinit_R",
         "marriedinit_1","marriedinit_NA",
         "heightinit_RX","heightinit_R",
         "weightinit_RX","weightinit_R",
         "log_vlplus1_init_RX","log_vlplus1_init_R","vl_0_1")

theX = c("new_cupd_1","log_cd4plus1") # LHS

theX_cov = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R") # RHS (from createind and createmis)
theX_cov_1 = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R","log_cd4preBplus1_RX","log_cd4preBplus1_R","new_cd4_cf_RX","new_cd4_cf_R") # RHS (from createind and createmis)

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

theA = "onarv_cf_1"

newdat = lagdat(newdat, c(theA, theX_cov, theY_ind)) 

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

newdat = newdat[,unique(c("myID","time",theY,theV,theX,theX_cov_1,theY_ind,theA,theA_lag1,theX_cov_lag1,theY_ind_lag1)),with=F]

save(newdat, file = "20-01-12_newdat.RData")

dim(newdat)# 678815     34
length(unique(newdat$myID))#76740

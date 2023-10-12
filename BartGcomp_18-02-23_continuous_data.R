
setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\my working datasets")

setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\my working datasets")


library(readstata13)
library(lubridate)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(data.table)
library(zoo)

#d1 = read.dta13(".\\Stata Process\\hivdemoggic2008_2016Stata.dta")
#d2 = read.dta13(".\\Stata Process\\hivvisitgic2008_2016Stata.dta")
#setDT(d1)
#setDT(d2)
#setkey(d1, patient_id)
#setkey(d2, patient_id)
#d3 = merge(d1,d2, all.y=TRUE)
#tmp = d3[dod <= enroldate,]
#length(unique(tmp$patient_id))

# merged data of demog and visits, done in STATA, with abnormal death eliminated
dat = read.dta13(".\\denom_all.dta")

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


### !!!!!!!!!!!!!!!!!!!!!!!!!!!
### update the transfer variable
# 1 for transfer within ampath
# 2 for transfer out of ampath
dt$transfer[dt$encounter_type =="AMPATH transfer form"] = 2

length(unique(dt$patient_id))#96045 the very raw

### SELECTION CRITERIA
# patients with only one adult initial
tmp = dt[, sum(encounter_type == "adult initial"), by=patient_id]
keepid_1init = tmp$patient_id[tmp$V1 == 1] # those IDs with only 1 adult init
dt1 = dt[patient_id %in% keepid_1init,]

dim(dt1)#1653119  
length(unique(dt1$patient_id))# 92880

### 1% of adult init record (578 ppl) have missing return visit date ---> DELETE THESE
tmp = dt1[encounter_type == "adult initial" & is.na(returnvisitdate)]
dt1 = dt1[!(patient_id %in% tmp$patient_id)]

dim(dt1)#1642075
length(unique(dt1$patient_id))#92299

# patients with only one adult initial in the first four records
tmp = dt1[, num[encounter_type == "adult initial"], by = 
            patient_id]
keepid_in4 = tmp$patient_id[tmp$V1 < 5] 
dt1 = dt1[patient_id %in% keepid_in4,]

dim(dt1)
length(unique(dt1$patient_id))#92050

###########################################################
# clean the global environment and keep only list / dataframe 
rm(list = ls.str(mode = 'character'))
rm(list = ls.str(mode = 'numeric'))


###########################################################
# IS PATIENT ID IN ORDER?
#is.unsorted(dt1$patient_id, na.rm = TRUE, strictly = FALSE)
# SORT BY PATIENT ID AND NUM TO REDUCE POSSIBLE MISTAKES IN THE FOLLOWING CLEANING PROCESS
dt1 = dt1[order(patient_id,num)]

### DEFINE NEW INDEX NEW_NUM, ADULT INITIAL RECORD AS 0, AND RECORDS BEFORE THAT AS NEGATIVE
dt1[, ai := which(encounter_type=="adult initial"), by = patient_id]
dt1$new_num = dt1$num - dt1$ai
dt1[, ai:= NULL]

###########################################################
### DEFINE RECORD_TYPE, INDICATOR OF ENCOUNTER, 1 WHEN EXISTS A RETURN VISIT DATE

### DEFINE MEASUREMENT POINTS AND ENCOUNTER PROCESS
# RECORD_TYPE =  1 AS ENCOUNTER, 0 AS NON-ENCOUNTER MEASUREMENT, -1 AS INVALID MEASUREMENT (NA IN WEIGHT, CD4, VL), 0.5 AS FORCED INIT ENCOUNTER (THOSE WITHOUT RVD)
dt1$record_type = 1*(!is.na(dt1$returnvisitdate))

#table(dt1$record_type)
#0       1 
#119015 1516620 

###########################################################

### CARRY FORWARD MEASUREMENTS TO ADULT INITIAL RECORD FROM NEG INDEX RECORDS, 
### IFF MEASUREMENT WAS MADE WITHIN 31 DAYS, FOR HEIGHT WEIGHT CD4 AND VL

# 1. extract patient_id and their records with new_num<=0. keep only those having new_num<0.
# 2. for each HWCV, get the carried forward value and the correponding date at which the value was measured
# 3. for day gap (adult ini date - measurement date of the carried forward value) > 30, set to NA for the corresponding HWCV at new_num==0
# 4. update to dt1, check if the nonmissing percentage increased (should be)

cfdata = function(data, tmp_before, colnum, gap){
  var = colnames(data)[colnum]
  
  # cf
  # mt for measurement time of the non-missing values of var
  tmp_before$mt = tmp_before$apptdate
  tmp_before$mt[is.na(tmp_before[,get(var)])] = NA 
  # var1 : carried forward values of var
  tmp_before[, var1 := na.locf(get(var),na.rm = F), by = patient_id]
  # mt1 : carried forward values of mt
  tmp_before[, mt1 := na.locf(mt,na.rm=F), by = patient_id] # mt1 tells you which date this carried forward record was measured
  # we will void those carried forwarded number when timediff > gap days
  timediff = as.numeric(tmp_before$apptdate - tmp_before$mt1) 
  tmp_before[new_num==0 & timediff > gap, var1 := NA]
  # update these cf meaurements to new_num==0 in data
  tmp_before1 = tmp_before[tmp_before$new_num==0] 
  data[(data$patient_id %in% tmp_before1$patient_id) & new_num==0, (var):= tmp_before1$var1]
  return(data[,get(var)])  
  
}



cbdata = function(data, tmp_after, colnum, gap){
  var = colnames(data)[colnum]
  # mt for measurement time of the non-missing values of var
  tmp_after$mt = tmp_after$apptdate
  tmp_after$mt[is.na(tmp_after[,get(var)])] = NA 
  # var1 : carried backward values of var
  tmp_after[, var1 := na.locf(get(var),na.rm = F,fromLast = TRUE), by = patient_id]
  # mt1 : carried backward values of mt
  tmp_after[, mt1 := na.locf(mt,na.rm=F, fromLast = TRUE), by = patient_id] 
  # we will void those carried forwarded number when timediff > gap days
  timediff = as.numeric(tmp_after$mt1 - tmp_after$apptdate) 
  tmp_after[new_num==0 & timediff > gap, var1 := NA]
  # update these cb meaurements to new_num==0 in data
  tmp_after1 = tmp_after[tmp_after$new_num==0] 
  data[(data$patient_id %in% tmp_after1$patient_id) & new_num==0, (var):= tmp_after1$var1]
  return(data[,get(var)])  
}


data = dt1
varlist = c("height" ,"weight" , "cd4" ,"viralload","whostage", "married" )
cfcbcols = which(colnames(data) %in% varlist)

# prepare for carry forward
tmp_before = data[new_num <= 0 ]
tmp_before[tmp_before == -99] = NA
tmp_before[, prev_num := c(NA,new_num[-.N]), by = patient_id ]
tmp_before = tmp_before[-which(tmp_before$new_num==0 & is.na(tmp_before$prev_num)),] # exclude from target those with no negative new_num

# carry forward all V and X
cfdt = lapply(cfcbcols, function(k) cfdata(data, tmp_before, k,gap=31))
vname1 = paste0(colnames(data)[cfcbcols],1)
data[, (vname1) := cfdt]

# prepare for carry backward

tmp_after = data[new_num >= 0 ]
tmp_after[tmp_after == -99] = NA
tmp_after[, post_num := c(new_num[-1],NA), by = patient_id ]
tmp_after = tmp_after[-which(tmp_after$new_num==0 & is.na(tmp_after$post_num)),] # exclude from target those with no records after adult init record

# only carry backward V
vname2 = c("height1", "whostage1",  "married1")
cfcbcols1 = which(colnames(data) %in% vname2)
#colnames(data)[cfcbcols1]

cbdt = lapply(cfcbcols1, function(k) cbdata(data, tmp_after, k,gap=31))
vname22 = colnames(data)[cfcbcols1]
data[, (vname22) := cbdt]



rm(list = ls.str(mode = 'numeric'))

# DEFINE INIT HEIGHT, WHOSTAGE, MARRIED, AGE

dt1[, heightinit := rep(height1[new_num==0],.N), by=patient_id]

dt1[, whostageinit := rep(whostage1[new_num==0],.N), by=patient_id]

dt1[, marriedinit := rep(married1[new_num==0],.N), by=patient_id]

dt1[, ageinit := rep(age[new_num==0],.N), by=patient_id]

###################################################
# delete records after dod (259 records)
dt1 = dt1[-which(!is.na(dod) & apptdate>dod),] #1635635 ->1635376   

# define LAST ENCOUNTER & LAST RECORD
# LENC SHOULD BE THE LAST ENCOUNTER (RECORD_TYPE==1) BUT NOT THE LAST RECORD!!!!
dt1[, whichlast := last(which(record_type==1)),by=patient_id] # index of the last encounter
dt1[, lenc := ifelse(!is.na(whichlast), c(rep(0,whichlast[.N]-1),1,rep(0,.N - whichlast[.N])), 0 ) ,by=patient_id]
dt1[, whichlast:= NULL]

dt1[, lrec := c(rep(0,.N-1),1),by = patient_id]
#tmp1 = dt1[patient_id==73548]
#tmp1 = dt1[patient_id==75099]
sum(dt1$lrec)#92050
sum(dt1$lenc)#92050
###########################################################
### CLEAN XFER RECORDS
# dd = dt1

# tranfer should be pushed to the last encounter and not conflict with death
#    (a)if both death and transfer happened, it's death
#    (b)set transfer to 2 at lenc == 1 if transfer happened before AND death!=1
#    (c)set transfer to 0 for non-last-encounter records (lenc!=1)  
#look = dt1[patient_id == 78349]
#look = dt1[patient_id == 74358]

dt1$transfer[is.na(dt1$transfer)] = 0

#(a)
dt1$transfer[dt1$death==1] = 0

dt1[,transfer1 := max(transfer,na.rm = T),by=patient_id] # 1: transfer within, 2: transfer out

#(b)
dt1$transfer[dt1$lenc == 1 & dt1$transfer1==2] = 2

#(c)
dt1$transfer[dt1$lenc != 1] = 0

###########################################################
### KEEP ONLY THOSE NEW_NUM >= 0 ---> define and use dt2 for data analysis 

dt2 = dt1[new_num>=0,]
length(unique(dt2$patient_id)) # 92050
###########################################################
# DEFINE TIME

dt2$apptdate1 = as.numeric(dt2$apptdate)
dt2$rvddate1 = as.numeric(dt2$returnvisitdate)

dt2[, apptime := apptdate1 - apptdate1[1], by=patient_id]
dt2[, rvdtime := rvddate1 - apptdate1[1], by=patient_id]

# DEFINE DATA CLOSURE DATE
dt2$closure = as.numeric(max(dt2$apptdate)) # 2016-8-24

# HOW MANY GRID POINT FOR EACH PERSON?
dt2[, length_in_care:= closure[new_num==0] - apptdate1[new_num==0],by=patient_id]

dt2[, apptdate1 := NULL]
dt2[, rvddate1 := NULL]
dt2[, closure := NULL]


### patients' behavior
dt2[, Rj := c(NA,rvdtime[-.N]),by=patient_id]
dt2[, Aj_Rj:= apptime - Rj] # gap between actual show up and the scheduled date

### clinicians' behavior
dt2[, Rj_Aj:= rvdtime - apptime] # gap between actual show up and the scheduled date

####################
save(dt2, file = "continuous_data.RData")
####################

intervals: [Rj,apptime]
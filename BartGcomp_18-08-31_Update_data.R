setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)
######################################################################################################################
### ZEUS script::

library(GcompBART)
library(allBART)
library(lubridate)# date processing
library(data.table)
library(zoo)# carry backward forward
library(foreign)#stata read in
library(nnet)# multinomial log-linear models
library(Hmisc)#rmultinom function

library("parallel")

source("Functions_new.R")
######################################################################################################################
######################################################################################################################

#dat = fread("./makeGridData/my_analysis_data.csv")
dat = fread("my_analysis_data.csv")

dat = dat[,c("patient_id","time","apptdate","onarv_cf","trt_switch","ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit","new_cupd","new_wupd","new_vupd","outcome","log_weight","log_cd4plus1","vl_0","log_vlplus1"),with=F]
length(unique(dat$patient_id))

theV = c("ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit")
theX = c("new_wupd","log_weight",
         "new_cupd","log_cd4plus1",
         "new_vupd","vl_0","log_vlplus1")


# outcome: 0 disengage, 1 engage, 2 transfer_outcome, 3 death_outcome
# RELEVEL THE OUTCOME TO 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
dat$outcome = dat$outcome+1

### need to check data processing
### delete the 16 ppl who were xferred out at time==200
# check.... something wrong with xfer definition
todelete = dat[time == 200 & outcome == 3, "patient_id"]
delloc = which(dat$patient_id %in% todelete$patient_id)
dat = dat[-delloc,]



### SETTING UP VARIABLES

theA = "onarv_cf"

theY = "outcome"

# create new ID
num = dat[, .N,by=patient_id]
dat$myID = rep(1:length(num[[2]]), times = num[[2]])

dat[, log_weight_init := log_weight[1],by = patient_id]
#dat[, log_cd4plus1_init := log_cd4plus1[1],by = patient_id]
dat[, log_vlplus1_init := log_vlplus1[1],by = patient_id]
dat[, vl_0_init := vl_0[1], by = patient_id]

### FACTORIZE CATEGORICAL VARIABLES WITH <= c 
newdat = factorizedata(dat,c = 10)
levels(newdat$TravelTime) =c("missing","below30min","over2h","1to2h","30to60min") # do this because symbols like "<" or space are not eligible for variable names in regressions

### create missing category (for factors)/ (RX,R) (for continuous variables)
newdat = createmis(newdat)

### create indicator variables for factors
newdat = createind(newdat, keep = c(theA,theY,theX)) # variables NOT in 'keep' discards the original factor variable when its indicators are generated

theV = c("ageinit",
         "male_1",
         "yrenrol_2009","yrenrol_2010","yrenrol_2011","yrenrol_2012","yrenrol_2013","yrenrol_2014", "yrenrol_2015", "yrenrol_2016",
         "TravelTime_below30min","TravelTime_over2h","TravelTime_1to2h","TravelTime_30to60min",
         "whostageinit_2","whostageinit_3","whostageinit_4","whostageinit_NA",
         "marriedinit_1","marriedinit_NA",
         "heightinit_RX","heightinit_R",
         "log_weight_init_RX","log_weight_init_R",
         "log_vlplus1_init_RX","log_vlplus1_init_R","vl_0_init_1")

#theV = c("ageinit",
#         "male_1",
#         "yrenrol_2009","yrenrol_2010","yrenrol_2011","yrenrol_2012","yrenrol_2013","yrenrol_2014", "yrenrol_2015", "yrenrol_2016",
#         "TravelTime_below30min","TravelTime_over2h","TravelTime_1to2h","TravelTime_30to60min",
#         "whostageinit_2","whostageinit_3","whostageinit_4","whostageinit_NA",
#         "marriedinit_1","marriedinit_NA",
#         "heightinit_RX","heightinit_R",
#         "log_cd4plus1_init_RX",  "log_cd4plus1_init_R",
#         "log_weight_init_RX","log_weight_init_R",
#         "log_vlplus1_init_RX","log_vlplus1_init_R","vl_0_init_1")

theX = c("new_cupd_1","log_cd4plus1") # LHS

theX_cov = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R") # RHS (from createind and createmis)

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

theA = "onarv_cf_1"

newdat = lagdat(newdat, c(theA, theX_cov, theY_ind)) 

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

newdat = newdat[,unique(c("myID","time",theY,theV,theX,theX_cov,theY_ind,theA,theA_lag1,theX_cov_lag1,theY_ind_lag1)),with=F]

save(newdat, file = "newdat_0831.RData")

######################################################################################################################
theA = "onarv_cf"

theY = "outcome"


theV = c("ageinit",
          "male_1",
          "yrenrol_2009","yrenrol_2010","yrenrol_2011","yrenrol_2012","yrenrol_2013","yrenrol_2014", "yrenrol_2015", "yrenrol_2016",
          "TravelTime_below30min","TravelTime_over2h","TravelTime_1to2h","TravelTime_30to60min",
          "whostageinit_2","whostageinit_3","whostageinit_4","whostageinit_NA",
          "marriedinit_1","marriedinit_NA",
          "heightinit_RX","heightinit_R",
          "log_weight_init_RX","log_weight_init_R",
          "log_vlplus1_init_RX","log_vlplus1_init_R","vl_0_init_1")


theX = c("new_cupd_1","log_cd4plus1") # LHS

theX_cov = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R") # RHS (from createind and createmis)

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

theA = "onarv_cf_1"

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

load("newdat_0831.RData")

######################################################################################################################

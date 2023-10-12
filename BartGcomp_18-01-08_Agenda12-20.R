# AGENDA 12-13 (SHARED GOOGLE DOC IN TOPIC2)


setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")
#setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")

### PACKAGES

library(lubridate)# date processing
library(data.table)
library(zoo)# carry backward forward
library(foreign)#stata read in
library(nnet)# multinomial log-linear models
library(Hmisc)#rmultinom function

library("parallel")
library("MCMCpack")# bayesian analyses
library("MNP")
#library("micEcon")

# plotting
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

dir = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART"
#dir = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART"
source(paste0(dir,"\\Functions_new.R"))

options(scipen = 10)

#setwd("./Topic2");source("Functions_new.R")
dat = fread("my_analysis_data.csv")
dat = dat[,c("patient_id","time","apptdate","onarv_cf","trt_switch","ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit","new_cupd","new_wupd","new_vupd","outcome","log_weight","log_cd4plus1","vl_0","log_vlplus1"),with=F]
length(unique(dat$patient_id))
dim(dat) # 772868
colnames(dat)

theV = c("ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit")
theX = c("new_wupd","log_weight",
         "new_cupd","log_cd4plus1",
         "new_vupd","vl_0","log_vlplus1")


# outcome: 0 disengage, 1 engage, 2 transfer_outcome, 3 death_outcome
# RELEVEL THE OUTCOME TO 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
dat$outcome = dat$outcome+1

#table(dat$outcome)
#1      2      3      4 
#297386 466298   2575   6609 

#table(dat$outcome)/772868
#1           2           3           4 
#0.384782395 0.603334593 0.003331746 0.008551266 
######################################################################################################################
### SETTING UP VARIABLES

theA = "onarv_cf"

theY = "outcome"

# create new ID
num = dat[, .N,by=patient_id]
dat$myID = rep(1:length(num[[2]]), times = num[[2]])

### FACTORIZE CATEGORICAL VARIABLES WITH <= c 
newdat = factorizedata(dat,c = 10)
levels(newdat$TravelTime) =c("missing","below30min","over2h","1to2h","30to60min") # do this because symbols like "<" or space are not eligible for variable names in regressions

### create missing category (for factors)/ (RX,R) (for continuous variables)
newdat = createmis(newdat)

### create indicator variables for factors
newdat = createind(newdat, keep = c(theA,theY,theX)) # variables NOT in 'keep' discards the original factor variable when its indicators are generated
# the R from createmis is a non-factor binary vector, hence does not have new indicator created
# however, binary varialbes in theA and theX are factors, so new non-factor binary variables are created correspondingly, i.e. new_cupd is factor, new_cupd_1 is non-factor



#printvar(colnames(newdat))

theV = c("ageinit",
         "male_1",
         "yrenrol_2009","yrenrol_2010","yrenrol_2011","yrenrol_2012","yrenrol_2013","yrenrol_2014", "yrenrol_2015", "yrenrol_2016",
         "TravelTime_below30min","TravelTime_over2h","TravelTime_1to2h","TravelTime_30to60min",
         "whostageinit_2","whostageinit_3","whostageinit_4","whostageinit_NA",
         "marriedinit_1","marriedinit_NA",
         "heightinit_RX","heightinit_R")

theX = c("new_cupd_1","log_cd4plus1") # LHS

theX_cov = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R") # RHS (from createind and createmis)

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

theA = "onarv_cf_1"


#all(c(theA,theX,theY,theV) %in% colnames(newdat))

newdat = lagdat(newdat, c(theA, theX_cov, theY_ind)) 

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

newdat = newdat[,unique(c("myID","time",theY,theV,theX,theX_cov,theY_ind,theA,theA_lag1,theX_cov_lag1,theY_ind_lag1)),with=F]

# coerce all variables into numeric for feeding into MCMCpack
#for(i in 1:ncol(newdat)){
#  if(is.factor(newdat[,get(colnames(newdat)[i])])){
#    print(colnames(newdat)[i])
#    newdat[,(colnames(newdat)[i]):= as.numeric(as.character(get(colnames(newdat)[i])))]
#  }
#}
######################################################################################################################
### (3)
#[S1 | V, A0]
#[S2 | V, A1, S1]
#[S3 | V, A2, S2]

### Model building

parell = TRUE # TRUE if on server, FALSE otherwise
ncore = 5 # NULL if parall = FALSE

maxt= 13 # for identifiability of parameters. No xfer or death in the observed data at t=15
mod =  vector("list", maxt) 
modRHS =  vector("list", maxt) 

setting = "test"
ndraw = 100 
nburn = 10
nthin = 3
basey = 1
nchain = 2
npost = floor((ndraw - nburn) / (nthin+1))* nchain

setting = 1
ndraw = 1100 +1
nburn = 100
nthin = 4 # number of discarded samples between two valid samples
basey = 1 # disengagement as reference level
nchain = 1 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain

setting = 2
ndraw = 1100 +1
nburn = 100
nthin = 4 # number of discarded samples between two valid samples
basey = 1 # disengagement as reference level
nchain = 5 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain


for(t in 1:maxt){
  fitdat = newdat[time==200*t,]
  if(t==1){
    RHS = c(theV,theA_lag1)
    RHS = fmlfilter(fitdat,RHS)
    modRHS[[t]] = RHS
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    mod[[t]] = mnp_chains(fml, fitdat, basey, ndraw, nburn, nthin, nchain, doparallel = parell, ncore)
    print(t)
  } else {
    RHS = c(theV,theA_lag1,theY_ind_lag1)
    RHS = fmlfilter(fitdat,RHS)
    modRHS[[t]] = RHS
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    mod[[t]] =  mnp_chains(fml, fitdat, basey, ndraw, nburn, nthin, nchain, doparallel = parell, ncore)
    print(t)
  }
}

save(ndraw,nburn,nthin,npost, basey,mod,file = paste0("MNPmod_setting",setting,".RData"))

#load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/mod.RData")
#load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/mod_setting1.RData")

####################################################

### Outcome simulation 

A = 1 #theA_lag1
nlatent = length(unique(newdat$outcome)) - 1
nsub = sum(newdat$time==0)
ylevels = as.numeric(levels(newdat$outcome))

Smat = vector("list", maxt) 



ptm <- proc.time()

for(t in 1:maxt){
  
  fitdat = newdat[time==200*t,]
  theVhere = theV
  theVhere = fmlfilter(fitdat,theVhere)
  Vdist = newdat[time==0,theVhere,with=F]
  
  dat.pred = data.frame(1, Vdist, A)
  colnames(dat.pred) = c("Intercept", theVhere, theA_lag1)
  
  if(t==1){
    
    Xmathere = getXmat(fixdat = dat.pred, simdat = NULL, simname = NULL, fmlRHS = c("Intercept", modRHS[[t]]))
    
    Smat[[t]] = MNPpredict(m = mod[[t]], Xmat = Xmathere, nsub, basey, ylevels, doparallel = parell)# basey is disengagement

    print(t)
    
  } else {
    
    # need to define theY_ind_lag1 from simulated values; indicator of engagement at t-1 (either eng or diseng)
    Xmathere = getXmat(fixdat = dat.pred, simdat = c((Smat[[t-1]] == 2)*1) , simname = theY_ind_lag1, fmlRHS = c("Intercept", modRHS[[t]]))

    Smat[[t]] = MNPpredict(m = mod[[t]], Xmat = Xmathere, nsub, basey, ylevels, doparallel = parell)# basey is disengagement
    
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
    
    print(t)
    
  }
  
}

proc.time() - ptm

save(Smat,file = paste0("MNPSmat_setting",setting,".RData"))


##############################################################

tmp = Smat
tmp = cbind(1,Smat)
tmp = as.factor(tmp)
levels(tmp) = c("diengaged","engaged","transferred out","died")
tmp = relevel(tmp, "engaged")

a = data.frame(time = rep(200*(0:maxt),each = nrow(Smat)), outcome = tmp)
counts <- table(a$outcome,a$time)

png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\BayesSimEngA",A,".png"),width = 600, height = 400)
barplot(counts, main=paste0("Simulated Engagement Status at a=",A,"\n No Time-varying Covariates"),
        xlab="Days since baseline", ylab="Number of observations", col=1:length(levels(a$outcome)))
legend("bottomleft", 
       legend = rownames(counts), 
       fill = 1:4)
dev.off()

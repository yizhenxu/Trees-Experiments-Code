#setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")
setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")

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

#dir = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART"
dir = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART"
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
# however, binary variables in theA and theX are factors, so new non-factor binary variables are created correspondingly, i.e. new_cupd is factor, new_cupd_1 is non-factor



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

#newdat = normalizeX(newdat)# this makes RX 0 assignment difficult. skip for now.

# coerce all variables into numeric for feeding into MCMCpack
#for(i in 1:ncol(newdat)){
#  if(is.factor(newdat[,get(colnames(newdat)[i])])){
#    print(colnames(newdat)[i])
#    newdat[,(colnames(newdat)[i]):= as.numeric(as.character(get(colnames(newdat)[i])))]
#  }
#}

######################################################################################################################
### (4)
#[L1 | V, A0, X0, L0]
#[X1 | V, A0, X0, L0, L1 = 1] # simulated Xt remains NA until an Lt == 1
#[S1 | V, A0, X1, L1]
#[L2 | V, A1, X1, L1, S1]
#[X2 | V, A1, X1, L1, S1, L2 = 1]
#[S2 | V, A1, S1, X2, L2]



### Model building

parell = FALSE # TRUE if on server, FALSE otherwise
ncore = 5 # NULL if parall = FALSE

maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15
#-------------------------------
setting = "test"
maxt= 3
ndraw = 100 
nburn = 10
nthin = 2
basey = 1
nchain = 2
npost = floor((ndraw - nburn) / (nthin+1))* nchain
#-------------------------------
setting = 2
maxt= 13
ndraw = 1100 +1
nburn = 100
nthin = 4 # number of discarded samples between two valid samples
basey = 1 # disengagement as reference level
nchain = 5 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain

#-------------------------------

setting = 3
maxt= 3
ndraw = 800 +1
nburn = 200
nthin = 1 # number of discarded samples between two valid samples
basey = 1 # disengagement as reference level
nchain = 1 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain
#-------------------------------

tmp = (length(theX)+1)*maxt
mod =  vector("list", tmp) 
modfml =  vector("list", tmp) 


# initialize index count
ct = 1


ptm <- proc.time()
for(t in 1:maxt){
  fitdat = newdat[time==200*t,]
  if(t==1){
    
    # [L1 | V, A0, X0, L0]

    RHS = c(theV, theA_lag1, theX_cov_lag1)
    # exclude lag X_R from models at t=1 because that is the same as lag L at baseline
    RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    mod[[ct]] = MCMC_chains(type = "logit",fml, fitdat, ndraw, nburn, nthin, nchain, doparallel = parell, ncore)

    ct = ct+1
    
    # [X1 | V, A0, X0, L0, L1 = 1] # simulated Xt remains NA until a first Lt == 1
    
    RHS = c(theV, theA_lag1, theX_cov_lag1)
    # exclude lag X_R from models at t=1 because that is the same as lag L at baseline
    RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
    RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
    fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    mod[[ct]] = MCMC_chains(type = "regress",fml, fitdat[get(theX[1])==1,], ndraw, nburn, nthin, nchain, doparallel = parell, ncore)
    
    ct = ct+1
    
    # [S1 | V, A0, X1, L1]
    
    RHS = c(theV,theA_lag1,theX_cov)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    mod[[ct]] = mnp_chains(fml, fitdat, basey, ndraw, nburn, nthin, nchain, doparallel = parell, ncore)
    
    ct = ct+1
    
    print(t)
  } else {
    
    # [Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    
    RHS = c(theV, theA_lag1, theX_cov_lag1,theY_ind_lag1)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    mod[[ct]] = MCMC_chains(type = "logit",fml, fitdat, ndraw, nburn, nthin, nchain, doparallel = parell, ncore)
    
    ct = ct+1
    
    # [Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}, Lt = 1]
    
    RHS = c(theV, theA_lag1, theX_cov_lag1, theY_ind_lag1)
    RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
    fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    mod[[ct]] = MCMC_chains(type = "regress",fml, fitdat[get(theX[1])==1,], ndraw, nburn, nthin, nchain, doparallel = parell, ncore)
    
    ct = ct+1
    
    # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]
    
    RHS = c(theV, theA_lag1, theX_cov, theY_ind_lag1)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    mod[[ct]] = mnp_chains(fml, fitdat, basey, ndraw, nburn, nthin, nchain, doparallel = parell, ncore)
    
    ct = ct+1
    print(t)
  }
}
proc.time() - ptm

a = rep(NA, length(mod))
for(k in 1:length(mod)){
  a[k] = nrow(mod[[k]])
}
a

l = which(a < 1000)[1]
keep = floor(l/3)*3
maxt = keep/3
mod = mod[1:keep]# exclude the time points with data too small for model fitting
save(ndraw,nburn,nthin,npost,nchain,basey,mod,modfml,file = paste0("0222BayesModwCD4_setting",setting,".RData"))

#
# test 76903
# setting 2 77035 --- 17 hours

#load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS RData/BayesModwCD4_settingtest.RData")
load("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData\\BayesModwCD4_setting2.RData")
#######################################################
### CHECK if any step in the loop was problematic
### WILL AUTOMATE THIS!


#[1] 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000
#[16] 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000
#[31] 1000 1000 1000 1000 1000 1000 1000    5 1000

# obviously the X_13 model needs refitting
which(a!=npost)

t = 13
ct = 38
fitdat = newdat[time==200*t,]
RHS = c(theV, theA_lag1, theX_cov_lag1, theY_ind_lag1)
RHS = fmlfilter(fitdat[!is.na(get(theX[2])),],RHS)
fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
modfml[[ct]] = fml
mod[[ct]] = MCMC_chains(type = "regress",fml, fitdat[!is.na(get(theX[2])),], ndraw, nburn, nthin, nchain, doparallel = parell, ncore)
## use those data with non-missing X, but not those with upd==1 --> these two concepts are the same at t==1, but different afterwards

save(ndraw,nburn,nthin,npost,nchain,basey,mod,modfml,file = paste0("2-13modwCD4_setting",setting,".RData"))

#load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS RData/2-13modwCD4_setting2.RData")
# file too big
####################################################

### Outcome simulation 

A = 0 #theA_lag1
nlatent = length(unique(newdat$outcome)) - 1
nsub = sum(newdat$time==0)
ylevels = as.numeric(levels(newdat$outcome))

Smat = vector("list", maxt) 
Xmat = vector("list", maxt)

for(i in 1:maxt){
  Xmat[[i]] = array(NA, dim = c(nsub, npost, length(theX)))
}

Vdist = cbind(1,newdat[time==0,theV,with=F])
Vdist = Vdist[rep(1:nsub,npost),]

ptm <- proc.time()

ct = 1

for(t in 1:maxt){
  
  if(t==1){
    
    print(ct)
    #[L1 | V, A0, X0, L0]
    
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, time = t, modhere = modfml[[ct]], simX = Xmat, policy_lag = A, dathere = newdat)
    
    # these 3 lines are done in MCMCpredict for others, we do it like this for only for CD4 related variables in the fist round of sim bc this is faster
    bX = tcrossprod(as.matrix(dat.pred), mod[[ct]])
    predp = logistic(bX)
    Xmat[[t]][,,1] = matrix(rbinom(length(predp),1,c(predp)), ncol = npost)
    
    ct = ct+1
    
    print(ct)
    #[X1 | V, A0, X0, L0] 
    #dat.pred = makeX(Xname = theX, Vname = theV, time = t, modhere = mod[[ct]], simX = Xmat, policy_lag = A, dathere = newdat)
    # same RHS as in modeling L1
    bX = tcrossprod(as.matrix(dat.pred), mod[[ct]][,1:(ncol(mod[[ct]])-1)])
    sig = sqrt(mod[[ct]][,ncol(mod[[ct]])])
    tmp = rnorm(c(bX))
    yhat = tmp *rep(sig, each = nsub) + c(bX)
    postY = matrix(yhat, ncol = npost)
    
    L = Xmat[[t]][,,1]
    mis_ind = which(L == 0) # people who have no CD4 update
    postY[mis_ind] = NA
    Xmat[[t]][,,2] = postY
    
    ct = ct+1
    
    print(ct)
    #[S1 | V, A0, X1, L1]
    
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, time = t, modhere = modfml[[ct]], simX = Xmat, policy_lag = A, dathere = newdat)
    Smat[[t]] = MNPpredict(m = mod[[ct]], Xmat = dat.pred, nsub, basey, ylevels, doparallel = parell)# basey is disengagement
    
    ct = ct+1

    print(paste0("t=",t))
    
  } else {
  
    print(ct)
    #[Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, yl = ylevels, time = t, modhere = modfml[[ct]], simX = Xmat, simY = Smat, policy_lag = A, dathere = newdat)
    Xmat[[t]][,,1] = MCMCpredict(type = "logit", m = mod[[ct]], nsub, Xhere = dat.pred)
    
    ct = ct+1
    
    print(ct)
    #[Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, yl = ylevels, time = t, modhere = modfml[[ct]], simX = Xmat, simY = Smat, policy_lag = A, dathere = newdat)
    postY = MCMCpredict(type = "regress", m = mod[[ct]], nsub, Xhere = dat.pred)
    
    L = Xmat[[t]][,,1]
    X_prev = Xmat[[t-1]][,,2]
    mis_ind = which(L == 0) # people who have no CD4 update keeps the old value
    postY[mis_ind] = X_prev[mis_ind]
    Xmat[[t]][,,2] = postY
    
    ct = ct+1
    
    print(ct)
    #[St | V, A_{t-1}, S_{t-1}, Xt, Lt]
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, yl = ylevels, time = t, modhere = modfml[[ct]], simX = Xmat, simY = Smat, policy_lag = A, dathere = newdat)
    Smat[[t]] = MNPpredict(m = mod[[ct]], Xmat = dat.pred, nsub, basey, ylevels, doparallel = parell)# basey is disengagement
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
    
    ct = ct+1
    
    print(paste0("t=",t))

    
  }
  
}

proc.time() - ptm

save(Xmat,Smat,file = paste0("0222BayesSmat_setting",setting,"_a=",A,".RData"))#177759
#save(Xmat,Smat,file = paste0("BayesSmat_setting",setting,".RData"))

#load("Smat_setting2.RData")# in ZEUS
##############################################################
##############################################################
tmp = unlist(Smat)
#tmp = Smat
tmp1 = c(rep(2,length(Smat[[1]])),tmp)
#tmp = as.factor(tmp) # THIS IS SLOW AS HELL
#levels(tmp) = c("diengaged","engaged","transferred out","died")
#tmp = relevel(tmp, "engaged")

a = data.frame(time = rep(200*(0:maxt),each = length(Smat[[1]])), outcome = tmp1)
#counts <- table(a$outcome,a$time)

longdat = setDT(a)
ldc = longdat[, .N, by = c("time","outcome")]
ldc = rbind(ldc,list(time = c(0,0,0),outcome = c(1,3,4),N=c(0,0,0) ))
ldc = ldc[order(time,outcome)]
counts = matrix(ldc$N,nrow = 4)

ctmp = counts[1,]
counts[1,] = counts[2,]
counts[2,] = ctmp
#counts[1,1] = counts[2,1];counts[2,1]=0

rownames(counts) = c("engaged","diengaged","transferred out","died")
counts = simplify2array(lapply(1:ncol(counts), function(x) counts[,x]/sum(counts[,x])))
colnames(counts) = 200*(0:maxt)


#sort, make table
save(counts,file = paste0("0222Bayescounts_setting",setting,"_a=",A,".RData"))

load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS RData/0222Bayescounts_setting2_a=0.RData")
load("I:\\Topic2\\0222Bayescounts_setting2_a=0.RData")

png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\BayesCD4SimEngA",A,"_infoPrior.png"),width = 600, height = 400)
png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\BayesCD4SimEngA",A,"_infoPrior.png"),width = 600, height = 400)

barplot(counts, main=paste0("Bayesian Simulated Engagement Status at a=",A,"\n Time-varying Covariates: CD4 info"),
        xlab="Days since baseline", ylab="Number of observations", col=1:4)

legend("bottomleft", 
       legend = rownames(counts), 
       fill = 1:4)
dev.off()

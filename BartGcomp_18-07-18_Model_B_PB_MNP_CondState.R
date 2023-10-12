#devtools::install_github("yizhenxu/GcompBART")
#remove.packages("GcompBART")

#setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")

# SHOULD HAVE PICKED basey AS THE OUTCOME MODE. USE basey = "1" starting second Smod at t = 10

setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)
######################################################################################################################
### ZEUS script::

sink(paste0("output_",screennum,".txt")) #write all output to txt -- source('script.R',echo=T) to write input too

library(GcompBART)
library(MNP)
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
######################################################################################################################
theA = "onarv_cf"

theY = "outcome"

theV = c("ageinit",
         "male_1",
         "yrenrol",
         "TravelTime_RX","TravelTime_R",
         "whostageinit_RX","whostageinit_R",
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

load("newdat.RData") # defined in 18-02-26 script
#load("./makeGridData/newdat.RData")

#OUTCOME: 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
newdat$outcome = as.numeric(newdat$outcome) #IMPORTANT! de-factorize outcome

######################################################################################################################
### optional: pick nsub people
seed = 1
nsub = 90000
nn = max(newdat$myID)
set.seed(seed)
id = sample(1:nn, nsub, replace = F)

###regular
newdat = newdat[myID %in% id,]
### 

fitdat = newdat[time == 200,]
table(fitdat$outcome)
fitdat = newdat[time == 400,]
table(fitdat$outcome)


### Separate data into training set and test set

ntr = 40000 # floor(nsub * 0.7) # sample size of training set
nte = 10000

set.seed(1)
samp = sample(id,nsub)

trid = samp[1:ntr]  
trid = sort(trid)
trd = copy(newdat[myID %in% trid,])

teid = samp[(ntr+1):(ntr+nte)]  
teid = sort(teid)
ted = copy(newdat[myID %in% teid,])

gcid = samp[(ntr+nte+1):nsub]  
gcid = sort(gcid)
gcd = copy(newdat[myID %in% gcid,])


fitdat = trd[time == 200,]
table(fitdat$outcome)
fitdat = trd[time == 400,]
table(fitdat$outcome)

fitdat = ted[time == 200,]
table(fitdat$outcome)
fitdat = ted[time == 400,]
table(fitdat$outcome)
###
##############################################################
### Build models based on trd (training data)
#parell = TRUE # TRUE if on server, FALSE otherwise
#ncore = 5 # NULL if parall = FALSE

#maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15

#-------------------------------
setting = 1
maxt= 10 # cannot fit for t = 11, sample size not fittable
nd = 2000
nb = 100
nt = 100
#basey = 1 # disengagement as reference level
#nchain = 5 # number of MCMC chains

Prior_cont = list(nu = 3, sigq = 0.9,
                  ntrees = nt,
                  kfac = 2,
                  pswap = 0.1, pbd = 0.5, pb = 0.25,
                  alpha = 0.95, beta = 2.0,
                  nc = 100, minobsnode = 10)
Prior_binary = list(ntrees = nt,
                    kfac = 2,
                    pswap = 0.1, pbd = 0.5, pb = 0.25,
                    alpha = 0.95, beta = 2.0,
                    nc = 100, minobsnode = 10)
Mcmc_cont = Mcmc_binary = list(burn=nb, ndraws = nd)
Prior_mult = function(p){
  return(list(nu = p-1+3, V = diag(p-1),
              ntrees = nt,
              kfac = 2,
              pswap = 0.1, pbd = 0.5, pb = 0.25,
              alpha = 0.95, beta = 2.0,
              nc = 100, minobsnode = 10))
}

Mcmc_mult = function(p){
  return(list(sigma0 = diag(p-1), burn = nb, ndraws = nd, nSigDr = 50))
}
#-------------------------------

ndraw = nd+nb 
nburn = nb
nthin = 0 # number of discarded samples between two valid samples
basey1 = "2" #basey for t = 1
basey = c(rep("2", 2*(maxt-1)-1), "1") #basey for t>1
nchain = 1 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain

##############################################################
### Build models based on trd (training data)
#parell = TRUE # TRUE if on server, FALSE otherwise
#ncore = 5 # NULL if parall = FALSE

#maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15

### (1) Bayesian with CD4 information


tmp = (length(theX)+1)*maxt
mod =  vector("list", tmp) # replace the third model at each t>1 with 2 models in Smod
modfml =  vector("list", tmp) 

Smod =  vector("list", 2*(maxt-1)) 


load("B_PB_MNP_screen46698_ct33_nsub90000.RData")

# initialize index count
ct = 1
sct = 1

ptm <- proc.time()
for(t in 1:maxt){
  fitdat = trd[time==200*t,]
  if(t==1){
    
    # [L1 | V, A0, X0, L0]
    
    ct = ct+1
    
    # [X1 | V, A0, X0, L0, L1 = 1] # simulated Xt remains NA until a first Lt == 1
    
    ct = ct+1
    
    # [S1 | V, A0, X1, L1]
    
    ct = ct+1
    
    print(t)
    
  } else {
    
    # [Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    
    ct = ct+1
    
    # [Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}, Lt = 1]
    
    ct = ct+1
    
    # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]

    #RHS = c(theV, theA_lag1, theX_cov, theY_ind_lag1)
    RHS = c(theV, theA_lag1, theX_cov)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    ### S_{t-1} = disengaged
    Smod[[sct]] = mnp_chains(as.formula(fml), fitdat[get(theY_ind_lag1)==0,], basey[sct], ndraw, nburn, nthin, nchain, doparallel = FALSE)
    sct = sct + 1
    ### S_{t-1} = engaged
    Smod[[sct]] = mnp_chains(as.formula(fml), fitdat[get(theY_ind_lag1)==1,], basey[sct], ndraw, nburn, nthin, nchain, doparallel = FALSE)
    sct = sct + 1
    
    ct = ct+1
    print(t)
  }
}
proc.time() - ptm

save(mod,Smod,modfml,file = paste0("CondSModel_screen",screennum,"_ct",ct,"_ntr",ntr,".RData"))


##################################################################################
##################################################################################
library(mailR)
a =  proc.time()-ptm
status = as.character(a[3])

sender <- "916800617qq@gmail.com"
recipients <- c("yizhen_xu@brown.edu")
send.mail(from = sender,
          to = recipients,
          subject = paste0(screennum," complete"),
          body = status,
          smtp = list(host.name = "smtp.gmail.com", port = 465,
                      user.name = "916800617qq@gmail.com",
                      passwd = "Xyz123456!", ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

sink()




#setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)
######################################################################################################################
### ZEUS script::

sink(paste0("output_",screennum,".txt")) #write all output to txt -- source('script.R',echo=T) to write input too

library(GcompBART)
library(allBART)
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

#OUTCOME: 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
newdat$outcome = as.numeric(newdat$outcome) #IMPORTANT! de-factorize outcome

######################################################################################################################
### optional: pick nsub people
seed = 1
nsub = 50000
nn = max(newdat$myID)
set.seed(seed)
id = sample(1:nn, nsub, replace = F)

###regular
newdat = newdat[myID %in% id,]
### 

### Separate data into training set and test set

ntr = 20000 # floor(nsub * 0.7) # sample size of training set
nte = 10000

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
ngc = length(unique(gcid))
##############################################################
### Build models based on trd (training data)

#maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15

#-------------------------------
setting = 1
maxt= 10 # cannot fit for t = 11, sample size not fittable for cond S
ndraw = 1000
nburn = 100
ntree = 100
ncut = 100
nsd = 50 # max repeated draw of Sigma for Mpbart
seed = 88
#basey = 2 # the outcome level used in RHS of model fitting (theY_ind)
basey1 = "2" #basey for t = 1
basey = c(rep("2", 2*(maxt-1)-1), "1") #basey for t>1
nmod = 3*maxt # number of models

Prior_cont = list(nu = 3, sigq = 0.9,
                  ntrees = ntree,
                  kfac = 2,
                  pswap = 0.1, pbd = 0.5, pb = 0.25,
                  alpha = 0.95, beta = 2.0,
                  nc = ncut, minobsnode = 10)
Prior_binary = list(ntrees = ntree,
                    kfac = 2,
                    pswap = 0.1, pbd = 0.5, pb = 0.25,
                    alpha = 0.95, beta = 2.0,
                    nc = ncut, minobsnode = 10)
Mcmc_cont = Mcmc_binary = list(burn=nburn, ndraws = ndraw)
Prior_mult = function(p){
  return(list(nu = p-1+3, V = diag(p-1),
              ntrees = ntree,
              kfac = 2,
              pswap = 0.1, pbd = 0.5, pb = 0.25,
              alpha = 0.95, beta = 2.0,
              nc = ncut, minobsnode = 10))
}

Mcmc_mult = function(p){
  return(list(sigma0 = diag(p-1), burn = nburn, ndraws = ndraw,nSigDr = nsd))
}
#-------------------------------

##############################################################
### Build models based on trd (training data)
tmp = (length(theX)+1)*maxt
mod =  vector("list", tmp) # replace the third model at each t>1 with 2 models in Smod
modfml =  vector("list", tmp) 
Smod =  vector("list", 2*(maxt-1)) 

runtime = rep(NA, nmod + 1)


Smat = vector("list", maxt)
Scov = vector("list", maxt)
Xmat = vector("list", maxt)

for(i in 1:maxt){
  Xmat[[i]] = array(NA, dim = c(ngc, ndraw, length(theX)))
}

Vdist = cbind(gcd[time==0,theV,with=F])
Vdist = Vdist[rep(1:ngc,ndraw),]


# initialize index count
ct = 1

ptm <- proc.time()
for(t in 1:maxt){
  fitdat = trd[time==200*t,]
  if(t==1){
    
    # [L1 | V, A0, X0, L0]
    
    ## trd
    RHS = c(theV, theA_lag1, theX_cov_lag1)
    # exclude lag X_R from models at t=1 because that is the same as lag L at baseline
    RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    runtime[ct] = proc.time()[3]
    tmp = model_bart(as.formula(fml), data = fitdat, type = "binary",
                     Prior = Prior_binary,
                     Mcmc = Mcmc_binary)
    
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod[[ct]] = tmp[-del]
    remove(tmp)
    
    ##########################
    ct = ct+1
    
    # [X1 | V, A0, X0, L0, L1 = 1] # simulated Xt remains NA until a first Lt == 1
    
    ## trd
    RHS = c(theV, theA_lag1, theX_cov_lag1)
    # exclude lag X_R from models at t=1 because that is the same as lag L at baseline
    RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
    RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
    fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    runtime[2] = proc.time()[3]
    tmp = model_bart(as.formula(fml), data = fitdat[get(theX[1])==1,], type = "continuous",
                     Prior = Prior_cont,
                     Mcmc = Mcmc_cont)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod[[ct]] = tmp[-del]
    remove(tmp)
    
    ##########################
    ct = ct+1
    
    # [S1 | V, A0, X1, L1]
    
    ## trd
    ytab = table(fitdat$outcome)
    ylevels = as.numeric(names(ytab)[ytab > 0])
    p = length(ylevels) # number of outcome categories
    
    RHS = c(theV,theA_lag1,theX_cov)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    runtime[ct] = proc.time()[3]
    tmp = model_bart(as.formula(fml), data = fitdat, type = "multinomial",
                     base = basey1,
                     Prior = Prior_mult(p),
                     Mcmc = Mcmc_mult(p),
                     correction = FALSE)
    
    del = which(names(tmp) %in% c("treefit","samp_y"))
    if(length(del)>0){
      mod[[ct]] = tmp[-del]
      remove(tmp)
    } else {
      mod[[ct]] = tmp
    }
    
    
    ##########################
    ct = ct+1
    
    print(t)
    
  } else {
    
    # [Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    
    ## trd
    RHS = c(theV, theA_lag1, theX_cov_lag1,theY_ind_lag1)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    runtime[ct] = proc.time()[3]
    tmp = model_bart(as.formula(fml), data = fitdat, type = "binary",
                     Prior = Prior_binary,
                     Mcmc = Mcmc_binary)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod[[ct]] = tmp[-del]
    remove(tmp)
    
    ##########################
    ct = ct+1
    
    # [Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}, Lt = 1]
    
    ## trd
    RHS = c(theV, theA_lag1, theX_cov_lag1, theY_ind_lag1)
    RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
    fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    runtime[ct] = proc.time()[3]
    tmp = model_bart(as.formula(fml), data = fitdat[get(theX[1])==1,], type = "continuous",
                     Prior = Prior_cont,
                     Mcmc = Mcmc_cont)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod[[ct]] = tmp[-del]
    remove(tmp)
    
    ##########################
    ct = ct+1
    
    # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]
    
    #RHS = c(theV, theA_lag1, theX_cov, theY_ind_lag1)
    RHS = c(theV, theA_lag1, theX_cov)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    ### S_{t-1} = disengaged
    ## trd
    ytab = table(fitdat[get(theY_ind_lag1)==0,]$outcome)
    ylevels = as.numeric(names(ytab)[ytab > 0])
    p = length(ylevels) # number of outcome categories
    
    if(p>2){
      tmp = model_bart(as.formula(fml), data =  fitdat[get(theY_ind_lag1)==0,], type = "multinomial",
                       base = basey[(t-1)*2-1],
                       Prior = Prior_mult(p),
                       Mcmc = Mcmc_mult(p),
                       correction = 0)
      del = which(names(tmp) %in% c("treefit","samp_y"))
      if(length(del)>0){
        Smod[[(t-1)*2-1]] = tmp[-del]
        remove(tmp)
      } else {
        Smod[[(t-1)*2-1]] = tmp
      }
      
      
    } else if(p==2){
      
      fitdathere = fitdat[get(theY_ind_lag1)==0,]
      Y = rep(0,nrow(fitdathere))
      Y[fitdathere$outcome == ylevels[2]] = 1
      fitdathere$outcome = Y
      tmp = model_bart(as.formula(fml), data = fitdathere, type = "binary",
                       Prior = Prior_binary,
                       Mcmc = Mcmc_binary)
      del = which(names(tmp) %in% c("treefit","samp_y"))
      Smod[[(t-1)*2-1]] = tmp[-del]
      remove(tmp)
      
    }
    
    ### S_{t-1} = engaged
    ytab = table(fitdat[get(theY_ind_lag1)==1,]$outcome)
    ylevels = as.numeric(names(ytab)[ytab > 0])
    p = length(ylevels) # number of outcome categories
    
    if(p>2){
      tmp = model_bart(as.formula(fml), data =  fitdat[get(theY_ind_lag1)==1,], type = "multinomial",
                       base = basey[(t-1)*2],
                       Prior = Prior_mult(p),
                       Mcmc = Mcmc_mult(p),
                       correction = 0)
      del = which(names(tmp) %in% c("treefit","samp_y"))
      if(length(del)>0){
        Smod[[(t-1)*2]] = tmp[-del]
        remove(tmp)
      } else {
        Smod[[(t-1)*2]] = tmp
      }
      
    } else if(p==2){
      
      fitdathere = fitdat[get(theY_ind_lag1)==1,]
      Y = rep(0,nrow(fitdathere))
      Y[fitdathere$outcome == ylevels[2]] = 1
      fitdathere$outcome = Y
      tmp = model_bart(as.formula(fml), data = fitdathere, type = "binary",
                       Prior = Prior_binary,
                       Mcmc = Mcmc_binary)
      del = which(names(tmp) %in% c("treefit","samp_y"))
      if(length(del)>0){
        Smod[[(t-1)*2]] = tmp[-del]
        remove(tmp)
      } else {
        Smod[[(t-1)*2]] = tmp
      }
      
    }
    
    ##########################
    ct = ct+1
    print(t)
  }
}
proc.time() - ptm

save(mod,Smod,modfml,file = paste0("CondSModel_screen",screennum,"_ct",ct,".RData"))


##################################################################################
##################################################################################
library(mailR)

sender <- "916800617qq@gmail.com"
recipients <- c("yizhen_xu@brown.edu")
send.mail(from = sender,
          to = recipients,
          subject = paste0(screennum," complete"),
          body = "status",
          smtp = list(host.name = "smtp.gmail.com", port = 465,
                      user.name = "916800617qq@gmail.com",
                      passwd = "Xyz123456!", ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

sink()



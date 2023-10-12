# screen 242855
library(lubridate)# date processing
library(data.table)
library(zoo)# carry backward forward
library(foreign)#stata read in
library(nnet)# multinomial log-linear models
library(Hmisc)#rmultinom function

library("parallel")
library("MCMCpack")# bayesian analyses
library("MNP")

setwd("./Topic2");source("Functions_new.R")

######################################################################################################################

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

newdat = lagdat(newdat, c(theA, theX_cov, theY_ind)) 

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

newdat = newdat[,unique(c("myID","time",theY,theV,theX,theX_cov,theY_ind,theA,theA_lag1,theX_cov_lag1,theY_ind_lag1)),with=F]

save(newdat, file = "newdat.RData")

######################################################################################################################
theA = "onarv_cf"

theY = "outcome"

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

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

load("newdat.RData")

######################################################################################################################

### Separate data into training set and test set

nsub = max(newdat$myID)
set.seed(1)
samp = sample(1:nsub)

trid = samp[1:floor(nsub * 0.7)]  
trid = sort(trid)
trd = copy(newdat[myID %in% trid,])

teid = samp[(floor(nsub * 0.7)+1):nsub]  
teid = sort(teid)
ted = copy(newdat[myID %in% teid,])

######################################################################################################################
### Build models based on trd (training data)
parell = TRUE # TRUE if on server, FALSE otherwise
ncore = 5 # NULL if parall = FALSE

maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15

#-------------------------------
setting = 2
maxt= 11
ndraw = 1100 +1
nburn = 100
nthin = 4 # number of discarded samples between two valid samples
basey = 1 # disengagement as reference level
nchain = 5 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain

#-------------------------------

##############################################
### (1) Bayesian with CD4 information


tmp = (length(theX)+1)*maxt
mod =  vector("list", tmp) 
modfml =  vector("list", tmp) 


# initialize index count
ct = 1


ptm <- proc.time()
for(t in 1:maxt){
  fitdat = trd[time==200*t,]
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

save(ndraw,nburn,nthin,npost,nchain,basey,mod,modfml,file = paste0("0226EvalBayesModwCD4_setting",setting,".RData"))
# 242888
##############################################
### (2) Bayesian without time varying information (MNP)

mod =  vector("list", maxt) 
modfml =  vector("list", maxt) 

for(t in 1:maxt){
  fitdat = trd[time==200*t,]
  if(t==1){
    RHS = c(theV,theA_lag1)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[t]] = fml
    mod[[t]] = mnp_chains(fml, fitdat, basey, ndraw, nburn, nthin, nchain, doparallel = parell, ncore)
    print(t)
  } else {
    RHS = c(theV,theA_lag1,theY_ind_lag1)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[t]] = fml
    mod[[t]] =  mnp_chains(fml, fitdat, basey, ndraw, nburn, nthin, nchain, doparallel = parell, ncore)
    print(t)
  }
}

save(ndraw,nburn,nthin,npost,nchain,basey,mod,modfml,file = paste0("0226EvalMNP_setting",setting,".RData"))
# 243516
##############################################

a = rep(NA, length(mod))
for(k in 1:length(mod)){
  a[k] = nrow(mod[[k]])
}
a


l = which(a < npost)[1]
t = ceiling(l/3)
ct = l

# rerun the model
##############################################
### (3) Frequentist with CD4 information


fitmod = function(y, fml, fitdat){
  
  a = length(unique(y))
  
  if(is.factor(y) & a > 2){
    m =  multinom(fml , data=fitdat) 
  } else{
    if(a == 2){
      m = glm(fml, family=binomial(link='logit'), data=fitdat) 
    } else {
      m = lm(fml, data=fitdat)
    }
  }
  return(m)
}

tmp = (length(theX)+1)*maxt
mod = vector("list", tmp) 
modfml =  vector("list", tmp) 

ct = 1
for(t in 1:maxt){
  
  fitdat = trd[time==200*t,]
  
  if(t==1){
    
    RHS = c(theV, theA_lag1, theX_cov_lag1)
    # exclude lag X_R from models at t=1 because that is the same as lag L at baseline
    RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
    RHS = fmlfilter(fitdat, RHS)
    fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml

    a = fitdat[,get(theX[1])]
    mod[[ct]] = fitmod(a, fml, fitdat) 
    ct = ct + 1
    
    RHS = c(theV, theA_lag1, theX_cov_lag1)
    # exclude lag X_R from models at t=1 because that is the same as lag L at baseline
    RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
    RHS = fmlfilter(fitdat[get(theX[1])==1,], RHS)
    fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    a = fitdat[get(theX[1])==1, get(theX[2])]
    mod[[ct]] = fitmod(a, fml, fitdat[get(theX[1])==1,])  # only fit for those current upd==1
    ct = ct + 1
    
    RHS = c(theV,theA_lag1,theX_cov)
    RHS = fmlfilter(fitdat, RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    mod[[ct]] =  multinom(fml , data=fitdat) 
    ct = ct + 1
    
    print(t)
  } else {
    
    RHS = c(theV,theA_lag1, theX_cov_lag1, theY_ind_lag1)
    RHS = fmlfilter(fitdat, RHS)
    fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    a = fitdat[,get(theX[1])]
    mod[[ct]] = fitmod(a, fml, fitdat) 
    ct = ct + 1
    
    RHS = c(theV,theA_lag1, theX_cov_lag1, theY_ind_lag1)
    RHS = fmlfilter(fitdat[get(theX[1])==1,], RHS)
    fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    a = fitdat[get(theX[1])==1, get(theX[2])]
    mod[[ct]] = fitmod(a, fml, fitdat[get(theX[1])==1,])
    ct = ct + 1
    
    RHS = c(theV,theA_lag1,theX_cov,theY_ind_lag1)
    RHS = fmlfilter(fitdat, RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    mod[[ct]] =  multinom(fml , data=fitdat) 
    ct = ct + 1
    
    print(t)
  }
}
##############################################
### (4) Frequentist without time varying information (Multinomial)


mod =  vector("list", maxt) 
modfml =  vector("list", maxt) 


for(t in 1:maxt){
  fitdat = trd[time==200*t,]
  if(t==1){
    
    RHS = c(theV,theA_lag1)
    RHS = fmlfilter(fitdat, RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[t]] = fml
  
    mod[[t]] =  multinom(fml , data = fitdat) 
    
    print(t)
  } else {
    RHS = c(theV,theA_lag1,theY_ind_lag1)
    RHS = fmlfilter(fitdat, RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[t]] = fml
 
    mod[[t]] =  multinom(fml , data = fitdat) 
    print(t)
  }
}


######################################################################################################################
### Evaluate Bayesian models based on test data

#load("I:\\Topic2\\newdat.RData")

#load("I:\\Topic2\\0222BayesModwCD4_setting2.RData")

#source("I:\\Topic2\\Functions_new.R")



### Evaluation --- per model
Eval_permod = function(mod, modfml, ted, times, types, intypes, maxt){
  
  lenmod = length(mod)
  nummod = length(types)
  
  E = vector("list", nummod)
  for(i in 1:nummod){
    if(types[i] == "mult"){
      E[[i]] = matrix(NA,5,maxt)
      rownames(E[[i]]) = c("y=1","y=2","y=3","y=4","overall")
      colnames(E[[i]]) = paste0("t=",1:maxt)
    } else {
      E[[i]] = rep(NA, maxt)
    }
  }
  
  ylist = vector("list",lenmod)
  
  PM = vector("list",lenmod)
  
  for(ct in 1:lenmod){
    
    t = times[ct]
    
    fv = fmlvars(modfml[[ct]])
 
    dat.pred = ted[time==200*t,fv$RHS[-1],with=F]
    dat.pred = cbind(1, dat.pred)
    colnames(dat.pred)[1] = fv$RHS[1]
    
    y =  unlist(ted[time==200*t,fv$LHS,with = F])
    
    ylist[[ct]] = y
    
    numy = length(unique(y))
    
    if (numy == 2){# logistic
      bX = tcrossprod(as.matrix(dat.pred), mod[[ct]])
      predp = logistic(bX)
      PM[[ct]] = matrix(rbinom(length(predp),1,c(predp)), ncol = npost)
      yM = matrix(rep(y, npost),ncol = npost)
      E[[intypes[ct]]][t] = mean(PM[[ct]] == yM)
      
    } else if (numy < 10) { # multinomial
      
      fitdat = trd[time==200*t,]
      ylevels = as.numeric(unique(fitdat$outcome))
      #ylevels = as.numeric(unique(y))
      latenty = sort(ylevels[-which(ylevels==basey)])
      nlatent = length(latenty)
      
      # simulate from multivariate normal (mu,Sigma)
      
      nhalfSig = nlatent*(nlatent - 1)/2 + nlatent
      sigmas = m[, (ncol(m) - nhalfSig + 1): ncol(m)]
      params = m[, 1:(ncol(m) - nhalfSig)]
      
      BbyY = getParamByY(params,nlatent)
 
      mu = lapply(1:nlatent, function(yl) tcrossprod(as.matrix(dat.pred), BbyY[[yl]]) )
      mu_array = simplify2array(mu) # 3 dims: nsub x npost x nlatent
      
      postY = lapply(1:npost, function(j) getYhat(j, nlatent, basey, latenty, mu_array, sigmas))
      PM[[ct]] = simplify2array(postY)
      remove(postY)
      
      yM = matrix(rep(y, npost),ncol = npost)
     
      for(k in 1:4){
        ind = which(yM==k)
        E[[intypes[ct]]][k,t] = mean(PM[[ct]][ind] == yM[ind])
      }
      E[[intypes[ct]]][5,t] = mean(PM[[ct]] == yM)
      
    } else { # continuous
      upd = unlist(ted[time==200*t, theX[1],with = F])
      ind = which(upd==1)
      y = y[ind]
      dat.pred = dat.pred[ind,]
      
      bX = tcrossprod(as.matrix(dat.pred), mod[[ct]][,1:(ncol(mod[[ct]])-1)])
      sig = sqrt(mod[[ct]][,ncol(mod[[ct]])])
      tmp = rnorm(c(bX))
      yhat = tmp *rep(sig, each = nrow(bX)) + c(bX)
      PM[[ct]] = matrix(yhat, ncol = npost)
      yM = matrix(rep(y, npost),ncol = npost)
      E[[intypes[ct]]][t] = mean((PM[[ct]] - yM)^2)
    }

    print(ct)
  }
  
  res = list(Pred = PM, Y = ylist, Accuracy = E)
  
  return(res)
}

# Bayesian G Comp
times = rep(1:maxt,each = 3) # time interval index for each ct model 
types = c("logit","cont","mult") # how many and which type should each E[[i]] be
intypes = rep(1:3, maxt) # which E[[i]] to go into

BGC_allres = Eval_permod(mod, modfml, ted, times, types, intypes, maxt)

save(BGC_allres, file = paste0("0226EvalBayesModwCD4_setting",setting,"_allres.RData"))

# MNP
times = c(1:maxt) # time interval index for each ct model 
types = "mult" # how many and which type should each E[[i]] be
intypes = rep(1, maxt) # which E[[i]] to go into

MNP_allres = Eval_permod(mod, modfml, ted, times, types, intypes, maxt)

save(MNP_allres, file = paste0("0226EvalMNP_setting",setting,"_allres.RData"))

######################################################################################################################
### Evaluate ordinary models based on test data


predmod = function(y, mod, fitdat){
  
  a = length(unique(y))
  
  if(is.factor(y) & a > 2){
    cp = predict(mod, newdata = fitdat, "probs")  
    yhat = rMultinom(cp, 1) 
  } else{
    if(a == 2){
      cp = predict(mod, newdata = fitdat, type="response")  
      yhat =  rbinom(n=nrow(fitdat), size=1, prob=cp)
    } else {
      yhat = predict(mod, newdata = fitdat) 
      sighat = summary(mod)$sigma
      epsilon = rnorm(nrow(fitdat),mean = 0, sd = sighat)
      yhat = yhat+epsilon
    }
  }
  return(yhat)
}



### Evaluation --- per model
Eval_permod_freq = function(mod, modfml, ted, times, types, intypes, maxt){
  
  lenmod = length(mod)
  nummod = length(types)
  
  E = vector("list", nummod)
  for(i in 1:nummod){
    if(types[i] == "mult"){
      E[[i]] = matrix(NA,5,maxt)
      rownames(E[[i]]) = c("y=1","y=2","y=3","y=4","overall")
      colnames(E[[i]]) = paste0("t=",1:maxt)
    } else {
      E[[i]] = rep(NA, maxt)
    }
  }
  

  for(ct in 1:lenmod){
    
    t = times[ct]
    
    fv = fmlvars(modfml[[ct]])
    
    dat.pred = ted[time==200*t,fv$RHS[-1],with=F]
    dat.pred = cbind(1, dat.pred)
    colnames(dat.pred)[1] = fv$RHS[1]
    
    y =  unlist(ted[time==200*t,fv$LHS,with = F])
       
    numy = length(unique(y))
    
    if (numy == 2){# logistic
      
      yhat = predmod(y, mod[[ct]], dat.pred)
      E[[intypes[ct]]][t] = mean(yhat == y)
      
    } else if (numy < 10) { # multinomial
      
      yhat = predmod(y, mod[[ct]], dat.pred)
      
      for(k in 1:4){
        ind = which(y==k)
        E[[intypes[ct]]][k,t] = mean(yhat[ind] == y[ind])
      }
      E[[intypes[ct]]][5,t] = mean(yhat == y)
      
    } else { # continuous
      upd = unlist(ted[time==200*t, theX[1],with = F])
      ind = which(upd==1)
      y = y[ind]
      dat.pred = dat.pred[ind,]
      
      yhat = predmod(y, mod[[ct]], dat.pred)
      E[[intypes[ct]]][t] = mean((yhat - y)^2)
    }
    
    print(ct)
  }

  
  return(E)
}

# Freq w CD4 info
times = rep(1:maxt,each = 3) # time interval index for each ct model 
types = c("logit","cont","mult") # how many and which type should each E[[i]] be
intypes = rep(1:3, maxt) # which E[[i]] to go into

FGC_allres = Eval_permod_freq(mod, modfml, ted, times, types, intypes, maxt)

save(FGC_allres, file = paste0("0226EvalFreqModwCD4_allres.RData"))

# Freq w.o. CD4
times = c(1:maxt) # time interval index for each ct model 
types = "mult" # how many and which type should each E[[i]] be
intypes = rep(1, maxt) # which E[[i]] to go into

FMN_allres = Eval_permod_freq(mod, modfml, ted, times, types, intypes, maxt)

save(FMN_allres, file = paste0("0226EvalFreqModOutcome_allres.RData"))

######################################################################################################################
# put all allres RData into one
setwd("./Topic2")
setting = 2

load(paste0("0226EvalBayesModwCD4_setting",setting,"_allres.RData"))
load(paste0("0226EvalMNP_setting",setting,"_allres.RData"))
load("0226EvalFreqModwCD4_allres.RData")
load("0226EvalFreqModOutcome_allres.RData")

BGC = BGC_allres$Accuracy
BMN = MNP_allres$Accuracy
FGC = FGC_allres
FMN = FMN_allres
save(BGC, BMN, FGC, FMN, file = "0226Eval_allres.RData")

######################################################################################################################
load("I:\\Topic2\\0226Eval_allres.RData")

BGC[[1]] - FGC[[1]]
BGC[[2]] - FGC[[2]]
BGC[[3]] - FGC[[3]]

BMN[[1]] - FMN[[1]]

png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\Accuract_w_CD4info.png"),width = 900, height = 650)

m <- matrix(c(1:8,9,9,9,9),nrow = 3,ncol = 4,byrow = TRUE)

layout(mat = m,heights = c(0.4,0.4,0.2))

par(mar = c(2,2.5,2,2))
#1
plot(1:maxt, BGC[[1]], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5)
title( "CD4 Update (Binary)", line = 0.4,cex.main = 1.5)
lines(1:maxt, FGC[[1]] ,col="red")
#2
plot(1:maxt, BGC[[2]], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5)
title("log(CD4+1) | Update"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FGC[[2]] ,col="red")
#3
plot(1:maxt, BGC[[3]][5,], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5)
title("Outcome"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FGC[[3]][5,] ,col="red")
#4
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#5
plot(1:maxt, BGC[[3]][1,], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5)
title("Engagement"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FGC[[3]][1,] ,col="red")
#6
plot(1:maxt, BGC[[3]][2,], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5)
title("Disengagement"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FGC[[3]][2,] ,col="red")
#7
plot(1:maxt, BGC[[3]][3,], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5, ylim = range(c(BGC[[3]][3,],FGC[[3]][3,])))
title("Xferred Out"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FGC[[3]][3,] ,col="red")
#8
plot(1:maxt, BGC[[3]][4,], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5, ylim = range(c(BGC[[3]][4,],FGC[[3]][4,])))
title("Death"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FGC[[3]][4,] ,col="red")
#9 
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,legend = c("Bayesian", "non-Bayesian"), 
       cex=1.5, fill = 1:2, horiz = TRUE)
#reset par
par(mar=c(0,0,0,0))
dev.off()

#############################################


png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\Accuract_MN.png"),width = 900, height = 650)

m <- matrix(c(1,2,2,2,3,4,5,6,7,7,7,7),nrow = 3,ncol = 4,byrow = TRUE)

layout(mat = m,heights = c(0.4,0.4,0.2))

par(mar = c(2,2.5,2,2))
#1
plot(1:maxt, BMN[[1]][5,], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5)
title("Outcome"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FMN[[1]][5,] ,col="red")
#2
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#3
plot(1:maxt, BMN[[1]][1,], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5)
title("Engagement"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FMN[[1]][1,] ,col="red")
#4
plot(1:maxt, BMN[[1]][2,], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5)
title("Disengagement"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FMN[[1]][2,] ,col="red")
#5
plot(1:maxt, BMN[[1]][3,], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5, ylim = range(c(BMN[[1]][3,],FMN[[1]][3,])))
title("Xferred Out"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FMN[[1]][3,] ,col="red")
#6
plot(1:maxt, BMN[[1]][4,], type = "l", 
     ylab = "Out-of-sample Mean Accuracy", xlab = "Time Index", cex.axis=1.5, ylim = range(c(BMN[[1]][4,],FMN[[1]][4,])))
title("Death"  , line = 0.4,cex.main = 1.5) 
lines(1:maxt, FMN[[1]][4,] ,col="red")
#9 
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,legend = c("Bayesian", "non-Bayesian"), 
       cex=1.5, fill = 1:2, horiz = TRUE)
#reset par
par(mar=c(0,0,0,0))
dev.off()


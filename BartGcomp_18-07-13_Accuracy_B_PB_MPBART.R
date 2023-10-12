#setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)

#20 trees
load("B_PB_MPBART_screen46844_ct15_ntr40000.RData")
load("B_PB_MPBART_screen46887_ct33_ntr40000.RData")

#LM
load("B_PB_MNP_screen46698_ct16_nsub90000.RData")
load("B_PB_MNP_screen46763_ct34_nsub90000.RData")

#50 trees
load("B_PB_MPBART_screen46873_ct33_ntr40000.RData")

#100 trees
load("B_PB_MPBART_screen46538_ct15_ntr40000.RData")
load("B_PB_MPBART_screen46495_ct33_ntr40000.RData")

for(i in 16:33){
  mod[[i]] = mod1[[i]]
  modfml[[i]] = modfml1[[i]]
}

#46844 BART-PBART-MPBART accuracy 20 trees 100 burn 2000 draws t=1,...,11
#46698 BART-PBART-MNP accuracy t=1,...,11

######################################################################################################################

library(GcompBART)

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
######################################################################################################################
### Evaluation by simulation outcomes --- per model on test set
# types and intypes give the evalulation output format

Eval_permod_bart = function(mod, modfml, ted, times, types, intypes, maxt){
  
  npost = mod[[1]]$ndraws
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
    
    fv = fmlvars(modfml[[ct]], itc = 0)
    
    dat.pred = ted[time==200*t,fv$RHS,with=F]
    y =  unlist(ted[time==200*t,fv$LHS,with = F])
    
    if(mod[[ct]]$type == "continuous"){
      upd = unlist(ted[time==200*t, theX[1],with = F])
      ind = which(upd==1)
      y = y[ind]
      dat.pred = dat.pred[ind,]
    }
    
    PM = predict_bart(mod[[ct]], dat.pred)$samp_y
    yM = matrix(rep(y, npost),ncol = npost)
    
    if(mod[[ct]]$type == "binary"){
      
      E[[intypes[ct]]][t] = mean(PM == yM)
      
    } else if(mod[[ct]]$type == "multinomial"){
      
      for(k in 1:4){
        ind = which(yM==k)
        E[[intypes[ct]]][k,t] = mean(PM[ind] == yM[ind])
      }
      E[[intypes[ct]]][5,t] = mean(PM == yM)
      
    } else {
      
      E[[intypes[ct]]][t] = mean((PM - yM)^2)
      
    }
    
    
    print(ct)
  }

  return(E)
}

# Bayesian G Comp
times = rep(1:maxt,each = 3) # time interval index for each ct model 
types = c("binary","cont","mult") # how many and which type should each E[[i]] be
intypes = rep(1:3, maxt) # which E[[i]] to go into

bart_allres = Eval_permod_bart(mod, modfml, ted, times, types, intypes, maxt)

###BART-PBART-MPBART
Eval_permod_bpbmnp = function(mod, modfml, ted, times, types, intypes, maxt){
  npost = mod[[1]]$ndraws
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
    
    fv = fmlvars(modfml[[ct]], itc = 0)
    y =  unlist(ted[time==200*t,fv$LHS,with = F])
    
    if(ct %% 3 == 0){
      
      fitdat = trd[time==200*t,]
      ytab = table(fitdat$outcome)
      ylevels = as.numeric(names(ytab)[ytab > 0])
      p = length(ylevels) # number of outcome categories
      m = mod[[ct]]
      
      if(p>2){#do MNP
        fv = fmlvars(modfml[[ct]])
        dat.pred = ted[time==200*t,fv$RHS[-1],with=F]
        dat.pred = cbind(1, dat.pred)
        colnames(dat.pred)[1] = fv$RHS[1]
        
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
        PM = simplify2array(postY)
        remove(postY)
        
      } else {#multinomial reduces to binary
        PM = predict_bart(mod[[ct]], dat.pred)$samp_y
        
        PM[which(PM==0)] = ylevels[1]
        PM[which(PM==1)] = ylevels[2]
      }
    } else {#cont or binary
      dat.pred = ted[time==200*t,fv$RHS,with=F]
      
      if(mod[[ct]]$type == "continuous"){ 
        upd = unlist(ted[time==200*t, theX[1],with = F])
        ind = which(upd==1)
        y = y[ind]
        dat.pred = dat.pred[ind,]
      }
      
      PM = predict_bart(mod[[ct]], dat.pred)$samp_y
      
    }
    
    
    yM = matrix(rep(y, npost),ncol = npost)
    
    if( types[intypes[ct]] == "binary"){
      
      E[[intypes[ct]]][t] = mean(PM == yM)
      
    } else if(types[intypes[ct]] == "mult"){
      
      for(k in 1:4){
        ind = which(yM==k)
        E[[intypes[ct]]][k,t] = mean(PM[ind] == yM[ind])
      }
      E[[intypes[ct]]][5,t] = mean(PM == yM)
      
    } else {
      
      E[[intypes[ct]]][t] = mean((PM - yM)^2)
      
    }
    
    
    print(ct)
  }
  #res = list(Pred = PM, Y = ylist, Accuracy = E)
  return(E)
}



# Bayesian G Comp
times = rep(1:maxt,each = 3) # time interval index for each ct model 
types = c("binary","cont","mult") # how many and which type should each E[[i]] be
intypes = rep(1:3, maxt) # which E[[i]] to go into

allres = Eval_permod_bpbmnp(mod, modfml, ted, times, types, intypes, maxt)

######################################################################################################################
### Evaluation by POSTERIOR MODE --- per model on test set
# types and intypes give the evalulation output format
yhat = function(vec){
  names(which.max(table(vec))) 
}

Mode_bart = function(mod, modfml, ted, times, types, intypes, maxt){
  
  npost = mod[[1]]$ndraws
  lenmod = length(mod)
  nummod = length(types)
  
  E = vector("list", nummod)
  
  for(i in 1:nummod){
    E[[i]] = rep(NA, maxt)
  }
  
  for(ct in 1:lenmod){
    
    t = times[ct]
    
    fv = fmlvars(modfml[[ct]], itc = 0)
    
    dat.pred = ted[time==200*t,fv$RHS,with=F]
    y =  unlist(ted[time==200*t,fv$LHS,with = F])
    
    if(mod[[ct]]$type == "continuous"){
      upd = unlist(ted[time==200*t, theX[1],with = F])
      ind = which(upd==1)
      y = y[ind]
      dat.pred = dat.pred[ind,]
    }
    
    PM = predict_bart(mod[[ct]], dat.pred)$samp_y
    PY = apply(PM, 1, yhat)
    
    if(mod[[ct]]$type == "binary"){
      
      E[[intypes[ct]]][t] = mean(PY == y)
      
    } else if(mod[[ct]]$type == "multinomial"){
      
      E[[intypes[ct]]][t] = mean(PY == y)
      
    } else {
      
      yM = matrix(rep(y, npost),ncol = npost)
      E[[intypes[ct]]][t] = mean((PM - yM)^2)
      
    }
    
    
    print(ct)
  }
  
  return(E)
}

#look at only part of the posterior, the last nshow posteriors
#nthin*nshow <= npost
#nthin: pick one sample out of every nthin samples
#nshow: pick nshow samples in total
Mode_bart_part = function(mod, modfml, ted, times, types, intypes, maxt, nshow, nthin, seed){
  
  set.seed(seed)
  npost = mod[[1]]$ndraws
  
  lenmod = length(mod)
  nummod = length(types)
  
  E = vector("list", nummod)
  E1 = vector("list", nummod)
  
  for(i in 1:nummod){
    E[[i]] = rep(NA, maxt)
    E1[[i]] = rep(NA, maxt)
  }
  
  for(ct in 1:lenmod){
    
    t = times[ct]
    
    fv = fmlvars(modfml[[ct]], itc = 0)
    
    dat.pred = ted[time==200*t,fv$RHS,with=F]
    y =  unlist(ted[time==200*t,fv$LHS,with = F])
    
    if(mod[[ct]]$type == "continuous"){
      upd = unlist(ted[time==200*t, theX[1],with = F])
      ind = which(upd==1)
      y = y[ind]
      dat.pred = dat.pred[ind,]
    }
    
    #PM = predict_bart(mod[[ct]], dat.pred)$samp_y[,(npost-nshow+1):npost]
    PM = predict_bart(mod[[ct]], dat.pred)$samp_y[,(npost-nshow*nthin+1)+(0:(nshow-1))*nthin]
    PY = apply(PM, 1, yhat)
    
    if(mod[[ct]]$type == "binary"){
      
      E[[intypes[ct]]][t] = mean(PY == y)
      
    } else if(mod[[ct]]$type == "multinomial"){
      
      E[[intypes[ct]]][t] = mean(PY == y)
      
    } else {
      
      yM = matrix(rep(y, nshow),ncol = nshow)
      E[[intypes[ct]]][t] = mean((PM - yM)^2)
      
    }
    
    
    yM = matrix(rep(y, nshow),ncol = nshow)
    
    if(mod[[ct]]$type == "binary"){
      
      E1[[intypes[ct]]][t] = mean(PM == yM)
      
    } else if(mod[[ct]]$type == "multinomial"){
      
      E1[[intypes[ct]]][t] = mean(PM == yM)
      
    } else {
      
      E1[[intypes[ct]]][t] = mean((PM - yM)^2)
      
    }
    
    print(ct)
  }
  
  return(list(Mode_E = E, Sim_E = E1))
}

# Bayesian G Comp
times = rep(1:maxt,each = 3) # time interval index for each ct model 
types = c("binary","cont","mult") # how many and which type should each E[[i]] be
intypes = rep(1:3, maxt) # which E[[i]] to go into

bart_mode = Mode_bart(mod, modfml, ted, times, types, intypes, maxt)

bart_mode_part = Mode_bart_part(mod, modfml, ted, times, types, intypes, maxt, 1000,2,1)

###BART-PBART-MPBART
yhat = function(vec){
  names(which.max(table(vec))) 
}

Mode_bpbmnp = function(mod, modfml, ted, times, types, intypes, maxt){
  npost = mod[[1]]$ndraws
  lenmod = length(mod)
  nummod = length(types)
  
  E = vector("list", nummod)
  
  for(i in 1:nummod){
    E[[i]] = rep(NA, maxt)
  }
  
  for(ct in 1:lenmod){
    
    t = times[ct]
    
    fv = fmlvars(modfml[[ct]], itc = 0)
    y =  unlist(ted[time==200*t,fv$LHS,with = F])
    
    if(ct %% 3 == 0){
      
      fitdat = trd[time==200*t,]
      ytab = table(fitdat$outcome)
      ylevels = as.numeric(names(ytab)[ytab > 0])
      p = length(ylevels) # number of outcome categories
      m = mod[[ct]]
      
      if(p>2){#do MNP
        fv = fmlvars(modfml[[ct]])
        dat.pred = ted[time==200*t,fv$RHS[-1],with=F]
        dat.pred = cbind(1, dat.pred)
        colnames(dat.pred)[1] = fv$RHS[1]
        
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
        PM = simplify2array(postY)
        remove(postY)
        
      } else {#multinomial reduces to binary
        PM = predict_bart(mod[[ct]], dat.pred)$samp_y
        
        PM[which(PM==0)] = ylevels[1]
        PM[which(PM==1)] = ylevels[2]
      }
    } else {#cont or binary
      dat.pred = ted[time==200*t,fv$RHS,with=F]
      
      if(mod[[ct]]$type == "continuous"){ 
        upd = unlist(ted[time==200*t, theX[1],with = F])
        ind = which(upd==1)
        y = y[ind]
        dat.pred = dat.pred[ind,]
      }
      
      PM = predict_bart(mod[[ct]], dat.pred)$samp_y
      
    }
    
    PY = apply(PM, 1, yhat)
    
    if( types[intypes[ct]] == "binary"){
      
      E[[intypes[ct]]][t] = mean(PY == y)
      
    } else if(types[intypes[ct]] == "mult"){
      
      E[[intypes[ct]]][t] = mean(PY == y)
      
    } else {
      
      yM = matrix(rep(y, npost),ncol = npost)
      E[[intypes[ct]]][t] = mean((PM - yM)^2)
      
    }
    
    
    print(ct)
  }
  #res = list(Pred = PM, Y = ylist, Accuracy = E)
  return(E)
}



# Bayesian G Comp
times = rep(1:maxt,each = 3) # time interval index for each ct model 
types = c("binary","cont","mult") # how many and which type should each E[[i]] be
intypes = rep(1:3, maxt) # which E[[i]] to go into

moderes = Mode_bpbmnp(mod, modfml, ted, times, types, intypes, maxt)

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

# no Scov and train_pred train_y for 10000 ppl and all
#save(mod, modfml, bart_allres, file = paste0("Bart_screen",screennum,"_ct",ct,"_nsub",nsub,".RData"))

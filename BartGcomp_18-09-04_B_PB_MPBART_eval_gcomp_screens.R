
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

A = 0

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
Jiao = FALSE
ylevels = c("1","2","3","4")

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
###BART-PBART-MPBART
yhat = function(vec){
  names(which.max(table(vec))) 
}

newAcc_bart = function(mod, Smod, modfml, ted, times, types, intypes, maxt, basey1, basey){
  npost = mod[[1]]$ndraws
  lenmod = length(times)
  nummod = length(types)
  
  E = vector("list", nummod)
  
  for(i in 1:nummod){
    if(types[i] == "mult"){
      E[[i]] = matrix(NA,8,maxt)
      rownames(E[[i]]) = c("y=1","y=2","y=3","y=4","average","mode","tree_average","tree_mode")
      colnames(E[[i]]) = paste0("t=",1:maxt)
    } else if(types[i] == "binary"){
      E[[i]] = matrix(NA,2,maxt)
      rownames(E[[i]]) = c("average","mode")
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
      d = ted[time==200*t,]
      
      if(t == 1){
        PB = predict_bart(mod[[ct]], d)
        PM = PB$samp_y
        PMT = PB$samp_treefit
      } else{
        PB = predict_bart(Smod[[(t-1)*2-1]], d)
        PM = PB$samp_y
        PMT = PB$samp_treefit
        
        PB1 = predict_bart(Smod[[(t-1)*2]], d)
        PM1 = PB1$samp_y
        PMT1 = PB1$samp_treefit
        
        ind = which(d[,theY_ind_lag1,with = F]==1)#check
        PM[ind,] = PM1[ind,]
        PMT[ind,] = PMT1[ind,]
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
      PY = apply(PM, 1, yhat)
      E[[intypes[ct]]][1,t] = mean(PM == yM)
      E[[intypes[ct]]][2,t] = mean(PY == y)
      
    } else if(types[intypes[ct]] == "mult"){
      
      for(k in 1:4){
        
        E[[intypes[ct]]][k,t] = mean(((PM != k) & (yM != k)) | (PM==yM))
        
      }
      E[[intypes[ct]]][5,t] = mean(PM == yM)
      
      PY = apply(PM, 1, yhat)
      E[[intypes[ct]]][6,t] = mean(PY == y)
      
      E[[intypes[ct]]][7,t] = mean(PMT == yM)
      
      PY = apply(PMT, 1, yhat)
      E[[intypes[ct]]][8,t] = mean(PY == y)
      
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

allres = newAcc_bart(mod, modfml, ted, times, types, intypes, maxt, basey1, basey)

##############################################################
### Build models based on trd (training data)

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

    ## gcd
    policy_lag = rep(A, ngc)
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, gcd)
    Xmat[[t]][,,1] = predict_bart(mod[[ct]],dat.pred)$samp_y
    
    ##########################
    ct = ct+1
    
    # [X1 | V, A0, X0, L0, L1 = 1] # simulated Xt remains NA until a first Lt == 1
    
    ## gcd
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, gcd)
    postY = predict_bart(mod[[ct]],dat.pred)$samp_y
    
    L = Xmat[[t]][,,1]
    mis_ind = which(L == 0) # people who have no CD4 update
    postY[mis_ind] = NA
    Xmat[[t]][,,2] = postY
    remove(postY)
    
    ##########################
    ct = ct+1
    
    # [S1 | V, A0, X1, L1]
    
    ## gcd
    policy_lag = rep(A, nrow(Vdist))
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, gcd)
    tmp = predict_bart(mod[[ct]],dat.pred, Gcomp = TRUE, mppred = TRUE)
    Smat[[t]] = tmp$samp_treefit
    
    ##########################
    ct = ct+1
    
    print(t)
    
  } else {
    
    # [Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    
    ## gcd Smat
    policy_lag = rep(A, nrow(Vdist))
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, gcd)
    Xmat[[t]][,,1] = predict_bart(mod[[ct]],dat.pred, Gcomp = TRUE)$samp_y
    
    ##########################
    ct = ct+1
    
    # [Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}, Lt = 1]
    
    ## gcd Smat
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, gcd)
    postY = predict_bart(mod[[ct]],dat.pred, Gcomp = TRUE)$samp_y
    
    L = Xmat[[t]][,,1]
    X_prev = Xmat[[t-1]][,,2]
    mis_ind = which(L == 0) # people who have no CD4 update keeps the old value
    postY[mis_ind] = X_prev[mis_ind]
    Xmat[[t]][,,2] = postY
    
    remove(postY)
    remove(L)
    remove(X_prev)
    remove(mis_ind)
    
    ##########################
    ct = ct+1
    
    # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]
    
    ### S_{t-1} = disengaged
    ## trd
    ytab = table(fitdat[get(theY_ind_lag1)==0,]$outcome)
    yl = as.numeric(names(ytab)[ytab > 0])
    p = length(yl) # number of outcome categories
    
    if(p>2){
      
      ## gcd
      dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, gcd)
      Smat[[t]] = predict_bart(Smod[[(t-1)*2-1]],dat.pred, Gcomp = TRUE, mppred = TRUE)$samp_treefit
      
    } else if(p==2){
      
      ## gcd
      dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, gcd)
      postY = predict_bart(Smod[[(t-1)*2-1]],dat.pred, Gcomp = TRUE)$samp_y
      Smat[[t]] = postY
      Smat[[t]][which(postY==0)] = yl[1]
      Smat[[t]][which(postY==1)] = yl[2]
      
      remove(postY)
      
    }
    
    ### S_{t-1} = engaged
    ytab = table(fitdat[get(theY_ind_lag1)==1,]$outcome)
    yl = as.numeric(names(ytab)[ytab > 0])
    p = length(yl) # number of outcome categories
    
    if(p>2){

      ## gcd
      tmp1 = predict_bart(Smod[[(t-1)*2]],dat.pred, Gcomp = TRUE, mppred = TRUE)$samp_treefit
      eng_ind = which(Smat[[t-1]] == 2)
      Smat[[t]][eng_ind] = tmp1[eng_ind]
      remove(tmp1)
      
      done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
      Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
      
    } else if(p==2){
      
      ## gcd
      postY = predict_bart(Smod[[(t-1)*2]],dat.pred, Gcomp = TRUE)$samp_y
      tmp1 = postY
      tmp1[which(postY==0)] = yl[1]
      tmp1[which(postY==1)] = yl[2]
      
      remove(postY)
      
      eng_ind = which(Smat[[t-1]] == 2)
      Smat[[t]][eng_ind] = tmp1[eng_ind]
      remove(tmp1)
      
      done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
      Smat[[t]][done_ind] = Smat[[t-1]][done_ind] 
    }

    ct = ct+1
    print(t)
  }
}
proc.time() - ptm

save(Xmat,Smat,allres,file = paste0("CondSModel_screen",screennum,"_ct",ct,"_ntr",ntr,".RData"))


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



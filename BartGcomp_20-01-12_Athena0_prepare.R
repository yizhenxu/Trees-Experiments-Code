#setwd("./Topic2")
library(GcompBART)
library(MNP)
library(lubridate)# date processing
library(data.table)
library(zoo)# carry backward forward
library(foreign)#stata read in
library(nnet)# multinomial log-linear models
library(Hmisc)#rmultinom function

source("Functions_new.R")
summ1 = function(M, s, op){# M: testn rows, nd cols
  vecp = apply(M, 2, function(x) mean(x==s))
  res = mean(vecp)
  
  if(op){
    svecp = sort(vecp)
    lp = length(svecp)
    ci = svecp[c(ceiling(lp*0.1), floor(lp*0.9))]
    res = c(res,ci)
  }
  
  return(res)
}

summ = function(M, s, op){# M: testn rows, nd cols
  vecp = apply(M, 2, function(x) mean(x==s))
  res = mean(vecp)
  
  if(op){
    svecp = sort(vecp)
    lp = length(svecp)
    ci = svecp[c(ceiling(lp*0.05), floor(lp*0.95))]
    res = c(res,ci)
  }
  
  return(res)
}

taby = function(y){
  return(unlist(lapply(1:4, function(s) mean(y==s))))
}

#get the outcome level of the largest frequency 
ymode = function(vec){
  names(which.max(table(vec))) 
}

accuracyfun = function(bmpfit,testd,trdy,tedy){
  # posterior prediction by the entire distribution
  # mode accuracy, mean accuracy
  res = rep(NA, 4); names(res) = c("Train Mean", "Train Mode", "Test Mean", "Test Mode")
  nd = bmpfit$ndraws
  # train prediction
  BTr = bmpfit$samp_y
  yTr = matrix(rep(trdy, nd),ncol = nd)
  PmodeTr = apply(BTr, 1, ymode)
  res[1] = mean(BTr == yTr)
  res[2] = mean(PmodeTr == trdy)
  # test prediction
  BTe = predict_bart(obj = bmpfit, newdata = testd)$samp_y
  yTe = matrix(rep(tedy, nd),ncol = nd)
  PmodeTe = apply(BTe, 1, ymode)
  res[3] = mean(BTe == yTe)
  res[4] = mean(PmodeTe == tedy)
  return(res)
}


summRL = function(M){# M: testn rows, nd cols
  vecp = apply(M, 2, function(x) mean(x))
  res = mean(vecp)
  
  svecp = sort(vecp)
  lp = length(svecp)
  ci = svecp[c(ceiling(lp*0.05), floor(lp*0.95))]
  res = c(res,ci)
  
  return(res)
}

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
         "weightinit_RX","weightinit_R",
         "log_vlplus1_init_RX","log_vlplus1_init_R","vl_0_1")


theX = c("new_cupd_1","log_cd4plus1") # LHS

theX_cov = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R") # RHS (from createind and createmis)
# for X1, use (X[0],X[-1])
theX_cov_1 = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R","log_cd4preBplus1_RX","log_cd4preBplus1_R") # RHS (from createind and createmis)

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

theA = "onarv_cf_1"

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")
theX_cov_1_lag1 = c(theX_cov_lag1,"log_cd4preBplus1_RX","log_cd4preBplus1_R")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

load("20-01-12_newdat.RData") # defined in BartGcomp_19-12-13_Reclean_AMPATH_data script

#OUTCOME: 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
newdat$outcome = as.numeric(newdat$outcome) #IMPORTANT! de-factorize outcome

######################################################################################################################
### optional: pick nsub people
seed = 1
nsub =   max(newdat$myID)

### Separate data into training set and test set

ntr = 50000 # floor(nsub * 0.7) # sample size of training set
#nte = 10000
ngc = 10000

set.seed(seed)
samp = sample(1:nsub, nsub, replace = F)

trid = samp[1:ntr]  
trid = sort(trid)
trd = copy(newdat[myID %in% trid,])

#teid = samp[(ntr+1):(ntr+nte)]  
teid = samp[(ntr+1):nsub]  
teid = sort(teid)
ted = copy(newdat[myID %in% teid,])

#gcid = samp[(ntr+nte+1):nsub]  
#gcid = sort(gcid)
#gcd = copy(newdat[myID %in% gcid,])
#gcd = copy(newdat)
#ngc = length(unique(gcd$myID))

gcid = sample(1:nsub, ngc, replace = F)
gcid = sort(gcid)
gcd = copy(newdat[myID %in% gcid,])
##############################################################
##############################################################
### Build models based on trd (training data)
#parell = TRUE # TRUE if on server, FALSE otherwise
#ncore = 5 # NULL if parall = FALSE

#maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15

#-------------------------------

#maxt= 10 # cannot fit for t = 11, sample size not fittable
if(setting == 0){nd = 2000;nb = 2000;nt = 50}
if(setting == 1){nb = 5000; nd = 5000; nt = 50}
if(setting == 2){nb = 5000; nd = 5000; nt = 100}
if(setting == 3){nb = 10000; nd = 5000; nt = 50}
if(setting == 4){nb = 10000; nd = 5000; nt = 100}
#
basey = 2
CatSeed = 10

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

Prior_mult = function(p, ntree){
  return(list(nu = p-1+1, V = diag(p-1),
              ntrees = ntree,
              kfac = 2,
              pswap = 0.1, pbd = 0.5, pb = 0.25,
              alpha = 0.95, beta = 2.0,
              nc = 100, minobsnode = 10))
}

Mcmc_mult = function(p,w0 = NULL,sig0 = NULL,nb, nd){
  res = list(sigma0 = diag(p-1), burn = nb, ndraws = nd, nSigDr = 50)  
  if(!is.null(w0)){
    res = append(res, list("w0" = w0), length(res))
  } 
  if(!is.null(sig0)){
    res$sigma0 = sig0
  }
  return(res)
}

##############################################################
##############################################################
if(0){
  
  #combine models
  maxt = 10
  mod1L =  vector("list", maxt) # binary models
  mod2L =  vector("list", maxt) # continuous models
  mod3_D =  vector("list", maxt-1) # categorical models for those previously disengaged
  mod3_E =  vector("list", maxt) # categorical models for those previously engaged
  
  modfml =  vector("list", (length(theX)+1)*maxt) 
  
  for(t in 1:maxt){
    
    ### Binary and Continuous
    load( paste0("ContBin20-01-12_t",t,"_setting",setting,".RData") ) 
    modfml[[3*(t-1)+1]] = fml1; mod1L[[t]] = mod1 # binary
    modfml[[3*(t-1)+2]] = fml2; mod2L[[t]] = mod2 # continuous
    remove(fml1, mod1, fml2, mod2)
    
    ### Categorical
    load(paste0("Cat20-01-12_t",t,"_setting",setting,".RData")) 
    modfml[[3*t]] = fml
    if(t==1){
      mod3_E[[1]] = Bmod2
    } else {
      mod3_D[[t-1]] = Bmod2[[1]]
      mod3_E[[t]] = Bmod2[[2]]
    }
    remove(fml, res, baseyl,  Bmod2)
  }#for t
  
  
  #
  
  tab = matrix(0,19,5)
  k=1
  for(t in 1:maxt){
    fitdat = trd[time==200*t,]
    if(t==1){
      
      ytab = table(fitdat$outcome)
      ind = as.numeric(names(ytab))
      tab[k,ind] = ytab
      tab[k, 5] = nrow(fitdat)
      k=k+1
      
    } else {
      
      for(prevS in c(0,1)){ #### S_{t-1} = 0 - disengaged, 1 - engaged
        
        d = fitdat[get(theY_ind_lag1) == prevS,]
        ytab = table(d$outcome)
        ind = as.numeric(names(ytab))
        tab[k,ind] = ytab
        tab[k, 5] = nrow(d)
        k=k+1
      } # prevS
      
    }
   
  }
  ################################################
  ### Data excerpt
  varlist = c("myID","time","outcome","onarv_cf_1","new_cupd_1", "log_cd4plus1_RX","log_cd4preBplus1_RX",
    "log_cd4preBplus1_R","ageinit",
    "male_1","yrenrol","TravelTime_RX","whostageinit_RX","whostageinit_R" ,
    "marriedinit_1","heightinit_RX","weightinit_RX","log_vlplus1_init_RX","vl_0_1")
  d = newdat[,varlist,with=F]
  ind = d[time == 0 &new_cupd_1==1 & log_cd4preBplus1_R==0,]
  see = d[myID %in% ind$myID,]
  
  rbind(d[myID == 3662,], d[myID==60050,])
  
  
}
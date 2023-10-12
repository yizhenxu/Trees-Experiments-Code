######################################################################################################################
for(i in 9:18) Smod[[i]] = Smod1[[i]]
for(i in 16:length(modfml)) modfml[i] = modfml1[i]
######################################################################################################################

#setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)
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

maxt = 10
basey1 = "2"
basey = c(rep("2", 2*(maxt-1)-1), "1") 

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


mnppred = function(d, m, fml, basey, fitdat){
  
  ytab = table(fitdat$outcome)
  ylevels = as.numeric(names(ytab)[ytab > 0])
  
  fv = fmlvars(fml)
  dat.pred = d[, fv$RHS[-1], with=F]
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
  return(PM)
}

yhat = function(vec){
  names(which.max(table(vec))) 
}


######################################################################################################################
###BART-PBART-MPBART
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

ptm <- proc.time()
allres = newAcc_bart(mod, modfml, ted, times, types, intypes, maxt, basey1, basey)

######################################################################################################################
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
#
sink()

#setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)
######################################################################################################################
### ZEUS script::

sink(paste0("output_",screennum,".txt")) #write all output to txt -- source('script.R',echo=T) to write input too

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

######################################################################################################################
### Model building

parell = TRUE # TRUE if on server, FALSE otherwise
ncore = 5 # NULL if parall = FALSE

maxt= 10 # for identifiability of parameters. No xfer or death in the observed data at t=15

#-------------------------------
setting = 2
maxt= 11
ndraw = 1100 
nburn = 100
nthin = 0 # number of discarded samples between two valid samples
#basey = 2 # disengagement as reference level
basey1 = "2" #basey for t = 1
basey = c(rep("2", 2*(maxt-1)-1), "1") #basey for t>1
nchain = 1 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain

#-------------------------------

tmp = (length(theX)+1)*maxt
mod =  vector("list", tmp) # replace the third model at each t>1 with 2 models in Smod
modfml =  vector("list", tmp) 
Smod =  vector("list", 2*(maxt-1)) 

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
    
    #RHS = c(theV, theA_lag1, theX_cov, theY_ind_lag1)
    RHS = c(theV, theA_lag1, theX_cov)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    ### S_{t-1} = disengaged
    Smod[[(t-1)*2-1]] = mnp_chains(as.formula(fml), fitdat[get(theY_ind_lag1)==0,], basey[(t-1)*2-1], ndraw, nburn, nthin, nchain, doparallel = FALSE)
    
    ### S_{t-1} = engaged
    Smod[[(t-1)*2]] = mnp_chains(as.formula(fml), fitdat[get(theY_ind_lag1)==1,], basey[(t-1)*2], ndraw, nburn, nthin, nchain, doparallel = FALSE)
    
    ct = ct+1
    print(t)
  }
}
proc.time() - ptm

save(mod,Smod,modfml,file = paste0("BayesModwCD4_screen",screennum,".RData"))

### Outcome simulation 

A = 0 #theA_lag1
nlatent = length(unique(gcd$outcome)) - 1
ylevels = as.numeric(levels(gcd$outcome))

Smat = vector("list", maxt) 
Xmat = vector("list", maxt)

for(i in 1:maxt){
  Xmat[[i]] = array(NA, dim = c(ngc, npost, length(theX)))
}

Vdist = cbind(1,gcd[time==0,theV,with=F])
Vdist = Vdist[rep(1:ngc,npost),]

ptm <- proc.time()

ct = 1

for(t in 1:maxt){
  
  if(t==1){
    
    print(ct)
    #[L1 | V, A0, X0, L0]
    
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, time = t, modhere = modfml[[ct]], simX = Xmat, policy_lag = A, dathere = gcd)
    
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
    sig = mod[[ct]][,ncol(mod[[ct]])]
    tmp = rnorm(c(bX))
    yhat = (tmp + c(bX))*sqrt(rep(sig, each = nsub))
    postY = matrix(yhat, ncol = npost)
    
    L = Xmat[[t]][,,1]
    mis_ind = which(L == 0) # people who have no CD4 update
    postY[mis_ind] = NA
    Xmat[[t]][,,2] = postY
    
    ct = ct+1
    
    print(ct)
    #[S1 | V, A0, X1, L1]
    
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, time = t, modhere = modfml[[ct]], simX = Xmat, policy_lag = A, dathere = gcd)
    Smat[[t]] = MNPpredict(m = mod[[ct]], Xmat = dat.pred, nsub, basey, ylevels, doparallel = parell)# basey is disengagement
    
    ct = ct+1
    
    print(paste0("t=",t))
    
  } else {
    
    print(ct)
    #[Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, yl = ylevels, time = t, modhere = modfml[[ct]], simX = Xmat, simY = Smat, policy_lag = A, dathere = gcd)
    Xmat[[t]][,,1] = MCMCpredict(type = "logit", m = mod[[ct]], nsub, Xhere = dat.pred)
    
    ct = ct+1
    
    print(ct)
    #[Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, yl = ylevels, time = t, modhere = modfml[[ct]], simX = Xmat, simY = Smat, policy_lag = A, dathere = gcd)
    postY = MCMCpredict(type = "regress", m = mod[[ct]], nsub, Xhere = dat.pred)
    
    L = Xmat[[t]][,,1]
    X_prev = Xmat[[t-1]][,,2]
    mis_ind = which(L == 0) # people who have no CD4 update keeps the old value
    postY[mis_ind] = X_prev[mis_ind]
    Xmat[[t]][,,2] = postY
    
    ct = ct+1
    
    print(ct)
    #[St | V, A_{t-1}, S_{t-1}, Xt, Lt]
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, yl = ylevels, time = t, modhere = modfml[[ct]], simX = Xmat, simY = Smat, policy_lag = A, dathere = gcd)
    
    fitdat = trd[time==200*t,]
    fitdat1 =  fitdat[get(theY_ind_lag1)==0,]
    ytab = table(fitdat1$outcome)
    ylevels = as.numeric(names(ytab)[ytab > 0])
    
    Smat[[t]] = MNPpredict(m = Smod[[(t-1)*2-1]], Xmat = dat.pred, nsub, basey[(t-1)*2-1], ylevels, doparallel = parell)# basey is disengagement
    
    fitdat1 =  fitdat[get(theY_ind_lag1)==1,]
    ytab = table(fitdat1$outcome)
    ylevels = as.numeric(names(ytab)[ytab > 0])
    tmp1 = MNPpredict(m = Smod[[(t-1)*2]], Xmat = dat.pred, nsub, basey[(t-1)*2], ylevels, doparallel = parell)# basey is disengagement
    
    eng_ind = which(Smat[[t-1]] == 2)
    Smat[[t]][eng_ind] = tmp1[eng_ind]
    remove(tmp1)
    
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
    
    ct = ct+1
    
    print(paste0("t=",t))
    
    
  }
  
}

proc.time() - ptm
#
save(Xmat,Smat,file = paste0("BayesSmat_screen",screennum,".RData"))

#load("Smat_setting2.RData")# in ZEUS
##############################################################
### Evaluation

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

newAcc_lmmnp = function(mod, Smod, modfml, ted, times, types, intypes, maxt, basey1, basey, npost){
  
  lenmod = length(mod)
  nummod = length(types)
  
  E = vector("list", nummod)
  
  for(i in 1:nummod){
    if(types[i] == "mult"){
      E[[i]] = matrix(NA,6,maxt)
      rownames(E[[i]]) = c("y=1","y=2","y=3","y=4","average","mode")
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
    
    dat.pred = ted[time==200*t,fv$RHS[-1],with=F]
    dat.pred = cbind(1, dat.pred)
    colnames(dat.pred)[1] = fv$RHS[1]
    
    if(ct %% 3 == 0){
      
      fitdat = trd[time==200*t,]
      d = ted[time==200*t,]
      
      if(t == 1){
        PM = mnppred(d, mod[[ct]], modfml[[ct]], basey1, fitdat)
      } else{
        PM = mnppred(d, Smod[[(t-1)*2-1]], modfml[[ct]], basey[(t-1)*2-1], fitdat[get(theY_ind_lag1)==0,])
        PM1 = mnppred(d, Smod[[(t-1)*2]], modfml[[ct]], basey[(t-1)*2], fitdat[get(theY_ind_lag1)==1,])
        ind = which(d[,theY_ind_lag1,with = F]==1)#check
        PM[ind,] = PM1[ind,] 
      }
      
      
    }
    
    
    if( types[intypes[ct]] == "binary"){
      
      yM = matrix(rep(y, npost),ncol = npost)

      bX = tcrossprod(as.matrix(dat.pred), mod[[ct]])
      predp = logistic(bX)
      PM = matrix(rbinom(length(predp),1,c(predp)), ncol = npost)
      
      PY = apply(PM, 1, yhat)
      E[[intypes[ct]]][1,t] = mean(PM == yM)
      E[[intypes[ct]]][2,t] = mean(PY == y)
      
    } else if(types[intypes[ct]] == "mult"){
      
      yM = matrix(rep(y, npost),ncol = npost)
      
      for(k in 1:4){
        
        E[[intypes[ct]]][k,t] = mean(((PM != k) & (yM != k)) | (PM==yM))
        
      }
      E[[intypes[ct]]][5,t] = mean(PM == yM)
      
      PY = apply(PM, 1, yhat)
      E[[intypes[ct]]][6,t] = mean(PY == y)
      
    } else if(types[intypes[ct]] == "continuous"){
      
      upd = unlist(ted[time==200*t, theX[1],with = F])
      ind = which(upd==1)
      y1 = y[ind]
      dat.pred1 = dat.pred[ind,]
      
      bX = tcrossprod(as.matrix(dat.pred1), mod[[ct]][,1:(ncol(mod[[ct]])-1)])
      sig = sqrt(mod[[ct]][,ncol(mod[[ct]])])
      tmp = rnorm(c(bX))
      yhat = tmp *rep(sig, each = nrow(bX)) + c(bX)
      PM = matrix(yhat, ncol = npost)
      
      yM = matrix(rep(y1, npost),ncol = npost)
      
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

allres = newAcc_lmmnp(mod, modfml, ted, times, types, intypes, maxt, basey1, basey, npoost)

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
save(counts,file = paste0("Bayescounts_setting",setting,".RData"))

load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS RData/Bayescounts_setting2.RData")


png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\BayesCD4SimEngA",A,".png"),width = 600, height = 400)
png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\BayesCD4SimEngA",A,".png"),width = 600, height = 400)

barplot(counts, main=paste0("Bayesian Simulated Engagement Status at a=",A,"\n Time-varying Covariates: CD4 info"),
        xlab="Days since baseline", ylab="Number of observations", col=1:4)legend("bottomleft", 
                                                                                  legend = rownames(counts), 
                                                                                  fill = 1:4)
dev.off()

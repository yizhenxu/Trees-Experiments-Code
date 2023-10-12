# Notes 
# Compare Accuracies for two models combined AND each model separately
#1. regular setting with nt=50 for all levels, try: (a) bylatent = 1, (b) bylatent = 0, (c) mympbartmod2
#2. run covariate importance test for bylatent = 1, pick covariates for small categories and rerun for nt=50
#3. number of trees, (a,b) with a for big levels, b for small levels, try: (a) c(100, 20), (b) c(50,20), (c) c(50,40)
#4. try pick covariates and diff nt

setting = 6

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


Prior_mult = function(p, nt){
  return(list(nu = p-1+3, V = diag(p-1),
              ntrees = nt,
              kfac = 2,
              pswap = 0.1, pbd = 0.5, pb = 0.25,
              alpha = 0.95, beta = 2.0,
              nc = 100, minobsnode = 10))
}

Mcmc_mult = function(p,w0 = NULL,sig0 = NULL){
  res = list(sigma0 = diag(p-1), burn = nb, ndraws = nd, nSigDr = 50)  
  if(!is.null(w0)){
    res = append(res, list("w0" = w0), length(res))
  } 
  if(!is.null(sig0)){
    res$sigma0 = sig0
  }
  return(res)
}


#assign big/small to categories as level-specific ntrees
#big when level frequency exceeds prop
if(0){
  mytrees = function(big, small, ytab, refy, prop){
    yt = as.numeric(ytab)
    yn = names(ytab)
    rind = which(yn==refy)
    yp = yt/sum(yt)
    res = rep(big, length(yt))
    res[yp < prop] = small
    res = res[-rind]
    return(res)
  }
}

mytrees = function(big, small, ytab, refy, num){
  yt = as.numeric(ytab)
  yn = names(ytab)
  rind = which(yn==refy)
  res = rep(big, length(yt))
  res[yt < num] = small
  res = res[-rind]
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
         "log_weight_init_RX","log_weight_init_R",
         "log_vlplus1_init_RX","log_vlplus1_init_R","vl_0_init_1")


theX = c("new_cupd_1","log_cd4plus1") # LHS

theX_cov = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R") # RHS (from createind and createmis)

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

theA = "onarv_cf_1"

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

load("newdat.RData") # defined in 18-06-22 script
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

### Separate data into training set and test set

ntr = 50000 # floor(nsub * 0.7) # sample size of training set
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
ngc = length(unique(gcid))
##############################################################
##############################################################
### Build models based on trd (training data)
#parell = TRUE # TRUE if on server, FALSE otherwise
#ncore = 5 # NULL if parall = FALSE

#maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15

#-------------------------------
#maxt= 10 # cannot fit for t = 11, sample size not fittable
nd = 2000
nb = 1000
nf = 500
blt = F #bytree F, bylatent T


nburn = nb
nthin = 3 # number of discarded samples between two valid samples
basey = "2" # engagement as reference level
nchain = 4 # number of MCMC chains
npost = nd 
ndraw = ceiling(nd/nchain) *(nthin+1) + nb

t = 9
fitdat = trd[time==200*t,]
print(dim(fitdat))

Bmod2 = vector("list", 2)
baseyl = rep(0,2)
ntL = vector("list", 2)

# [St | V, A_{t-1}, Xt, Lt, S_{t-1}]

RHS = c(theV, theA_lag1, theX_cov)
RHS = fmlfilter(fitdat,RHS)

fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))

for(prevS in c(0,1)){ #### S_{t-1} = 0 - disengaged, 1 - engaged
  
  ind = prevS+1
  
  ## trd
  d = fitdat[get(theY_ind_lag1) == prevS,]
  ytab = table(d$outcome)
  print(ytab)
  
  ylevels = as.numeric(names(ytab)[ytab > 0])
  p = length(ylevels) # number of outcome categories
  
  yprop = ytab/sum(ytab)
  print(rbind(ytab,yprop))
  
  Yname = theY
  Xname = vector("list", p-1)
  for(i in 1:(p-1)){
    Xname[[i]] = RHS
  }
  
  mod =  tryCatch(mnp_chains(fml, d, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE),
                  error = function(e) NULL )
  if(length(mod)==0){#ref eng fail, try diseng
    mod =  tryCatch(mnp_chains(fml, d, "1", ndraw, nburn, nthin, nchain, doparallel = FALSE),
                    error = function(e) NULL )
    if(length(mod)==0){
      baseyl[ind] = NA
      doit = 0
    } else {
      baseyl[ind] = "1"
      doit = 1
    }
  } else {
    baseyl[ind] = basey
    doit = 1
  }
  
  if(doit){
    ntL[[ind]] = mytrees(big=200,small=50,ytab,refy=baseyl[ind],num=200)
    
    ###
    fv = fmlvars(fml);
    dat.pred = d[, fv$RHS[-1], with=F]
    dat.pred = cbind(1, dat.pred)
    colnames(dat.pred)[1] = fv$RHS[1]
    
    latenty = sort(ylevels[-which(ylevels==basey)])
    nlatent = length(latenty)
    
    # simulate from multivariate normal (mu,Sigma)
    
    nhalfSig = nlatent*(nlatent - 1)/2 + nlatent
    sigmas = mod[, (ncol(mod) - nhalfSig + 1): ncol(mod)]
    params = mod[, 1:(ncol(mod) - nhalfSig)]
    
    BbyY = getParamByY(params,nlatent)
    
    mu = lapply(1:nlatent, function(yl) tcrossprod(as.matrix(dat.pred), BbyY[[yl]]) )
    mu_array = simplify2array(mu) # 3 dims: nsub x npost x nlatent
    
    m = mu_array[,npost,]
    sig = getSymat(sigmas[npost,], nlatent)
    w = mvrnorm(nrow(m), rep(0, nlatent), sig) + m
    
    ## diagnostic for all 
    Bmod2[[ind]] = model_bart(Yname = Yname, Xname = Xname, data = d, type = "multinomial",
                              base = baseyl[ind],
                              Prior = Prior_mult(p = p, nt = ntL[[ind]]),
                              Mcmc = Mcmc_mult(p = p, w0 = w, sig0 = sig),
                              correction = FALSE,
                              fitMNP = nf, bylatent = blt)
    
    
  } #doit
} # prevS

save(fml,baseyl,ntL,Yname,Xname, Bmod2, file = paste0("Mod3Run",setting,".RData"))


##############################################################
##############################################################


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

res = matrix(NA, ncol = 4, nrow = 4)

# the first dim: y, mean post, 5% post, 95% post
# the second dim: eng, diseng, xfer, death
d = ted[time==200*t,]
y =  unlist(ted[time==200*t,Yname,with = F])

PB = predict_bart(Bmod2[[1]], d)
PM = PB$samp_y
#PMT = PB$samp_treefit

PB1 = predict_bart(Bmod2[[2]], d)
PM1 = PB1$samp_y
#PMT1 = PB1$samp_treefit

ind = which(d[,theY_ind_lag1,with = F]==1)#check
PM[ind,] = PM1[ind,]
#PMT[ind,] = PMT1[ind,]

res[2:4,] = matrix(unlist(lapply(1:4, function(s) summ(PM,s,1))), ncol=4)

res[1, ] = taby(y)


save(res, file = paste0("res",setting,".RData"))

##################################################

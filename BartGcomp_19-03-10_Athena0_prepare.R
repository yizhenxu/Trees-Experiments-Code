# No pre-fit by MNP
# bytree back to the MPBART version
# burn 1000, draw 4000, 5000 in total

#setwd("./Topic2")
library(GcompBART)
library(MNP)
library(lubridate)# date processing
library(data.table)
library(zoo)# carry backward forward
library(foreign)#stata read in
library(nnet)# multinomial log-linear models
library(Hmisc)#rmultinom function
library(mlogit)

library("parallel")

source("Functions_new.R")

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


getInitMNP = function(s, fml, d, basey, ndraw, nburn, nthin, nchain, npost){
  set.seed(s)
  mod =  tryCatch(mnp_chains(fml, d, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE),
                  error = function(e) NULL )
  if(length(mod)==0){
    set.seed(s)
    mod =  tryCatch(mnp_chains(fml, d, "1", ndraw, nburn, nthin, nchain, doparallel = FALSE),
                    error = function(e) NULL )
    if(length(mod)==0){
      yref = NA
      doit = 0
    } else {
      yref = "1"
      doit = 1
    }
  } else {
    yref = basey
    doit = 1
  }
  w = NA; sig = NA;
  if(doit){
    
    ###
    fv = fmlvars(fml);
    dat.pred = d[, fv$RHS[-1], with=F]
    dat.pred = cbind(1, dat.pred)
    colnames(dat.pred)[1] = fv$RHS[1]
    
    latenty = sort(ylevels[-which(ylevels==yref)])
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
  }
    
  return(list(doit = doit, yref = yref, w = w, sig = sig))
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

maxt= 10 # cannot fit for t = 11, sample size not fittable
nd = 4000
nb = 1000
nt = 50
nf = 0
blt = F #bytree F, bylatent T
CatSeed = 888888

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
  return(list(nu = p-1+3, V = diag(p-1),
              ntrees = ntree,
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
#-------------------------------
# not doing fitMNP for this round
nburn = nb
nthin = 3 # number of discarded samples between two valid samples
basey = "2" # engagement as reference level
nchain = 4 # number of MCMC chains
npost = nd 
ndraw = ceiling(nd/nchain) *(nthin+1) + nb

##############################################################
##############################################################
### summary plot of newdat outcome
if(0){
  a = newdat$outcome
  levels(a) = c("Disengaged","Engaged","Transferred out","Died")
  a = relevel(a, "Engaged")
  counts <- table(a,newdat$time)
  barplot(counts, main="Observed Outcome Status",
          xlab="Days since baseline", ylab="Number of observations", col=1:length(levels(a)),
          legend = rownames(counts))
  
  library(ggplot2)
  cPalette <- c("#0072B2","#E69F00","#FF3333","#339966")
  pathto = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft3\\Accuracy Comparison - Diag Plots\\19-03-05\\Plots"
  
  tab = data.frame(Status = factor( rep(rownames(counts), each = ncol(counts)), levels = c("Died","Transferred out","Disengaged","Engaged") ),
                   X = rep(as.numeric(colnames(counts)), nrow(counts)),
                   Y = c(t(counts)))
  ggplot(data = tab, aes(x=X, y=Y, fill=Status)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values = c("Died"=cPalette[1],"Transferred out"=cPalette[2],"Disengaged"=cPalette[3],"Engaged"=cPalette[4])) +
    ggtitle("Observed Outcome Status")+theme(plot.title = element_text(hjust = 0.5))+
    xlab("Days since baseline") + ylab("Number of observations")
  ggsave(paste0(pathto,"\\ObaEng.png"),width = 5.5, height = 3.5)
}
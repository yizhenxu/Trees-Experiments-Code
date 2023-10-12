#devtools::install_github("yizhenxu/GcompBART")
#remove.packages("GcompBART")

#setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")

#setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)

t = 10
#sc = c(54060 ,54098,54081,54110,54143,54125)
#sc = c(105607 ,
#       105434 ,
#       105635 ,
#       105582 ,
#       105663 ,
#       105649 ,
#       117991 ,
#       117976 ,
#       117958)
screennum = sc[t-1]

t=8
sc = c(12071 ,
       11755 ,
       12052 ,
       11902 ,
       12085 ,
       12099 ,
       11888 ,
       12195)

screennum = sc[t]

setwd("./Topic2");source("script_54060.R",echo=T)
source("script_12071.R",echo=T)
####################################################################################################################
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

### Separate data into training set and test set

ntr = 50000 # floor(nsub * 0.7) # sample size of training set
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
##############################################################
### Build models based on trd (training data)
#parell = TRUE # TRUE if on server, FALSE otherwise
#ncore = 5 # NULL if parall = FALSE

#maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15

#-------------------------------
setting = 1
maxt= 11
nd = 2000
nb = 100
nt = 50
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

#ndraw = nd+nb 
#nburn = nb
#nthin = 0 # number of discarded samples between two valid samples
#basey = "2" # disengagement as reference level
#nchain = 1 # number of MCMC chains
#npost = floor((ndraw - nburn) / (nthin+1)) * nchain


nburn = nb
nthin = 3 # number of discarded samples between two valid samples
basey = "2" # engagement as reference level
nchain = 4 # number of MCMC chains
npost = nd 
ndraw = ceiling(nd/nchain) *(nthin+1) + nb
nf = 1000
##############################################################
##############################################################

t=1
fitdat = trd[time==200*t,]
# [S1 | V, A0, X1, L1]

ytab = table(fitdat$outcome)
ylevels = as.numeric(names(ytab)[ytab > 0])
p = length(ylevels) # number of outcome categories

RHS = c(theV,theA_lag1,theX_cov)
RHS = fmlfilter(fitdat,RHS)
fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))

mod = mnp_chains(fml, fitdat, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE)
###
fv = fmlvars(fml); d = fitdat
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

#save(mod,w,m,sig, file = "tmp.RData")

###################################################
## trd

## diagnostic for all 
Bmod2 = model_bart(as.formula(fml), data = fitdat, type = "multinomial2",
                   base = basey,
                   Prior = Prior_mult(p),
                   Mcmc = Mcmc_mult(p),
                   correction = FALSE,
                   w = w, sig = sig, fitMNP = nf)
###################################################

Bmod1 = model_bart(as.formula(fml), data = fitdat, type = "multinomial2",
                   base = basey,
                   Prior = Prior_mult(p),
                   Mcmc = Mcmc_mult(p),
                   correction = FALSE,
                   w = w, sig = sig, fitMNP = 1000)





##

Bmod0 = model_bart(as.formula(fml), data = fitdat, type = "multinomial",
                  base = basey,
                  Prior = Prior_mult(p),
                  Mcmc = Mcmc_mult(p),
                  correction = FALSE)


##

load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS CODE/Bmod0.RData")
load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS CODE/Bmod.RData")
load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS CODE/mod.RData")
par(mfrow=c(2,3))
j=1
plot(Bmod$Tree_Depth[j,],type="l",main=paste0("latent ",j))
lines(Bmod0$Tree_Depth[j,],col="red")

plot(Bmod$Tree_Num_Nodes[j,],type="l",main=paste0("latent ",j))
lines(Bmod0$Tree_Num_Nodes[j,],col="red")

plot(Bmod$Tree_Num_Leaves[j,],type="l",main=paste0("latent ",j))
lines(Bmod0$Tree_Num_Leaves[j,],col="red")

j=2
plot(Bmod$Tree_Depth[j,],type="l",main=paste0("latent ",j))
lines(Bmod0$Tree_Depth[j,],col="red")

plot(Bmod$Tree_Num_Nodes[j,],type="l",main=paste0("latent ",j))
lines(Bmod0$Tree_Num_Nodes[j,],col="red")

plot(Bmod$Tree_Num_Leaves[j,],type="l",main=paste0("latent ",j))
lines(Bmod0$Tree_Num_Leaves[j,],col="red")

#MPBART sig
par(mfrow=c(1,2))
plot(Bmod$sigmasample[(0:1999)*4+1],type="l",ylim = range(Bmod$sigmasample))
for(j in c(2,4)){
  lines(Bmod$sigmasample[(0:1999)*4+j],col=j)
}

plot(Bmod0$sigmasample[(0:1999)*4+1],type="l",ylim = range(Bmod0$sigmasample))
for(j in c(2,4)){
  lines(Bmod0$sigmasample[(0:1999)*4+j],col=j)
}



#MNP beta
a = 4;b = 10
par(mfrow = c(a,b))
for(i in 1:(a*b))
plot(mod[,i])

#MNP sig
a = 1;b = 3
par(mfrow = c(a,b))
for(i in 1:(a*b))
  plot(mod[,42+i])



##############################################################
##############################################################

t=2

source("script_54060.R",echo=T)
fitdat = trd[time==200*t,]
Bmod2 = vector("list",2)

# [St | V, A_{t-1}, Xt, Lt, S_{t-1}]

#RHS = c(theV, theA_lag1, theX_cov, theY_ind_lag1)
RHS = c(theV, theA_lag1, theX_cov)
RHS = fmlfilter(fitdat,RHS)
fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))

### S_{t-1} = disengaged
## trd
d = fitdat[get(theY_ind_lag1)==0,]
ytab = table(d$outcome)
ylevels = as.numeric(names(ytab)[ytab > 0])
p = length(ylevels) # number of outcome categories

mod = mnp_chains(fml, d, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE)
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
Bmod2[[1]] = model_bart(as.formula(fml), data = d, type = "multinomial2",
                   base = basey,
                   Prior = Prior_mult(p),
                   Mcmc = Mcmc_mult(p),
                   correction = FALSE,
                   w = w, sig = sig, fitMNP = nf)
###################################################

### S_{t-1} = engaged
## trd
d = fitdat[get(theY_ind_lag1)==1,]
ytab = table(d$outcome)
ylevels = as.numeric(names(ytab)[ytab > 0])
p = length(ylevels) # number of outcome categories

mod = mnp_chains(fml, d, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE)
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
Bmod2[[2]] = model_bart(as.formula(fml), data = d, type = "multinomial2",
                        base = basey,
                        Prior = Prior_mult(p),
                        Mcmc = Mcmc_mult(p),
                        correction = FALSE,
                        w = w, sig = sig, fitMNP = nf)
###################################################


#png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot\\Continuous_",t,".png"),width = 600, height = 400)
Bmod2$ndraws = ncol(Bmod2$Tree_Depth)
Bmod2$type = "multinomial"
png(paste0("Multinomial_S1_",t,".png"),width = 600, height = 400)
DiagPlot(Bmod2,1)
dev.off()

png(paste0("Multinomial_S1_",t,"_sig.png"),width = 600, height = 400)
plot(Bmod2$sigmasample)
dev.off()

# t>1
# previously disengaged
Bmod2[[1]]$ndraws = ncol(Bmod2[[1]]$Tree_Depth)
Bmod2[[1]]$type = "multinomial"
png(paste0("Multinomial_S0_",t,".png"),width = 600, height = 400)
DiagPlot(Bmod2[[1]],1)
dev.off()

png(paste0("Multinomial_S0_",t,"_sig.png"),width = 600, height = 400)
plot(Bmod2[[1]]$sigmasample)
dev.off()

# previously engaged
Bmod2[[2]]$ndraws = ncol(Bmod2[[2]]$Tree_Depth)
Bmod2[[2]]$type = "multinomial"
png(paste0("Multinomial_S1_",t,".png"),width = 600, height = 400)
DiagPlot(Bmod2[[2]],1)
dev.off()

png(paste0("Multinomial_S1_",t,"_sig.png"),width = 600, height = 400)
plot(Bmod2[[2]]$sigmasample)
dev.off()


library(lubridate)# date processing
library(data.table)
library(zoo)# carry backward forward
library(foreign)#stata read in
library(nnet)# multinomial log-linear models
library(Hmisc)#rmultinom function

library("parallel")

setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
#setwd("./Topic2");
source("Functions_new.R")
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

######################################################################################################################

#load("newdat.RData") # defined in 18-02-26 script
load("./makeGridData/newdat.RData")

newdat$outcome = as.numeric(newdat$outcome) #IMPORTANT! de-factorize outcome


### optional: pick 1000 people
nn = 1000
set.seed(1)
nsub = sum(newdat$time==0)
id = sample(1:nsub, nn, replace = F)
newdat = newdat[myID %in% id,]
###


#-------------------------------
setting = "only for this script"
maxt= 1
ndraw = 1000
nburn = 100
ntree = 50
ncut = 100
nsd = 50 # max repeated draw of Sigma for Mpbart
seed = 88
basey = 2 # the outcome level used in RHS of model fitting (theY_ind)
nmod = 3 # number of models

#-------------------------------

# initialization

set.seed(seed)

#runtime = rep(NA, nmod + 1)

modfml =  vector("list", nmod)
mod = vector("list", nmod)

Smat = vector("list", maxt)
Scov = vector("list", maxt)
Xmat = vector("list", maxt)

train_pred = vector("list", nmod)
train_y = vector("list", nmod)

for(i in 1:maxt){
  Xmat[[i]] = array(NA, dim = c(nsub, ndraw, length(theX)))
}

Vdist = cbind(newdat[time==0,theV,with=F])
Vdist = Vdist[rep(1:nsub,ndraw),]

##########################
### Model building

# not this -- [S1 | V, A0, X1, L1]
# use this -- [S1 | V, A0, X0, L0]
nsub = sum(newdat$time==0)
t = 1
fitdat = newdat[time==200*t,]

ytab = table(fitdat$outcome)
ylevels = as.numeric(names(ytab)[ytab > 0])
p = length(ylevels) # number of outcome categories


RHS = c(theV, theA_lag1, theX_cov_lag1)
# exclude lag X_R from models at t=1 because that is the same as lag L at baseline
RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
RHS = fmlfilter(fitdat,RHS)
fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))

policy_lag = rep(A, nsub)
dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, fml, Xmat, Smat, policy_lag, newdat)

mod[[1]] = model_bart(as.formula(fml), data = fitdat, type = "multinomial",
                      base = basey,
                      Prior = list(nu = p-1+3, V = diag(p-1),
                                   ntrees = ntree,
                                   kfac = 2,
                                   pswap = 0.1, pbd = 0.5, pb = 0.25,
                                   alpha = 0.95, beta = 2.0,
                                   nc = ncut, minobsnode = 10),
                      Mcmc = list(sigma0 = diag(p-1), burn = nburn, ndraws = ndraw,
                                  nSigDr = nsd))
png(filename = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-06-07_demo_GcompBART\\IncPlot.png", width = 1000, height = 600) 
op = par()$mar
par( mar=c(10, op[-1]) ) 
DiagPlot(mod[[1]], 0)
par(mar = op)
dev.off()

png(filename = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-06-07_demo_GcompBART\\DiagPlot.png", width = 700, height =900) 
op = par()$mar
par( mar=c(op[1],7, op[3:4]) ) 
par(cex.lab=2, cex.axis = 1.3, cex.main = 2)
DiagPlot(mod[[1]], 1,byrow=FALSE)
par(mar = op)
dev.off()

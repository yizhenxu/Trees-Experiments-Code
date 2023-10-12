#!/bin/bash
#SBATCH -J APacc
#SBATCH --time 72:00:00
#SBATCH --output=acc.out
#SBATCH --error=acc.err
#SBATCH --mail-type=all
#SBATCH --mail-user=xuyizhen00@gmail.com

#Rscript BartGcomp_19-11-26_mbart_AMPATHacc.R

## run the current R script

# sbatch mbash0.sh

rm(list=ls()) 

gc()

#devtools::install_github("yizhenxu/GcompBART")
#library(GcompBART)
library(data.table)
library(BART)

#ulimit::memory_limit(15000)

#library(unix)
#rlimit_all(); rlimit_as()
# set a soft limit
#lim = rlimit_as(2e10);print(lim) # 20G
# reset 
#rlimit_as(cur = lim$max)
setwd("./MPBART")


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
         "weightinit_RX","weightinit_R")

theX_cov = c("new_cupd_1","new_cd4_cf_RX","new_cd4_cf_R") # RHS (from createind and createmis)

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

theA = "onarv_cf_1"

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

load("20-01-12_newdat.RData") # defined in BartGcomp_19-12-13_Reclean_AMPATH_data script

#OUTCOME: 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
newdat$outcome = as.numeric(newdat$outcome) #IMPORTANT! de-factorize outcome

######################################################################################################################
#####################################
seed = 1
nsub =   max(newdat$myID)

### Separate data into training set and test set

ntr = 50000 # floor(nsub * 0.7) # sample size of training set

set.seed(seed)
samp = sample(1:nsub, nsub, replace = F)

trid = samp[1:ntr]  
trid = sort(trid)
trd = copy(newdat[myID %in% trid,])

#teid = samp[(ntr+1):(ntr+nte)]  
teid = samp[(ntr+1):nsub]  
teid = sort(teid)
ted = copy(newdat[myID %in% teid,])

### Restrict to t = 1
t = 1

trd = trd[time==200*t,]

ted = ted[time==200*t,]

# rename outcome to y
yind = which(colnames(trd)==theY)
colnames(trd)[yind] = "y"
yind = which(colnames(ted)==theY)
colnames(ted)[yind] = "y"
theY = "y"

print(dim(trd))
print(dim(ted))

# [S1 | V, A0, X1, L1]
fmlfilter = function(fitdat, RHS){
  todelete  = c()
  for(i in 1:length(RHS)){
    a = fitdat[,get(RHS[i])]
    if(length(unique(a))==1) todelete = c(todelete,i)
  }
  if(length(todelete)>0) RHS = RHS[-todelete]
  return(RHS)
}

RHS = c(theV,theA_lag1,theX_cov)
RHS = fmlfilter(trd,RHS)
fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))

#####################################
### Build models based on trd (training data)
#parell = TRUE # TRUE if on server, FALSE otherwise
#ncore = 5 # NULL if parall = FALSE

#maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15

#-------------------------------

#maxt= 10 # cannot fit for t = 11, sample size not fittable
nd = 10000 # posterior draws
nb = 10000 # burn in
nt = 100 # number of trees
p=length(unique(trd$y)) # number of outcome levels

CatSeed = 888888


#####################################

source("BartGcomp_19-11-22_Functions.R")

sfn = paste0("CK_AMPATHacc_b",nb,"d",nd,"_AMPATH.RData")

set.seed(CatSeed)

xtrain = as.matrix(trd[,RHS,with=F])
ytrain = as.matrix(trd[,theY,with=F])
cdat = rbind(trd,ted)
xtest  = as.matrix(cdat[,RHS,with=F])

post <- mbart(x.train = xtrain, y.train = ytrain, x.test = xtest, ndpost = nd, nskip = nb, keepevery = 1, ntree=nt)

tmp = accuracyfun_mbart(post, cdat$y, nrow(xtrain), nrow(xtest))

save(tmp, file = sfn)




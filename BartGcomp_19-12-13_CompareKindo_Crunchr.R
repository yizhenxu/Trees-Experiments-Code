#devtools::install_github("yizhenxu/GcompBART")
library(GcompBART)
library(data.table)
setwd("/home/idies/workspace/Storage/yxu143/persistent/MPBART")
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
nd = 30000 # posterior draws
nb = 50000 # burn in
nt = 100 # number of trees
p=length(unique(trd$y)) # number of outcome levels

CatSeed = 888888


#####################################

source("BartGcomp_19-11-22_Functions.R")

for(sim in 1:3){
  if(sim ==  1){KD = F; nuC = 1; Vm = diag(2); by = "1"; da = TRUE}
  if(sim ==  2){KD = F; nuC = 1; Vm = diag(2); by = "1"; da = FALSE}
  if(sim ==  3){KD = T; nuC = 1; Vm = diag(2); by = "1"}
  
  if(sim ==  4){KD = F; nuC = 1; Vm = diag(2); by = "2"; da = TRUE}
  if(sim ==  5){KD = F; nuC = 1; Vm = diag(2); by = "2"; da = FALSE}
  if(sim ==  6){KD = T; nuC = 1; Vm = diag(2); by = "2"}
  
  if(sim ==  7){KD = F; nuC = 1; Vm = diag(2); by = "4"; da = TRUE}
  if(sim ==  8){KD = F; nuC = 1; Vm = diag(2); by = "4"; da = FALSE}
  if(sim ==  9){KD = T; nuC = 1; Vm = diag(2); by = "4"}
  
  
  set.seed(CatSeed)
  
  # do_alpha2_prior = FALSE valid only when KD = F then both mympbart0 (Alg2) and mympbart1 (Alg1) reduce to Alg 3
  # TRUE: update alpha2 in both Step 1 & 3; FALSE: update alpha2 only in Step 3 (Alg3)
  
  if(KD == F){
    if(da){
      bmpfit0 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                           base = by,
                           Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                           Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                           correction = FALSE, Kindo = KD, do_alpha2_prior = da)
      
      save(bmpfit0, file = paste0("CK_AMPATH_b",nb,"d",nd,"s",sim,".RData"))
      
      acc0 = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
      save(acc0, file = paste0("CK_AMPATHacc_b",nb,"d",nd,"s",sim,".RData"))
      
    } else {
      bmpfit2 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                           base = by,
                           Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                           Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                           correction = FALSE, Kindo = KD, do_alpha2_prior = da)
      
      save(bmpfit2, file = paste0("CK_AMPATH_b",nb,"d",nd,"s",sim,".RData"))
      
      acc2 = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)
      save(acc2, file = paste0("CK_AMPATHacc_b",nb,"d",nd,"s",sim,".RData"))
    }
    
  }
  
  if(KD == T){
    bmpfit1 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                         base = by,
                         Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                         Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                         correction = FALSE, Kindo = KD)
    
    save(bmpfit1, file = paste0("CK_AMPATH_b",nb,"d",nd,"s",sim,".RData"))
    
    acc1 = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
    save(acc1, file = paste0("CK_AMPATHacc_b",nb,"d",nd,"s",sim,".RData"))
    
    
  }
  
}

#############################################################################

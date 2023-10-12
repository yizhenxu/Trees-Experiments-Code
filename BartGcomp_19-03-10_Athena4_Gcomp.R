rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
A = as.numeric(slurm_arrayid)
# theA_lag1 A = 0,1,2; A==2 means A = 1{CD4<=350}
print(paste0("A is ",A))

if(0){
  #!/bin/bash
  #SBATCH -J ParaGcomp
  #SBATCH -n 1
  #SBATCH -N 1
  #SBATCH --time 72:00:00
  #SBATCH --output=job_%A_%a.out
  #SBATCH --error=job_%A_%a.err
  #SBATCH --array=0-2
  #SBATCH --mail-type=all
  #SBATCH --mail-user=yizhen_xu@brown.edu
  
  #Rscript BartGcomp_19-03-10_Athena4_Gcomp.R
}

# TP plot
# Gcomp barplot
# diffCI

source("BartGcomp_19-03-10_Athena0_prepare.R")

load("Mod19-03-10.RData")

thinning = 3; ndthin = nd/(thinning + 1)
##############################################################
### Outcome simulation 

nlatent = length(unique(gcd$outcome)) - 1
nsub = sum(gcd$time==0)
ytab = table(gcd$outcome)
ylevels = as.numeric(names(ytab)[ytab > 0])

Smat = vector("list", maxt) 
Xmat = vector("list", maxt)

for(i in 1:maxt){
  Xmat[[i]] = array(NA, dim = c(nsub, ndthin, length(theX)))
}

Vdist = cbind(1,gcd[time==0,theV,with=F])
Vdist = Vdist[rep(1:nsub,ndthin),]

ptm <- proc.time()

for(t in 1:maxt){
  
  if(t==1){
    
    if(!(A %in% c(0,1))){
      #theA_lag1 as a function of theX_lag1 
      I0 = gcd[time==0, "new_cupd_1", with=F]$new_cupd_1
      C0 = gcd[time==0, "log_cd4plus1_RX", with=F]$log_cd4plus1_RX
      policy_lag = 1*(I0==1 & (C0<= log(351)))    #rep(A, nsub)
    } else {
      policy_lag = rep(A, nsub)
    }
    
    #[L1 | V, A0, X0, L0]
    
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, c(1,2,4), t, modfml[[3*(t-1)+1]], Xmat, Smat, policy_lag, gcd)
    
    Xmat[[t]][,,1] = predict_bart(mod1[[t]], thinning, dat.pred)$samp_y
    
    #[X1 | V, A0, X0, L0] 
    #dat.pred = makeX(Xname = theX, Vname = theV, time = t, modhere = mod[[ct]], simX = Xmat, policy_lag = A, dathere = gcd)
    # same RHS as in modeling L1
    postY = predict_bart(mod2[[t]], thinning, dat.pred)$samp_y
    
    L = Xmat[[t]][,,1]
    mis_ind = which(L == 0) # people who have no CD4 update
    postY[mis_ind] = NA
    Xmat[[t]][,,2] = postY
    remove(postY)
    
    #[S1 | V, A0, X1, L1]
    
    if(!(A %in% c(0,1))){
      policy_lag = policy_lag[rep(1:nsub,ndthin)]
    } else {
      policy_lag = rep(A, nrow(Vdist))
    }
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, c(1,2,4), t, modfml[[3*t]], Xmat, Smat, policy_lag, gcd)
    
    Smat[[t]] = predict_bart(mod3_E[[t]], thinning,dat.pred, Gcomp = TRUE)$samp_y
    
    print(paste0("t=",t))
    
  } else {
    
    if(!(A %in% c(0,1))){
      # define A_{t-1}
      I0 = Xmat[[t-1]][,,1]
      C0 = Xmat[[t-1]][,,2]
      dim(I0) = NULL
      dim(C0) = NULL
      policy_lag = 1*(I0==1 & (C0<= log(351)))    #rep(A, nrow(Vdist))
    }
    
    #[Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[3*(t-1)+1]], Xmat, Smat, policy_lag, gcd)
    Xmat[[t]][,,1] = predict_bart(mod1[[t]], thinning,dat.pred, Gcomp = TRUE)$samp_y
    
    #[Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[3*(t-1)+2]], Xmat, Smat, policy_lag, gcd)
    postY = predict_bart(mod2[[t]], thinning,dat.pred, Gcomp = TRUE)$samp_y
    
    L = Xmat[[t]][,,1]
    X_prev = Xmat[[t-1]][,,2]
    mis_ind = which(L == 0) # people who have no CD4 update keeps the old value
    postY[mis_ind] = X_prev[mis_ind]
    Xmat[[t]][,,2] = postY
    
    remove(postY)
    
    #[St | V, A_{t-1}, S_{t-1}, Xt, Lt]
    
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[3*t]], Xmat, Smat, policy_lag, gcd)
    
    Smat[[t]] = predict_bart(mod3_D[[t-1]], thinning,dat.pred, Gcomp = TRUE)$samp_y
    
    tmp1 = predict_bart(mod3_E[[t]], thinning,dat.pred, Gcomp = TRUE)$samp_y
    
    eng_ind = which(Smat[[t-1]] == 2)
    Smat[[t]][eng_ind] = tmp1[eng_ind]
    remove(tmp1)
    
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
    
    print(paste0("t=",t))
    
    
  }
  
}

proc.time() - ptm

save(Xmat,Smat,file = paste0("Gcomp_A",A,".RData"))

####################################################################
### TP plot
library(data.table)

tmp = unlist(Smat)

m = length(Smat[[1]])
tmp1 = c(rep(2,m),tmp[1:((maxt-1)*m)])

longdat = data.table(time = rep(200*(1:maxt),each = length(Smat[[1]])), outcome_lag = tmp1, outcome = tmp)
ldc = longdat[, .N, by = c("time", "outcome_lag", "outcome")]
ldc_denom = longdat[, .N, by = c("time","outcome_lag")]

ldc = ldc[-which(ldc$outcome_lag %in% c(3,4)),]
ldc = ldc[order(time, outcome_lag, outcome)]
ldc_denom = ldc_denom[-which(ldc_denom$outcome_lag %in% c(3,4)),]
ldc_denom = ldc_denom[order(time, outcome_lag)]

eldc = ldc[outcome_lag == 2,]
eldc = rbind(eldc,list(time = 200, outcome_lag = 2, outcome = 3, N = 0 ))
eldc = eldc[order(time, outcome, outcome_lag)]

dldc = ldc[outcome_lag == 1,]

eldc_denom = ldc_denom[outcome_lag == 2,]
dldc_denom = ldc_denom[outcome_lag == 1,]

engtp = matrix(eldc$N/rep(eldc_denom$N, each = 4), nrow = 4)
rownames(engtp) = c("diengaged","engaged","transferred out","died")
colnames(engtp) = unique(eldc$time)
dengtp = matrix(dldc$N/rep(dldc_denom$N, each = 3), nrow = 3)
rownames(dengtp) = c("diengaged","engaged","died")
colnames(dengtp) = unique(dldc$time)

save(engtp,dengtp,file = paste0("TP_A",A,".RData"))

####################################################################
###Gcomp bartplot

#tmp = unlist(Smat)
tmp1 = c(rep(2, m),tmp)

longdat = data.table(time = rep(200*(0:maxt),each = m), outcome = tmp1)
ldc = longdat[, .N, by = c("time","outcome")]

tt = sum(ldc$time==200)
if(tt == 3){
  ldc = rbind(ldc,list(time = c(0,0,0,200),outcome = c(1,3,4,3),N=c(0,0,0,0) ))
} else {
  ldc = rbind(ldc,list(time = c(0,0,0),outcome = c(1,3,4),N=c(0,0,0) ))
}

ldc = ldc[order(time,outcome)]
counts = matrix(ldc$N,nrow = 4)

ctmp = counts[1,]
counts[1,] = counts[2,]
counts[2,] = ctmp

rownames(counts) = c("engaged","diengaged","transferred out","died")
counts = simplify2array(lapply(1:ncol(counts), function(x) counts[,x]/sum(counts[,x])))
colnames(counts) = 200*(0:maxt)

save(counts,file = paste0("Counts_A",A,".RData"))

####################################################################
### Prepare for diffCI
summ = function(vec){
  return(unlist(lapply(1:4, function(x) mean(vec==x))))
}
engci = function(i, Smat){
  res = apply(Smat[[i]],2, summ)
  return(res)
}

mat = lapply(1:maxt, function(i) engci(i, Smat)) #in A = 1 and A = 1{CD4 < 350}  

save(mat,file = paste0("CImat_A",A,".RData"))


#!/bin/bash
#SBATCH -J A2t5
#SBATCH --time 72:00:00
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.err
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

#Rscript BartGcomp_20-01-12_Athena4_Gcomp_thin1000_2.R

rm(list=ls()) 

A=2
# theA_lag1 A = 0,1,2; A==2 means A = 1{CD4<=350}
print(paste0("A is ",A))


setting = 4; maxt = 10; 
source("BartGcomp_20-01-12_Athena0_prepare.R")
skip = 4; ndthin = floor(nd/(skip+1))

X0 = gcd[time ==0 ,log_cd4plus1_RX]
X0[gcd[time ==0,log_cd4plus1_R] == 1] = NA

X0Mat = X0[rep(1:nsub,ndthin)]
##############################################################
### Outcome simulation 

nlatent = length(unique(gcd$outcome)) - 1
ytab = table(gcd$outcome)
ylevels = as.numeric(names(ytab)[ytab > 0])

load("Gcomp_A2_t4_setting4_20-01-12.RData")#t5#

Vdist = cbind(1,gcd[time==0,theV,with=F])
Vdist = Vdist[rep(1:nsub,ndthin),]

I0 = gcd[time==0, log_cd4preBplus1_R] #missing indicator
C0 = gcd[time==0, log_cd4preBplus1_RX]
policy_lag = 1*(I0==0 & (C0<= log(351)))    #rep(A, nsub)
policy_lag = policy_lag[rep(1:nsub,ndthin)]
for(t in 2:4){#t5#
  A_prev = policy_lag
  I0 = Smat[[t-1]]
  if(t==2){
    C0 = X0Mat
    remove(X0Mat); gc()
  } else {
    C0 = Xmat[[t-2]][,,2]
  }
  dim(I0) = NULL
  dim(C0) = NULL
  policy_lag = 1*(A_prev== 1) + 1*(A_prev==0 & I0==2 & !is.na(C0) & (C0<= log(351)))  
}

for(t in 5:maxt){#t5#
  
  if(t==1){
    
    #[A0 | X(-1)]
    if(!(A %in% c(0,1))){
      #theA_lag1 as a function of theX_lag1 
      I0 = gcd[time==0, log_cd4preBplus1_R] #missing indicator
      C0 = gcd[time==0, log_cd4preBplus1_RX]
      policy_lag = 1*(I0==0 & (C0<= log(351)))    #rep(A, nsub)
    } else {
      policy_lag = rep(A, nsub)
    }
    
    load( paste0("ContBin20-01-12_t",t,"_setting",setting,".RData") )  
    
    #[R1 | V, A0, X0, X(-1)]
    
    theV0 = c(theV,"log_cd4preBplus1_RX","log_cd4preBplus1_R")
    dat.pred = makeX_bart(theX, theX_cov, theV0, Vdist, c(1,2,4), t,fml1, Xmat, Smat, policy_lag, gcd)
    
    Xmat[[t]][,,1] = predict_bart(mod1,skip, dat.pred)$samp_y
    
    # [L1 | V, A0, X0, X(-1), R1 = 1] 
    #dat.pred = makeX(Xname = theX, Vname = theV, time = t, modhere = mod[[ct]], simX = Xmat, policy_lag = A, dathere = gcd)
    # same RHS as in modeling L1
    postY = predict_bart(mod2,skip, dat.pred)$samp_y
    
    remove(fml1, fml2, mod1, mod2); gc()
    
    R = Xmat[[t]][,,1]
    # X(1) = X(0) for people who have no CD4 update at time==0 (days (0,12])
    mis_ind = which(R == 0) 
    postY[mis_ind] = X0Mat[mis_ind] 
    Xmat[[t]][,,2] = postY
    
    remove(postY)
    #if(A!=2){remove(X0Mat) ; gc()} 
    
    #[S1 | V, A0, R1, L1]
    
    load(paste0("Cat20-01-12_t",t,"_setting",setting,".RData")) 
    
    if(!(A %in% c(0,1))){
      policy_lag = policy_lag[rep(1:nsub,ndthin)]
    } else {
      policy_lag = rep(A, nrow(Vdist))
    }
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, c(1,2,4), t, fml, Xmat, Smat, policy_lag, gcd)
    
    Smat[[t]] = predict_bart(Bmod2, skip,dat.pred, Gcomp = TRUE)$samp_y
    
    remove(fml, res, baseyl,  Bmod2); gc()
    
    print(paste0("t=",t))
    
  } else {# t>1
    
    #[A(t-1) | A(t-2), L(t-2), S(t-1)]
    if(!(A %in% c(0,1))){
      
      A_prev = policy_lag
      I0 = Smat[[t-1]]
      if(t==2){
        C0 = X0Mat
        remove(X0Mat); gc()
      } else {
        C0 = Xmat[[t-2]][,,2]
      }
      dim(I0) = NULL
      dim(C0) = NULL
      policy_lag = 1*(A_prev== 1) + 1*(A_prev==0 & I0==2 & !is.na(C0) & (C0<= log(351)))    #rep(A, nrow(Vdist))
      
    }
    
    load( paste0("ContBin20-01-12_t",t,"_setting",setting,".RData") ) 
    
    #[Rt | V, A_{t-1}, R_{t-1}, L_{t-1}, S_{t-1}]
    
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, fml1, Xmat, Smat, policy_lag, gcd)
    Xmat[[t]][,,1] = predict_bart(mod1, skip,dat.pred, Gcomp = TRUE)$samp_y
    
    #[Lt | V, A_{t-1}, R_{t-1}, L_{t-1}, S_{t-1}]
    
    #dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, fml2, Xmat, Smat, policy_lag, gcd)
    postY = predict_bart(mod2, skip,dat.pred, Gcomp = TRUE)$samp_y
    
    remove(fml1, fml2, mod1, mod2); gc()
    
    R = Xmat[[t]][,,1]
    X_prev = Xmat[[t-1]][,,2]
    mis_ind = which(R == 0) 
    postY[mis_ind] = X_prev[mis_ind]
    Xmat[[t]][,,2] = postY
    
    remove(postY)
    
    #[St | V, A_{t-1}, S_{t-1}, Rt, Lt]
    
    load(paste0("Cat20-01-12_t",t,"_setting",setting,".RData")) 
    
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, fml, Xmat, Smat, policy_lag, gcd)
    
    Smat[[t]] = predict_bart(Bmod2[[1]], skip, dat.pred, Gcomp = TRUE)$samp_y
    
    tmp1 = predict_bart(Bmod2[[2]], skip, dat.pred, Gcomp = TRUE)$samp_y
    
    remove(fml, res, baseyl,  Bmod2); gc()
    
    eng_ind = which(Smat[[t-1]] == 2)
    Smat[[t]][eng_ind] = tmp1[eng_ind]
    remove(tmp1)
    
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
    
    print(paste0("t=",t))
    
    
  }#if t
  
  
  save(Xmat,Smat,file = paste0("Gcomp_A",A,"_t",t,"_setting",setting,"_20-01-12.RData"))
  
}#for t

save(Xmat,Smat,file = paste0("Gcomp_A",A,"_setting",setting,"_20-01-12.RData"))

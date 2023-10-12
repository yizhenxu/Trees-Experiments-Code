args <- commandArgs(trailingOnly=TRUE)
A = eval(parse(text=args[1]))

# theA_lag1 A = 0,1,2; A==2 means A = 1{CD4<=350}

source("BartGcomp_18-10-27_Athena0_prepare.R")

load("Mod18-10-27.RData")

##############################################################
### Outcome simulation 

nlatent = length(unique(gcd$outcome)) - 1
nsub = sum(gcd$time==0)
ytab = table(gcd$outcome)
ylevels = as.numeric(names(ytab)[ytab > 0])

Smat = vector("list", maxt) 
Xmat = vector("list", maxt)

for(i in 1:maxt){
  Xmat[[i]] = array(NA, dim = c(nsub, nd, length(theX)))
}

Vdist = cbind(1,gcd[time==0,theV,with=F])
Vdist = Vdist[rep(1:nsub,nd),]

#ptm <- proc.time()

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
    
    Xmat[[t]][,,1] = predict_bart(mod1[[t]],dat.pred)$samp_y
    
    #[X1 | V, A0, X0, L0] 
    #dat.pred = makeX(Xname = theX, Vname = theV, time = t, modhere = mod[[ct]], simX = Xmat, policy_lag = A, dathere = gcd)
    # same RHS as in modeling L1
    postY = predict_bart(mod2[[t]],dat.pred)$samp_y
    
    L = Xmat[[t]][,,1]
    mis_ind = which(L == 0) # people who have no CD4 update
    postY[mis_ind] = NA
    Xmat[[t]][,,2] = postY
    remove(postY)
    
    #[S1 | V, A0, X1, L1]
    
    if(!(A %in% c(0,1))){
      policy_lag = policy_lag[rep(1:nsub,nd)]
    } else {
      policy_lag = rep(A, nrow(Vdist))
    }
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, c(1,2,4), t, modfml[[3*t]], Xmat, Smat, policy_lag, gcd)
    
    Smat[[t]] = predict_bart(mod3_E[[t]],dat.pred, Gcomp = TRUE)$samp_y

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
    Xmat[[t]][,,1] = predict_bart(mod1[[t]],dat.pred, Gcomp = TRUE)$samp_y
    
    #[Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[3*(t-1)+2]], Xmat, Smat, policy_lag, gcd)
    postY = predict_bart(mod2[[t]],dat.pred, Gcomp = TRUE)$samp_y
    
    L = Xmat[[t]][,,1]
    X_prev = Xmat[[t-1]][,,2]
    mis_ind = which(L == 0) # people who have no CD4 update keeps the old value
    postY[mis_ind] = X_prev[mis_ind]
    Xmat[[t]][,,2] = postY
    
    remove(postY)
    
    #[St | V, A_{t-1}, S_{t-1}, Xt, Lt]
    
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[3*t]], Xmat, Smat, policy_lag, gcd)
    
    Smat[[t]] = predict_bart(mod3_D[[t-1]],dat.pred, Gcomp = TRUE)$samp_y
    
    tmp1 = predict_bart(mod3_E[[t]],dat.pred, Gcomp = TRUE)$samp_y
    
    eng_ind = which(Smat[[t-1]] == 2)
    Smat[[t]][eng_ind] = tmp1[eng_ind]
    remove(tmp1)
    
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
    
    print(paste0("t=",t))
    
    
  }
  
}

#proc.time() - ptm
#
save(Xmat,Smat,file = paste0("Gcomp_A",A,".RData"))

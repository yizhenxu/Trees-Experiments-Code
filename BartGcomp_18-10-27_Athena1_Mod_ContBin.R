sink(paste0("output_",screennum,".txt")) #write all output to txt -- source('script.R',echo=T) to write input too

source("BartGcomp_18-10-27_Athena0_prepare.R")

mod1 =  vector("list", maxt) # binary models
mod2 =  vector("list", maxt) # continuous models
mod3_D =  vector("list", maxt-1) # categorical models for those previously disengaged
mod3_E =  vector("list", maxt) # categorical models for those previously engaged

modfml =  vector("list", (length(theX)+1)*maxt) 

ptm <- proc.time()
for(t in 1:maxt){
  
  fitdat = trd[time==200*t,]
  
  if(t==1){
    
    # [L1 | V, A0, X0, L0]
    
    ## trd
    RHS = c(theV, theA_lag1, theX_cov_lag1)
    # exclude lag X_R from models at t=1 because that is the same as lag L at baseline
    RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
    modfml[[3*(t-1)+1]] = fml
    
    tmp = model_bart(as.formula(fml), data = fitdat, type = "binary",
                     Prior = Prior_binary,
                     Mcmc = Mcmc_binary)
    
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod1[[t]] = tmp[-del]
    remove(tmp)
    
    # [X1 | V, A0, X0, L0, L1 = 1] # simulated Xt remains NA until a first Lt == 1
    
    ## trd
    RHS = c(theV, theA_lag1, theX_cov_lag1)
    # exclude lag X_R from models at t=1 because that is the same as lag L at baseline
    RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
    RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
    fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
    modfml[[3*(t-1)+2]] = fml
    
    tmp = model_bart(as.formula(fml), data = fitdat[get(theX[1])==1,], type = "continuous",
                     Prior = Prior_cont,
                     Mcmc = Mcmc_cont)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod2[[t]] = tmp[-del]
    remove(tmp)
    
    print(t)
    
  } else {
    
    # [Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    
    ## trd
    RHS = c(theV, theA_lag1, theX_cov_lag1,theY_ind_lag1)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
    modfml[[3*(t-1)+1]] = fml
    
    tmp = model_bart(as.formula(fml), data = fitdat, type = "binary",
                     Prior = Prior_binary,
                     Mcmc = Mcmc_binary)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod1[[t]] = tmp[-del]
    remove(tmp)
    
    # [Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}, Lt = 1]
    
    ## trd
    RHS = c(theV, theA_lag1, theX_cov_lag1, theY_ind_lag1)
    RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
    fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
    modfml[[3*(t-1)+2]] = fml
    
    tmp = model_bart(as.formula(fml), data = fitdat[get(theX[1])==1,], type = "continuous",
                     Prior = Prior_cont,
                     Mcmc = Mcmc_cont)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod2[[t]] = tmp[-del]
    remove(tmp)
    
    print(t)
  }
}
proc.time() - ptm

save(mod1,mod2,modfml,file = "ContBin.RData")


##################################################################################
##################################################################################

sink()
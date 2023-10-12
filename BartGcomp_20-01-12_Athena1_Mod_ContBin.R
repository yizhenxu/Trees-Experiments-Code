#!/bin/bash
#SBATCH -J CB
#SBATCH --time 72:00:00
#SBATCH --array=1-10
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.err
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

#Rscript BartGcomp_20-01-12_Athena1_Mod_ContBin.R

rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer

t = as.numeric(slurm_arrayid)#1:10
setting = 4 

source("BartGcomp_20-01-12_Athena0_prepare.R")

fitdat = trd[time==200*t,]
  
if(t==1){
  
  #[R1 | V, A0, X0, X(-1)]
  
  ## trd
  RHS = c(theV, theA_lag1, theX_cov_1_lag1) ### USING X[-1]!!!!!!!!
  RHS = fmlfilter(fitdat,RHS)
  fml1 = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
  
  tmp = model_bart(as.formula(fml1), data = fitdat, type = "binary",
                   Prior = Prior_binary,
                   Mcmc = Mcmc_binary)
  
  del = which(names(tmp) %in% c("treefit","samp_y"))
  mod1 = tmp[-del]
  remove(tmp)
  
  # [L1 | V, A0, X0, X(-1), R1 = 1] # simulated Xt remains NA until a first Rt == 1
  
  ## trd
  RHS = c(theV, theA_lag1, theX_cov_1_lag1)
  RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
  fml2 = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
  
  tmp = model_bart(as.formula(fml2), data = fitdat[get(theX[1])==1,], type = "continuous",
                   Prior = Prior_cont,
                   Mcmc = Mcmc_cont)
  
  del = which(names(tmp) %in% c("treefit","samp_y"))
  mod2 = tmp[-del]
  remove(tmp)
  
} else {
  
  # [Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
  
  ## trd
  RHS = c(theV, theA_lag1, theX_cov_lag1,theY_ind_lag1)
  RHS = fmlfilter(fitdat,RHS)
  fml1 = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
  
  tmp = model_bart(as.formula(fml1), data = fitdat, type = "binary",
                   Prior = Prior_binary,
                   Mcmc = Mcmc_binary)
  
  del = which(names(tmp) %in% c("treefit","samp_y"))
  mod1 = tmp[-del]
  remove(tmp)
  
  # [Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}, Lt = 1]
  
  ## trd
  RHS = c(theV, theA_lag1, theX_cov_lag1, theY_ind_lag1)
  RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
  fml2 = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
  
  tmp = model_bart(as.formula(fml2), data = fitdat[get(theX[1])==1,], type = "continuous",
                   Prior = Prior_cont,
                   Mcmc = Mcmc_cont)
  
  del = which(names(tmp) %in% c("treefit","samp_y"))
  mod2 = tmp[-del]
  remove(tmp)
  
}  
  
save(mod1,mod2,fml1,fml2,file = paste0("ContBin20-01-12_t",t,"_setting",setting,".RData"))


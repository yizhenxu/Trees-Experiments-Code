#!/bin/bash
#SBATCH -J ampath
#SBATCH --time 72:00:00
#SBATCH --array=7-9
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.err
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

#Rscript BartGcomp_19-12-13_AMPATH_AthenaSep.R



# SEPARATE KINDO=T AND KINDO=F INTO PARALLEL RUNS

rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
sim = as.numeric(slurm_arrayid)#1:9

library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")
datseed = 99 # seed for simulating data in Setting.R
source("BartGcomp_19-12-13_AMPATH_setting.R") #uniform prior for correlation

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
    acc = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
    
  } else {
    bmpfit2 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                         base = by,
                         Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                         Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                         correction = FALSE, Kindo = KD, do_alpha2_prior = da)
    
    save(bmpfit2, file = paste0("CK_AMPATH_b",nb,"d",nd,"s",sim,".RData"))
    acc = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)
  }
  
}

if(KD == T){
  bmpfit1 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                       base = by,
                       Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                       Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                       correction = FALSE, Kindo = KD)
  
  save(bmpfit1, file = paste0("CK_AMPATH_b",nb,"d",nd,"s",sim,".RData"))
  acc = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
  
}

print(acc)


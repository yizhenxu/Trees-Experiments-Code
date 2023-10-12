# SEPARATE KINDO=T AND KINDO=F INTO PARALLEL RUNS

rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
sim = as.numeric(slurm_arrayid)#1:20

if(0){ # do sims to do Alg 3
  ct = as.numeric(slurm_arrayid)#1:10
  s = seq(1,19,2); sim = s[ct]
}

library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")
datseed = 99 # seed for simulating data in Setting.R
source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation

set.seed(CatSeed)

# do_alpha2_prior = FALSE valid only when KD = F then both mympbart0 (Alg2) and mympbart1 (Alg1) reduce to Alg 3
# TRUE: update alpha2 in both Step 1 & 3; FALSE: update alpha2 only in Step 3

if(KD == F){
  bmpfit0 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                       base = by,
                       Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                       Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                       correction = FALSE, Kindo = KD, do_alpha2_prior = TRUE)
  
  save(bmpfit0, file = paste0("CK11-26_b",nb,"d",nd,"s",sim,".RData"))
  
  bmpfit2 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                       base = by,
                       Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                       Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                       correction = FALSE, Kindo = KD, do_alpha2_prior = FALSE)
  
  save(bmpfit2, file = paste0("CK11-26_b",nb,"d",nd,"s",sim,"mg",mg,".RData"))

}

if(KD == T){
  bmpfit1 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                       base = by,
                       Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                       Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                       correction = FALSE, Kindo = KD)

	save(bmpfit1, file = paste0("CK11-26_b",nb,"d",nd,"s",sim,".RData"))

}




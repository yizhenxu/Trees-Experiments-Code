rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer, 6 comparisons, 50 round of simulations, id = 0:(6*50-1) = 0:299
id = as.numeric(slurm_arrayid) # seed for simulating data in Setting.R
#0,...,49; 50,...,99;100,...,149;150,...,199;...,299

library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")

Ns = 50 # number of data simulations for each setting
# sim from 9 to 20 # reference level = 3
sim = floor(id/Ns)*2+9 #9 11 13 15 17 19
datseed = id - floor(id/Ns)*Ns #0:49; seed for simulating data in Setting.R
source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation

destfile = paste0("CK12-10-2019_sim",sim,"_ct",datseed,"_alg3.txt")

if(!file.exists(destfile)){
  bmpfit2 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                       base = by,
                       Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                       Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                       correction = FALSE, Kindo = FALSE, do_alpha2_prior = FALSE)
  c = bmpfit2$sigmasample; c = matrix(c, nrow = 4); cm = apply(c,1,mean); cs = apply(c,1,sd)
  acc = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)
  res = c(cm, cs, acc)
  write(res,file =  paste0("CK12-10-2019_sim",sim,"_ct",datseed,"_alg3.txt"), sep="\n")
  
}

if(!file.exists(destfile)){
  resm = c(sim, datseed)
  write(resm,file =  "check3.txt",append=T, sep="\n")
}




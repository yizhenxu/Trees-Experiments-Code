rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
sim = as.numeric(slurm_arrayid)


library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")

#sim = 1 setting1 (more balanced), 2 setting2 (extreme unbalanced)
print(paste0("Simulation Setting ",sim))

source(paste0("BartGcomp_19-11-22_Setting",sim,".R"))

ylev = 1:3
bmpfit0 = vector("list",length(ylev))
bmpfit1 = vector("list",length(ylev))

tab0 = matrix(NA, length(ylev), 4)
tab1 = matrix(NA, length(ylev), 4)
colnames(tab0) = colnames(tab1) = c("Train Mean", "Train Mode", "Test Mean", "Test Mode") 

for(j in 1:length(ylev)){
  set.seed(CatSeed)
  
  bmpfit0[[j]] = model_bart(as.formula(fml), data = trd, type = "multinomial",
                            base = ylev[j],
                            Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                            Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                            correction = FALSE, Kindo = FALSE)
  bmpfit1[[j]] = model_bart(as.formula(fml), data = trd, type = "multinomial",
                            base = ylev[j],
                            Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                            Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                            correction = FALSE, Kindo = TRUE)
  
  tab0[j,] = accuracyfun(bmpfit = bmpfit0[[j]],testd = ted,trdy = trd$y,tedy = ted$y)
  tab1[j,] = accuracyfun(bmpfit = bmpfit1[[j]],testd = ted,trdy = trd$y,tedy = ted$y)
  
}

save(bmpfit0,bmpfit1, file = paste0("CK11-23-2019_b",nb,"d",nd,"_s",sim,"_Ref.RData"))

# summary of Sigma matrix
sigtab = noquote(cbind(SigSumm(bmpfit1[[3]]),SigSumm(bmpfit0[[3]])))
sigtab

# aggregated accuracy
acctab = cbind(tab1, tab0)
acctab

save(sigtab,acctab, file = paste0("CK11-23-2019_b",nb,"d",nd,"s",sim,"_tab.RData"))



rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
#slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
#sim = as.numeric(slurm_arrayid)


library(GcompBART)
source("BartGcomp_19-04-17_CompareKindo_Functions.R")

sim = 6 # reference level demonstration on setting 2 (smallest lev & ref = 3) and setting 6 (smallest lev 1, ref 3) 
print(paste0("Simulation Setting ",sim))

source(paste0("BartGcomp_19-04-17_CompareKindo_Setting",sim,".R"))

ylev = 1:3
bmpfit0 = vector("list",length(ylev))
bmpfit1 = vector("list",length(ylev))
tab = vector("list",length(ylev))

load(paste0("CompareKindo_",sim,"_Ref.RData"))
load(paste0("CompareKindo_",sim,"_tab_Ref.RData"))

for(j in c(1)){
  set.seed(CatSeed)
  
  bmpfit0[[j]] = model_bart(as.formula(fml), data = trd, type = "multinomial",
                       base = ylev[j],
                       Prior = Prior_mult(p = p, ntree = nt),
                       Mcmc = Mcmc_mult(p = p),
                       correction = FALSE, Kindo = FALSE)
  bmpfit1[[j]] = model_bart(as.formula(fml), data = trd, type = "multinomial",
                       base = ylev[j],
                       Prior = Prior_mult(p = p, ntree = nt),
                       Mcmc = Mcmc_mult(p = p),
                       correction = FALSE, Kindo = TRUE)
  
  pdf(paste0("TDPlot",sim,"_",j,".pdf"),width = 7, height = 7)
  pm1 = p-1
  par(mfrow=c(2,pm1))
  TDPlot(bmpfit1[[j]])
  TDPlot(bmpfit0[[j]])
  dev.off()
  
  pdf(paste0("SigPlot",sim,"_",j,".pdf"),width = 7, height = 3.5)
  par(mfrow=c(1,2))
  a = bmpfit0[[j]]$sigmasample; plot(a)
  a = bmpfit1[[j]]$sigmasample; plot(a)
  dev.off()
  
  a = accuracyfun(bmpfit = bmpfit0[[j]],testd = ted,trdy = trd$y,tedy = ted$y)
  b = accuracyfun(bmpfit = bmpfit1[[j]],testd = ted,trdy = trd$y,tedy = ted$y)
  
  tab[[j]] = rbind(b,a); colnames(tab[[j]]) = names(a); rownames(tab[[j]]) = c("Kindo","NewAlg") 

}
save(bmpfit0,bmpfit1, file = paste0("CompareKindo_",sim,"_Ref.RData"))

tab

save(tab, file = paste0("CompareKindo_",sim,"_tab_Ref.RData"))

rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
sim = as.numeric(slurm_arrayid)

#devtools::install_github("yizhenxu/GcompBART")
library(GcompBART)
source("BartGcomp_19-04-17_CompareKindo_Functions.R")

estModel = 1
accucalc = 1

#sim = 1 setting1, 2 setting2, 3 AMPATH grid data t=1, 4 UCI Contraceptive Method Choice, 5 UCI Diabetes Readmission
print(paste0("Simulation Setting ",sim))

source(paste0("BartGcomp_19-04-17_CompareKindo_Setting",sim,".R"))


by = ymode(trd$y)
#by = 3

print(paste0("Simulation Setting ",sim))

source(paste0("BartGcomp_19-04-17_CompareKindo_Setting",sim,".R"))

if(estModel){
  set.seed(CatSeed)
  bmpfit0 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                       base = by,
                       Prior = Prior_mult(p = p, ntree = nt),
                       Mcmc = Mcmc_mult(p = p),
                       correction = FALSE, Kindo = FALSE)
  bmpfit1 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                       base = by,
                       Prior = Prior_mult(p = p, ntree = nt),
                       Mcmc = Mcmc_mult(p = p),
                       correction = FALSE, Kindo = TRUE)
  save(bmpfit0,bmpfit1, file = paste0("CompareKindo_",sim,"_by",by,".RData"))
} else {
  load(paste0("CompareKindo_",sim,".RData"))
}

pdf(paste0("TDPlot",sim,"_by",by,".pdf"),width = 7, height = 7)
pm1 = p-1
par(mfrow=c(2,pm1))
TDPlot(bmpfit1)
TDPlot(bmpfit0)
dev.off()

pdf(paste0("SigPlot",sim,"_by",by,".pdf"),width = 7, height = 3.5)
par(mfrow=c(1,2))
a = bmpfit0$sigmasample; plot(a)
a = bmpfit1$sigmasample; plot(a)
dev.off()

noquote(cbind(SigSumm(bmpfit1),SigSumm(bmpfit0)))

if(accucalc){
  a = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
  b = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
  
  tab = rbind(b,a); colnames(tab) = names(a); rownames(tab) = c("Kindo","NewAlg") 
  
  tab
  
  table(trd$y)
  table(ted$y)
  
  save(tab, file = paste0("CompareKindo_",sim,"_by",by,"_tab.RData"))
}
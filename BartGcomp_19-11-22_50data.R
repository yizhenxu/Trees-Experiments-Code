rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer, ct ranges 1:50
ct = as.numeric(slurm_arrayid)

setting = 1

by = 3 # base y (reference level) is 3 from data generation

source("BartGcomp_19-11-22_50data_settings.R") 
# ct used as seed for data generation
# setting used in this sourced setting file too

library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")


set.seed(CatSeed)

bmpfit0 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                          base = by,
                          Prior = Prior_mult(p = p, nu_C = nuC, ntree = nt),
                          Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                          correction = FALSE, Kindo = FALSE)
bmpfit1 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                          base = by,
                          Prior = Prior_mult(p = p, nu_C = nuC, ntree = nt),
                          Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                          correction = FALSE, Kindo = TRUE)


a = bmpfit0$sigmasample; a = matrix(a, nrow = 4); am = apply(a,1,mean); as = apply(a,1,sd)
b = bmpfit1$sigmasample; b = matrix(b, nrow = 4); bm = apply(b,1,mean); bs = apply(b,1,sd)
resm = c(ct, bm, am)
ress = c(ct, bs, as)
write(resm,file =  paste0("CK11-23-2019_b",nb,"d",nd,"s",setting,"_SigM.txt"),append=T, sep="\n")
write(ress,file =  paste0("CK11-23-2019_b",nb,"d",nd,"s",setting,"_SigS.txt"),append=T, sep="\n")

if(0){
  pathfrom = "C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/CompareKindo/nb2000nd3000/"
  all = read.table(paste0(pathfrom,"CK11-23-2019_b",nb,"d",nd,"s",setting,"_Sig.txt"),sep="\n")
  all = read.table(paste0(pathfrom,"CK11-23-2019_b",nb,"d",nd,"s",setting,"_acc.txt"),sep="\n")
  all = matrix(all$V1,nrow = 9)
  allm = apply(all[2:9,],1,mean)
  #nd = ncol(all)
  sdm = apply(all[2:9,],1,sd)
  #nct = 50
  #CI = apply(all[2:9,],1,function(x) sort(x)[c(ceiling(nct*0.05), floor(nct*0.95))])
  # 9 rows: ct, sig 11 12 21 22 for alg1, sig 11 12 21 22 for alg2 
}

acc = rep(NA,9) # "Train Mean", "Train Mode", "Test Mean", "Test Mode" for alg 1, alg 2
acc[1] = ct
acc[2:5] = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
acc[6:9] = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
write(acc,file =  paste0("CK11-23-2019_b",nb,"d",nd,"s",setting,"_acc.txt"),append=T, sep="\n")

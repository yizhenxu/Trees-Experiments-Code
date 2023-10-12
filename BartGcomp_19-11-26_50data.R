#!/bin/bash
#SBATCH -J 50data1
#SBATCH --time 72:00:00
#SBATCH --array=[0-299]%10
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

#Rscript BartGcomp_19-11-26_50data_1.R



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

#set.seed(CatSeed)
acc = rep(NA,13) # "Train Mean", "Train Mode", "Test Mean", "Test Mode" for alg 1, alg 2, alg3
acc[1] = datseed

bmpfit0 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                     base = by,
                     Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                     Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                     correction = FALSE, Kindo = FALSE)
a = bmpfit0$sigmasample; a = matrix(a, nrow = 4); am = apply(a,1,mean); as = apply(a,1,sd)
acc[6:9] = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
remove(a); remove(bmpfit0)

bmpfit1 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                     base = by,
                     Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                     Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                     correction = FALSE, Kindo = TRUE)
b = bmpfit1$sigmasample; b = matrix(b, nrow = 4); bm = apply(b,1,mean); bs = apply(b,1,sd)
acc[2:5] = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
remove(b); remove(bmpfit1)

bmpfit2 = model_bart(as.formula(fml), data = trd, type = "multinomial",
                     base = by,
                     Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                     Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                     correction = FALSE, Kindo = FALSE, do_alpha2_prior = FALSE)
c = bmpfit2$sigmasample; c = matrix(c, nrow = 4); cm = apply(c,1,mean); cs = apply(c,1,sd)
acc[10:13] = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)
remove(c); remove(bmpfit2)

resm = c(datseed, bm, am, cm)
ress = c(datseed, bs, as, cs)

res = cbind(resm, ress, acc)

save(res, file = paste0("CK12-10-2019_sim",sim,"_ct",datseed,".RData"))
#write(resm,file =  paste0("CK12-10-2019_sim",sim,"_SigM.txt"),append=T, sep="\n")
#write(ress,file =  paste0("CK12-10-2019_sim",sim,"_SigS.txt"),append=T, sep="\n")
#write(acc,file =  paste0("CK12-10-2019_sim",sim,"_acc.txt"),append=T, sep="\n")

if(0){
  pathfrom = "C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/CompareKindo/out/"
  pathfrom = "G:/My Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/CompareKindo/out/"
  
  # Setting 1: sim 9 13 17
  # Setting 2: sim 11 15 19
  
  L1 = L2 = vector("list",3)
  names(L1) = names(L2) = c("SigM", "SigS", "ACC")
  for(i in 1:3){
    if(i<3){
      L1[[i]] = matrix(NA, 3, 6); rownames(L1[[i]]) = c("unif","neg","pos"); colnames(L1[[i]]) = rep(paste0("A",1:3),2)
      L2[[i]] = matrix(NA, 3, 6); rownames(L2[[i]]) = c("unif","neg","pos"); colnames(L2[[i]]) = rep(paste0("A",1:3),2)
    } else {
      L1[[i]] = matrix(NA, 6, 6); rownames(L1[[i]]) = rep(c("unif","neg","pos"),2); colnames(L1[[i]]) = rep(paste0("A",1:3),2)
      L2[[i]] = matrix(NA, 6, 6); rownames(L2[[i]]) = rep(c("unif","neg","pos"),2); colnames(L2[[i]]) = rep(paste0("A",1:3),2)
    }
  }
  
  # Setting 1

  summtab = function(s,L){
    for(i in 1:3){
      sim = s[i]
      for(alg in 1:3){
        a = c()
        for(datseed in 0:49){
          b = read.table(paste0("CK12-10-2019_sim",sim,"_ct",datseed,"_alg",alg,".txt"),sep="\n")
          a = cbind(a, b$V1[c(1,2,5,6,9:12)])
        }
        am = apply(a, 1, mean); as = apply(a, 1, sd)
        am = format(round(am, 3),nsmall=3) ; as = format(round(as, 3),nsmall=3)
        L[[1]][i,alg] = paste0(am[1]," (",as[1],")")
        L[[1]][i,alg+3] = paste0(am[2]," (",as[2],")")
        L[[2]][i,alg] = paste0(am[3]," (",as[3],")")
        L[[2]][i,alg+3] = paste0(am[4]," (",as[4],")")
        L[[3]][i,alg] = paste0(am[5]," (",as[5],")")
        L[[3]][i,alg+3] = paste0(am[7]," (",as[7],")")
        L[[3]][i+3,alg] = paste0(am[6]," (",as[6],")")
        L[[3]][i+3,alg+3] = paste0(am[8]," (",as[8],")")
        
      }#alg
    }#sim i
    return(L)
  }#function
  L1 =summtab(c(9,13,17),L1)
  L2 =summtab(c(11,15,19),L2)
  
  noquote(L1)
  noquote(L2)
  
  for(sim in seq(9,19,2)){
    for(datseed in 0:49){
      read.table(paste0("CK12-10-2019_sim",sim,"_ct",datseed,"_alg2.txt"),sep="\n")
    }
    
    M = read.table(paste0(pathfrom,"CK12-10-2019_sim",sim,"_SigM.txt"),sep="\n")
    S = read.table(paste0(pathfrom,"CK12-10-2019_sim",sim,"_SigS.txt"),sep="\n")
    A = read.table(paste0(pathfrom,"CK12-10-2019_sim",sim,"_acc.txt"),sep="\n")
    M = matrix(M$V1,nrow = 13)
    S = matrix(S$V1,nrow = 13)
    A = matrix(A$V1,nrow = 13)
    mM = apply(M[2:9,],1,mean);sdM = apply(M[2:9,],1,sd)
    resM = matrix(paste0(round(mM,2), " (",round(sdM,2),")"), nrow = 2 )
    mS = apply(S[2:9,],1,mean);sdS = apply(S[2:9,],1,sd)
    resS = matrix(paste0(round(mS,2), " (",round(sdS,2),")"), nrow = 2 )
    mA = apply(A[2:9,],1,mean);sdA = apply(A[2:9,],1,sd)
    resA = matrix(paste0(round(mA,2), " (",round(sdA,2),")"), nrow = 2 )
    print(sim)
    print(noquote(resM));print(noquote(resS));print(noquote(resA))
  }
  
  #ndatseed = 50
  #CI = apply(all[2:9,],1,function(x) sort(x)[c(ceiling(nct*0.05), floor(nct*0.95))])
  # 9 rows: ct, sig 11 12 21 22 for alg1, sig 11 12 21 22 for alg2 
}
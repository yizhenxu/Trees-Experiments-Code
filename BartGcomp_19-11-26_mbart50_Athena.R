#!/bin/bash
#SBATCH -J acc50
#SBATCH --time 96:00:00
#SBATCH --array=1-100
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.err
#SBATCH --mail-type=all
#SBATCH --mail-user=xuyizhen00@gmail.com

#Rscript BartGcomp_19-11-26_mbart50_Athena.R

# data simulated and accuracy summarized under settings 1 and 2

## run the current R script

# sbatch mbash2.sh

rm(list=ls()) 

gc()
#memory.limit()
# set max memory usage is 10G
#memory.size(max=10000)

# coerce the value to a numeric
simidx = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # 1:100

book = matrix(NA, ncol = 2, nrow = 100)
colnames(book) = c("sim", "datseed")
book[,1] = rep(c(1,3) ,each = 50); book[,2] = rep(0:49, 2);

sim = book[simidx, 1]

datseed = book[simidx, 2] # seed for simulating data in Setting.R
source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation

sfn = paste0("./res/CK11-26_b",nb,"d",nd,"s",sim,"_ct",datseed,"_mbart.RData")

if(!file.exists(sfn)){
  library(BART)
  library(GcompBART)
  source("BartGcomp_19-11-22_Functions.R")
  
  set.seed(CatSeed)
  
  post <- mbart(x.train = trd[,1:6], y.train = trd[,7], x.test = newdat[,1:6], ndpost = nd, nskip = nb, keepevery = 1, ntree=nt)
  
  tmp = accuracyfun_mbart(post, newdat$y, ntr, n)
  
  save(tmp, file = sfn)
}


### calculate ratio of each category averaged across replicates
{
  
  summ = function(x){
    x = x[!is.na(x)]
    sx = sort(x)
    n = length(x)
    r = paste0(round(mean(x),2), " (",round(sx[ceiling(n*0.025)],2),", ",round(sx[floor(n*0.975)],2),")")
    return(r)
  }
  
  
  book = matrix(NA, ncol = 2, nrow = 100)
  colnames(book) = c("sim", "datseed")
  book[,1] = rep(c(1,3) ,each = 50); book[,2] = rep(0:49, 2);
  
  trdtab = c()
  tedtab = c()
  
  for(simidx in 1:100){
    sim = book[simidx, 1]
    
    datseed = book[simidx, 2] # seed for simulating data in Setting.R
    source("/Volumes/GoogleDrive/My Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/BartGcomp_19-11-26_Setting.R") #uniform prior for correlation
    
    trdtab = rbind(trdtab, table(trd$y))
    
    tedtab = rbind(tedtab, table(ted$y))
    
  }
  
  trdtab1 = t(apply(trdtab, 1, function(x) x/sum(x)))
  noquote(apply(trdtab1[1:50,],2, summ)) #0.45 (0.44, 0.46) 0.25 (0.24, 0.26)  0.3 (0.29, 0.31) 
  noquote(apply(trdtab1[51:100,],2, summ))#0.32 (0.3, 0.33) 0.65 (0.64, 0.66) 0.03 (0.03, 0.04) 
  round(apply(trdtab1[1:50,],2,mean),2)# 0.45 0.25 0.30 
  round(apply(trdtab1[51:100,],2,mean),2)# 0.32 0.65 0.03 
  
  
  tedtab1 = t(apply(tedtab, 1, function(x) x/sum(x)))
  noquote(apply(trdtab1[1:50,],2, summ)) #0.45 (0.44, 0.46) 0.25 (0.24, 0.26)  0.3 (0.29, 0.31) 
  noquote(apply(trdtab1[51:100,],2, summ))#0.32 (0.3, 0.33) 0.65 (0.64, 0.66) 0.03 (0.03, 0.04) 
  round(apply(tedtab1[1:50,],2,mean),2)# 0.45 0.25 0.30  
  round(apply(tedtab1[51:100,],2,mean),2)# 0.32 0.65 0.03 
  
  
}


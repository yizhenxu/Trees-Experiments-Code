#!/bin/bash
#SBATCH -J acc
#SBATCH --time 72:00:00
#SBATCH --array=1-2
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.err
#SBATCH --mail-type=all
#SBATCH --mail-user=xuyizhen00@gmail.com

#Rscript BartGcomp_19-11-26_mbart_Athena.R

# data simulated and accuracy summarized under settings 1 and 2

## run the current R script

# sbatch mbash1.sh

rm(list=ls()) 

gc()
#memory.limit()
# set max memory usage is 10G
#memory.size(max=10000)

# coerce the value to a numeric
simidx = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # 1,2

### Table 1 Accuracy for setting 1 & uniform prior
###     Alg2/Alg3   KINDO
### Ref 1 sim = 1, sim = 2  
### Ref 2 sim = 5, sim = 6
### Ref 3 sim = 9, sim = 10

### Table 2 Accuracy for setting 2 & uniform prior
### Ref 1 sim = 3, sim = 4  
### Ref 2 sim = 7, sim = 8
### Ref 3 sim =11, sim = 12

v = c(1,3)
sim = v[simidx]

library(BART)
library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")
datseed = 99 # seed for simulating data in Setting.R
source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation

set.seed(CatSeed)

post <- mbart(x.train = trd[,1:6], y.train = trd[,7], x.test = newdat[,1:6], ndpost = nd, nskip = nb, keepevery = 1, ntree=nt)

tmp = accuracyfun_mbart(post, newdat$y, ntr, n)

save( tmp, file = paste0("./res/CK11-26_b",nb,"d",nd,"s",sim,"_mbart.RData"))


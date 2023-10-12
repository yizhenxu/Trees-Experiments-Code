#!/bin/bash
#SBATCH -J Acc
#SBATCH --time 72:00:00
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.err
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

#Rscript BartGcomp_19-12-13_AMPATH_Acc.R

rm(list=ls()) 


library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")

sim = 1 # by = 1
datseed = 99 # seed for simulating data in Setting.R
source("BartGcomp_19-12-13_AMPATH_setting.R") #uniform prior for correlation
# load bmpfit0,bmpfit1,bmpfit2
load(paste0("CK_AMPATH_b",nb,"d",nd,"s",sim+0,".RData"))
load(paste0("CK_AMPATH_b",nb,"d",nd,"s",sim+1,".RData"))
load(paste0("CK_AMPATH_b",nb,"d",nd,"s",sim+2,".RData"))

Accby1 =  vector("list", 3)
Accby1[[1]] = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
Accby1[[2]] = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
Accby1[[3]] = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)
remove(bmpfit0, bmpfit1, bmpfit2)

sim = 4 # by = 2
datseed = 99 # seed for simulating data in Setting.R
source("BartGcomp_19-12-13_AMPATH_setting.R") #uniform prior for correlation
# load bmpfit0,bmpfit1,bmpfit2
load(paste0("CK_AMPATH_b",nb,"d",nd,"s",sim+0,".RData"))
load(paste0("CK_AMPATH_b",nb,"d",nd,"s",sim+1,".RData"))
load(paste0("CK_AMPATH_b",nb,"d",nd,"s",sim+2,".RData"))

Accby2 =  vector("list", 3)
Accby2[[1]] = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
Accby2[[2]] = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
Accby2[[3]] = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)

save(Accby1, Accby2, file = "AccAMPATH.RData")

#########################################################################
load("AccAMPATH.RData")

sim = 7 # by = 4
datseed = 99 # seed for simulating data in Setting.R
source("BartGcomp_19-12-13_AMPATH_setting.R") #uniform prior for correlation
# load bmpfit0,bmpfit1,bmpfit2
load(paste0("CK_AMPATH_b",nb,"d",nd,"s",sim+0,".RData"))
load(paste0("CK_AMPATH_b",nb,"d",nd,"s",sim+1,".RData"))
load(paste0("CK_AMPATH_b",nb,"d",nd,"s",sim+2,".RData"))

Accby3 =  vector("list", 3)
Accby3[[1]] = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
Accby3[[2]] = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
Accby3[[3]] = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)

save(Accby1, Accby2, Accby3, file = "AccAMPATH.RData")
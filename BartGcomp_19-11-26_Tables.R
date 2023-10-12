rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
cp = as.numeric(slurm_arrayid) #1:6

library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")
datseed = 99 # seed for simulating data in Setting.R

if(cp == 1){
  ### Table 1 Accuracy for setting 1 & uniform prior
  ###     Alg2/Alg3   KINDO
  ### Ref 1 sim = 1, sim = 2  
  ### Ref 2 sim = 5, sim = 6
  ### Ref 3 sim = 9, sim = 10
  
  loc = matrix(c(1,5,9,2,6,10),nrow=3)
  
  Tab1 = matrix(NA, 6, 6); rownames(Tab1) = colnames(Tab1) = rep(c("1","2","3"),2); 
  # accuracyfun output : "Train Mean", "Train Mode", "Test Mean", "Test Mode"
  # Tab1 output: first 3 rows - Mean, last 3 rows - Mode; first 3 cols - Train, last 3 cols - Test
  
  sim = loc[1,1] # not using KD nuC Vm from variable sim in plots.R
  source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation
  
  for(ref in 1:3){
    sim0 = loc[ref,1]; sim1 = loc[ref,2]
    # load bmpfit0, bmpfit1, bmpfit2
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim1,".RData"))
    alg1 = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab1[ref, 1] = alg1[1]; Tab1[ref, 4] = alg1[3]; Tab1[ref+3, 1] = alg1[2]; Tab1[ref+3, 4] = alg1[4]
    remove(alg1);remove(bmpfit1)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,".RData"))
    alg2 = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab1[ref, 2] = alg2[1]; Tab1[ref, 5] = alg2[3]; Tab1[ref+3, 2] = alg2[2]; Tab1[ref+3, 5] = alg2[4]
    remove(alg2);remove(bmpfit0)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,"mg",FALSE,".RData"))
    alg3 = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab1[ref, 3] = alg3[1]; Tab1[ref, 6] = alg3[3]; Tab1[ref+3, 3] = alg3[2]; Tab1[ref+3, 6] = alg3[4]
    remove(alg3);remove(bmpfit2)
    
  }
  
  remove(trd);remove(ted)
  
  cat("\n Table 1 \n\n")
  
  print(round(Tab1,2))
  
  save(Tab1, file = paste0("CK12-11_Table",cp,".RData"))
  
}

if(cp == 2){
  ### Table 2 Accuracy for setting 2 & uniform prior
  ### Ref 1 sim = 3, sim = 4  
  ### Ref 2 sim = 7, sim = 8
  ### Ref 3 sim =11, sim = 12
  
  loc = matrix(c(3,7,11,4,8,12),nrow=3)
  Tab2 = matrix(NA, 6, 6); rownames(Tab2) = colnames(Tab2) = rep(c("1","2","3"),2); 
  
  sim = loc[1,1] # not using KD nuC Vm from variable sim in plots.R
  source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation
  
  for(ref in 1:3){
    sim0 = loc[ref,1]; sim1 = loc[ref,2]
    # load bmpfit0, bmpfit1, bmpfit2
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim1,".RData"))
    alg1 = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab2[ref, 1] = alg1[1]; Tab2[ref, 4] = alg1[3]; Tab2[ref+3, 1] = alg1[2]; Tab2[ref+3, 4] = alg1[4]
    remove(alg1);remove(bmpfit1)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,".RData"))
    alg2 = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab2[ref, 2] = alg2[1]; Tab2[ref, 5] = alg2[3]; Tab2[ref+3, 2] = alg2[2]; Tab2[ref+3, 5] = alg2[4]
    remove(alg2);remove(bmpfit0)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,"mg",FALSE,".RData"))
    alg3 = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab2[ref, 3] = alg3[1]; Tab2[ref, 6] = alg3[3]; Tab2[ref+3, 3] = alg3[2]; Tab2[ref+3, 6] = alg3[4]
    remove(alg3);remove(bmpfit2)
    
  }
  
  remove(trd);remove(ted)
  
  cat("\n Table 2 \n\n")
  
  print(round(Tab2,2))
  
  save(Tab2, file = paste0("CK12-11_Table",cp,".RData"))
}

if(cp == 3){
  ### Table 3 Setting 1 Ref 3
  ### Prior uniform sim = 9 , sim = 10
  ### Prior negativ sim = 13, sim = 14
  ### Prior positiv sim = 17, sim = 18
  
  loc = matrix(c(9,13,17,10,14,18),nrow=3)
  Tab3 = matrix(NA, 6, 6)
  rownames(Tab3) = rep(c("unif","neg","pos"),2); colnames(Tab3) = rep(c("1","2","3"),2); 
  
  sim = loc[1,1] # not using KD nuC Vm from variable sim in plots.R
  source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation
  
  for(ref in 1:3){
    sim0 = loc[ref,1]; sim1 = loc[ref,2]
    # load bmpfit0, bmpfit1, bmpfit2
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim1,".RData"))
    alg1 = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab3[ref, 1] = alg1[1]; Tab3[ref, 4] = alg1[3]; Tab3[ref+3, 1] = alg1[2]; Tab3[ref+3, 4] = alg1[4]
    remove(alg1);remove(bmpfit1)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,".RData"))
    alg2 = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab3[ref, 2] = alg2[1]; Tab3[ref, 5] = alg2[3]; Tab3[ref+3, 2] = alg2[2]; Tab3[ref+3, 5] = alg2[4]
    remove(alg2);remove(bmpfit0)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,"mg",FALSE,".RData"))
    alg3 = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab3[ref, 3] = alg3[1]; Tab3[ref, 6] = alg3[3]; Tab3[ref+3, 3] = alg3[2]; Tab3[ref+3, 6] = alg3[4]
    remove(alg3);remove(bmpfit2)
    
  }
  remove(trd);remove(ted)
  
  cat("\n Table 3 \n\n")
  
  print(round(Tab3,2))
  
  save(Tab3, file = paste0("CK12-11_Table",cp,".RData"))
}

if(cp == 4){
  
  ### Table 4 Setting 1 Ref 3
  ### Prior uniform sim = 9 , sim = 10
  ### Prior negativ sim = 13, sim = 14
  ### Prior positiv sim = 17, sim = 18
  
  loc = matrix(c(9,13,17,10,14,18),nrow=3)
  
  Tab4 = matrix(NA, 3, 6); rownames(Tab4) = c("unif","neg","pos"); colnames(Tab4) = rep(c("sig11", "sig12"), 3)
  
  sim = loc[1,1] # not using KD nuC Vm from variable sim in plots.R
  source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation
  
  for(ref in 1:3){
    sim0 = loc[ref,1]; sim1 = loc[ref,2]
    # load bmpfit0, bmpfit1, bmpfit2
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim1,".RData"))
    tmp = SigSumm(bmpfit1); Tab4[ref,1:2] = c(tmp[1,1],tmp[1,2])
    remove(bmpfit1)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,".RData"))
    tmp = SigSumm(bmpfit0); Tab4[ref,3:4] = c(tmp[1,1],tmp[1,2])
    remove(bmpfit0)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,"mg",FALSE,".RData"))
    tmp = SigSumm(bmpfit2); Tab4[ref,5:6] = c(tmp[1,1],tmp[1,2])
    remove(bmpfit2)
    
    Tab4 = noquote(Tab4)
  }
  remove(trd);remove(ted)
  
  cat("\n Table 4 \n\n")
  
  print(Tab4)
  
  save(Tab4, file = paste0("CK12-11_Table",cp,".RData"))
}

if(cp == 5){
  ### Table 5 Setting 2 Ref 3
  ### Prior uniform sim = 11, sim = 12
  ### Prior negativ sim = 15, sim = 16
  ### Prior positiv sim = 19, sim = 20
  
  loc = matrix(c(11,15,19,12,16,20),nrow=3)
  
  Tab5 = matrix(NA, 6, 6)
  rownames(Tab5) = rep(c("unif","neg","pos"),2); colnames(Tab5) = rep(c("1","2","3"),2); 
  
  sim = loc[1,1] # not using KD nuC Vm from variable sim in plots.R
  source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation
  
  for(ref in 1:3){
    sim0 = loc[ref,1]; sim1 = loc[ref,2]
    # load bmpfit0, bmpfit1, bmpfit2
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim1,".RData"))
    alg1 = accuracyfun(bmpfit = bmpfit1,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab5[ref, 1] = alg1[1]; Tab5[ref, 4] = alg1[3]; Tab5[ref+3, 1] = alg1[2]; Tab5[ref+3, 4] = alg1[4]
    remove(alg1);remove(bmpfit1)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,".RData"))
    alg2 = accuracyfun(bmpfit = bmpfit0,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab5[ref, 2] = alg2[1]; Tab5[ref, 5] = alg2[3]; Tab5[ref+3, 2] = alg2[2]; Tab5[ref+3, 5] = alg2[4]
    remove(alg2);remove(bmpfit0)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,"mg",FALSE,".RData"))
    alg3 = accuracyfun(bmpfit = bmpfit2,testd = ted,trdy = trd$y,tedy = ted$y)
    Tab5[ref, 3] = alg3[1]; Tab5[ref, 6] = alg3[3]; Tab5[ref+3, 3] = alg3[2]; Tab5[ref+3, 6] = alg3[4]
    remove(alg3);remove(bmpfit2)
    
  }
  remove(trd);remove(ted)
  
  cat("\n Table 5 \n\n")
  
  print(round(Tab5,2),)
  
  save(Tab5, file = paste0("CK12-11_Table",cp,".RData"))
}

if(cp == 6){
  
  ### Table 6 Setting 2 Ref 3
  ### Prior uniform sim = 11, sim = 12
  ### Prior negativ sim = 15, sim = 16
  ### Prior positiv sim = 19, sim = 20
  
  loc = matrix(c(11,15,19,12,16,20),nrow=3)
  
  Tab6 = matrix(NA, 3, 6); rownames(Tab6) = c("unif","neg","pos"); colnames(Tab6) = rep(c("sig11", "sig12"), 3)
  
  sim = loc[1,1] # not using KD nuC Vm from variable sim in plots.R
  source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation
  
  for(ref in 1:3){
    sim0 = loc[ref,1]; sim1 = loc[ref,2]
    # load bmpfit0, bmpfit1, bmpfit2
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim1,".RData"))
    tmp = SigSumm(bmpfit1); Tab6[ref,1:2] = c(tmp[1,1],tmp[1,2])
    remove(bmpfit1)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,".RData"))
    tmp = SigSumm(bmpfit0); Tab6[ref,3:4] = c(tmp[1,1],tmp[1,2])
    remove(bmpfit0)
    
    load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,"mg",FALSE,".RData"))
    tmp = SigSumm(bmpfit2); Tab6[ref,5:6] = c(tmp[1,1],tmp[1,2])
    remove(bmpfit2)
    
    Tab6 = noquote(Tab6)
  }
  remove(trd);remove(ted)
  
  cat("\n Table 6 \n\n")
  
  print(Tab6)
  
  save(Tab6, file = paste0("CK12-11_Table",cp,".RData"))
}

if(cp == 7){# join all tables
  for(k in 1:6){
    load( paste0("CK12-11_Table",k,".RData"))
  }
  save(Tab1,Tab2,Tab3,Tab4,Tab5,Tab6, file = "CK12-11_Table.RData")
  noquote(rep(paste0("Alg",1:3),2))
  print(round(Tab1,2))
  print(round(Tab2,2))
  print(round(Tab3,2))
  print(Tab4)
  print(round(Tab5,2))
  print(Tab6)

}
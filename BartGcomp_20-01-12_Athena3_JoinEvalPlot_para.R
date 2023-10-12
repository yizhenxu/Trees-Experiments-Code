#!/bin/bash
#SBATCH -J ACC
#SBATCH --time 72:00:00
#SBATCH --array=1-10
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.err
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

#Rscript BartGcomp_20-01-12_Athena3_JoinEvalPlot_para.R

rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer

t = as.numeric(slurm_arrayid)#1:10
setting = 4 ; maxt = 10

##############################################################

### Accuracy - Categorical Models 

### Categorical

source("BartGcomp_20-01-12_Athena0_prepare.R")

td = ted[time==200*t,]

ty =  unlist(ted[time==200*t, theY, with = F])

acc = matrix(NA, ncol=4, nrow=4)
mmacc = rep(NA, 2)

load(paste0("Cat20-01-12_t",t,"_setting",setting,".RData")) 

if(t==1){
  
  # [S1 | V, A0, X1, L1]
  PB = predict_bart(obj = Bmod2,  newdata = td)
  PM = PB$samp_y
  acc[1,] = taby(ty)
  acc[2:4, ] = matrix(unlist(lapply(1:4, function(s) summ(PM,s,1))), ncol=4) 
  
} else {
  
  for(prevS in c(0,1)){ #### S_{t-1} = 0 - disengaged, 1 - engaged
    
    ind = prevS+1
    
    ## ted
    if(prevS == 0){
      
      PB0 = predict_bart(obj = Bmod2[[1]],  newdata = td)
      PM0 = PB0$samp_y
      
    } else { #prevS == 1
      
      acc[1,] = taby(ty)
      
      #Condition 1: observed prop of xfer and death within 95% posterior CI
      PB1 = predict_bart(obj = Bmod2[[2]],  newdata = td)
      PM1 = PB1$samp_y
      
      tind = which(td[,theY_ind_lag1,with = F]==1)
      PM = PM0; PM[tind,] = PM1[tind,]
      
      # the first dim: y, mean post, 5% post, 95% post
      # the second dim:  diseng, eng, xfer, death
      acc[2:4, ] = matrix(unlist(lapply(1:4, function(s) summ(PM,s,1))), ncol=4) 
      
      a = acc[1,]-acc[3,]; b = acc[4,]-acc[1,]
      cat(paste0("\n Accuracy status ",all(c(a[3:4],b[3:4])>=0),"\n"))
      
    }## prevS if
  } # prevS for
  
}#t>1


# mean, mode
yTe = matrix(rep(ty, nd),ncol = nd)
Pmode = apply(PM, 1, ymode)
mmacc[1] = mean(PM == yTe)
mmacc[2] = mean(Pmode == ty)

#remove(fml, baseyl,  Bmod2)

save(acc, mmacc, file = paste0("resCat_20-01-12_t",t,".RData"))


##########################################################################
##########################################################################
### Accuracy - Continuous and Binary Models 

print(t)

d = ted[time==200*t,]

load( paste0("ContBin20-01-12_t",t,"_setting",setting,".RData") ) 

rlacc = matrix(NA, 4, 2)

# Binary Model

fv = fmlvars(fml1, itc = 0)
y =  unlist(ted[time==200*t,fv$LHS,with = F])

dat.pred = ted[time==200*t,fv$RHS,with=F]

PM = predict_bart(mod1,0, dat.pred)$samp_y

rlacc[2:4, 1] = summRL(PM)
rlacc[1, 1] = mean(y)

# Continuous Model

fv = fmlvars(fml2, itc = 0)
y =  unlist(ted[time==200*t,fv$LHS,with = F])

dat.pred = ted[time==200*t,fv$RHS,with=F]

upd = unlist(ted[time==200*t, theX[1],with = F])
ind = which(upd==1)
y = y[ind]
dat.pred = dat.pred[ind,]
PM = predict_bart(mod2,0, dat.pred)$samp_y

rlacc[2:4, 2] = summRL(PM)
rlacc[1, 2] = mean(y)

#remove(fml1, fml2, mod1, mod2)

save(rlacc, file = paste0("resRL_20-01-12_t",t,".RData"))

##########################################################################
if(0){ 
  maxt = 10
  # the 1st dim: y, mean post, 5% post, 95% post
  # the 2nd dim: time from 1 to maxt
  # the 3rd dim: eng, diseng, xfer, death
  resCat = array(NA, dim = c(4,maxt,4)) 
  
  # row 1: posterior mean
  # row 2: posterior mode
  # cols: time from 1 to maxt
  resAcc = matrix(NA, 2, maxt)
  rownames(resAcc) = c("Test Mean", "Test Mode")
  
  for(t in 1:maxt){
    load( paste0("resCat_20-01-12_t",t,".RData"))
    resCat[, t, ] = acc
    resAcc[,t] = mmacc
    remove(acc, mmacc)
  }
  
  save(resCat, resAcc, file = "resCat_20-01-12.RData")
  
  # the first dim: y, mean post, 5% post, 95% post
  # the third dim: R, L
  resRL = array(NA, dim = c(4,maxt,2)) 
  
  for(t in 1:maxt){
    load(paste0("resRL_20-01-12_t",t,".RData"))
    resRL[, t, ] = rlacc
    remove(rlacc)
  }
  
  save(resRL, file = "resRL_20-01-12.RData")
  
  ############################################
  
  load("resRL_20-01-12.RData")
  load("resCat_20-01-12.RData")
  maxt = 10
  
  ### Plot Outcome Model Accuracy
  
  pdS = as.data.frame(matrix(NA, nrow = 4*maxt, ncol = 6))
  
  colnames(pdS) = c("Outcome", "Time", "obs", "Posteriors", "p5", "p95")
  pdS[,1] = rep(c("Disengaged","Engaged","Transferred","Died"), each = maxt)
  pdS[,2] = rep(1:maxt, 4)
  pdS[,3:6] = rbind(t(resCat[,,1]), t(resCat[,,2]), t(resCat[,,3]), t(resCat[,,4]) ) 
  pdS$Outcome = factor(pdS$Outcome, levels = c("Disengaged","Engaged","Transferred","Died"))
  pdS = pdS[-which(pdS[,1]=="Transferred" & pdS[,2]==1),]
  
  library("ggplot2")
  
  pS=ggplot(pdS, aes(Time)) + 
    geom_line(aes(y=Posteriors),colour="red") + 
    geom_point(aes(y=obs), colour="blue") +
    geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.5)+  facet_wrap(~ Outcome, scales="free_y")
  
  pdf("AccuracyS_20-01-12.pdf",width = 6, height = 4)
  #postscript(file="Accuracy19-03-10.eps",onefile=FALSE,width = 6, height = 4,paper = "special") 
  plot(pS)
  dev.off() 
  
  ### Plot CD4 Model Accuracy
  pd = as.data.frame(matrix(NA, nrow = 2*maxt, ncol = 6))
  
  pd[,1] = rep(c("Rt","Lt"), each = maxt)
  pd[,2] = rep(1:maxt, 2)
  pd[,3:6] = rbind(t(resRL[,,1]),t(resRL[,,2])) 
  colnames(pd) = c("Outcome", "Time", "obs", "Posteriors", "p5", "p95")
  pd$Outcome = factor(pd$Outcome)
  levels(pd$Outcome) = c("Continuous","Binary")
  
  library("ggplot2")
  
  p=ggplot(pd, aes(Time)) + 
    geom_line(aes(y=Posteriors),colour="red") + 
    geom_point(aes(y=obs), colour="blue") +
    geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.5)+  facet_wrap(~ Outcome, scales="free_y")+ 
    theme(strip.text.x = element_text(size = 12))
  
  #postscript(file="AccuracyRL19-03-10.eps",onefile=FALSE,width = 7, height = 3.5,paper = "special") 
  #pathto = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft3\\Accuracy Comparison - Diag Plots\\19-03-05\\Plots"
  #png(filename = paste0(pathto,"\\AccuracyRL19-03-10.png"), width =500, height = 200)
  #pdf("AccuracyRL_20-01-12.pdf",width = 7, height = 3.5)
  #plot(p)
  #dev.off()
  
  pdALL = rbind(pd,pdS)
  pdALL$Time = pdALL$Time * 200
  pALL=ggplot(pdALL, aes(Time)) + 
    geom_line(aes(y=Posteriors),colour="red") + 
    geom_point(aes(y=obs), colour="blue") +
    geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.5)+  facet_wrap(~ Outcome, ncol=2, scales="free_y")+ 
    theme(strip.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 10, angle = 90),
          axis.title.x = element_text(size = 10))
  
  pdf("Accuracy_20-01-12.pdf",width = 7, height = 7)
  #postscript(file="Accuracy19-03-10.eps",onefile=FALSE,width = 6, height = 4,paper = "special") 
  plot(pALL)
  dev.off() 
}

rm(list=ls()) 

# 0 if comparing Algorithm 1 & 2, 1 if comparing Algorithm 1 & 2 & 3
type = 1

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
cp = as.numeric(slurm_arrayid) #1:10

library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")

# load parameters based on index cp: compare sim (1,2), (3,4), (5,6), (7,8),...,(19,20)
sim0 = cp*2 - 1; sim1 = cp*2
sim = sim0 # not using KD nuC Vm from variable sim in plots.R
datseed = 99 # seed for simulating data in Setting.R
source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation



# load bmpfit0,bmpfit1,bmpfit2
load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,".RData"))
load(paste0("CK11-26_b",nb,"d",nd,"s",sim1,".RData"))
load(paste0("CK11-26_b",nb,"d",nd,"s",sim0,"mg",FALSE,".RData"))

#ylev = 1:3
pm1 = bmpfit0$ndim

### plot choices
#1: 4x4 whole diagnostic plots burn+npost. row: latents for ALg1, Alg2; col: percent acceptance, num nodes, num leaves, tree depth 
#2: 4x4 whole diagnostic plots npost.
#3: 2x2 Tree depth burn+npost.
#4: 2x2 Tree depth npost.
#5: 1x2 Sigma npost.
#6: 1x1 Sigma correlation histogram.

ch = c(2,4,5,6)
#ch = 5

if(type == 0){
  if(1 %in% ch){#1: 4x4 whole diagnostic plots burn+npost
    pdf(paste0("PLOT1_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 14, height = pm1*3.5)
    par(mfrow=c(pm1,4))
    ALLPlot2(bmpfit0, bmpfit1, by, all=T)
    dev.off()
  }
  
  if(2 %in% ch){#2: pm1 x 4 whole diagnostic plots npost
    pdf(paste0("PLOT2_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 14, height = pm1*3.5)
    par(mfrow=c(pm1,4))
    ALLPlot2(bmpfit0, bmpfit1, by, all=F)
    dev.off()
  }
  
  if(3 %in% ch){#3: 1 x pm1 Tree depth burn+npost
    pdf(paste0("PLOT3_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 7, height = 3.5)
    pm1 = p-1
    par(mfrow=c(1,pm1))
    TDPlot2(bmpfit0,bmpfit1, by, all=T)
    dev.off()
  }
  
  if(4 %in% ch){#4: 1 x pm1 Tree depth npost
    pdf(paste0("PLOT4_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 7, height = 3.5)
    pm1 = p-1
    par(mfrow=c(1,pm1))
    TDPlot2(bmpfit0,bmpfit1, by, all=F)
    dev.off()
  }
  
  if(5 %in% ch){#5: 1x2 Sigma (sig11, sig12) npost; sig22 = C - sig11
    pdf(paste0("PLOT5_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 7, height = 3.5)
    ns = 2 #num of entries in Sigma = (pm1^2 -pm1)/2 + pm1m = 3
    ind = (1:nd - 1)*pm1*pm1 # index +1 for sig11, +2 for sig12, +4 for sig22
    par(mfrow=c(1,ns))
    a = bmpfit0$sigmasample ; b = bmpfit1$sigmasample
    for(s in c(1,2)){
      scplot2(1:nd, a[ind + s], b[ind + s], 
              xlab = "MCMC Iteration", ylab = ifelse(s==1, "Sig11", "Sig12"), 
              main = "")      
    }
    
    dev.off()
  }
  
  if(6 %in% ch){#6: 1x1 Sigma correlation histogram
    pdf(paste0("PLOT6_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 4, height = 4)
    par(mfrow=c(1,1))
    ind = (1:nd - 1)*pm1*pm1 +2
    
    a = bmpfit1$sigmasample_prior; #prior kindo
    b = bmpfit1$sigmasample; #posterior kindo
    a0 = bmpfit0$sigmasample_prior; #prior my method
    b0 = bmpfit0$sigmasample; #posterior my method
    
    p1 <- hist(c(a[ind],a0[ind]), plot=FALSE, breaks=30); p1$counts=p1$counts/sum(p1$counts)
    p2 <- hist(b[ind], plot=FALSE, breaks=ceiling((max(b[ind])-min(b[ind]))/0.033)); p2$counts=p2$counts/sum(p2$counts)     
    p20 <- hist(b0[ind], plot=FALSE, breaks=ceiling((max(b0[ind])-min(b0[ind]))/0.033)); p20$counts=p20$counts/sum(p20$counts)     
    
    plot( p1, col=rgb(0,0,1,1/4), main = "", cex.main=0.8,
          ylim=range(c(p1$counts,p2$counts,p20$counts)), xlab = "Correlation", ylab="Density")  # first histogram
    plot( p20, col='#00000088',  add=T)  # third, my posterior
    plot( p2, col='#FF000088',  add=T)  # second, kindo posterior
    dev.off()
  }  
} 


if(type == 1){ # comparing all 3 Algorithms
  if(1 %in% ch){#1: 4x4 whole diagnostic plots burn+npost
    pdf(paste0("PLOT1_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 14, height = pm1*3.5)
    par(mfrow=c(pm1,4))
    ALLPlot3(bmpfit0, bmpfit1, bmpfit2, by, all=T)
    dev.off()
  }
  
  if(2 %in% ch){#2: pm1 x 4 whole diagnostic plots npost
    pdf(paste0("PLOT2_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 14, height = pm1*3.5)
    par(mfrow=c(pm1,4))
    ALLPlot3(bmpfit0, bmpfit1, bmpfit2, by, all=F)
    dev.off()
  }
  
  if(3 %in% ch){#3: 1 x pm1 Tree depth burn+npost
    pdf(paste0("PLOT3_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 7, height = 3.5)
    pm1 = p-1
    par(mfrow=c(1,pm1))
    TDPlot3(bmpfit0,bmpfit1, bmpfit2, by, all=T)
    dev.off()
  }
  
  if(4 %in% ch){#4: 1 x pm1 Tree depth npost
    pdf(paste0("PLOT4_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 7, height = 3.5)
    pm1 = p-1
    par(mfrow=c(1,pm1))
    TDPlot3(bmpfit0,bmpfit1, bmpfit2, by, all=F)
    dev.off()
  }
  
  if(5 %in% ch){#5: 1x2 Sigma (sig11, sig12) npost; sig22 = C - sig11
    pdf(paste0("PLOT5_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 7, height = 3.5)
    ns = 2 #num of entries in Sigma = (pm1^2 -pm1)/2 + pm1m = 3
    ind = (1:nd - 1)*pm1*pm1 # index +1 for sig11, +2 for sig12, +4 for sig22
    par(mfrow=c(1,ns))
    a = bmpfit0$sigmasample ; b = bmpfit1$sigmasample; cc = bmpfit2$sigmasample
    for(s in c(1,2)){
      scplot3(1:nd, a[ind + s], b[ind + s], cc[ind + s], 
              xlab = "MCMC Iteration", ylab = ifelse(s==1, "Sig11", "Sig12"), 
              main = "")      
    }
    
    dev.off()
  }
  
  if(6 %in% ch){#6: 1x1 Sigma correlation histogram
    pdf(paste0("PLOT6_setting",sim0,"x",sim1,"_lev",by,".pdf"),width = 4, height = 4)
    par(mfrow=c(1,1))
    ind = (1:nd - 1)*pm1*pm1 +2
    
    a = bmpfit1$sigmasample_prior; #prior kindo
    b = bmpfit1$sigmasample; #posterior kindo
    a0 = bmpfit0$sigmasample_prior; #prior alg 2
    b0 = bmpfit0$sigmasample; #posterior alg 2
    a3 = bmpfit2$sigmasample_prior; #prior alg 3
    b3 = bmpfit2$sigmasample; #posterior alg 3
    
    p1 <- hist(c(a[ind],a0[ind],a3[ind]), plot=FALSE, breaks=30); p1$counts=p1$counts/sum(p1$counts)
    p2 <- hist(b[ind], plot=FALSE, breaks=ceiling((max(b[ind])-min(b[ind]))/0.033)); p2$counts=p2$counts/sum(p2$counts)     
    p20 <- hist(b0[ind], plot=FALSE, breaks=ceiling((max(b0[ind])-min(b0[ind]))/0.033)); p20$counts=p20$counts/sum(p20$counts)     
    p23 <- hist(b3[ind], plot=FALSE, breaks=ceiling((max(b3[ind])-min(b3[ind]))/0.033)); p23$counts=p23$counts/sum(p23$counts)     
    
    plot( p1, col=rgb(0,0,1,1/4), main = "", cex.main=0.8,
          ylim=range(c(p1$counts,p2$counts,p20$counts,p23$counts)), xlab = "Correlation", ylab="Density")  # first histogram
    plot( p20, col='#00000088',  add=T)  # alg2 posterior black
    plot( p2, col='#FF000088',  add=T)  # alg1 kindo posterior red
    plot( p23, col='#0080FF88', add=T)  # alg3 posterior blue
    dev.off()
  }  
} 

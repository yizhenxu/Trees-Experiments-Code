rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
sim = as.numeric(slurm_arrayid)

library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")


# load parameters based on setting index sim
source(paste0("BartGcomp_19-11-22_Setting",sim,".R"))

# load bmpfit0,bmpfit1
#load(paste0("CK11-23-2019_setting",sim,"_Ref.RData")) 
load(paste0("CK11-23-2019_b",nb,"d",nd,"_s",sim,"_Ref.RData"))


### plot choices
#1: 4x4 whole diagnostic plots burn+npost. row: latents for ALg1, Alg2; col: percent acceptance, num nodes, num leaves, tree depth 
#2: 4x4 whole diagnostic plots npost.
#3: 2x2 Tree depth burn+npost.
#4: 2x2 Tree depth npost.
#5: 1x2 Sigma npost.
#6: 1x1 Sigma correlation histogram.

ch = 1:6

ylev = 1:3
pm1 = p-1

for(j in 1:length(ylev)){ # j is the ref level
  
  if(1 %in% ch){#1: 4x4 whole diagnostic plots burn+npost
    pdf(paste0("PLOT1_setting",sim,"_lev",j,".pdf"),width = 14, height = 14)
    par(mfrow=c(2*pm1,4))
    ALLPlot(bmpfit1[[j]], j, all=T)
    ALLPlot(bmpfit0[[j]], j, all=T)
    dev.off()
  }
  
  if(2 %in% ch){#2: 4x4 whole diagnostic plots npost
    pdf(paste0("PLOT2_setting",sim,"_lev",j,".pdf"),width = 14, height = 14)
    pm1 = p-1
    par(mfrow=c(2*pm1,4))
    ALLPlot(bmpfit1[[j]], j, all=F)
    ALLPlot(bmpfit0[[j]], j, all=F)
    dev.off()
  }
  
  if(3 %in% ch){#3: 2x2 Tree depth burn+npost
    pdf(paste0("PLOT3_setting",sim,"_lev",j,".pdf"),width = 7, height = 7)
    pm1 = p-1
    par(mfrow=c(2,pm1))
    TDPlot(bmpfit1[[j]], j, all=T)
    TDPlot(bmpfit0[[j]], j, all=T)
    dev.off()
  }
  
  if(4 %in% ch){#4: 2x2 Tree depth npost
    pdf(paste0("PLOT4_setting",sim,"_lev",j,".pdf"),width = 7, height = 7)
    pm1 = p-1
    par(mfrow=c(2,pm1))
    TDPlot(bmpfit1[[j]], j, all=F)
    TDPlot(bmpfit0[[j]], j, all=F)
    dev.off()
  }
  
  if(5 %in% ch){#5: 2x3 Sigma (sig11, sig12, sig22) npost
    pdf(paste0("PLOT5_setting",sim,"_lev",j,".pdf"),width = 10.5, height = 7)
    ns = 3 #num of entries in Sigma = (pm1^2 -pm1)/2 + pm1
    ind = (1:nd - 1)*pm1*pm1 # index +1 for sig11, +2 for sig12, +4 for sig22
    par(mfrow=c(2,ns))
    a = bmpfit1[[j]]$sigmasample ; b = bmpfit0[[j]]$sigmasample
    for(s in c(1,2,4)){
      scatter.smooth(1:nd, a[ind + s], lpars = list(col = "red", lwd = 3, lty = 3),
                     xlab = "MCMC Iteration", ylab = ifelse(s==1, "Sig11", ifelse(s==2,"Sig12", "Sig22")), 
                     main = "")      
    }
    for(s in c(1,2,4)){
      scatter.smooth(1:nd, b[ind + s], lpars = list(col = "red", lwd = 3, lty = 3),
                     xlab = "MCMC Iteration", ylab = ifelse(s==1, "Sig11", ifelse(s==2,"Sig12", "Sig22")), 
                     main = "" )      
    }
    
    dev.off()
  }
  
  if(6 %in% ch){#6: 1x1 Sigma correlation histogram
    pdf(paste0("PLOT6_setting",sim,"_lev",j,".pdf"),width = 4, height = 4)
    par(mfrow=c(1,1))
    ind = (1:nd - 1)*pm1*pm1 +2
    
    a = bmpfit1[[j]]$sigmasample_prior; #prior kindo
    b = bmpfit1[[j]]$sigmasample; #posterior kindo
    a0 = bmpfit0[[j]]$sigmasample_prior; #prior my method
    b0 = bmpfit0[[j]]$sigmasample; #posterior my method
    
    p1 <- hist(c(a[ind],a0[ind]), plot=FALSE, breaks=30); p1$counts=p1$counts/sum(p1$counts)
    p2 <- hist(b[ind], plot=FALSE, breaks=20); p2$counts=p2$counts/sum(p2$counts)     
    p20 <- hist(b0[ind], plot=FALSE, breaks=20); p20$counts=p20$counts/sum(p20$counts)     
    
    plot( p1, col=rgb(0,0,1,1/4), main = "Prior and Posteriors (Alg1 Red, Alg2 Green)", cex.main=0.8,
          ylim=range(c(p1$counts,p2$counts,p20$counts)), xlab = "Correlation", ylab="Density")  # first histogram
    plot( p2, col=rgb(1,0,0,0.5),  add=T)  # second, kindo posterior
    plot( p20, col=rgb(0,1,0,0.5),  add=T)  # third, my posterior
    
    dev.off()
  }
  
}
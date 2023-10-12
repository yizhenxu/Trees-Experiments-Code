### On my laptop
#devtools::install_github("yizhenxu/GcompBART")
library(GcompBART)
#source("BartGcomp_19-03-02_Athena0_prepare.R")

####################################################################
load("Mod19-03-05.RData")
setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft3\\Accuracy Comparison - Diag Plots\\19-03-05\\Convergence Plots")
# mod1, mod2, mod3_D, mod3_E

if(0){
  for(t in 1:length(mod1)){
    Bmod = mod1[[t]]
    png(paste0("R_",t,".png"),width = 600, height = 400)
    DiagPlot(Bmod,1)
    dev.off()
  }
  
  for(t in 1:length(mod2)){
    Bmod = mod2[[t]]
    png(paste0("L_",t,".png"),width = 600, height = 400)
    DiagPlot(Bmod,1)
    dev.off()
  }
  
  for(t in 1:length(mod3_D)){
    Bmod = mod3_D[[t]]
    nd = Bmod$ndraws
    Bmod$ndraws = ncol(Bmod$Tree_Depth)
    Bmod$type = "multinomial"
    png(paste0("S0_",t+1,".png"),width = 600, height = 400)
    DiagPlot(Bmod,1)
    dev.off()
    png(paste0("S0_",t+1,"_sig.png"),width = 600, height = 400)
    plot(rep(1:nd,each = (Bmod$ndim)^2), Bmod$sigmasample)
    dev.off()
  }
  
  for(t in 1:length(mod3_E)){
    Bmod = mod3_E[[t]]
    nd = Bmod$ndraws
    Bmod$ndraws = ncol(Bmod$Tree_Depth)
    Bmod$type = "multinomial"
    png(paste0("S1_",t,".png"),width = 600, height = 400)
    DiagPlot(Bmod,1)
    dev.off()
    png(paste0("S1_",t,"_sig.png"),width = 600, height = 400)
    plot(rep(1:nd,each = (Bmod$ndim)^2), Bmod$sigmasample)
    dev.off()
  }
  
}
# plot only tree depth

TDPlot <- function(bfit, byrow = TRUE){
  nd = bfit$ndraws
  pm1 = bfit$ndim
  par(mfrow=c(1,pm1))
  for(i in 1:pm1){
    a = bfit$Tree_Depth
    scatter.smooth(1:nd, a[i,(bfit$burn + bfit$fitMNP + 1):(ncol(a))], lpars =
                     list(col = "red", lwd = 3, lty = 3)
                   , xlab = "MCMC Iteration", ylab = "Average Tree Depth")
  }
  par(mfrow=c(1,1))
}


for(t in 1:length(mod3_E)){
  Bmod = mod3_E[[t]]
  png(paste0("S1_",t,"_TD_FixBytree.png"),width = 600, height = 200)
  TDPlot(Bmod)
  dev.off()
  
}

for(t in 1:length(mod3_D)){
  Bmod = mod3_D[[t]]
  nd = Bmod$ndraws
  Bmod$ndraws = ncol(Bmod$Tree_Depth)
  Bmod$type = "multinomial"
  png(paste0("S0_",t+1,"_TD.png"),width = 600, height = 200)
  TDPlot(Bmod)
  dev.off()
  
}
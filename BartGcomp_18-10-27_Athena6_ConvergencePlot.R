#source("BartGcomp_18-10-27_Athena0_prepare.R")

####################################################################
load("Mod18-10-27.RData")

# mod1, mod2, mod3_D, mod3_E

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

# plot only tree depth

TDPlot <- function(bfit, byrow = TRUE){
  nd = bfit$ndraws
  pm1 = bfit$ndim
  par(mfrow=c(1,pm1))
  for(i in 1:pm1){
    scatter.smooth(1:nd, bfit$Tree_Depth[i,], lpars =
                     list(col = "red", lwd = 3, lty = 3)
                   , xlab = "MCMC Iteration", ylab = "Average Tree Depth")
  }
  par(mfrow=c(1,1))
}
  
        
for(t in 1:length(mod3_E)){
  Bmod = mod3_E[[t]]
  nd = Bmod$ndraws
  Bmod$ndraws = ncol(Bmod$Tree_Depth)
  Bmod$type = "multinomial"
  png(paste0("S1_",t,"_TD.png"),width = 600, height = 200)
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
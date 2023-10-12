source("BartGcomp_19-03-10_Athena0_prepare.R")

thin = 3

ACFPlot <- function(bfit, nthin){
  nd = bfit$ndraws
  pm1 = bfit$ndim
  a = bfit$Tree_Depth
  par(mfrow=c(1,pm1))
  for(i in 1:pm1){
    b = a[i,(bfit$burn + bfit$fitMNP + 1):(ncol(a))]
    samp = a[seq(1,length(b),nthin+1)]
    acf(ts(samp, names = "Tree Depth"), lag.max = 50)
  }
  #par(mfrow=c(1,1))
}

for(t in 1:maxt){
  load(paste0("Cat19-03-10_",t,".RData"))
  
  if(t==1){
    pdf(paste0("S1_",t,"_ACF.pdf"),width = 7, height = 3.5)
    #postscript(file=paste0("S1_",t,"_ACF.eps"),onefile=FALSE,width = 7, height = 3.5,paper = "special") 
    ACFPlot(Bmod2,thin)
    dev.off()
    
  } else {
    pdf(paste0("S0_",t,"_ACF.pdf"),width = 7, height = 3.5)
    #postscript(file=paste0("S0_",t,"_ACF.eps"),onefile=FALSE,width = 7, height = 3.5,paper = "special") 
    ACFPlot(Bmod2[[1]],thin)
    dev.off()
    pdf(paste0("S1_",t,"_ACF.pdf"),width = 6, height = 2)
    #postscript(file=paste0("S1_",t,"_ACF.eps"),onefile=FALSE,width = 6, height = 2,paper = "special") 
    ACFPlot(Bmod2[[2]],thin)
    dev.off()
    
  }
}
  
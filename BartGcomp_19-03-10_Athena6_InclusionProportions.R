bfit = Bmod2
  pm1 = bfit$ndim
  

  setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\Poster and Presentation\\2019Committee Meeting")
  
  png(filename="VarImp.png", width = 600,height = 400)
  par(mfrow=c(1,pm1))
  par(mar=c(11,4,4,4))
  for(i in 1:pm1){
    tmp = bfit$Inclusion_Proportions[[i]]
    barplot(sort(tmp,decreasing = T), ylab = "Inclusion Proportion", las=2, main=paste0("latent ",i))
  }
  
  par(mfrow=c(1,1))
  dev.off()
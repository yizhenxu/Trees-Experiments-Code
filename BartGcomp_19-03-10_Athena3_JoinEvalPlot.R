# Join models at different t
# Accuracy matrix
# Convergence plot
# Accuracy plot

source("BartGcomp_19-03-10_Athena0_prepare.R")

load("ContBin.RData")

mod3_D =  vector("list", maxt-1) # categorical models for those previously disengaged
mod3_E =  vector("list", maxt) # categorical models for those previously engaged
baseylist = vector("list", maxt)
ntlist = vector("list", maxt)
seedlist = vector("list", maxt)

resCat = array(NA, dim = c(4,maxt,4)) # Accuracy - Categorical Models 
# the first dim: y, mean post, 5% post, 95% post
# the second dim: time from 1 to maxt
# the third dim: eng, diseng, xfer, death


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
  #par(mfrow=c(1,1))
}

### start
for(t in 1:maxt){
  load(paste0("Cat19-03-10_",t,".RData"))
  
  modfml[[3*t]] = fml
  
  resCat[, t, ] = res
  baseylist[[t]] = baseyl
  ntlist[[t]] = ntL
  seedlist[[t]] = seedvec
  
  if(t==1){
    mod3_E[[1]] = Bmod2
    #pdf(paste0("S1_",t,"_TD.pdf"),width = 7, height = 3.5)
    postscript(file=paste0("S1_",t,"_TD.eps"),onefile=FALSE,width = 7, height = 3.5,paper = "special") 
    TDPlot(Bmod2)
    dev.off()
    
  } else {
    mod3_D[[t-1]] = Bmod2[[1]]
    mod3_E[[t]] = Bmod2[[2]]
    remove(baseyl)
    
    #pdf(paste0("S0_",t,"_TD.pdf"),width = 7, height = 3.5)
    postscript(file=paste0("S0_",t,"_TD.eps"),onefile=FALSE,width = 7, height = 3.5,paper = "special") 
    TDPlot(Bmod2[[1]])
    dev.off()
    #pdf(paste0("S1_",t,"_TD.pdf"),width = 6, height = 2)
    postscript(file=paste0("S1_",t,"_TD.eps"),onefile=FALSE,width = 6, height = 2,paper = "special") 
    TDPlot(Bmod2[[2]])
    dev.off()
    
  }
  
  remove(Bmod2)
  remove(fml)
  remove(res)
  remove(ntL)
  remove(seedvec)
}

print(seedlist)
print(baseylist)
resCat

save(mod1, mod2, mod3_D, mod3_E, baseylist, ntlist, seedlist, modfml, file = "Mod19-03-10.RData")
save(resCat, file = "resCat.RData")

pd = as.data.frame(matrix(NA, nrow = 4*maxt, ncol = 6))

colnames(pd) = c("Outcome", "Time", "obs", "Posteriors", "p5", "p95")
pd[,1] = rep(c("Engaged","Disengaged","Transferred","Died"), each = maxt)
pd[,2] = rep(1:maxt, 4)
pd[,3:6] = rbind(t(resCat[,,1]), t(resCat[,,2]), t(resCat[,,3]), t(resCat[,,4]) ) 
pd$Outcome = factor(pd$Outcome, levels = c("Engaged","Disengaged","Transferred","Died"))
pd = pd[-which(pd[,1]=="Transferred" & pd[,2]==1),]

library("ggplot2")

p=ggplot(pd, aes(Time)) + 
  geom_line(aes(y=Posteriors),colour="red") + 
  geom_point(aes(y=obs), colour="blue") +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.5)+  facet_wrap(~ Outcome, scales="free_y")

#pdf("Accuracy19-03-10.pdf",width = 6, height = 4)
postscript(file="Accuracy19-03-10.eps",onefile=FALSE,width = 6, height = 4,paper = "special") 
plot(p)
dev.off()

##########################################################################
### Accuracy for Continuous and Binary Models

summRL = function(M){# M: testn rows, nd cols
  vecp = apply(M, 2, function(x) mean(x))
  res = mean(vecp)
  
  svecp = sort(vecp)
  lp = length(svecp)
  ci = svecp[c(ceiling(lp*0.05), floor(lp*0.95))]
  res = c(res,ci)
  
  return(res)
}

resRL = array(NA, dim = c(4,maxt,2)) 
# the first dim: y, mean post, 5% post, 95% post
# the third dim: R, L

## R
for(t in 1:maxt){
  
  print(t)
  
  d = ted[time==200*t,]
  
  fv = fmlvars(modfml[[3*(t-1)+1]], itc = 0)
  y =  unlist(ted[time==200*t,fv$LHS,with = F])
  
  dat.pred = ted[time==200*t,fv$RHS,with=F]
  
  PM = predict_bart(mod1[[t]],0, dat.pred)$samp_y
  
  resRL[2:4, t, 1] = summRL(PM)
  
  resRL[1, t, 1] = mean(y)
  
}

## L
for(t in 1:maxt){
  
  print(t)
  
  d = ted[time==200*t,]
  
  fv = fmlvars(modfml[[3*(t-1)+2]], itc = 0)
  y =  unlist(ted[time==200*t,fv$LHS,with = F])
  
  dat.pred = ted[time==200*t,fv$RHS,with=F]
  
  upd = unlist(ted[time==200*t, theX[1],with = F])
  ind = which(upd==1)
  y = y[ind]
  dat.pred = dat.pred[ind,]
  PM = predict_bart(mod2[[t]],0, dat.pred)$samp_y
  
  resRL[2:4, t, 2] = summRL(PM)
  
  resRL[1, t, 2] = mean(y)
  
}

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

postscript(file="AccuracyRL19-03-10.eps",onefile=FALSE,width = 7, height = 3.5,paper = "special") 
#pathto = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft3\\Accuracy Comparison - Diag Plots\\19-03-05\\Plots"
#png(filename = paste0(pathto,"\\AccuracyRL19-03-10.png"), width =500, height = 200)
plot(p)
dev.off()
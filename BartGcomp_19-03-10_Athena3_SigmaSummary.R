
source("BartGcomp_19-03-10_Athena0_prepare.R")

# plot only tree depth

SigSumm <- function(bfit){
  nd = bfit$ndraws
  pm1 = bfit$ndim
  sig = bfit$sigmasample
  ss = matrix(sig, nrow = pm1^2)
  sm = matrix(apply(ss,1,mean), nrow = pm1)
  
  tmp = apply(ss,1,function(x) sort(x)[c(ceiling(nd*0.05), floor(nd*0.95))])
  s05 = matrix(tmp[1,], nrow = pm1)
  s95 = matrix(tmp[2,], nrow = pm1)
  
  sm = round(sm,3);s05 = round(s05,3);s95 = round(s95,3)
  res = matrix(NA, pm1, pm1)
  for(i in 1:pm1){
    for(j in 1:pm1){
      res[i,j] = paste0(sm[i,j], " (",s05[i,j],", ",s95[i,j],")")
    }
  }
  return(res)
}


### start
SigList = vector("list",2*maxt-1)
ind = 1
for(t in 1:maxt){
  load(paste0("Cat19-03-10_",t,".RData"))
  
  if(t==1){
    SigList[[ind]] = noquote(SigSumm(Bmod2))
    ind = ind + 1
  } else {
    SigList[[ind]] = noquote(SigSumm(Bmod2[[1]]))
    ind = ind + 1
    SigList[[ind]] = noquote(SigSumm(Bmod2[[2]]))
    ind = ind + 1
  }
  
}

save(SigList, file = "SigList.RData")


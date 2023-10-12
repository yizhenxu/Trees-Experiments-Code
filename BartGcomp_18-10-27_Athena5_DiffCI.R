
source("BartGcomp_18-10-27_Athena0_prepare.R")

####################################################################
summ = function(vec){
  return(unlist(lapply(1:4, function(x) mean(vec==x))))
}
engci = function(i, Smat){
  res = apply(Smat[[i]],2, summ)
  return(res)
}
#m = vector("list",3);m[[1]] = matrix(sample(1:3,120,replace=T),20,6);m[[2]] = matrix(sample(1:3,120,replace=T),20,6);m[[3]] = matrix(sample(1:3,120,replace=T),20,6)
#lapply(1:3, function(i) engci(i, m))

#A = 1
load("Gcomp_A1.RData")

tmp1 = lapply(1:maxt, function(i) engci(i, Smat)) #in A = 1 and A = 1{CD4 < 350}  

remove(Smat)
remove(Xmat)

#A = 1{CD4 < 350}  
load("Gcomp_A2.RData")

tmp2 = lapply(1:maxt, function(i) engci(i, Smat)) #in A = 1 and A = 1{CD4 < 350}  

#tab
nd = ncol(tmp1[[1]])
tab = vector("list",4) #one for each outcome category
  
for(k in 1:4){
  tab[[k]] = matrix(NA, 3, maxt)
  
  for(i in 1:maxt){
    dif = tmp1[[i]][k,] - tmp2[[i]][k,]
    tab[[k]][1,i] = mean(dif)
    a = sort(dif)
    tab[[k]][2,i] = a[floor(nd*0.025)]
    tab[[k]][3,i] = a[ceiling(nd*0.975)]
    
  }
}


save(tab, file = "tab.RData")

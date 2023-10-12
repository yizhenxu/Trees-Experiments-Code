args <- commandArgs(trailingOnly=TRUE)
A = eval(parse(text=args[1]))

# theA_lag1 A = 0,1,2; A==2 means A = 1{CD4<=350}

source("BartGcomp_18-10-27_Athena0_prepare.R")

load(paste0("Gcomp_A",A,".RData"))
####################################################################

library(data.table)

tmp = unlist(Smat)

tmp1 = c(rep(2,length(Smat[[1]])),tmp)

longdat = data.table(time = rep(200*(0:maxt),each = length(Smat[[1]])), outcome = tmp1)
ldc = longdat[, .N, by = c("time","outcome")]

tt = sum(ldc$time==200)
if(tt == 3){
  ldc = rbind(ldc,list(time = c(0,0,0,200),outcome = c(1,3,4,3),N=c(0,0,0,0) ))
} else {
  ldc = rbind(ldc,list(time = c(0,0,0),outcome = c(1,3,4),N=c(0,0,0) ))
}

ldc = ldc[order(time,outcome)]
counts = matrix(ldc$N,nrow = 4)

ctmp = counts[1,]
counts[1,] = counts[2,]
counts[2,] = ctmp
#counts[1,1] = counts[2,1];counts[2,1]=0

rownames(counts) = c("engaged","diengaged","transferred out","died")
counts = simplify2array(lapply(1:ncol(counts), function(x) counts[,x]/sum(counts[,x])))
colnames(counts) = 200*(0:maxt)

save(counts,file = paste0("Counts_A",A,".RData"))

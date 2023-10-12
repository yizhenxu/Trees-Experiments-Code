args <- commandArgs(trailingOnly=TRUE)
A = eval(parse(text=args[1]))

# theA_lag1 A = 0,1,2; A==2 means A = 1{CD4<=350}

source("BartGcomp_18-10-27_Athena0_prepare.R")

load(paste0("Gcomp_A",A,".RData"))
####################################################################

library(data.table)


tmp = unlist(Smat)

m = length(Smat[[1]])
tmp1 = c(rep(2,m),tmp[1:((maxt-1)*m)])

longdat = data.table(time = rep(200*(1:maxt),each = length(Smat[[1]])), outcome_lag = tmp1, outcome = tmp)
ldc = longdat[, .N, by = c("time", "outcome_lag", "outcome")]
ldc_denom = longdat[, .N, by = c("time","outcome_lag")]

ldc = ldc[-which(ldc$outcome_lag %in% c(3,4)),]
ldc = ldc[order(time, outcome_lag, outcome)]
ldc_denom = ldc_denom[-which(ldc_denom$outcome_lag %in% c(3,4)),]
ldc_denom = ldc_denom[order(time, outcome_lag)]

eldc = ldc[outcome_lag == 2,]
eldc = rbind(eldc,list(time = 200, outcome_lag = 2, outcome = 3, N = 0 ))
eldc = eldc[order(time, outcome, outcome_lag)]

dldc = ldc[outcome_lag == 1,]

eldc_denom = ldc_denom[outcome_lag == 2,]
dldc_denom = ldc_denom[outcome_lag == 1,]

engtp = matrix(eldc$N/rep(eldc_denom$N, each = 4), nrow = 4)
rownames(engtp) = c("diengaged","engaged","transferred out","died")
colnames(engtp) = unique(eldc$time)
dengtp = matrix(dldc$N/rep(dldc_denom$N, each = 3), nrow = 3)
rownames(dengtp) = c("diengaged","engaged","died")
colnames(dengtp) = unique(dldc$time)

save(engtp,dengtp,file = paste0("TP_A",A,".RData"))

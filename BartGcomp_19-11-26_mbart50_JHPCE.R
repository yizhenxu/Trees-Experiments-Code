# data simulated and accuracy summarized under settings 1 and 2

## run the current R script
# cd ./MPBART
# qsub mbash2.sh


rm(list=ls()) 

gc()
#memory.limit()
# set max memory usage is 10G
#memory.size(max=10000)


# coerce the value to an integer
simidx = as.numeric(Sys.getenv('SGE_TASK_ID'))# 1:100



book = matrix(NA, ncol = 2, nrow = 100)
colnames(book) = c("sim", "datseed")
book[,1] = rep(c(1,3) ,each = 50); book[,2] = rep(0:49, 2);

sim = book[simidx, 1]

library(BART)
library(GcompBART)
source("BartGcomp_19-11-22_Functions.R")
datseed = book[simidx, 2] # seed for simulating data in Setting.R
source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation
set.seed(CatSeed)

post <- mbart(x.train = trd[,1:6], y.train = trd[,7], x.test = newdat[,1:6], ndpost = nd, nskip = nb, keepevery = 1, ntree=nt)

tmp = accuracyfun_mbart(post, newdat$y, ntr, n)

save(tmp, file = paste0("./res/CK11-26_b",nb,"d",nd,"s",sim,"_ct",datseed,"_mbart.RData"))


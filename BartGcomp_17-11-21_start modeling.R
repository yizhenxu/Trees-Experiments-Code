setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")

setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")
library(BayesTree)
library("mlogit") # for function mFormula in myGcomp2
devtools::install_github("yizhenxu/myGcomp2", args="--no-multiarch")
library(myGcomp2)

library(readstata13)
library(lubridate)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(data.table)
library(zoo)

dat = fread("my_analysis_data1.csv")
length(unique(dat$patient_id))
dim(dat)
colnames(dat)

### select a subset of data to test on 
idset = dat[dat$time==0,]$patient_id

set.seed(1)
m = 5000 # 0 when using all data, m is the number of individuals used for analysis

if(m == 0){
  ind = 1:length(idset)
}else{
  ind = sample(1:length(idset),m,replace = F)
}

dat.sub = dat[dat$patient_id %in% idset[ind],] 
dim(dat.sub)

# OUTCOME: 0 disengage, 1 engage, 2 transfer_outcome, 3 death_outcome
# for modling log cd4, weight, and vl, only use those with missing indicator = 0
# V: ageinit, male, yrenrol, TravelTime, whostageinit, marriedinit, heightinit
# X: miss_weight, log_weight, miss_cd4, log_cd4plus1, miss_vl, vl_0, log_vlplus1



table(dat.sub[time==200,outcome])
sum(dat$miss_vl==0 & dat$time==0)





burn = 100;ndraws = 1000 ;ntree = 200;beta = 2

tmp = dat.sub[time==200*1, patient_id]

V = dat.sub[dat.sub$time == 0 & patient_id%in%tmp,c("ageinit", "male", "yrenrol", "TravelTime", "whostageinit", "marriedinit", "heightinit")] 
S1 = dat.sub$outcome[dat.sub$time==200]
A0 = dat.sub[dat.sub$time == 0 & patient_id%in%tmp,]$onarv_cf # on ARV status at enrollment


traind = cbind(A0,V,S1)

a0 = 1
testd = cbind(a0,V,S1) ; colnames(testd) = colnames(traind)


setpara = function(p){
  Priors = list(nu=p+2, V=(p+2)*diag(p-1), ntrees=ntree, kfac=3.0, pbd=1.0, pb=0.5 , beta = beta, alpha = 0.95, nc = 100, priorindep = FALSE, minobsnode = 10)
  Mcmcs = list(sigma0 = diag(p - 1), keep = 1, burn = burn, ndraws = ndraws, keep_sigma_draws=FALSE)
  return(list(Priors, Mcmcs))
}


l = length(unique(S1))


modmpbart = mpbart1(as.factor(S1) ~  1 | A0 + ageinit + male ,
                    train.data =  traind ,
                    test.data =  testd,
                    base = NULL,
                    varying = NULL,
                    sep = NULL,
                    Prior = setpara(l)[[1]],
                    Mcmc = setpara(l)[[2]],
                    seedvalue = 99)

dim(modmpbart$vec_class_test)
table(c(modmpbart$vec_class_test))
aa = modmpbart$vec_class_test[800:1000,]
tmp.s1.p = apply(aa,2,function(x) c(mean(x==1),mean(x==2),mean(x==3)))
apply(t(tmp.s1.p),2,mean)

head(modmpbart$class_prob_test)
apply(modmpbart$class_prob_test,2,mean)
apply(modmpbart$class_prob_train,2,mean)


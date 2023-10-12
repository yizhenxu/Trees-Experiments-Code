#dir.create("./Topic2")
setwd("./Topic2")

#devtools::install_github("yizhenxu/BARTGcompv1")
#remove.packages("BARTGcompv1")
library(BARTGcompv1)
library(BayesTree)
library(lubridate)
library(data.table)
library(zoo)

dat = fread("my_analysis_data.csv")

colnames(dat)

dat$cd4c = ifelse(is.na(dat$log_cd4plus1)==1, 2, ifelse(dat$log_cd4plus1<log(350+1), 0, 1))
#0      1      2
#335704 307105 123652

#outcome: 0 disengage 1 engage 2 xfer 3 death
#0      1      2      3
#294415 462699   2700   6647

dat[, cd4c_1 := c(NA,cd4c[-.N]), by = patient_id]
dat[, onarvcf_1 := c(NA,onarv_cf[-.N]), by = patient_id]



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

dat.sub$m0 = 1*(dat.sub$marriedinit == 0 & !is.na(dat.sub$marriedinit ))
dat.sub$m1 = 1*(dat.sub$marriedinit == 1 & !is.na(dat.sub$marriedinit ))
dat.sub$m2 = 1*(is.na(dat.sub$marriedinit ))

for(i in 2008:2016){
  dat.sub[, paste0("y",i) := 1*(dat.sub$yrenrol == i)]
}



t = 1 
tmpdat = dat.sub[time==200*t, ]

# Model X(t)
#V = dat.sub[dat.sub$time == 0 & patient_id%in%tmpdat$patient_id,c("ageinit", "male", "yrenrol", "TravelTime", "whostageinit", "marriedinit", "heightinit")] 

V = dat.sub[dat.sub$time == 0 & patient_id%in%tmpdat$patient_id,c("ageinit", "male", paste0("y",2008:2016), paste0("m",0:2))] 
Aprev = tmpdat$onarvcf_1 # on ARV status at enrollment
Xprev = tmpdat$cd4c_1
Xt = tmpdat$cd4c
table(Xt)
St = tmpdat$outcome
table(St)

a0 = 1
traind = cbind(Aprev,V,Xt)
testd = cbind(a0,V,Xt) ; colnames(testd) = colnames(traind)

l = length(unique(Xt))



burn = 100;ndraws = 1000 ;ntree = 200;beta = 2

setpara = function(p){
  Priors = list(nu=p+2, V=(p+2)*diag(p-1), ntrees=ntree, kfac=3.0, pbd=1.0, pb=0.5 , beta = beta, alpha = 0.95, nc = 100, priorindep = FALSE, minobsnode = 10)
  Mcmcs = list(sigma0 = diag(p - 1), keep = 1, burn = burn, ndraws = ndraws, keep_sigma_draws=FALSE)
  return(list(Priors, Mcmcs))
}



#print(c(paste0("+ y", 2008:2016),paste0("+ m",0:2)) ,quote=F)
ptm = proc.time()
modmpbart = mpbart1(as.factor(Xt) ~  1 | Aprev + ageinit + male+ y2008 + y2009 + y2010 + y2011 +
                      y2012 + y2013 + y2014 + y2015 + y2016 + m0    + m1    + m2   ,
                    train.data =  traind ,
                    test.data =  testd,
                    base = NULL,
                    varying = NULL,
                    sep = NULL,
                    Prior = setpara(l)[[1]],
                    Mcmc = setpara(l)[[2]],
                    seedvalue = 99)
proc.time()-ptm




dim(modmpbart$vec_class_test)
prob_test = apply(modmpbart$vec_class_test,2,function(x) c(mean(x==1),mean(x==2),mean(x==3)))




table(c(modmpbart$vec_class_test))
aa = modmpbart$vec_class_test[800:1000,]
tmp.s1.p = apply(aa,2,function(x) c(mean(x==1),mean(x==2),mean(x==3)))
apply(t(tmp.s1.p),2,mean)

head(modmpbart$class_prob_test)
apply(modmpbart$class_prob_test,2,mean)
apply(modmpbart$class_prob_train,2,mean)




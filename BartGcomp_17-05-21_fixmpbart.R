#devtools::install_github("yizhenxu/myGcomp2")
devtools::install_github("yizhenxu/myGcomp2", args="--no-multiarch")
library(myGcomp2)


library(data.table)
options(scipen = 10)
#dat<- read.csv("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\my_dat200days_from_enrol.csv", header=T, sep=",")
dat<- read.csv("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\my_dat200days_from_enrol.csv", header=T, sep=",")
#dat<- read.csv("my_dat200days_from_enrol.csv", header=T, sep=",")
head(dat)
dat$Eng_1 = as.numeric(dat$state_1 == 1)
dat$Dis_1 = as.numeric(dat$state_1 == 2)

dat$arv_E = dat$Eng_1*dat$onarvcf
dat$arv_D = dat$Dis_1*dat$onarvcf

# <350:0, >=350: 1, missing: 2
dat$cd4c<-ifelse(is.na(dat$cd4)==1, 2, ifelse(dat$cd4<350, 0, 1))

#lagged cd4 category
# create lagged onarv
dat = as.data.table(dat)
dat[, cd4c_1 := c(NA,cd4c[-.N]), by = patient_id]
dat[, cd4_1 := c(NA,cd4[-.N]), by = patient_id]
dat[, onarvcf_1 := c(NA,onarvcf[-.N]), by = patient_id]
dat[, maxobsgrid := max(gridtime), by = patient_id]
setDF(dat)

#baseline age category : age at enrollment is available
dat$bage = as.factor(ifelse(dat$ageatenroll<35, 0, 1))


#reduce data size
dat.sub = subset(dat, select=c("patient_id", "gridtime", "male", "bage", "ageatenroll",
                               "cd4", "cd4_1","cd4c", "cd4c_1", "onarvcf", "onarvcf_1",
                               "Eng_1", "Dis_1", "arv_E", "arv_D",
                               "state", "state_1") )
#assign new subject id : 1 to n
dat.sub$bage = as.numeric(dat.sub$bage)
dat.sub$idx=as.numeric(factor(dat.sub$patient_id))
dat.sub$cd4c_11 = 1*(dat.sub$cd4c_1==1)
dat.sub$cd4c_10 = 1*(dat.sub$cd4c_1==0)
dat.sub$cd4c1 = 1*(dat.sub$cd4c==1)
dat.sub$cd4c0 = 1*(dat.sub$cd4c==0)
dat.sub$state_1 = dat.sub$state_1-1 # prev state  state_1 in {0 engage, 1 disengage}
dat.sub$state_1 = 1 - dat.sub$state_1 # prev state  state_1 in {1 engage, 0 disengage}
# state in {1,2,3,4,5}
dat.sub$cd4c = dat.sub$cd4c+1 # cd4c category as {1,2,3}
n = length(unique(dat.sub$idx))
n
#######################################################################
sum(dat.sub$gridtime==0 & !is.na(dat.sub$cd4)) #31535
temp = dat.sub[dat.sub$gridtime==0,]
ind = which(!is.na(temp$cd4))
idset = temp[ind,]$patient_id
set.seed(1)

m = 5000 # 0 when using all data, m is the number of individuals used for analysis

if(m == 0){
  ind = 1:length(idset)
}else{
  ind = sample(1:length(idset),m,replace = F)
}



dat.sub = dat.sub[dat.sub$patient_id %in% idset[ind],] #31535 patients, dim 232252x22
dim(dat.sub)

# take out the one transfer out at grid one for now
dat1 = dat.sub[dat.sub$gridtime==1,]
table(dat1$state)
pid.s3 = dat1$patient_id[dat1$state==3]
ind.discard = which(dat.sub$patient_id %in% pid.s3)
if(length(ind.discard)!=0) dat.sub = dat.sub[-ind.discard,]

n = length(unique(dat.sub$idx))
n


burn = 100;ndraws = 1000 ;ntree = 200;beta = 2

V = dat.sub[dat.sub$gridtime == 0,c("male","ageatenroll")] #baseline gender age
X0 = dat.sub$cd4[dat.sub$gridtime==0]
A0 = dat.sub[dat.sub$gridtime == 0,]$onarvcf # on ARV status at enrollment
traind = cbind(A0,V,X0)

a0 = 1
testd = cbind(a0,V,X0) ; colnames(testd) = colnames(traind)


setpara = function(p){
  Priors = list(nu=p+2, V=(p+2)*diag(p-1), ntrees=ntree, kfac=3.0, pbd=1.0, pb=0.5 , beta = beta, alpha = 0.95, nc = 100, priorindep = FALSE, minobsnode = 10)
  Mcmcs = list(sigma0 = diag(p - 1), keep = 1, burn = burn, ndraws = ndraws, keep_sigma_draws=FALSE)
  return(list(Priors, Mcmcs))
}


S1 = dat.sub$state[dat.sub$gridtime==1]
l = length(unique(S1))
traind = cbind(S1,A0,V,X0)

modmpbart = mpbart1(as.factor(S1) ~  1 | A0 + ageatenroll + male + X0,
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

library(Hmisc)
s1.mpbart.p  = as.numeric(rMultinom(modmpbart$class_prob_test,1))
table(s1.mpbart.p)/length(s1.mpbart.p)

#########################################
install.packages("MNP")
library("MNP")
modmnp = mnp(as.factor(S1) ~  A0 + ageatenroll + male + X0,
                    data =  traind,n.draws=1000,base="5")
s1.mnp.p = predict(modmnp,testd)
s1.mnp.pp = apply(s1.mnp.p$y,2,function(x) c(mean(x==0),mean(x==1),mean(x==2)))
head(t(s1.mnp.pp))
apply(s1.mnp.pp,1,mean)
# results from mnp looks all good
apply(modmnp$param[,11:13],2,mean)

#########################################
S1 = dat.sub$state[dat.sub$gridtime==1]
table(S1[A0==1])/sum(A0)
table(S1)/length(S1)
#########################################

library(Hmisc)
library(nnet)
modmult = multinom(S1 ~ A0 + ageatenroll + male + X0 , data=traind)
cp = predict(modmult, newdata = testd, "probs")
apply(cp,2,mean)

set.seed(8)
s1.mult =  rMultinom(cp, 1)

table(s1.mult)/length(s1.mult)

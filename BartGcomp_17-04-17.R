########################################
########################################
########################################
### LOG
### 1. test passed for function bart1:
###    subTest (N) : number of subjects in test data
###    ndraws (K) : number of posterior draws to keep
###    V : baseline variables
###    E1star (E1*) : first round of posterior draws, K draws for each subject
###    -------------------------
###    E1*
###    V1 -> e1*(11) ... e1*(1K)
###    .         .    .     .
###    .         .    .     .
###    .         .    .     .
###    VN -> e1*(N1) ... e1*(NK)
###    --------------------------
###    E2star (E2*) : second round of posterior draws, 1 draw for each e1*(nk)
###    ---------------------------------------------------
###    E2*
###    [e1*(11),V1] -> e2*(11)  ...  [e1*(1K),V1] -> e2*(1K)
###          .            .      .         .           .
###          .            .      .         .           .
###          .            .      .         .           .
###    [e1*(N1),VN] -> e2*(N1)  ...  [e1*(NK),VN] -> e2*(NK)
###    ---------------------------------------------------
###    where e2*(1k) ... e2*(Nk) are generated from the same kth posterior model
###    draw from bart at time 2
###    and posterior model draws are from the same chain
###    x.test input for bart1 should have the form
###    ------------
###    [e1*(11),V1]
###          .
###          .
###          .
###    [e1*(N1),VN]
###          .
###          .
###          .
###    [e1*(1K),V1]
###          .
###          .
###          .
###    [e1*(NK),VN]
###    ------------
###    x.test should have N x K rows
###    and the yhat.test output of bart1 :
###    --------------------
###    e2*(11) ... e2*(N1)
###      .      .    .
###      .      .    .
###      .      .    .
###    e2*(1K) ... e2*(NK)
###    --------------------
###    row number represents posterior draws, col number represents test subjects
###    this form is to correspond to the output from bart
###    the output needs to pnorm to get phat(Y=1) for binary outcomes
###
######################################
###
### 2. try doing the poster analysis for mpbart using multinomial states
###
###    (a) Does the results at grid 1 (200 days) match the proportion of the real data?
###
########################################
########################################
########################################

#setwd("C:\\Users\\yxu\\Google Drive\\Treatment Cascade - Stat Med\\Data\\GIC2016Aug\\Statadata")

library(data.table)

options(scipen = 10)
dat<- read.csv("C:\\Users\\yxu\\Google Drive\\Treatment Cascade - Stat Med\\Data\\GIC2016Aug\\Statadata\\dat200days_from_enrol.csv", header=T, sep=",")
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
dat[, onarvcf_1.t := c(NA,onarvcf[-.N]), by = patient_id]
setDF(dat)



#baseline age category : age at enrollment is available
dat$bage = as.factor(ifelse(dat$ageatenroll<35, 0, 1))

#reduce data size
dat.sub = subset(dat, select=c("patient_id", "gridtime", "male", "bage", "ageatenroll",
                               "cd4", "cd4_1","cd4c", "cd4c_1", "onarvcf", "onarvcf_1.t",
                               "Eng_1", "Dis_1", "arv_E", "arv_D",
                               "state", "state_1") )
#assign new subject id : 1 to n
colnames(dat.sub)[11] = "onarvcf_1"
dat.sub$bage = as.numeric(dat.sub$bage)
dat.sub$idx=as.numeric(factor(dat.sub$patient_id))
dat.sub$cd4c_11 = 1*(dat.sub$cd4c_1==1)
dat.sub$cd4c_10 = 1*(dat.sub$cd4c_1==0)
dat.sub$cd4c1 = 1*(dat.sub$cd4c==1)
dat.sub$cd4c0 = 1*(dat.sub$cd4c==0)
dat.sub$state_1 = dat.sub$state_1-1 # prev state  state_1 in {0 engage, 1 disengage}
# state in {1,2,3,4,5}
dat.sub$cd4c = dat.sub$cd4c+1 # cd4c category as {1,2,3}
n = length(unique(dat.sub$idx))
n
#######################################################################

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

n = length(unique(dat.sub$idx))
n

out = matrix(NA,4,4)
#colnames(out) = c("Perc_onARV","Perc_Eng","Perc_Eng|onARV1","Perc_Eng|onARV0","Perc_Eng|prevonARV1","Perc_Eng|prevonARV0")
colnames(out) = c("Perc_onARV","Perc_Eng","Perc_Eng|onARV1","Perc_Eng|onARV0")
rownames(out) = c(0,200,400,600)
for(k in 0:3){
  dat1 = dat.sub[dat.sub$gridtime==k ,]
  out[k+1,1] = mean(dat1$onarvcf)
  out[k+1,2] = mean(dat1$state==1)
  out[k+1,3] = mean(dat1[dat1$onarvcf==1,]$state==1)
  out[k+1,4] = mean(dat1[dat1$onarvcf==0,]$state==1)
  #out[k+1,5] = mean(dat1[dat1$onarvcf_1==1,]$state==1)
  #out[k+1,6] = mean(dat1[dat1$onarvcf_1==0,]$state==1)
}
out


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
# 0. replicate poster analysis for binary outcome

bdat = dat.sub

# settings and prior parameters for state predictions
burn = 100;ndraws = 1000 ;ntree = 200;beta = 2

## set up list of conditional models and fit the models to data
Scd4 <- list()
Cd4 <- list()


devtools::install_github("yizhenxu/myGcomp")
library(myGcomp)
library(Hmisc)  #rmultinom function
library(mpbart)
library(BayesTree)

########################################
#V = dat.sub[dat.sub$gridtime == 0,c("cd4","male","bage")] #baseline cd4 gender agecat
V = dat.sub[dat.sub$gridtime == 0,c("cd4","male","ageatenroll")] #baseline cd4 gender age
Y = (dat.sub[dat.sub$gridtime == 1,]$state==1)*1 # engagement indicator at 200 days
A = dat.sub[dat.sub$gridtime == 1,]$onarvcf # on ARV status of 0-200 days
a = 1
traind = cbind(A,V)
testd = cbind(a,V) ; colnames(testd) = colnames(traind)
Cd4[[1]] = bart(x.train = traind , y.train = Y, x.test = testd,
                power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)

#hist(apply(pnorm(Cd4[[1]]$yhat.test),1,mean))
#mean(apply(pnorm(Cd4[[1]]$yhat.test),1,mean))
set.seed(1)
tmp1 = c(Cd4[[1]]$yhat.test) #c(ndraws of ID 1, ... , ndraws of ID n)
Estar= rbinom(n*ndraws, 1, pnorm(tmp1))
E1star = matrix(Estar, nrow = ndraws, ncol = n)
post1 = apply(E1star,1,mean)
#hist(post1,breaks =15)
#mean(post1)
E1star = t(E1star)# E1*|A0=1,V    nrowTest x ndpost


#######################################
tmp2 = cbind(1:n,V) ; colnames(tmp2)[1] = "ID"
tmp2 = as.data.table(tmp2)
tmp2 = tmp2[rep(tmp2$ID, each = 2),]
tmp2 = setDF(tmp2)[,2:4]
testd1 = cbind(rep(0:1,n), tmp2) # (E1*, V)

ind = which(dat.sub[dat.sub$gridtime == 2,]$patient_id %in%dat.sub[dat.sub$gridtime == 1,]$patient_id)
V = V[ind,]

E1 = (dat.sub[dat.sub$gridtime == 2,]$state_1==0)*1 # 1 engage, 0 o.w.
Y = (dat.sub[dat.sub$gridtime == 2,]$state==1)*1
A = dat.sub[dat.sub$gridtime == 2,]$onarvcf
a = 1
traind = cbind(A,E1,V)
testd = cbind(a, testd1) ; colnames(testd) = colnames(traind)

Cd4[[2]] = bart(x.train = traind , y.train = Y, x.test = testd,
                power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)

tmp11 = c(Cd4[[2]]$yhat.test) #c(2*ndraws of ID 1, ... , 2*ndraws of ID n)
Estar= rbinom(2*n*ndraws, 1, pnorm(tmp11))
tmp11 = t(matrix(Estar, nrow = ndraws, ncol = 2*n))# E1*|A0=1,V

foo = function(x){
  aa = tmp11[c(2*x-1,2*x),] # first row if E1*=0, second row if E1*=1
  a1 = which(E1star[x,]==1)
  #a0 = which(E1star[x,]==0)
  res = aa[1,]
  res[a1] = aa[2,a1]
  return(res)
}
E2star = matrix(unlist(mapply(foo ,1:n)),nrow = ndraws,ncol=n)
post2 = apply(E2star,1,mean)
#hist(post2)


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
# 1. test for function bart1
tmp3 = cbind(1:nrow(V),V) ; colnames(tmp3)[1] = "ID"
tmp3 = as.data.table(tmp3)
tmp3 = tmp3[rep(tmp3$ID, ndraws),]
tmp3 = setDF(tmp3)[,2:ncol(tmp3)]

E1star.2 = E1star[ind,]
testd1 = cbind(c(E1star.2),tmp3)
testd = cbind(a, testd1) ; colnames(testd) = colnames(traind)

save(traind, Y, testd, beta, ntree, burn, ndraws,file = "C:\\Users\\yxu\\Desktop\\temp.RData")

# need to change bart1 such that
# kXTest = XTest[(nrowTest*(ndraws-1)+1):(nrowTest*ndraws),]
# and output as nrowTest x ndraws

load("C:\\Users\\yxu\\Desktop\\temp.RData")

mod1= bart1(x.train = traind , y.train = Y, x.test = testd,
                power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)

post2.bart1 = apply(pnorm(mod1$yhat.test),1,mean)
hist(post2.bart1)
# testing passssss!!! hurray
# should work well for continuous BART too

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
# 2. poster analysis for multinomial response
rm(list=ls())

# run the data

#bdat = dat.sub

# settings and prior parameters for state predictions
burn = 100;ndraws = 1000 ;ntree = 200;beta = 2

## set up list of conditional models and fit the models to data
Scd4 <- list()
Cd4 <- list()
setpara = function(p){
  Priors = list(nu=p-1, V= 0.5*diag(p - 1), ntrees=ntree, kfac=3.0, pbd=1.0, pb=0.5 , beta = beta, alpha = 0.95, nc = 100, priorindep = FALSE, minobsnode = 10)
  Mcmcs = list(sigma0 = diag(p - 1), keep = 1, burn = burn, ndraws = ndraws, keep_sigma_draws=FALSE)
  return(list(Priors, Mcmcs))
}

library(mpbart)

V = dat.sub[dat.sub$gridtime == 0,c("cd4","male","ageatenroll")] #baseline cd4 gender age
Y = dat.sub[dat.sub$gridtime == 1,]$state # status at 200 days
table(Y)/length(Y)
#1      2      5
#0.8686 0.0904 0.0410
A = dat.sub[dat.sub$gridtime == 1,]$onarvcf # on ARV status of 0-200 days
a = 1
traind = cbind(Y,A,V)
l = length(unique(traind[,1]))
testd = cbind(a,V) ; colnames(testd) = colnames(traind)[-1]

mod.mpb = mpbart(as.factor(Y) ~  1 | A + cd4 + male + ageatenroll,
                   train.data =  traind ,
                   test.data =  testd,
                   base = NULL,
                   varying = NULL,
                   sep = NULL,
                   Prior = setpara(l)[[1]],
                   Mcmc = setpara(l)[[2]],
                   seedvalue = 99)
tmp = mod.mpb$class_prob_test
apply(tmp,2,mean)
#1         2         5
#0.7768693 0.1152456 0.1078851

p = 3
x.train = traind[,-1]
y.train = Y
x.test = testd

mod.mpb1 = rmpbart1(x.train, y.train, x.test,
                   Prior = setpara(l)[[1]],
                   Mcmc = setpara(l)[[2]],
                   seedvalue = 99)
tmp1 = mod.mpb1$class_prob_test
apply(tmp1,2,mean)


mod.mpb0 = rmpbart(x.train, y.train, x.test,
                    Prior = setpara(l)[[1]],
                    Mcmc = setpara(l)[[2]],
                    seedvalue = 99)
tmp0 = mod.mpb0$class_prob_test
apply(tmp0,2,mean)

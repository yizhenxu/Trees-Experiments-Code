# 1. for presentation, present results for S2
# 2. because we are considering only treatment regime effect on patient retention, 
#    but we don't look at effect on mortality, so we should only include those
#    who will not die or transfer out within 400 days/ 2 time grids
# 3. by the problem set up, there is no LTFU at t = 0, 1, 2
#############################################################################
###    NOTES FOR BART1
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
library(myGcomp2)
#library(myGcomp)
#library(BayesTree)
library(data.table)

options(scipen = 10)
#dat<- read.csv("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\my_dat200days_from_enrol.csv", header=T, sep=",")
dat<- read.csv("my_dat200days_from_enrol.csv", header=T, sep=",")
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
#######################################################################
#######################################################################
#######################################################################
#######################################################################

# S1 | A0, X0, V

# 

burn = 200;ndraws = 1500 ;ntree = 200;beta = 2

V = dat.sub[dat.sub$gridtime == 0,c("male","ageatenroll")] #baseline gender age
X0 = dat.sub$cd4[dat.sub$gridtime==0]
A0 = dat.sub[dat.sub$gridtime == 0,]$onarvcf # on ARV status at enrollment
traind = cbind(A0,V,X0)

a0 = 1
testd = cbind(a0,V,X0) ; colnames(testd) = colnames(traind)

#set.seed(1)
set.seed(100)

S1_DEATH = (dat.sub$state[dat.sub$gridtime==1]==5)*1

modbart_DEATH = bart(x.train = traind , y.train = S1_DEATH , x.test = testd,
                power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)

ind = which(dat.sub$state[dat.sub$gridtime==1]==5)
traind1 = traind[-ind,]
S1_ENG = (dat.sub$state[dat.sub$gridtime==1]==1)*1
S1_ENG = S1_ENG[-ind]
modbart_ENG = bart(x.train = traind1 , y.train = S1_ENG , x.test = testd,
                   power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)


#########################################################
load("I:\\05-09_compare_bart_500draws_200burn.RData")
# 50647 save.image("05-09_compare_bart_500draws_100burn.RData")
# 173893 save.image("05-09_compare_bart_500draws_200burn.RData")
#173893 ndraws = 500 save.image("05-08_compare_bart.RData")
#174119 ndraws = 200 .............................._200draws....
#setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
#load("05-08_compare_bart_200draws.RData")
#load("05-08_compare_bart.RData")
#########################################################

s1_death = pnorm(modbart_DEATH$yhat.test)
s1_death.pp = matrix(s1_death,nrow=ndraws) # posterior Phat(S1=1)
tt = c(s1_death)
set.seed(1)
s1_death.p = matrix(rbinom(length(tt),1,tt),nrow=ndraws) # posterior draws 0-1

s1_eng = pnorm(modbart_ENG$yhat.test)
s1_eng.pp = matrix(s1_eng,nrow=ndraws) # posterior Phat(S1=1)
tt = c(s1_eng)
set.seed(8)
s1_eng.p = matrix(rbinom(length(tt),1,tt),nrow=ndraws) # posterior draws 0-1

tmp = c(s1_eng.p)
tmp[c(s1_death.p)==1] = 5
tmp[c(s1_eng.p)==0] = 2
s1.p = matrix(tmp, nrow = ndraws)

res.bart = apply(s1.p, 1, function(x) c(sum(x==1),sum(x==2),sum(x==5))/length(x))

res.bart = t(res.bart)
colnames(res.bart) = c("1","2","5")
apply(res.bart,2,mean)

S1 = dat.sub$state[dat.sub$gridtime==1]
table(S1[A0==1])/sum(A0)
table(S1)/length(S1)

S1.alive = S1[S1!=5]
table(S1.alive)/length(S1.alive)
apply(s1_eng,1,mean)


#######################################################
library(Hmisc)  #rmultinom function
library("mpbart")
library("mlogit") # for function mFormula in myGcomp2
devtools::install_github("yizhenxu/myGcomp2")
library("myGcomp2")
setpara1 = function(p){
  Priors = list(nu=p-1, V= 0.5*diag(p - 1), ntrees=ntree, kfac=3.0, pbd=1.0, pb=0.5 , beta = beta, alpha = 0.95, nc = 100, priorindep = FALSE, minobsnode = 10)
  Mcmcs = list(sigma0 = diag(p - 1), keep = 1, burn = burn, ndraws = ndraws, keep_sigma_draws=FALSE)
  return(list(Priors, Mcmcs))
}


S1 = dat.sub$state[dat.sub$gridtime==1]
l = length(unique(S1))
traind = cbind(S1,A0,V,X0)

modmpbart1 = mpbart1(as.factor(S1) ~  1 | A0 + ageatenroll + male + X0,
                   train.data =  traind , 
                   test.data =  testd,
                   base = NULL, 
                   varying = NULL,
                   sep = NULL,
                   Prior = NULL, 
                   Mcmc = NULL, 
                   seedvalue = 1)
       

s1.mpbart.p  = as.numeric(rMultinom(modmpbart1$class_prob_test,1))
table(s1.mpbart.p)/length(s1.mpbart.p)


dim(modmpbart1$vec_class_test)
table(c(modmpbart1$vec_class_test))/length(c(modmpbart1$vec_class_test))

aa = modmpbart1$vec_class_test
tmp.s1.p = apply(aa,2,function(x) c(mean(x==1),mean(x==2),mean(x==3)))
apply(t(tmp.s1.p),2,mean)

table(S1[A0==1])/sum(S1[A0==1])
#######################################################
setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
save.image("05-24_compare_mpbart_fixed_1000draws.RData")

save.image("05-08_compare_mpbart_200draws.RData")
load("05-08_compare_mpbart_200draws.RData")
#######################################################
#########################################
#install.packages("MNP")
library("MNP")
modmnp = mnp(as.factor(S1) ~  A0 + ageatenroll + male + X0,
             data =  traind,n.draws=1000)
s1.mnp.p = predict(modmnp,testd)
s1.mnp.pp = apply(s1.mnp.p$y,2,function(x) c(mean(x==0),mean(x==1),mean(x==2)))
head(t(s1.mnp.pp))
apply(s1.mnp.pp,1,mean)
# results from mnp looks all good
apply(modmnp$param[,11:13],2,mean)
#1:1        1:2        2:2 
#0.75455434 0.02378457 1.24544566 

#2:2        2:5        5:5 
#1.19486759 0.03530897 0.80513241 

save.image("05-24_compare_fixed_mpbart_mnp_1000draws.RData")
#########################################
S1 = dat.sub$state[dat.sub$gridtime==1]
table(S1[A0==1])/sum(A0)
table(S1)/length(S1)
#########################################


library(Hmisc) 
library(nnet)
modmult = multinom(S1 ~ A0 + ageatenroll + male + X0 , data=traind) 
cp = predict(modmult, newdata = testd, "probs") 
set.seed(8)
s1.mult =  rMultinom(cp, 1)  

table(s1.mult)/length(s1.mult)


getwd()
save.image("05-10_compare_all.RData")

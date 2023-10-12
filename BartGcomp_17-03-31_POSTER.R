rm(list=ls())
setwd("C:\\Users\\yxu\\Google Drive\\Treatment Cascade\\Data\\GIC2016Aug\\Statadata")
setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Treatment Cascade\\Data\\GIC2016Aug\\Statadata")


library(foreign)
library(nnet)
library(Hmisc)  #rmultinom function
library(ggplot2)


setwd("C:\\Users\\yxu\\Google Drive\\Treatment Cascade\\Data\\GIC2016Aug\\Statadata")

library(data.table)

options(scipen = 10)
dat<- read.csv("dat200days_from_enrol.csv", header=T, sep=",")
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

m = 0 # 0 when using all data, m is the number of individuals used for analysis

if(m == 0){
  ind = 1:length(idset)
}else{
  ind = sample(1:length(idset),m,replace = F)
}



dat.sub = dat.sub[dat.sub$patient_id %in% idset[ind],] #31535 patients, dim 232252x22
dim(dat.sub)
#sum(is.na(dat.sub[dat.sub$gridtime==0,]$cd4))

#idxset = unique(dat.sub$patient_id)
#dat.sub = dat.sub[dat.sub$patient_id%in% idxset[1:10000],]
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
library(Hmisc)  #rmultinom function
library(mpbart)
library(BayesTree)



# odd col: cd4. even col: state.
#Process_Data = matrix(NA, nrow = length(unique(dat.sub$patient_id)), ncol = 2*14)


bdat = dat.sub

# settings and prior parameters for state predictions
burn = 100;ndraws = 1000 ;ntree = 200;beta = 2

## set up list of conditional models and fit the models to data
Scd4 <- list()
Cd4 <- list()

ttt = proc.time()
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

set.seed(1)
tmp1 = c(Cd4[[1]]$yhat.test) #c(ndraws of ID 1, ... , ndraws of ID n)
Estar= rbinom(n*ndraws, 1, pnorm(tmp1))
E1star = matrix(Estar, nrow = ndraws, ncol = n)
post1 = apply(E1star,1,mean)
#hist(post1)

E1star = t(E1star)# E1*|A0=1,V


#######################################
tmp2 = cbind(1:n,rep(2,n),V) ; colnames(tmp2)[1:2] = c("ID","ndraw")
tmp2 = as.data.table(tmp2)
tmp2 = tmp2[rep(tmp2$ID, each = 2),]
tmp2 = setDF(tmp2)[,3:5]
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
  aa = tmp11[c(2*x-1,2*x),]
  a1 = which(E1star[x,]==1)
  a0 = which(E1star[x,]==0)
  res = aa[1,]
  res[a1] = aa[2,a1]
  return(res)
}
E2star = matrix(unlist(mapply(foo ,1:n)),nrow = ndraws,ncol=n)
post2 = apply(E2star,1,mean)
#hist(post2)
proc.time() - ttt

save.image(paste0("m",m,"ndraws",ndraws,"a",1,"age",".RData"))




##
#219048 m0ndraws1000a1 agecat
#222099 m0ndraws1000a1 age




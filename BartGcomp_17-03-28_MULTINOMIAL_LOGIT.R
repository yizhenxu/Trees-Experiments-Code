# do it for 3 time grids
# use continuous cd4
# use bartMachine
# ask to install rJava
# make the log and explain what caused that same prediction across time
# I did not generate label by predicted prob
# compare the class_prob_test of mpbart and multinom
# class_prob_test of mpbart are nearly deterministic

library(foreign)
library(nnet)
library(Hmisc)  #rmultinom function
library(ggplot2)
options(scipen = 10)

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
dat.sub$state_1 = dat.sub$state_1-1 # prev state  state_1 in {0, 1}
# state in {1,2,3,4,5}
dat.sub$cd4c = dat.sub$cd4c+1 # cd4c category as {1,2,3}
n = length(unique(dat.sub$idx))

#idxset = unique(dat.sub$patient_id)
#dat.sub = dat.sub[dat.sub$patient_id%in% idxset[1:1000],]
#n = length(unique(dat.sub$idx))

bdat = dat.sub
attach(bdat)
male0 = male[gridtime==0]
age0 = bage[gridtime==0]


## set up list of conditional models and fit the models to data
Scd4        = list()
Cd4        = list()

# time=1
# regenerate cd4 at time 1 using baseline information
# because cd4 is a confounder, we assume that X(0) -> A(0) -> X (1) and so forth
Cd4[[1]]       = multinom(cd4c ~ male + bage + onarvcf_1 + cd4c_1 , data=bdat[gridtime==1,])  #state_1==1 
Scd4[[1]]       = multinom(state ~ male + bage + cd4c + onarvcf, data=bdat[gridtime==1,]) 

#time=2
# cd4 at time t: now we can generate (X(t) | S(t-1), X(t-1), A(t-1), V)
for (i in 2:14){
  Cd4[[i]]   = multinom(cd4c ~ state_1 + male + bage + onarvcf_1 + cd4c_1, data=bdat[gridtime==i,])  
  
  Scd4[[i]]       = multinom(state ~ -1 + Eng_1 + arv_E + Dis_1 + arv_D 
                             + male + bage + cd4c + onarvcf, data=bdat[gridtime==i,])  
}

##### 2. G-comp under treat if CD4<350
cd4c.sim.t350  = matrix(0, nrow=n, ncol=(max(gridtime)-1))
Ssim.cd4.t350 = matrix(0, nrow=n, ncol=(max(gridtime)-1)) 

#at time=1
dat.cd4<-data.frame(cd4c_1 = cd4c_1[gridtime==1], male=male0, bage = age0, 
                    onarvcf_1 = ifelse(cd4c_1[gridtime==1]==0,1,0)  )
cp<-predict(Cd4[[1]], newdata = dat.cd4, "probs")  
cd4c.sim.t350[,1] <- rMultinom(cp, 1)  
dat.new<-data.frame(male=male0, bage = age0, cd4c = as.numeric(cd4c.sim.t350[,1]), 
                    onarvcf =ifelse(dat.cd4$onarvcf_1==1,1, ifelse(cd4c.sim.t350[,1]==0,1,0))  )
pp<-predict(Scd4[[1]], newdata = dat.new, "probs")  
Ssim.cd4.t350[,1]<-as.numeric(rMultinom(pp,1))

#time 2 - 14
#for (i in 2:14){
i=2
  dat.cd4<-data.frame(state_1 = Ssim.cd4.t350[,i-1], male=male0, bage = age0, 
                      cd4c_1 = cd4c.sim.t350[,i-1], onarvcf_1.t = dat.new$onarvcf   )
  cp<-predict(Cd4[[i]], newdata = dat.cd4, "probs")  
  cd4c.sim.t350[,i] <- rMultinom(cp, 1)   
  cd4c.sim.t350[,i][Ssim.cd4.t350[,i-1] %in% c(3,4,5)]<-cd4c.sim.t350[,i-1][Ssim.cd4.t350[,i-1] %in% c(3,4,5)]
  
  Eng_1.sim = ifelse(Ssim.cd4.t350[,i-1]==1,1,0)
  onarvcf.sim = ifelse(dat.new$onarvcf==1,1, ifelse(cd4c.sim.t350[,i-1]==0,1,0)) 
  arv_E.sim = Eng_1.sim*onarvcf.sim
  Dis_1.sim= ifelse(Ssim.cd4.t350[,i-1]==2,1,0)
  arv_D.sim = Dis_1.sim*onarvcf.sim
  
  dat.new<-data.frame(Eng_1=Eng_1.sim, arv_E = arv_E.sim,  Dis_1 = Dis_1.sim, arv_D  = arv_D.sim, 
                      male=male0, bage = age0, cd4c = cd4c.sim.t350[,i], 
                      onarvcf = onarvcf.sim)
  pp<-predict(Scd4[[i]], newdata = dat.new, "probs")  
  Ssim.cd4.t350[,i]<-as.numeric(rMultinom(pp,1))
  Ssim.cd4.t350[,i][Ssim.cd4.t350[,i-1] %in% c(3,4,5)]<- Ssim.cd4.t350[,i-1][Ssim.cd4.t350[,i-1] %in% c(3,4,5)]
#}

for (state in 1:5){  
  StateP.treat.cd4350[state,,sample] = apply(Ssim.cd4.t350==state, 2, mean)
}  


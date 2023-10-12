rm(list=ls())
gc()
options(java.parameters = "-UseGCOverheadLimit")
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

idxset = unique(dat.sub$patient_id)
dat.sub = dat.sub[dat.sub$patient_id%in% idxset[15001:25000],]
n = length(unique(dat.sub$idx))
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

library(Hmisc)  #rmultinom function
library(mpbart)
library(bartMachine)

set_bart_machine_num_cores(4)


# odd col: cd4. even col: state.
Process_Data = matrix(NA, nrow = length(unique(dat.sub$patient_id)), ncol = 2*14)
# insert a = 1
a = 1
bdat = dat.sub

# settings and prior parameters for state predictions
burn = 100;ndraws = 200 ;ntree = 200;beta = 2
setpara = function(p){
  Priors = list(nu=p-1, V= 0.5*diag(p - 1), ntrees=ntree, kfac=3.0, pbd=1.0, pb=0.5 , beta = beta, alpha = 0.95, nc = 100, priorindep = FALSE, minobsnode = 10)
  Mcmcs = list(sigma0 = diag(p - 1), keep = 1, burn = burn, ndraws = ndraws, keep_sigma_draws=FALSE)
  return(list(Priors, Mcmcs))
}

## set up list of conditional models and fit the models to data
Scd4 <- list()
Cd4 <- list()

tttt = proc.time()
# time point 1
# cd4 at time 1
traind = as.data.frame(bdat[bdat$gridtime==1,c("cd4","onarvcf_1", "cd4_1", "male", "bage")])
l = length(unique(traind[,1]))
testd = as.data.frame(bdat[bdat$gridtime==1,c("onarvcf_1", "cd4_1", "male", "bage")])
testd$onarvcf_1 = a

traind = traind[!is.na(traind[,1]),]
Cd4[[1]] = bartMachine(X = traind[,-1], y = traind[,1],num_trees = ntree,num_burn_in = burn,num_iterations_after_burn_in = ndraws,use_missing_data = T)
#Cd4[[1]]
Process_Data[,1] = predict(Cd4[[1]],testd)


# state at time 1
traind = as.data.frame(bdat[bdat$gridtime==1,c("state","onarvcf", "cd4", "male", "bage")])
l = length(unique(traind[,1]))
testd = as.data.frame(bdat[bdat$gridtime==1,c("onarvcf", "cd4", "male", "bage")])
testd$onarvcf = a
testd$cd4 = Process_Data[,1]

#Scd4[[1]] = rmpbart( x.train = traind[,-1] , y.train = as.factor(traind[,1]),
#                     x.test = testd, Prior = setpara(l)[[1]], Mcmc = setpara(l)[[2]], seedvalue = 99)

Scd4[[1]] = mpbart(as.factor(state) ~  1 | onarvcf + cd4 + male + bage,
                   train.data =  traind , 
                   test.data =  testd,
                   base = NULL, 
                   varying = NULL,
                   sep = NULL,
                   Prior = setpara(l)[[1]], 
                   Mcmc = setpara(l)[[2]], 
                   seedvalue = 99)
tmp = Scd4[[1]]$class_prob_test
ind = which(!is.na(colnames(tmp)))
tmp1 = tmp[,ind]/apply(tmp[,ind],1,sum)
Process_Data[,2] = as.numeric(rMultinom(tmp1,1)) #1-5


# time point 2 - 14
for(i in 2:14){
  
  # cd4 at the grid
  traind = as.data.frame(bdat[bdat$gridtime==i,c("cd4c","state_1","onarvcf_1", "cd4_1", "male", "bage")])
  traind = traind[!is.na(traind[,1]),]
  l = length(unique(traind[,1]))
  testd = as.data.frame(bdat[bdat$gridtime==1,c("state_1","onarvcf_1", "cd4_1", "male", "bage")])
  testd$state_1 =  as.numeric(Process_Data[,2*(i-1)])-1 # previous state predictions 1-5, minus 1 to get 0-4
  testd$onarvcf_1 = a  
  testd$cd4_1 = Process_Data[,2*(i-1)-1]
  idt = which(!is.na(testd$state_1) & testd$state_1 <=1) 
  testd = testd[idt,]
  if(l==2){
    Cd4[[i]] = bart(x.train = traind[,-1] , y.train = traind[,1]-1, x.test = testd,
                    power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)
    Process_Data[idt,2*i-1] =1*(apply(Cd4[[i]]$yhat.test,2,mean)>0)+1
  }
  if(l>2){
    
    Cd4[[1]] = bartMachine(X = traind[,-1], y = traind[,1],num_trees = ntree,num_burn_in = burn,num_iterations_after_burn_in = ndraws,use_missing_data = T)
    #Cd4[[1]]
    Process_Data[idt,2*i-1] = predict(Cd4[[i]],testd)
    
   
  }
  
  
  
  
  # state at the grid
  traind = as.data.frame(bdat[bdat$gridtime==i,c("state","state_1","onarvcf", "cd4c0", "cd4c1", "male", "bage")])
  l = length(unique(traind[,1]))
  testd = as.data.frame(bdat[bdat$gridtime==1,c("state_1","onarvcf", "cd4c0", "cd4c1", "male", "bage")])
  testd$state_1 =  as.numeric(Process_Data[,2*(i-1)])-1 # previous state predictions 1-5, minus 1 to get 0-4
  testd$onarvcf = a
  testd$cd4c0 = 1*(Process_Data[,2*i-1]==1) # current cd4c predictions 1-3
  testd$cd4c1 = 1*(Process_Data[,2*i-1]==2)
  idt = which(!is.na(testd$state_1) & testd$state_1 <=1)
  testd = testd[idt,]
  if(l==2){
    Scd4[[i]] = bart(x.train = traind[,-1] , y.train = traind[,1]-1, x.test = testd,
                     power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)
    Process_Data[idt,2*i] =1*(apply(Scd4[[i]]$yhat.test,2,mean)>0)+1
  }
  if(l>2){
    #Scd4[[i]] = rmpbart( x.train = traind[,-1] , y.train = as.factor(traind[,1]),
    #                     x.test = testd, Prior = setpara(l)[[1]], Mcmc = setpara(l)[[2]], seedvalue = 99)
    
    Scd4[[i]] = mpbart(as.factor(state) ~  1 | state_1 + onarvcf + cd4c0 + cd4c1 + male + bage,
                       train.data =  traind , 
                       test.data =  testd,
                       base = NULL, 
                       varying = NULL,
                       sep = NULL,
                       Prior = setpara(l)[[1]], 
                       Mcmc = setpara(l)[[2]], 
                       seedvalue = 99)
    Process_Data[idt,2*i] = as.numeric(rMultinom(Scd4[[i]]$class_prob_test,1)) #1-5
  }  
  
  
  
}

proc.time()-tttt

setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
save.image("a1ntree200beta2.RData")







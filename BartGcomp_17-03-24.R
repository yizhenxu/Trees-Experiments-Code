rm(list=ls())
setwd("C:\\Users\\yxu\\Google Drive\\Treatment Cascade - StatMed\\Data\\GIC2016Aug\\Statadata")

setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Treatment Cascade - StatMed\\Data\\GIC2016Aug\\Statadata")
library(foreign)
library(nnet)
library(Hmisc)  #rmultinom function
library(ggplot2)


setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Treatment Cascade\\Data\\GIC2016Aug\\Statadata")

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
#dat.sub = dat.sub[dat.sub$patient_id%in% idxset[1:10000],]
#n = length(unique(dat.sub$idx))

temp = dat.sub[dat.sub$gridtime==1 & dat.sub$onarvcf==0,]$state
table(temp)/length(temp)
#######################################################################

library(mpbart)
library(bartMachine)
bdat = dat.sub
p = 3
Prior1 = list(nu=p+2, V= diag(p - 1), ntrees=20, kfac=2.0, pbd=1.0, pb=0.5 , beta = 2.0, alpha = 0.95, nc = 100, priorindep = FALSE, minobsnode = 10)
Mcmc1 = list(sigma0 = diag(p - 1), keep = 1, burn = 100, ndraws = 200, keep_sigma_draws=FALSE)
traind = testd = as.data.frame(bdat[bdat$gridtime==1,])


timec = proc.time()

#multinom(cd4c ~ male + bage + onarvcf_1 + cd4c_1 , data=bdat[bdat$gridtime==1,], trace = F)
out1cd4 <- mpbart(as.factor(cd4c) ~  1 | onarvcf_1 + cd4c_11 +cd4c_12 + male + bage,
              train.data =  traind , 
              test.data =  testd,
              base = NULL, 
              varying = NULL,
              sep = NULL,
              Prior = Prior1, 
              Mcmc = Mcmc1, 
              seedvalue = 99)
proc.time() - timec

head(out1cd4$class_prob_train)
head(out1cd4$class_prob_test)
head(out1cd4$predicted_class_train)
head(out1cd4$predicted_class_test)

dim(out1cd4$class_prob_train)
dim(out1cd4$class_prob_test)
length(out1cd4$predicted_class_train)
length(out1cd4$predicted_class_test)

t1 = apply(out1cd4$class_prob_train,1,which.max)
t2 = apply(out1cd4$class_prob_test,1,which.max)

head(cbind(t1,out1cd4$predicted_class_train))
all(t1==out1cd4$predicted_class_train)

head(cbind(t2,out1cd4$predicted_class_test))
all(t2==out1cd4$predicted_class_test)

table(traind$cd4c, out1cd4$predicted_class_train)
table(testd$cd4c, out1cd4$predicted_class_test)
test_err <- sum(testd$cd4c !=
                  out1cd4$predicted_class_test)/length(testd$cd4c)
cat("test error :", test_err )



bdat = dat.sub
p = 3
Prior1 = list(nu=p+2, V= diag(p - 1), ntrees=20, kfac=2.0, pbd=1.0, pb=0.5 , beta = 2.0, alpha = 0.95, nc = 100, priorindep = FALSE, minobsnode = 10)
Mcmc1 = list(sigma0 = diag(p - 1), keep = 1, burn = 500, ndraws = 1, keep_sigma_draws=FALSE)
traind = testd = as.data.frame(bdat[bdat$gridtime==1,c("cd4c","onarvcf_1", "cd4c_11","cd4c_12", "male", "bage")])


timec = proc.time()

#multinom(cd4c ~ male + bage + onarvcf_1 + cd4c_1 , data=bdat[bdat$gridtime==1,], trace = F)
rout1cd4 <- rmpbart( x.train = traind[,-1] , y.train = as.factor(traind[,1]),
                     x.test = testd[,-1], Prior = Prior1,  Mcmc = Mcmc1, seedvalue = 99)

proc.time() - timec

table(traind[,1], rout1cd4$predicted_class_train )
table(testd[,1], rout1cd4$predicted_class_test)

test_err =  sum(testd[,1] !=
                  rout1cd4$predicted_class_test)/length(testd[,1])
cat("test error :", test_err )


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
### USE MPBART INSTEAD OF RMPBART
library(Hmisc)  #rmultinom function
library(mpbart)
#library(BayesTree)
#library(bartMachine)
#library(BayesTree)
#jgc <- function()
#{
#  .jcall("java/lang/System", method = "gc")
#} 
#options(java.parameters = "-Xmx16g")
#set_bart_machine_num_cores(4)
#jgc()
#c1 = bartMachine(X = traind[,-1], y = traind[,1], num_burn_in = 300, num_iterations_after_burn_in = 1,use_missing_data = TRUE)
#Process_Data[,1] = predict(c1, testd)
#c1 = bart(x.train = traind[,-1], y.train = traind[,1], x.test = testd,
#          nskip = 300, ndpost = 1)
### Notes
### Fit continuous cd4 by bartMachine only on those non-missing cd4 data
### in cd4 models, outcomes are complete, predictor cd4_1 may have missingness
### Fit multinomial state by mpbart

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
traind = as.data.frame(bdat[bdat$gridtime==1,c("cd4c","onarvcf_1", "cd4c_10", "cd4c_11", "male", "bage")])
l = length(unique(traind[,1]))
testd = as.data.frame(bdat[bdat$gridtime==1,c("onarvcf_1", "cd4c_10", "cd4c_11", "male", "bage")])
testd$onarvcf_1 = a

#Cd4[[1]] = rmpbart( x.train = traind[,-1] , y.train = as.factor(traind[,1]),
#                    x.test = testd, Prior = setpara(l)[[1]], Mcmc = setpara(l)[[2]], seedvalue = 99)
Cd4[[1]] = mpbart(as.factor(cd4c) ~  1 | onarvcf_1 + cd4c_10 +cd4c_11 + male + bage,
                  train.data =  traind , 
                  test.data =  testd,
                  base = NULL, 
                  varying = NULL,
                  sep = NULL,
                  Prior = setpara(l)[[1]], 
                  Mcmc = setpara(l)[[2]], 
                  seedvalue = 99)
Process_Data[,1] = as.numeric(rMultinom(Cd4[[1]]$class_prob_test,1)) #1-3


# state at time 1
traind = as.data.frame(bdat[bdat$gridtime==1,c("state","onarvcf", "cd4c0", "cd4c1", "male", "bage")])
l = length(unique(traind[,1]))
testd = as.data.frame(bdat[bdat$gridtime==1,c("onarvcf", "cd4c0", "cd4c1", "male", "bage")])
testd$onarvcf = a
testd$cd4c0 = 1*(Cd4[[1]]$predicted_class_test==1)
testd$cd4c1 = 1*(Cd4[[1]]$predicted_class_test==2)

#Scd4[[1]] = rmpbart( x.train = traind[,-1] , y.train = as.factor(traind[,1]),
#                     x.test = testd, Prior = setpara(l)[[1]], Mcmc = setpara(l)[[2]], seedvalue = 99)

Scd4[[1]] = mpbart(as.factor(state) ~  1 | onarvcf + cd4c0 + cd4c1 + male + bage,
                  train.data =  traind , 
                  test.data =  testd,
                  base = NULL, 
                  varying = NULL,
                  sep = NULL,
                  Prior = setpara(l)[[1]], 
                  Mcmc = setpara(l)[[2]], 
                  seedvalue = 99)
Process_Data[,2] = as.numeric(rMultinom(Scd4[[1]]$class_prob_test,1)) #1-5

# time point 2 - 14
for(i in 2:14){
  
  # cd4 at the grid
  traind = as.data.frame(bdat[bdat$gridtime==i,c("cd4c","state_1","onarvcf_1", "cd4c_10", "cd4c_11", "male", "bage")])
  l = length(unique(traind[,1]))
  testd = as.data.frame(bdat[bdat$gridtime==1,c("state_1","onarvcf_1", "cd4c_10", "cd4c_11", "male", "bage")])
  testd$state_1 =  as.numeric(Process_Data[,2*(i-1)])-1 # previous state predictions 1-5, minus 1 to get 0-4
  testd$onarvcf_1 = a  
  testd$cd4c_10 = 1*(Process_Data[,2*(i-1)-1]==1) # previous cd4c predictions 1-3
  testd$cd4c_11 = 1*(Process_Data[,2*(i-1)-1]==2)
  # only predict for those 0 (engage) and 1 (disengage) 
  # idt is index wrt ppl of the entire dataset
  idt = which(!is.na(testd$state_1) & testd$state_1 <=1) 
  testd = testd[idt,]
  if(l==2){
    Cd4[[i]] = bart(x.train = traind[,-1] , y.train = traind[,1]-1, x.test = testd,
                    power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)
    Process_Data[idt,2*i-1] =1*(apply(Cd4[[i]]$yhat.test,2,mean)>0)+1
  }
  if(l>2){
    #Cd4[[i]] = rmpbart( x.train = traind[,-1] , y.train = as.factor(traind[,1]),
    #                    x.test = testd, Prior = setpara(l)[[1]], Mcmc = setpara(l)[[2]], seedvalue = 99)
    
    Cd4[[i]] = mpbart(as.factor(cd4c) ~  1 | state_1 + onarvcf_1 + cd4c_10 +cd4c_11 + male + bage,
                      train.data =  traind , 
                      test.data =  testd,
                      base = NULL, 
                      varying = NULL,
                      sep = NULL,
                      Prior = setpara(l)[[1]], 
                      Mcmc = setpara(l)[[2]], 
                      seedvalue = 99)
    Process_Data[idt,2*i-1] = as.numeric(rMultinom(Cd4[[i]]$class_prob_test,1)) #1-3
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



for(i in 1:26){
  rr = table(Process_Data[,i])
  print(rr)
}

#a = 1 
#  all ppl --  118495
# burn = 100;ndraws = 200 ;ntree = 5;beta = 3
#  all ppl --  118536
# burn = 100;ndraws = 200 ;ntree = 200;beta = 3
#  all ppl --   118938
# burn = 100;ndraws = 200 ;ntree = 200;beta = 2


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
### USE MPBART INSTEAD OF RMPBART

library(Hmisc)  #rmultinom function
library(mpbart)


# odd col: cd4. even col: state.
Process_Data = matrix(NA, nrow = length(unique(dat.sub$patient_id)), ncol = 2*14)
# trt when cd4<350
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
traind = as.data.frame(bdat[bdat$gridtime==1,c("cd4c","onarvcf_1", "cd4c_10", "cd4c_11", "male", "bage")])
l = length(unique(traind[,1]))
testd = as.data.frame(bdat[bdat$gridtime==1,c("onarvcf_1", "cd4c_10", "cd4c_11", "male", "bage")])
testd$onarvcf_1 = testd$cd4c_10

#Cd4[[1]] = rmpbart( x.train = traind[,-1] , y.train = as.factor(traind[,1]),
#                    x.test = testd, Prior = setpara(l)[[1]], Mcmc = setpara(l)[[2]], seedvalue = 99)
Cd4[[1]] = mpbart(as.factor(cd4c) ~  1 | onarvcf_1 + cd4c_10 +cd4c_11 + male + bage,
                  train.data =  traind , 
                  test.data =  testd,
                  base = NULL, 
                  varying = NULL,
                  sep = NULL,
                  Prior = setpara(l)[[1]], 
                  Mcmc = setpara(l)[[2]], 
                  seedvalue = 99)
Process_Data[,1] = as.numeric(rMultinom(Cd4[[1]]$class_prob_test,1)) #1-3


# state at time 1
traind = as.data.frame(bdat[bdat$gridtime==1,c("state","onarvcf", "cd4c0", "cd4c1", "male", "bage")])
l = length(unique(traind[,1]))
testd = as.data.frame(bdat[bdat$gridtime==1,c("onarvcf", "cd4c0", "cd4c1", "male", "bage")])
testd$cd4c0 = 1*(Cd4[[1]]$predicted_class_test==1)
testd$cd4c1 = 1*(Cd4[[1]]$predicted_class_test==2)
testd$onarvcf = testd$cd4c0

#Scd4[[1]] = rmpbart( x.train = traind[,-1] , y.train = as.factor(traind[,1]),
#                     x.test = testd, Prior = setpara(l)[[1]], Mcmc = setpara(l)[[2]], seedvalue = 99)

Scd4[[1]] = mpbart(as.factor(state) ~  1 | onarvcf + cd4c0 + cd4c1 + male + bage,
                   train.data =  traind , 
                   test.data =  testd,
                   base = NULL, 
                   varying = NULL,
                   sep = NULL,
                   Prior = setpara(l)[[1]], 
                   Mcmc = setpara(l)[[2]], 
                   seedvalue = 99)
Process_Data[,2] = as.numeric(rMultinom(Scd4[[1]]$class_prob_test,1)) #1-5

# time point 2 - 14
for(i in 2:14){
  
  # cd4 at the grid
  traind = as.data.frame(bdat[bdat$gridtime==i,c("cd4c","state_1","onarvcf_1", "cd4c_10", "cd4c_11", "male", "bage")])
  l = length(unique(traind[,1]))
  testd = as.data.frame(bdat[bdat$gridtime==1,c("state_1","onarvcf_1", "cd4c_10", "cd4c_11", "male", "bage")])
  testd$state_1 =  as.numeric(Process_Data[,2*(i-1)])-1 # previous state predictions 1-5, minus 1 to get 0-4
  testd$cd4c_10 = 1*(Process_Data[,2*(i-1)-1]==1) # previous cd4c predictions 1-3
  testd$cd4c_11 = 1*(Process_Data[,2*(i-1)-1]==2)
  testd$onarvcf_1 = testd$cd4c_10  
  
  # only predict for those 0 (engage) and 1 (disengage) 
  # idt is index wrt ppl of the entire dataset
  idt = which(!is.na(testd$state_1) & testd$state_1 <=1) 
  testd = testd[idt,]
  if(l==2){
    Cd4[[i]] = bart(x.train = traind[,-1] , y.train = traind[,1]-1, x.test = testd,
                    power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)
    Process_Data[idt,2*i-1] =1*(apply(Cd4[[i]]$yhat.test,2,mean)>0)+1
  }
  if(l>2){
    #Cd4[[i]] = rmpbart( x.train = traind[,-1] , y.train = as.factor(traind[,1]),
    #                    x.test = testd, Prior = setpara(l)[[1]], Mcmc = setpara(l)[[2]], seedvalue = 99)
    
    Cd4[[i]] = mpbart(as.factor(cd4c) ~  1 | state_1 + onarvcf_1 + cd4c_10 +cd4c_11 + male + bage,
                      train.data =  traind , 
                      test.data =  testd,
                      base = NULL, 
                      varying = NULL,
                      sep = NULL,
                      Prior = setpara(l)[[1]], 
                      Mcmc = setpara(l)[[2]], 
                      seedvalue = 99)
    Process_Data[idt,2*i-1] = as.numeric(rMultinom(Cd4[[i]]$class_prob_test,1)) #1-3
  }
  
  
  
  
  # state at the grid
  traind = as.data.frame(bdat[bdat$gridtime==i,c("state","state_1","onarvcf", "cd4c0", "cd4c1", "male", "bage")])
  l = length(unique(traind[,1]))
  testd = as.data.frame(bdat[bdat$gridtime==1,c("state_1","onarvcf", "cd4c0", "cd4c1", "male", "bage")])
  testd$state_1 =  as.numeric(Process_Data[,2*(i-1)])-1 # previous state predictions 1-5, minus 1 to get 0-4
  testd$cd4c0 = 1*(Process_Data[,2*i-1]==1) # current cd4c predictions 1-3
  testd$cd4c1 = 1*(Process_Data[,2*i-1]==2)
  testd$onarvcf = testd$cd4c0
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
save.image("acd4350ntree200beta2.RData")




for(i in 1:26){
  rr = table(Process_Data[,i])
  print(rr)
}

#a = cd4<350 
#  all ppl --  118581
# burn = 100;ndraws = 200 ;ntree = 5;beta = 3
#  all ppl --   118880
# burn = 100;ndraws = 200 ;ntree = 200;beta = 3

#  all ppl --   118959
# burn = 100;ndraws = 200 ;ntree = 200;beta = 2

acd4350ntree5beta3.RData
acd4350ntree200beta3.RData
acd4350ntree200beta2.RData


########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################



setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
load("acd4350ntree200beta2.RData")


########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
library(zoo)
### make barplot
load("temp.RData")

Process_Data.t_350 = t(apply(Process_Data_350, 1, function(x) na.locf(x,na.rm=F)))
Process_Data.t = t(apply(Process_Data, 1, function(x) na.locf(x,na.rm=F)))

m=9
Process_Data1 = mapply(function(x) Process_Data.t[,2*x],1:m)
Process_Data1_350 = mapply(function(x) Process_Data.t_350[,2*x],1:m)
counteach = function(x){
  res = c()
  for(i in 1:5){
    res[i] = mean(x==i,na.rm=T)
  }
  return(res)
}
ta = mapply(function(x) counteach(Process_Data1[,x]), 1:m)
ta350 = mapply(function(x) counteach(Process_Data1_350[,x]), 1:m)

Probability_1 = c(ta)
Probability_350 = c(ta350)
Time = 200*(1:m)
Time = rep(Time, each =5)
State =  factor(rep(c("engage", "disengage","transfer", "LTFU", "death"),m),levels = c("death", "LTFU", "transfer", "disengage","engage"))
dat = data.frame(Time,State,Probability_1,Probability_350)

png(filename = "BARPLOTa1ntree200beta2.png",width = 550, height = 480)
p1 = ggplot(dat, aes(x=Time,y=Probability_1,fill=State)) + geom_bar(stat="identity") +
  ggtitle("Model Based State Probabilities") +
  labs(x="Days since enrollment",y="(Marginal) State Probability")
dev.off()


png(filename = "BARPLOTacd4350ntree200beta2.png",width = 550, height = 480)
ggplot(dat, aes(x=Time,y=Probability_350,fill=State)) + geom_bar(stat="identity") +
  ggtitle("Model Based State Probabilities") +
  labs(x="Days since enrollment",y="(Marginal) State Probability")
dev.off()

save(dat,file = "plotdata.RData")
##################################################################
library(ggplot2)
library(grid)
library(gridExtra)
p1 = ggplot(dat, aes(x=Time,y=Probability_1,fill=State)) + geom_bar(stat="identity") +
  ggtitle("Model Based State Probabilities onARV = 1") +
  labs(x="Days since enrollment",y="(Marginal) State Probability")
p2 = ggplot(dat, aes(x=Time,y=Probability_350,fill=State)) + geom_bar(stat="identity") +
  ggtitle("Model Based State Probabilities onARV = 1{CD4<350}") +
  labs(x="Days since enrollment",y="(Marginal) State Probability")
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1)
blank<-grid.rect(gp=gpar(col="white"))
t1 = arrangeGrob(p1 + theme(legend.position="none"),ncol=1,nrow=1)
t2 = arrangeGrob(p2 + theme(legend.position="none"),ncol=1,nrow=1)

png(filename = "BARPLOT.png",width = 1000, height = 500)
g4 = grid.arrange(t1,t2,
                  arrangeGrob(blank,mylegend,blank,nrow=1,widths=c(0.01,0.68,0.01)), 
                  nrow=1,ncol = 3, widths = c(1,1,0.3))
dev.off()

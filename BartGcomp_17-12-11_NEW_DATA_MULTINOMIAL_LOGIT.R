### IMPORTANT NOTE ON X SIMULATION: (different from the Young paper)
# When doing simulation for covariate, we simulate xxx_RX and xxx_upd. 
# Impute xxx_RX with 0 AND set miss_xxx to 1 if xxx_upd == 0 SINCE t=0 
# After the first xxx_upd == 1 in follow-up, miss_xxx is always 0
# miss_xxx is derivable from xxx_upd, so miss_xxx is only used as covariates and not part of the simulation

### THIS VERSION OF SIMULATION
# only categorical variables are considered FOR X here
# use multinomial regression
# missingness is considered as a category and miss_xxx is not used in analysis

#screen 228341
setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")
setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")

### PACKAGES

library(lubridate)# date processing
library(data.table)
library(zoo)# carry backward forward
library(foreign)#stata read in
library(nnet)# multinomial log-linear models
library(Hmisc)#rmultinom function
# plotting
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)


options(scipen = 10)

#setwd("./Topic2")
dat = fread("my_analysis_data.csv")
length(unique(dat$patient_id))
dim(dat)
colnames(dat)

# outcome: 0 disengage, 1 engage, 2 transfer_outcome, 3 death_outcome
######################################################################################################################
### SKIP THIS SESSION
# define variables from the previous grid: X_1, S_1, A_1
tmp = colnames(dat)
for(i in 1:ncol(dat)){
  if(i==1) tmp1 = paste0('\"',tmp[1],'\"')
  if(i>1) tmp1 =   paste0(tmp1, ",", '\"',tmp[i],'\"')
}
cat(tmp1)

theV = c("ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit")

theX = c("new_wupd","log_weight","log_weight_RX","miss_weight",
         "new_cupd","log_cd4plus1","log_cd4plus1_RX","miss_cd4",
         "new_vupd","vl_0","log_vlplus1","log_vlplus1_RX","miss_vl")

theY = "outcome"
######################################################################################################################
### DEFINE CATEGORICAL VARIABLES TO BE USED FOR ANALYSIS
getCAT = function(x,xcut,missind = NA){ # x1 < u <= x2 is category 2
  xcut = c(xcut, Inf)
  xc = unlist(lapply(x , function(u) which(u<=xcut)[1]))
  xc[is.na(xc)] = 0
  if(length(missind)>1) xc[missind==1] = 0
  return(xc)
}

dat$Clogweight = getCAT(dat$log_weight, log(c(50,55,60,65,75,100)))
dat$Clogcd4p1 = getCAT(dat$log_cd4plus1 , log(c(200,350,500))+1)
dat$Clog_vlp1 = getCAT(dat$log_vlplus1, log(c(1000,2600, 32000))+1) #Siberry 2012 vl to predict who stage 3,4

### DEFINE MISSING FOR WHOSTAGEINIT, MARRIEDINIT, AND HEIGHTINIT
# CATEGORIZE HEIGHTINIT
dat$Cheight = getCAT(dat$heightinit, c(150,155,160,165,170,175,180))
# ADD NA TO CATEGORY WHOSTAGE AND MARRIAGE
dat$whostageinit[is.na(dat$whostageinit)] = 0
dat$marriedinit[is.na(dat$marriedinit)] = 2 # because it already has value 0

theA = "onarv_cf"
theV = c("ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","Cheight")
theX = c("new_wupd","Clogweight",
         "new_cupd","Clogcd4p1",
         "new_vupd","vl_0","Clog_vlp1")
theY = "outcome"

theA_1 = paste0(theA,"_1")
theX_1 = paste0(theX,"_1")
theY_1 = paste0(theY,"_1")


### DEFINE LAGGED THEX AND THEY
# USE DATA.TABLE
lagfun = function(data,vname){
  vname1 = paste0(vname,"_1")
  data[, (vname1) := c(NA,get(vname)[-.N]), by = patient_id]
  return(data[,get(vname1)])
}

lagX = lapply(c(theX,theY,theA), function(vname) lagfun(dat,vname)  )
vname1 = paste0(c(theX,theY,theA),"_1")
dat[,(vname1) := lagX]


######################################################################################################################
### SET UP LIST OF CONDITIONAL MODELS AND FIT THE MODELS TO DATA

#order of X: "new_wupd"->"Clogweight"->"new_cupd"->"Clogcd4p1"->"new_vupd"->"vl_0"->"Clog_vlp1"
#order: X0->X1->S1->X2->S2->X3->S3->...->S15
#conditional models: X(t+1) ~ [S(t),X(t),A(t)], S(t+1) ~ [S(t),X(t+1),A(t)]
#in details, Xj(t+1) ~ [S(t),X(t),X_{1 to j-1}(t+1),A(t)]

### FACTORIZE CATEGORICAL VARIABLES WITH <= 10 LEVELS 
for(j in 1:ncol(dat)){
  vname = colnames(dat)[j]
  tmp = length(unique(dat[,get(vname)]))
  if(tmp <= 10) dat[,(vname) :=as.factor(get(vname))] 
}


maxt = 10
mod =  vector("list", (length(theX)+1)*maxt) 


cnt = 1
for(t in 1:maxt){
  print(t)
  
  if(t==1){
    
    for(p in 1:length(theX)){
      print(p)
      
      if(p == 1){
        fml = paste0(theX[1]," ~ ",paste(c(theV,theA_1,theX_1),collapse = " + "))
        mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
      } else {
        fml = paste0(theX[p]," ~ ",paste(c(theV,theA_1,theX_1,theX[1:(p-1)]),collapse = " + "))
        mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
      }
      
      cnt = cnt+1
      
    }#p
    
    fml = paste0(theY," ~ ",paste(c(theV,theA_1,theX),collapse = " + "))
    mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
    cnt = cnt+1
    
  } else {
    
    for(p in 1:length(theX)){
      print(p)
      
      if(p == 1){
        fml = paste0(theX[1]," ~ ",paste(c(theV,theA_1,theY_1,theX_1),collapse = " + "))
        mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
      } else {
        fml = paste0(theX[p]," ~ ",paste(c(theV,theA_1,theY_1,theX_1,theX[1:(p-1)]),collapse = " + "))
        mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
      }
      cnt = cnt+1
    }#p
    
    fml = paste0(theY," ~ ",paste(c(theV,theA_1,theY_1,theX),collapse = " + "))
    mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
    
  }
  cnt = cnt+1
}#t


getwd()
save(mod, file = "BartGcomp_17-12-11_mod.RData")
######################################################################################################################
### G comp

# create new ID
num = dat[, .N,by=patient_id]
dat$myID = rep(1:length(num[[2]]), times = num[[2]])


d = new
# create indicator variables
createind = function(d){
  newd = d
  for(j in 1:ncol(d)){
  
    if(is.factor(d[,get(colnames(d)[j])])){
      tmp = class.ind(d[,get(colnames(d)[j])])
      newd[,(colnames(d)[j]) := NULL]
      newd[, (paste0(colnames(d)[j], names(tmp)[-1])) := tmp[,-1]]
    }#if  
    
  }#j
  return(newd)
}


tmp = dat[time==0,]
res = matrix(NA, ncol = (length(theX)+1)*maxt, nrow = nrow(tmp))
rs = 1:nrow(tmp)

cnt = 1
for(t in 1:maxt){
  if(t==1){
    
    for(p in 1:length(theX)){
      print(p)
      
      if(p == 1){
        fml = paste0(theX[1]," ~ ",paste(c(theV,theA_1,theX_1),collapse = " + "))
        mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
        
        orname = c(theV,theA_1,theX_1)
        vname = c("myID",theV,theA,theX)
        new = dat[time==200*(t-1), vname, with=F]
        colnames(new) = orname
        
        
        
        newd = createind()
        prob = predict(mod[[1]], newdata = new, "probs")  
        res[,1] = rMultinom(prob, 1)  
        
          
      } else {
        fml = paste0(theX[p]," ~ ",paste(c(theV,theA_1,theX_1,theX[1:(p-1)]),collapse = " + "))
        mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
      }
      
      cnt = cnt+1
      
    }#p
    
    fml = paste0(theY," ~ ",paste(c(theV,theA_1,theX),collapse = " + "))
    mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
    cnt = cnt+1
    
  } else {
    
    for(p in 1:length(theX)){
      print(p)
      
      if(p == 1){
        fml = paste0(theX[1]," ~ ",paste(c(theV,theA_1,theY_1,theX_1),collapse = " + "))
        mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
      } else {
        fml = paste0(theX[p]," ~ ",paste(c(theV,theA_1,theY_1,theX_1,theX[1:(p-1)]),collapse = " + "))
        mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
      }
      cnt = cnt+1
    }#p
    
    fml = paste0(theY," ~ ",paste(c(theV,theA_1,theY_1,theX),collapse = " + "))
    mod[[cnt]] =  multinom(fml , data=dat[time==200*t,]) 
    
  }
  cnt = cnt+1 
  
   
}#t


#at time=1
dat.cd4<-data.frame(cd4c_1 = cd4c_1[gridtime==1], male=male0, bage = age0, 
                    onarvcf_1 = ifelse(cd4c_1[gridtime==1]==0,1,0)  )
cp<-predict(Cd4[[1]], newdata = dat.cd4, "probs")  
cd4c.sim.t350[,1] <- rMultinom(cp, 1)  
dat.new<-data.frame(male=male0, bage = age0, cd4c = as.numeric(cd4c.sim.t350[,1]), 
                    onarvcf =ifelse(dat.cd4$onarvcf_1==1,1, ifelse(cd4c.sim.t350[,1]==0,1,0))  )
pp<-predict(Scd4[[1]], newdata = dat.new, "probs")  
Ssim.cd4.t350[,1]<-as.numeric(rMultinom(pp,1))


######################################################################################################################
### BAR HIST OF OBSERVED OUTCOME

png("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\ObsEng.png")
a = dat$outcome
levels(a) =c("diengaged","engaged","transferred out","died")
a = relevel(a, "engaged")
counts <- table(a,dat$time)
barplot(counts, main="Observed Engagement Status",
        xlab="Days since baseline", ylab="Number of observations", col=1:length(levels(a)),
        legend = rownames(counts))
dev.off()
  
######################################################################################################################
### 2. G-comp under treat if CD4<350
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



######################################################################################################################

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


######################################################################################################################
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

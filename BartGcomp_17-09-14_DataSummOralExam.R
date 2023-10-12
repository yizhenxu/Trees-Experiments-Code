library(data.table)
options(scipen = 10)
#dat<- read.csv("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\my_dat200days_from_enrol.csv", header=T, sep=",")
dat<- read.csv("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\my_dat200days_from_enrol.csv", header=T, sep=",")

head(dat)
# state=1: engaged, 2: disengaged, 3: transferout, 4: ltfu, 5: death

dat$Eng_1 = as.numeric(dat$state_1 == 1)
dat$Dis_1 = as.numeric(dat$state_1 == 2)

dat$arv_E = dat$Eng_1*dat$onarvcf
dat$arv_D = dat$Dis_1*dat$onarvcf

# <350:0, >=350: 1, missing: 2
dat$cd4c<-ifelse(is.na(dat$cd4)==1, 2, ifelse(dat$cd4<350, 0, 1))

# create variables by patient
dat = as.data.table(dat)
dat[, cd4c_1 := c(NA,cd4c[-.N]), by = patient_id] #lagged CD4 category 
dat[, cd4_1 := c(NA,cd4[-.N]), by = patient_id] #lagged CD4
dat[, onarvcf_1 := c(NA,onarvcf[-.N]), by = patient_id] #lagged onARV
dat[, maxobsgrid := max(gridtime), by = patient_id] #max number of grids
setDF(dat)

#baseline age category (cutoff 35): age at enrollment is available
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

m = 0 # 0 when using all data, m is the number of individuals used for analysis

if(m == 0){
  ind = 1:length(idset)
}else{
  ind = sample(1:length(idset),m,replace = F)
}

dat.sub = dat.sub[dat.sub$patient_id %in% idset[ind],] #31535 patients, dim 232252x22
dim(dat.sub)
#######################################################################
# Summary of Observed Data

summtab = matrix(NA,2,8)
colnames(summtab) = c("j","tj","nj","Eng","Dis","Xfer","Death","onART")
rownames(summtab) = c("1","2")
summtab[,1] = c(1,2)
summtab[,2] = 200*summtab[,1]

fillsummtab = function(j){
  tempdat = dat.sub[dat.sub$gridtime==j,]
  temptab = table(tempdat$state)
  proptab = prop.table(temptab)
  res = rep(NA,6)
  res[1] = sum(temptab)
  res[2:5] = proptab[which(names(proptab)%in% c(1,2,3,5))]
  res[6] = prop.table(table(tempdat$onarvcf))[2]
  res[2:6] = res[2:6] * 100
  res = round(res,1)
  return(res)
}

summtab[1,3:8] = fillsummtab(1)
summtab[2,3:8] = fillsummtab(2)
summtab
#######################################################################
# Observed Transition Probabilities from t1 to t2

statenames = c("Engaged","Disengaged","Xfer out","Death")
probtab = matrix(0,4,4)
colnames(probtab) = rownames(probtab) = statenames

tempdat = dat.sub[dat.sub$gridtime==2,]
fill = table(tempdat$state_1,tempdat$state) #state_1: {1 engage, 0 disengage}

probtab[1,] = fill[2,]/sum(fill[2,])
probtab[2,] = fill[1,]/sum(fill[1,])
probtab[3,3] = probtab[4,4] = 1
round(probtab,2)
#######################################################################
# Run MPBART and BART for preliminary results?
# Reopen RStudio

devtools::install_github("yizhenxu/myGcomp2", args="--no-multiarch")
library(myGcomp2)

#
#
#
#####RELOAD DATA#####
#
#
#
tempdat = dat.sub[dat.sub$gridtime==2,]
ttt=table(tempdat$state[tempdat$onarvcf==1])
ttt=table(tempdat$state[tempdat$onarvcf_1==1])
prop.table(ttt)

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


library(myGcomp)
library(BayesTree)
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

m = 0 # 0 when using all data, m is the number of individuals used for analysis

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
summary(dat.sub$cd4[dat.sub$gridtime==0])
sum(is.na(dat.sub$cd4[dat.sub$gridtime==0])) #0
dat1 = dat.sub[dat.sub$gridtime==1,]
table(dat1$state) # no state==3
#######################################################################
#######################################################################
#######################################################################
#######################################################################
ttt = proc.time()
# settings and prior parameters for state predictions
burn = 100;ndraws = 500 ;ntree = 200;beta = 2

# set up list of conditional models and fit the models to data
Scd4 <- list()
Cd4 <- list()
Dcd4 <- list()

# use bart from BayesTree for the first CD4 conditional posterior draw
### X1 | V, X0, A0
V = dat.sub[dat.sub$gridtime == 0,c("male","ageatenroll")] #baseline gender age
X0 = dat.sub$cd4[dat.sub$gridtime==0]
A0 = dat.sub[dat.sub$gridtime == 0,]$onarvcf # on ARV status at enrollment
X1 = dat.sub$cd4[dat.sub$gridtime==1]
ind = !is.na(X1)
traind = cbind(A0,V,X0)
traind = traind[ind,]

a0 = 1
testd = cbind(a0,V,X0) ; colnames(testd) = colnames(traind)

#set.seed(1)
set.seed(100)
Cd4[[1]] = bart(x.train = traind , y.train = X1, x.test = testd,
                power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)

mu.x1.p = c(Cd4[[1]]$yhat.test)
sig.x1.p = rep(Cd4[[1]]$sigma,ncol(Cd4[[1]]$yhat.test))
x1.p = rnorm(n = length(mu.x1.p), mean = mu.x1.p, sd = sig.x1.p)

x1.p[x1.p<0] = 0 # cd4 has to be non-negative
x1.p = matrix(x1.p, nrow = ndraws, ncol = nrow(testd))

#######################################
post1 = apply(x1.p,2,mean) # cd4 posterior mean for each person
hist(post1)
hist(X1) # observed cd4 for each person

#######################################
# use bart1 from myGcomp

# testd is for everyone who survived at the previous grid. We only look at grid 1 now, and
# the previous grid is enrollment, hence testd contains everyone.

### Simulate death status at grid 1
### D1 | A0, X1, V

# make test data for D1
tmp3 = cbind(1:nrow(V),V) ; colnames(tmp3)[1] = "ID"
tmp3 = as.data.table(tmp3)
tmp3 = tmp3[rep(tmp3$ID, ndraws),]
tmp3 = setDF(tmp3)[,2:ncol(tmp3)]

testd = cbind(a0, c(t(x1.p)),tmp3)

# make train data for D1
D1 = (dat.sub[dat.sub$gridtime == 1,]$state==5)*1 # death indicator at 200 days

ind = which(dat.sub[dat.sub$gridtime == 1,]$patient_id %in%dat.sub[dat.sub$gridtime == 0,]$patient_id)

V = dat.sub[dat.sub$gridtime == 1,c("male","ageatenroll")]
A0 = dat.sub$onarvcf_1[dat.sub$gridtime == 1]
X1 = dat.sub$cd4[dat.sub$gridtime==1]
traind = cbind(A0, X1, V)
colnames(testd) = colnames(traind)

Dcd4[[1]]= bart1(x.train = traind , y.train = D1, x.test = testd,
                 power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)

d1 = pnorm(Dcd4[[1]]$yhat.test)
tt = c(d1)
set.seed(99)
d1.p = matrix(rbinom(length(tt),1,tt),nrow=ndraws)

#mean(d1.p) # 0.031683 overall mortality does decrease under a0 = 1
#hist(apply(d1.p,1,mean)) # posterior distribution of overall death rate
#abline(v = mean(D1)) # observed death rate
#######################################
# use bart1 from myGcomp

### Simulate engagement status at grid 1, among those death = 0
### S1 | A0, X1, V, D1 = 0

# make test data for S1, exclude D1 == 1
# testd remain unchanged, but TAKE OUT D1==1 AT THE VERY END OF THE SIMULATION!!! 
# will eventually have to change the myGcomp package and take out death on the way
# with test simulated death as an additional input to the bart1
#####################################

# make train data for S1, exclude state == 5 (death) from model training
dat1 = dat.sub[dat.sub$gridtime==1,]
ind.discard = which(dat1$state==5) # death indicator at 200 days

S1 = (dat1$state==1)*1 # engagement indicator at 200 days
S1 = S1[-ind.discard]
traind = traind[-ind.discard,]


Scd4[[1]]= bart1(x.train = traind , y.train = S1, x.test = testd,
                 power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)

s1 = pnorm(Scd4[[1]]$yhat.test)
s1.pp = matrix(s1,nrow=ndraws) # posterior Phat(S1=1)
tt = c(s1)
set.seed(8)
s1.p = matrix(rbinom(length(tt),1,tt),nrow=ndraws) # posterior draws of s1
# s1.pp is used for plotting posterior Phat, s1.p is used for the next simulation step

###############################################
### PLOT FOR THOSE SIMULATED ALIVE
#s1.p.d0 = s1.pp
#s1.p.d0[d1.p==1] = NA
#mean(s1.p,na.rm=T)  #0.90059
#mean(s1.p.d0,na.rm=T) #0.9249853

### posterior distribution of engagement rate AMONG ALIVE SUBJECTS!!!
#hist(apply(s1.p.d0,1,function(x) mean(x,na.rm=T))) 
#abline(v = mean(S1)) # observed engagement rate, AMONG ALIVE SUBJECTS!!! NOT OVERALL ENGAGEMENT
#mean(S1[traind$A0==1]) #0.8966667
#mean(S1[traind$A0==0]) #0.9070322

#######################################################
#dat1 = dat.sub[dat.sub$gridtime == 1,]
#S1.all = (dat1$state==1)*1
#D1.all = (dat1$state==5)*1
#ind = which(dat.sub[dat.sub$gridtime == 1,]$patient_id %in%dat.sub[dat.sub$gridtime == 0,]$patient_id)
#A0.all = dat.sub[dat.sub$gridtime == 0,]$onarvcf
#A0.all = A0.all[ind]
#mean(S1.all)
#mean(S1.all[A0.all==1])#0.8747967
#mean(S1.all[A0.all==0])#0.8677309
#table(A0.all,D1.all)
#################
# When death is included, engagenment among treated is higher then that when death is excluded
# this make sense, because people who were about to die visit clinic more eagerly and mostly likely on treatment
# the lower engagement in treated to untreated shows that, if a person is not about to die, 
# then generally treated engagement is worse than untreated.
#######################################################


### X2 | V, X1, A1, S1
V = dat.sub[dat.sub$gridtime == 2,c("male","ageatenroll")] #baseline gender age
X1 = dat.sub$cd4_1[dat.sub$gridtime==2]
A1 = dat.sub$onarvcf_1[dat.sub$gridtime == 2] # on ARV status at enrollment
S1 = dat.sub$state_1[dat.sub$gridtime==2]

X2 = dat.sub$cd4[dat.sub$gridtime==2]
traind = cbind(A1,V,X1,S1)

a1 = 1

# make test data (a1, V, x1.p, s1.p) 
V0 = dat.sub[dat.sub$gridtime == 0,c("male","ageatenroll")]
tmp = cbind(1:nrow(V0),V0) ; colnames(tmp)[1] = "ID"
tmp = as.data.table(tmp)
tmp = tmp[rep(tmp$ID, ndraws),]
tmp = setDF(tmp)[,2:ncol(tmp)]

testd = cbind(a1, tmp, c(t(x1.p)),c(t(s1.p)))
 colnames(testd) = colnames(traind)

Cd4[[2]] = bart1(x.train = traind , y.train = X2, x.test = testd,
                power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)

mu.x2.p = c(Cd4[[2]]$yhat.test)
sig.x2.p = rep(Cd4[[2]]$sigma,ncol(Cd4[[2]]$yhat.test))
set.seed(8)
x2.p = rnorm(n = length(mu.x2.p), mean = mu.x2.p, sd = sig.x2.p)

x2.p[x2.p<0] = 0 # cd4 has to be non-negative
x2.p = matrix(x2.p, nrow = ndraws)


#post2 = apply(x2.p,2,mean) # cd4 posterior mean for each person

###############################################
### PLOT FOR THOSE SIMULATED ALIVE
hist(post2)
hist(X2) # observed cd4 for each person
x2.p.d0 = x2.p
x2.p.d0[d1.p==1] = NA
mean(x2.p)  #357.3257
mean(x2.p.d0,na.rm=T) # 363
mean(X2) # 342.1735
### posterior distribution of engagement rate AMONG ALIVE SUBJECTS!!!
hist(apply(x2.p.d0,1,function(x) mean(x,na.rm=T))) 
abline(v = mean(X2)) # observed engagement rate, AMONG ALIVE SUBJECTS!!! NOT OVERALL ENGAGEMENT
mean(X2[traind$A1==1]) #407.897
mean(X2[traind$A1==0]) #294.6435

pid0 = dat.sub$patient_id[dat.sub$gridtime==0 & dat.sub$onarvcf==1]
pid1 = dat.sub$patient_id[dat.sub$gridtime==1 & dat.sub$onarvcf==1]
a0a1 = pid1[which(pid1 %in% pid0)]
pid2 = dat.sub$patient_id[dat.sub$gridtime==2]
pid2a0a1 = pid2[which(pid2 %in% a0a1)]
X2a0a1 = dat.sub$cd4[dat.sub$gridtime==2 & dat.sub$patient_id %in% pid2a0a1]
mean(X2a0a1) # mean of X2 for those a0==1 and a1==1
# 412.2759
###############################################



### S2 | V, X2, A1, S1
V = dat.sub[dat.sub$gridtime == 2,c("male","ageatenroll")] #baseline gender age
A1 = dat.sub$onarvcf_1[dat.sub$gridtime == 2] # on ARV status at enrollment
S1 = dat.sub$state_1[dat.sub$gridtime==2]
X2 = dat.sub$cd4[dat.sub$gridtime==2]

S2 = (dat.sub$state[dat.sub$gridtime==2]==1)*1
traind = cbind(A1,V,X2,S1)

a1 = 1

# make test data (a1, V, x2.p, s1.p) 
V0 = dat.sub[dat.sub$gridtime == 0,c("male","ageatenroll")]
tmp = cbind(1:nrow(V0),V0) ; colnames(tmp)[1] = "ID"
tmp = as.data.table(tmp)
tmp = tmp[rep(tmp$ID, ndraws),]
tmp = setDF(tmp)[,2:ncol(tmp)]

testd = cbind(a1, tmp, c(t(x2.p)),c(t(s1.p)))
colnames(testd) = colnames(traind)

set.seed(100)
Scd4[[2]] = bart1(x.train = traind , y.train = S2, x.test = testd,
                power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)


s2 = pnorm(Scd4[[2]]$yhat.test)
s2.pp = matrix(s2,nrow=ndraws) # posterior Phat(S1=1)
tt = c(s2)
set.seed(8)
s2.p = matrix(rbinom(length(tt),1,tt),nrow=ndraws) # posterior draws of s1
# s1.pp is used for plotting posterior Phat, s1.p is used for the next simulation step
#################################################
### PLOT FOR THOSE SIMULATED ALIVE
s2.p.d0 = s2.pp
s2.p.d0[d1.p==1] = NA
mean(s2.pp,na.rm=T)  # 0.7318083
mean(s2.p.d0,na.rm=T) #0.7341732

posts2 = apply(s2.p.d0,1,function(x) mean(x,na.rm=T))
mean(posts2) # 0.7341433
summary(posts2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.6825  0.7265  0.7358  0.7341  0.7440  0.7644

mean(S2) # 0.668263

### posterior distribution of engagement rate AMONG ALIVE SUBJECTS!!!
hist(posts2) 
abline(v = mean(S2)) # observed engagement rate, AMONG ALIVE SUBJECTS!!! NOT OVERALL ENGAGEMENT
mean(S2[A1==1]) # 0.7749437
mean(S2[A1==0]) # 0.5207469

pid0 = dat.sub$patient_id[dat.sub$gridtime==0 & dat.sub$onarvcf==1]
pid1 = dat.sub$patient_id[dat.sub$gridtime==1 & dat.sub$onarvcf==1]
a0a1 = pid1[which(pid1 %in% pid0)]
pid2 = dat.sub$patient_id[dat.sub$gridtime==2]
pid2a0a1 = pid2[which(pid2 %in% a0a1)]
S2a0a1 = 1*(dat.sub$state[dat.sub$gridtime==2 & dat.sub$patient_id %in% pid2a0a1]==1)
S2a0a1 = (dat.sub$state[dat.sub$gridtime==2]==1)*1
S2 = (dat.sub$state[dat.sub$gridtime==2]==1)*1
mean(S2a0a1,na.rm=T) # mean of X2 for those a0==1 and a1==1
# 0.6642729


#################################################
proc.time() - ttt

save.image("2017-04-28_Gcomp_s2_a=1_500draws.RData")

##
################################################
# do in new window
load("I:\\2017-04-29_Gcomp_s2_a=01_200draws_1.RData")

s2.p.d0 = s2.pp
s2.p.d0[d1.p==1] = NA
posts2 = apply(s2.p.d0,1,function(x) mean(x,na.rm=T))
hist(posts2,breaks = 15,xlim = c(0.66,max(posts2)))
abline(v = mean(S2)) 
mean(posts2) #  0.7330075
summary(posts2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.7169  0.7297  0.7331  0.7333  0.7369  0.7492 

mean(S2) # 0.668263


mean(s2.pp,na.rm=T)  # 0.7311191
mean(s2.p.d0,na.rm=T) #0.73332


posts2.01 = posts2
save(posts2.01, file = "I:\\2017-04-29_Gcomp_s2_a=01_200draws_1_posts2.RData")
##################################################
load("I:\\2017-04-28_Gcomp_s2_a=1_500draws.RData")
load("I:\\2017-04-29_Gcomp_s2_a=01_200draws_1_posts2.RData") # for posts2.01
# histogram show D1 is lower than observed for a0 = 1
mean(d1.p) # 0.03689478 overall mortality does decrease under a0 = 1
mean(D1) #0.0403057
hist(apply(d1.p,1,mean)) # posterior distribution of overall death rate
abline(v = mean(D1)) # observed death rate
postd1 = apply(d1.p,1,mean)


# histogram show S2 is higher than observed for (a0, a1) = (1, 1)
s2.p.d0 = s2.pp
s2.p.d0[d1.p==1] = NA
posts2 = apply(s2.p.d0,1,function(x) mean(x,na.rm=T))
hist(posts2,breaks = 15,xlim = c(0.66,max(posts2)))
abline(v = mean(S2)) 
mean(posts2) #  0.7333174
summary(posts2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.7169  0.7297  0.7331  0.7333  0.7369  0.7492 

mean(S2) # 0.668263


mean(s2.pp,na.rm=T)  # 0.7311191
mean(s2.p.d0,na.rm=T) #0.73332


out = matrix(NA,4,6)
colnames(out) = c("Perc_onARV","Perc_Eng","Perc_Eng|onARV1","Perc_Eng|onARV0","Perc_Eng|prevonARV1","Perc_Eng|prevonARV0")
#colnames(out) = c("Perc_onARV","Perc_Eng","Perc_Eng|onARV1","Perc_Eng|onARV0")
rownames(out) = c(0,200,400,600)
for(k in 0:3){
  dat1 = dat.sub[dat.sub$gridtime==k ,]
  out[k+1,1] = mean(dat1$onarvcf)
  out[k+1,2] = mean(dat1$state==1)
  out[k+1,3] = mean(dat1[dat1$onarvcf==1,]$state==1)
  out[k+1,4] = mean(dat1[dat1$onarvcf==0,]$state==1)
  out[k+1,5] = mean(dat1[dat1$onarvcf_1==1,]$state==1)
  out[k+1,6] = mean(dat1[dat1$onarvcf_1==0,]$state==1)
}
out



post = data.frame(Death_Day200 = postd1, Engage_Day400 = posts2)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
p1 = ggplot(data=post, aes(Death_Day200)) + 
  geom_histogram(breaks=seq(min(postd1), max(postd1), by =0.0008), 
                 col="red", 
                 aes(y = (..count..)/sum(..count..),fill=..count..)) +
  scale_fill_gradient("Count", low = "green", high = "red")+ 
  scale_y_continuous(labels=percent)+ theme(legend.position="none") +
  xlab('Posterior Death Rate (1 ~ 200 Days)') +
  ylab('Percentage')+ geom_vline(xintercept = mean(D1))+ 
  annotate("text", size = 7,x = 0.03, y = 0.1, label = c("a0 = 1"))+
  theme(text = element_text(size=18))

mean(postd1>=mean(D1))

p2 = ggplot(data=post, aes(Engage_Day400)) + 
  geom_histogram(breaks=seq(min(posts2), max(posts2), by =0.002), 
                 aes(y = (..count..)/sum(..count..),fill=..count..)) +
  scale_fill_gradient("Count", low = "green", high = "red")+ 
  scale_y_continuous(labels=percent)+ theme(legend.position="none") +
  xlab('Posterior Engagement Rate (201 ~ 400 Days)') +
  ylab('Percentage') + geom_vline(xintercept = mean(S2))+ 
  annotate("text", size =7,x = 0.7, y = 0.12, label = c('a[1] == 1'),parse=TRUE)+
  theme(text = element_text(size=18))


p3 = ggplot(data=data.frame(posts2.01), aes(posts2.01)) + 
  geom_histogram(breaks=seq(min(posts2.01), max(posts2.01), by =0.002), 
                 aes(y = (..count..)/sum(..count..),fill=..count..)) +
  scale_fill_gradient("Count", low = "green", high = "red")+ 
  scale_y_continuous(labels=percent)+ theme(legend.position="none") +
  xlab('Posterior Engagement Rate (201 ~ 400 Days)') +
  ylab('Percentage') + geom_vline(xintercept = mean(S2))+ 
  annotate("text", size = 7,x = 0.72, y = 0.2, label = c("(a0,a1) = (0, 1)"))+
  theme(text = element_text(size=18))

png("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Poster and Presentation\\2017Presentation\\S2a11.png" ,width = 550, height = 400)
plot(p2)
dev.off()

png("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Poster and Presentation\\2017Presentation\\S2a01.png" ,width = 550, height = 400)
plot(p3)
dev.off()

png("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Poster and Presentation\\2017Presentation\\D1a1.png" ,width = 550, height = 400)
plot(p1)
dev.off()


carrots = data.frame(Probability = posts2[1:200])
cukes = data.frame(Probability = posts2.01)
carrots$Regime <- 'a0 = 1 a1 = 1'
cukes$Regime <- 'a0 = 0 a1 = 1'

#and combine into your new data frame vegLengths
vegLengths <- rbind(carrots, cukes)

#aa = mean(S1[A==0])
#bb = mean(S1[A==1])

library(ggplot2)
#now make your lovely plot
ggplot(vegLengths, aes(Probability, fill = Regime)) + 
  geom_histogram(alpha = 0.6,binwidth = 0.0015, bins=30,aes(y = ..density..), position = 'identity')+
  #geom_density(alpha=0.5)  +
  #scale_y_continuous(breaks=seq(0.837,0.877,0.01))+
  geom_vline(xintercept = mean((dat.sub[dat.sub$gridtime == 2,]$state==1)*1 ))

ggplot(vegLengths, aes(Probability, fill = Regime)) + 
  geom_histogram(alpha = 0.6,binwidth = 0.0015, bins=30, position = 'identity', 
                 aes(y = (..count..)/sum(..count..)))+
  #geom_density(alpha=0.5)  +
  #scale_y_continuous(breaks=seq(0.837,0.877,0.01))+
  geom_vline(xintercept = mean((dat.sub[dat.sub$gridtime == 1,]$state==1)*1 ))



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


# open project myGcomp, and load myGcomp locally
# devtools::install_github("yizhenxu/myGcomp")
library(myGcomp)
#devtools::load_all()
#devtools::document()
#build and reload in build panel

library("data.table")

options(scipen = 10)
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

n = length(unique(dat.sub$idx))
n


#######################################################################
summary(dat.sub$cd4[dat.sub$gridtime==0])
sum(is.na(dat.sub$cd4[dat.sub$gridtime==0])) #0
#######################################################################
#######################################################################
#######################################################################
#######################################################################
ttt = proc.time()
# settings and prior parameters for state predictions
burn = 100;ndraws = 500 ;ntree = 200;beta = 2

## set up list of conditional models and fit the models to data
Scd4 <- list()
Cd4 <- list()

library(BayesTree)
# use bart from BayesTree for the first CD4 conditional posterior draw

V = dat.sub[dat.sub$gridtime == 0,c("male","ageatenroll")] #baseline gender age
X0 = dat.sub$cd4[dat.sub$gridtime==0]
A = dat.sub[dat.sub$gridtime == 0,]$onarvcf # on ARV status at enrollment
X1 = dat.sub$cd4[dat.sub$gridtime==1]
ind = !is.na(X1)
a = 0
traind = cbind(A,V,X0)
traind = traind[ind,]
testd = cbind(a,V,X0) ; colnames(testd) = colnames(traind)

#set.seed(1)
set.seed(100)
Cd4[[1]] = bart(x.train = traind , y.train = X1, x.test = testd,
                power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)

#dim(Cd4[[1]]$yhat.test)
#summary(Cd4[[1]]$yhat.test)
#summary(X1)
#length(Cd4[[1]]$sigma)
#hist(Cd4[[1]]$sigma)

mu.x1.p = c(Cd4[[1]]$yhat.test)
sig.x1.p = rep(Cd4[[1]]$sigma,ncol(Cd4[[1]]$yhat.test))

#set.seed(10)
x1.p = rnorm(n = length(mu.x1.p), mean = mu.x1.p, sd = sig.x1.p)

x1.p[x1.p<0] = 0 # cd4 has to be non-negative
x1.p = matrix(x1.p, nrow = ndraws, ncol = nrow(testd))
#hist(x1.p)
#hist(X1)

#tmp1 = c(Cd4[[1]]$yhat.test) #c(ndraws of ID 1, ... , ndraws of ID n)
#Estar= rbinom(n*ndraws, 1, pnorm(tmp1))
#E1star = matrix(Estar, nrow = ndraws, ncol = n)
E1star = x1.p
post1 = apply(E1star,2,mean)
#hist(post1,breaks =15)
#hist(X1)
#mean(post1)
E1star = t(E1star)# E1*|A0=1,V    nrowTest x ndraws

#######################################
# use bart1 from myGcomp

# testd is for everyone who survived at the previous grid. We only look at grid 1 now, and
# the previous grid is enrollment, hence testd contains everyone.
tmp3 = cbind(1:nrow(V),V) ; colnames(tmp3)[1] = "ID"
tmp3 = as.data.table(tmp3)
tmp3 = tmp3[rep(tmp3$ID, ndraws),]
tmp3 = setDF(tmp3)[,2:ncol(tmp3)]

testd1 = cbind(c(E1star),tmp3)
testd = cbind(a, testd1)

ind = which(dat.sub[dat.sub$gridtime == 1,]$patient_id %in%dat.sub[dat.sub$gridtime == 0,]$patient_id)
V = V[ind,]

S1 = (dat.sub[dat.sub$gridtime == 1,]$state==1)*1 # engagement indicator at 200 days
A = dat.sub[dat.sub$gridtime == 0,]$onarvcf # on ARV status at enrollment
A = A[ind]
X1 = dat.sub$cd4[dat.sub$gridtime==1]

a = 0
traind = cbind(A,X1,V)
colnames(testd) = colnames(traind)

Scd4[[1]]= bart1(x.train = traind , y.train = S1, x.test = testd,
                 power = beta, ntree = ntree, ndpost = ndraws, nskip = burn)

post2.bart1 = apply(pnorm(Scd4[[1]]$yhat.test),1,mean)
#hist(post2.bart1)
#mean(post2.bart1)


save.image(file = "2017-04-26_Gcomp_s1_a=0_500draws.RData")


###########################################
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

proc.time() - ttt

###########################################
# run in another new R window 
save.image(file = "2017-04-26_Gcomp_s1_cd4le350_500draws.RData")

load("I:\\2017-04-26_Gcomp_s1_cd4le350_500draws.RData")

mean(post2.bart1)
est.a_cd4le350 = post2.bart1
hist(est.a_cd4le350,breaks=13)
summary(est.a_cd4le350)

save(est.a_cd4le350, file = "I:\\tmp1.RData")

load("I:\\2017-04-26_Gcomp_s1_a=1_500draws.RData")

mean(post2.bart1)
est.a1 = post2.bart1
hist(est.a1,breaks=13)
summary(est.a1)
save(est.a1, file = "I:\\tmp2.RData")

#########################################
load("I:\\2017-04-26_Gcomp_s1_a=0_500draws.RData")

hist(post2.bart1,breaks=13)
mean(post2.bart1)


tt = c(pnorm(Scd4[[1]]$yhat.test))
set.seed(99)
kk = apply(matrix(rbinom(length(tt),1,tt),nrow=ndraws),1,mean)
tmp = sort(kk)[ceiling(length(kk)*0.05)]
tmp

est.a0 = post2.bart1
hist(est.a0,breaks=13)
summary(est.a0)

load("I:\\tmp1.RData")
load("I:\\tmp2.RData")
hist(est.a_cd4le350)
hist(est.a1,breaks=13)

load("I:\\2017presentationPlots.RData")
carrots = data.frame(Probability = est.a_cd4le350)
cukes = data.frame(Probability = est.a0)
caddy = data.frame(Probability = est.a1)
#Now, combine your two dataframes into one.  First make a new column in each.
carrots$Regime <- 'a0 = 1{X0<=350}'
cukes$Regime <- 'a0 = 0'
caddy$Regime <- 'a0 = 1'
#and combine into your new data frame vegLengths
vegLengths <- rbind(carrots, cukes,caddy)

#aa = mean(S1[A==0])
#bb = mean(S1[A==1])

library(ggplot2)
#now make your lovely plot
p0 = ggplot(vegLengths, aes(Probability, fill = Regime)) + 
  #geom_histogram(alpha = 0.6,binwidth = 0.0015, bins=30,aes(y = ..density..), position = 'identity')+
  geom_density(alpha=0.5)  +
  #scale_y_continuous(breaks=seq(0.837,0.877,0.01))+
  geom_vline(xintercept = mean((dat.sub[dat.sub$gridtime == 1,]$state==1)*1 ))+
  xlab('Posterior Engagement Rate (1 ~ 200 Days)')+
  theme(text = element_text(size=18))

png("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Poster and Presentation\\2017Presentation\\S1compare.png" ,width = 550, height = 400)
plot(p0)
dev.off()
#save(est.a_cd4le350,est.a0,est.a1,file = "I:\\2017presentationPlots.RData")


p0a1 = ggplot(data = data.frame(est.a1), aes(est.a1)) + 
  #geom_histogram(alpha = 0.6,binwidth = 0.0015, bins=30,aes(y = ..density..), position = 'identity')+
  geom_density(alpha=0.5)  +
  #scale_y_continuous(breaks=seq(0.837,0.877,0.01))+
  geom_vline(xintercept = mean((dat.sub[dat.sub$gridtime == 1,]$state==1)*1 ))+
  xlab('Posterior Engagement Rate (1 ~ 200 Days)')+
  theme(text = element_text(size=18))

p0a1 = ggplot(data=data.frame(est.a1), aes(est.a1)) + 
  geom_histogram(breaks=seq(min(est.a1), max(est.a1), by =0.0015), 
                 aes(y = (..count..)/sum(..count..),fill=..count..)) +
  scale_fill_gradient("Count", low = "green", high = "red")+ 
  scale_y_continuous(labels=percent)+ theme(legend.position="none") +
  xlab('Posterior Engagement Rate (1 ~ 200 Days)') +
  ylab('Percentage') + geom_vline(xintercept = mean((dat.sub[dat.sub$gridtime == 1,]$state==1)*1 ))+ 
  #annotate("text", size = 7,x = 0.72, y = 0.2, label = c("(a0,a1) = (0, 1)"))+
  theme(text = element_text(size=18))

png("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Poster and Presentation\\2017Presentation\\S1a1.png" ,width = 550, height = 400)
plot(p0a1)
dev.off()


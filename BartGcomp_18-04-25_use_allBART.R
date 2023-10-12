#devtools::install_github("yizhenxu/allBART")
#remove.packages("allBART")

library(allBART)

library(lubridate)# date processing
library(data.table)
library(zoo)# carry backward forward
library(foreign)#stata read in
library(nnet)# multinomial log-linear models
library(Hmisc)#rmultinom function

library("parallel")

setwd("./Topic2");source("Functions_new.R")
#setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
######################################################################################################################

theA = "onarv_cf"

theY = "outcome"

theV = c("ageinit",
         "male_1",
         "yrenrol_2009","yrenrol_2010","yrenrol_2011","yrenrol_2012","yrenrol_2013","yrenrol_2014", "yrenrol_2015", "yrenrol_2016",
         "TravelTime_below30min","TravelTime_over2h","TravelTime_1to2h","TravelTime_30to60min",
         "whostageinit_2","whostageinit_3","whostageinit_4","whostageinit_NA",
         "marriedinit_1","marriedinit_NA",
         "heightinit_RX","heightinit_R")

theX = c("new_cupd_1","log_cd4plus1") # LHS

theX_cov = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R") # RHS (from createind and createmis)

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

theA = "onarv_cf_1"

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

######################################################################################################################

load("newdat.RData") # defined in 18-02-26 script
#load("./makeGridData/newdat.RData")

newdat$outcome = as.numeric(newdat$outcome) #IMPORTANT! de-factorize outcome

screennum =95590
#88574 88587 88755 88773 95447 95466 95590
source(paste0("script_",screennum,".R"),echo=T)
######################################################################################################################
### ZEUS script::

# Dyn simulation
A = 0 #theA_lag1

### optional: pick 10000 people
set.seed(1)
nsub = sum(newdat$time==0)
id = sample(1:nsub, 10000, replace = F)
newdat = newdat[myID %in% id,]
###

sink(paste0("output_",screennum,".txt")) #write all output to txt -- source('script.R',echo=T) to write input too

##############################################################

### Model building

nsub = sum(newdat$time==0)

#-------------------------------
setting = 1
maxt= 11
ndraw = 1000
nburn = 100
ntree = 100
ncut = 100
nsd = 50 # max repeated draw of Sigma for Mpbart
seed = 88
basey = 2 # the outcome level used in RHS of model fitting (theY_ind)
nmod = 33 # number of models

#-------------------------------

# initialization

set.seed(seed)

runtime = rep(NA, nmod + 1)

modfml =  vector("list", nmod)
mod = vector("list", nmod)

Smat = vector("list", maxt)
Scov = vector("list", maxt)
Xmat = vector("list", maxt)

train_pred = vector("list", nmod)
train_y = vector("list", nmod)

for(i in 1:maxt){
  Xmat[[i]] = array(NA, dim = c(nsub, ndraw, length(theX)))
}

Vdist = cbind(newdat[time==0,theV,with=F])
Vdist = Vdist[rep(1:nsub,ndraw),]

##########################

t = 1
fitdat = newdat[time==200*t,]

##########################
ct = 1

# [L1 | V, A0, X0, L0]

RHS = c(theV, theA_lag1, theX_cov_lag1)
# exclude lag X_R from models at t=1 because that is the same as lag L at baseline
RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
RHS = fmlfilter(fitdat,RHS)
fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
modfml[[ct]] = fml

policy_lag = rep(A, nsub)
dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)

runtime[1] = proc.time()[3]
mod[[ct]] = PBART_call(as.formula(fml), data = fitdat, test.data = dat.pred,
                   Prior = list(ntrees = ntree,
                                kfac = 2,
                                pswap = 0.1, pbd = 0.5, pb = 0.25,
                                alpha = 0.95, beta = 2.0,
                                nc = ncut, minobsnode = 10),
                   Mcmc = list(burn=nburn, ndraws = ndraw))
Xmat[[t]][,,1] = mod[[ct]]$samp_test

train_y[[ct]] = setDF(fitdat[, theX[1], with = F])
train_pred[[ct]] = mod[[ct]]$samp_train
##########################
ct = ct+1

# [X1 | V, A0, X0, L0, L1 = 1] # simulated Xt remains NA until a first Lt == 1

RHS = c(theV, theA_lag1, theX_cov_lag1)
# exclude lag X_R from models at t=1 because that is the same as lag L at baseline
RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
modfml[[ct]] = fml

policy_lag = rep(A, nsub)
dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)

runtime[2] = proc.time()[3]

mod[[ct]] = BART_call(as.formula(fml), data = fitdat[get(theX[1])==1,], test.data = dat.pred,
                 Prior = list(nu = 3, sigq = 0.9,
                              ntrees = ntree,
                              kfac = 2,
                              pswap = 0.1, pbd = 0.5, pb = 0.25,
                              alpha = 0.95, beta = 2.0,
                              nc = ncut, minobsnode = 10),
                 Mcmc = list(burn = nburn, ndraws = ndraw))

postY =  mod[[ct]]$samp_test
mis_ind = which(Xmat[[t]][,,1] == 0) # people who have no CD4 update
postY[mis_ind] = NA
Xmat[[t]][,,2] = postY
remove(postY)

train_y[[ct]] = setDF(fitdat[get(theX[1])==1, theX[2], with = F])
train_pred[[ct]] = mod[[ct]]$samp_train
##########################
ct = ct+1

# [S1 | V, A0, X1, L1]

ytab = table(fitdat$outcome)
ylevels = as.numeric(names(ytab)[ytab > 0])
p = length(ylevels) # number of outcome categories

RHS = c(theV,theA_lag1,theX_cov)
RHS = fmlfilter(fitdat,RHS)
fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
modfml[[ct]] = fml

policy_lag = rep(A, nrow(Vdist))
dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)

runtime[3] = proc.time()[3]

mod[[ct]] = DynMPBART_call(as.formula(fml), data = fitdat, test.data = dat.pred,
                        base = basey,
                     Prior = list(nu = p-1+3, V = diag(p-1),
                                  ntrees = ntree,
                                  kfac = 2,
                                  pswap = 0.1, pbd = 0.5, pb = 0.25,
                                  alpha = 0.95, beta = 2.0,
                                  nc = ncut, minobsnode = 10),
                     Mcmc = list(sigma0 = diag(p-1), burn = nburn, ndraws = ndraw,
                                 nSigDr = nsd, keep_sigma_draws=T))

Smat[[t]] = mod[[ct]]$samp_test
Scov[[t]] = mod[[ct]]$sigmasample
train_y[[ct]] = setDF(fitdat[, theY, with = F])
train_pred[[ct]] = mod[[ct]]$samp_train

print(paste0("t=",t))
##################################################################################
##################################################################################
for(t in 2:maxt){
  fitdat = newdat[time==200*t,]

  ##########################
  ct = ct + 1

  # [Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]

  RHS = c(theV, theA_lag1, theX_cov_lag1,theY_ind_lag1)
  RHS = fmlfilter(fitdat,RHS)
  fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
  modfml[[ct]] = fml

  policy_lag = rep(A, nrow(Vdist))
  dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)

  runtime[ct] = proc.time()[3]
  mod[[ct]] = DynPBART_call(as.formula(fml), data = fitdat, test.data = dat.pred,
                         Prior = list(ntrees = ntree,
                                      kfac = 2,
                                      pswap = 0.1, pbd = 0.5, pb = 0.25,
                                      alpha = 0.95, beta = 2.0,
                                      nc = ncut, minobsnode = 10),
                         Mcmc = list(burn=nburn, ndraws = ndraw))
  Xmat[[t]][,,1] = mod[[ct]]$samp_test
  train_y[[ct]] = setDF(fitdat[, theX[1], with = F])
  train_pred[[ct]] = mod[[ct]]$samp_train
  ##########################
  ct = ct+1

  # [Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}, Lt = 1]

  RHS = c(theV, theA_lag1, theX_cov_lag1, theY_ind_lag1)
  RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
  fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
  modfml[[ct]] = fml

  policy_lag = rep(A, nrow(Vdist))
  dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)

  runtime[ct] = proc.time()[3]

  mod[[ct]] = DynBART_call(as.formula(fml), data = fitdat[get(theX[1])==1,], test.data = dat.pred,
                        Prior = list(nu = 3, sigq = 0.9,
                                     ntrees = ntree,
                                     kfac = 2,
                                     pswap = 0.1, pbd = 0.5, pb = 0.25,
                                     alpha = 0.95, beta = 2.0,
                                     nc = ncut, minobsnode = 10),
                        Mcmc = list(burn = nburn, ndraws = ndraw))

  postY =  mod[[ct]]$samp_test
  L = Xmat[[t]][,,1]
  X_prev = Xmat[[t-1]][,,2]
  mis_ind = which(L == 0) # people who have no CD4 update keeps the old value
  postY[mis_ind] = X_prev[mis_ind]
  Xmat[[t]][,,2] = postY
  # to save some space for R...
  remove(postY)
  remove(L)
  remove(X_prev)
  remove(mis_ind)

  train_y[[ct]] = setDF(fitdat[get(theX[1])==1, theX[2], with = F])
  train_pred[[ct]] = mod[[ct]]$samp_train
  ##########################
  ct = ct+1

  # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]

  ytab = table(fitdat$outcome)
  ylevels = as.numeric(names(ytab)[ytab > 0])
  p = length(ylevels) # number of outcome categories

  RHS = c(theV, theA_lag1, theX_cov, theY_ind_lag1)
  RHS = fmlfilter(fitdat,RHS)
  fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
  modfml[[ct]] = fml

  policy_lag = rep(A, nrow(Vdist))
  dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)

  runtime[ct] = proc.time()[3]

  if(p>2){
    mod[[ct]] = DynMPBART_call(as.formula(fml), data = fitdat, test.data = dat.pred,
                               base = basey,
                               Prior = list(nu = p-1+3, V = diag(p-1),
                                            ntrees = ntree,
                                            kfac = 2,
                                            pswap = 0.1, pbd = 0.5, pb = 0.25,
                                            alpha = 0.95, beta = 2.0,
                                            nc = ncut, minobsnode = 10),
                               Mcmc = list(sigma0 = diag(p-1), burn = nburn, ndraws = ndraw,
                                           nSigDr = nsd, keep_sigma_draws=T))

    Smat[[t]] = mod[[ct]]$samp_test
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]

  } else if(p==2){

    Y = rep(0,nrow(fitdat))
    Y[fitdat$outcome == ylevels[2]] = 1
    fitdat$outcome = Y
    mod[[ct]] = DynPBART_call(as.formula(fml), data = fitdat, test.data = dat.pred,
                              Prior = list(ntrees = ntree,
                                           kfac = 2,
                                           pswap = 0.1, pbd = 0.5, pb = 0.25,
                                           alpha = 0.95, beta = 2.0,
                                           nc = ncut, minobsnode = 10),
                              Mcmc = list(burn=nburn, ndraws = ndraw))

    postY = mod[[ct]]$samp_test
    postY[which(postY==0)] = ylevels[1]
    postY[which(postY==1)] = ylevels[2]
    Smat[[t]] = postY
    remove(postY)
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
    
  }
  Scov[[t]] = mod[[ct]]$sigmasample
  train_y[[ct]] = setDF(fitdat[, theY, with = F])
  train_pred[[ct]] = mod[[ct]]$samp_train
  
  print(paste0("t=",t))
}# end for t

##################################################################################
##################################################################################
runtime[ct+1] = proc.time()[3]

a = runtime[-1] - runtime[-nmod] #runtime for each model
a = a[!is.na(a)]
status = as.character(paste0(" | ", a, collapse = ""))


library(mailR)
sender <- "916800617qq@gmail.com"
recipients <- c("yizhen_xu@brown.edu")
send.mail(from = sender,
          to = recipients,
          subject = paste0(screennum," complete"),
          body = status,
          smtp = list(host.name = "smtp.gmail.com", port = 465,
                      user.name = "916800617qq@gmail.com",
                      passwd = "Xyz123456!", ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

sink()

# no Scov and train_pred train_y for 10000 ppl and all
save(newdat, modfml, Xmat,Smat,Scov,train_y, train_pred, file = paste0("Bart_screen",screennum,"_ct",ct,"_nsub",nsub,".RData"))

##################################################################################

dir = "C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS CODE/"
fn = "Bart_screen88574_ct33_nsub10000.RData"
A = 0
fn = "Bart_screen88587_ct33_nsub10000.RData"
A = 1
load(paste0(dir,fn))


#rm(list = ls())
load("Bart_screen88574_ct33_nsub20000.RData")
load("Bart_screen88587_ct33_nsub20000.RData")
load("Bart_screen95590_ct33_nsub20000.RData")
# save all the Smat's: Smat350 Smat1 Smat0 into tmp.RData 
##################################################################################
maxt = 11

tmp = unlist(Smat)
#tmp = Smat
tmp1 = c(rep(2,length(Smat[[1]])),tmp)
a = data.frame(time = rep(200*(0:maxt),each = length(Smat[[1]])), outcome = tmp1)
#counts <- table(a$outcome,a$time)

longdat = setDT(a)
ldc = longdat[, .N, by = c("time","outcome")]
ldc = rbind(ldc,list(time = c(0,0,0,200),outcome = c(1,3,4,3),N=c(0,0,0,0) ))
#ldc = rbind(ldc,list(time = c(0,0,0),outcome = c(1,3,4),N=c(0,0,0) ))
ldc = ldc[order(time,outcome)]
counts = matrix(ldc$N,nrow = 4)

ctmp = counts[1,]
counts[1,] = counts[2,]
counts[2,] = ctmp
#counts[1,1] = counts[2,1];counts[2,1]=0

rownames(counts) = c("engaged","diengaged","transferred out","died")
counts = simplify2array(lapply(1:ncol(counts), function(x) counts[,x]/sum(counts[,x])))
colnames(counts) = 200*(0:maxt)

save(counts, file = "tmp1.RData")
load("I:\\Topic2\\tmp1.RData")

A = "CD4ls350"
png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\BARTCD4SimEngA",A,".png"),width = 600, height = 400)
png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\BARTCD4SimEngA",A,".png"),width = 600, height = 400)

png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\MNPCD4SimEngA",A,".png"),width = 600, height = 400)
barplot(counts, main=paste0("Bayesian Simulated Engagement Status at a=",A,"\n Time-varying Covariates: CD4 info"),
        xlab="Days since baseline", ylab="Number of observations", col=1:4)
legend("bottomleft", 
       legend = rownames(counts), 
       fill = 1:4)
dev.off()



##################################################################################
# save all the Smat's: Smat350 Smat1 Smat0 into tmp.RData 
#extract their [[5]] (the 1000's day) overall engagement posteriors into tmp1.RData
load("I:\\Topic2\\tmp.RData")
a = Smat350[[5]]
b = Smat1[[5]]
c = Smat0[[5]]

#A = 0
h1 = apply(1*(c==2),2,mean)
#A = 1
h2 = apply(1*(b==2),2,mean)
#A = 1{cd4<350}
h3 = apply(1*(a==2),2,mean)

hd = rbind(cbind(h1, "a=0"), cbind(h2, "a=1"), cbind(h3 ,"a=1(cd4<350)"))
colnames(hd) = c("Posterior Engagement Rate", "Policy")
hd = as.data.frame(hd)
save(hd, file = "tmp1.RData")
load("I:\\Topic2\\tmp1.RData")

hd$EngRate = as.numeric(as.character(hd[,1]))

#library(devtools)
#install_github("easyGgplot2", "kassambara")
library(easyGgplot2)
#http://www.sthda.com/english/wiki/ggplot2-histogram-easy-histogram-graph-with-ggplot2-r-package


# Histogram plots with mean lines
plot = ggplot2.histogram(data=hd, xName='EngRate',
                  groupName='Policy', legendPosition="top",
                  alpha=0.5, addDensity=TRUE,
                  addMeanLine=TRUE, meanLineColor="white", meanLineSize=1.5)
print(plot)
ggplot2.customize(plot, mainTitle="Posterior Distribution of Engagement Rate at Day 1000",
                  xtitle="Overall Engagement Rate", ytitle="Histogram")

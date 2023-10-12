# AGENDA 12-13 (SHARED GOOGLE DOC IN TOPIC2)


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

dir = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART"
#dir = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART"
source(paste0(dir,"\\Functions.R"))

options(scipen = 10)

#setwd("./Topic2")
dat = fread("my_analysis_data.csv")
dat = dat[,c("patient_id","time","apptdate","onarv_cf","trt_switch","ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit","new_cupd","new_wupd","new_vupd","outcome","log_weight","log_cd4plus1","vl_0","log_vlplus1"),with=F]
length(unique(dat$patient_id))
dim(dat)
colnames(dat)

theV = c("ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit")
theX = c("new_wupd","log_weight",
         "new_cupd","log_cd4plus1",
         "new_vupd","vl_0","log_vlplus1")


# outcome: 0 disengage, 1 engage, 2 transfer_outcome, 3 death_outcome
# RELEVEL THE OUTCOME TO 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
dat$outcome = dat$outcome+1
######################################################################################################################
### (1) PART 1

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
### (1) PART 2

# transition prob (IN DRAFT 8)

tab = table(dat$outcome_1,dat$outcome)
tab[1,] = tab[1,]/sum(tab[1,])
tab[2,] = tab[2,]/sum(tab[2,])

prop = table(dat$outcome)
prop = prop/sum(prop)
tp = matrix(0,5,4)
rownames(tp) = c("overall","diengaged","engaged","transferred out","died")
colnames(tp) = c("diengaged","engaged","transferred out","died")
tp[1,] = prop
tp[2:3,] = tab
tp[4,3] = tp[5,4] = 1
round(tp,3)


# transition prob over time

tmp = unique(dat$time)[-1]

engtp = matrix(NA,length(tmp),4)
colnames(engtp) = c("diengaged","engaged","transferred out","died")
rownames(engtp) = tmp
dengtp = matrix(NA,length(tmp)-1,4)
colnames(dengtp) = c("diengaged","engaged","transferred out","died")
rownames(dengtp) = tmp[-1]

erestab = vector("list", length = length(tmp))
drestab = vector("list", length = length(tmp))


for(t in 1:length(tmp)){
  
  td = dat[time == tmp[t] & outcome_1=="engaged",]
  tab = table(td$outcome_1,td$outcome)
  
  ptab = lapply(1:4, function(k) tab[k,]/sum(tab[k,]))
  erestab[[t]]  = t(simplify2array(ptab))
  
  td = dat[time == tmp[t] & outcome_1=="disengaged",]
  tab = table(td$outcome_1,td$outcome)
  
  ptab = lapply(1:4, function(k) tab[k,]/sum(tab[k,]))
  drestab[[t]]  = t(simplify2array(ptab))
  
  engtp[t,] = erestab[[t]][2,]
  
  if(t>1) dengtp[t-1,] = drestab[[t]][1,]
  
}

dd = engtp
plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                     group = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
p1 = ggplot(plotdat, aes(time, transition_prob))+ geom_line(aes(group = group,col=group),size=1.2)+
  labs(title = "Transition from engagement")+ labs(y = "Transition Probabilities")+
  theme(plot.title = element_text(size=22),text = element_text(size=20),
        axis.text.x = element_text(size = 20))


dd = dengtp
plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                     group = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
p2 = ggplot(plotdat, aes(time, transition_prob))+ geom_line(aes(group = group,col=group),size=1.2)+
  labs(title = "Transition from disengagement")+ labs(y = "Transition Probabilities")+
  theme(plot.title = element_text(size=22),text = element_text(size=20),
        axis.text.x = element_text(size = 20))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1)

png(filename = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft8\\TPtime.png", width =800, height = 300)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=1,widths=c(10, 2))
dev.off()


######################################################################################################################
### SETTING UP VARIABLES

theA = "onarv_cf"

theY = "outcome"

# create new ID
num = dat[, .N,by=patient_id]
dat$myID = rep(1:length(num[[2]]), times = num[[2]])

### FACTORIZE CATEGORICAL VARIABLES WITH <= c 
newdat = factorizedata(dat,c = 10)
levels(newdat$TravelTime) =c("missing","below30min","over2h","1to2h","30to60min")

### create missing category (for factors)/ (RX,R) (for continuous variables)
newdat = createmis(newdat)

### create indicator variables
newdat = createind(newdat, keep = c(theA,theY,theX)) # variables NOT in keep discards the original factor variable when its indicators are generated

# create (RX,R) for variables with missingness


#printvar(colnames(newdat))

theV = c("ageinit",
         "male_1",
         "yrenrol_2009","yrenrol_2010","yrenrol_2011","yrenrol_2012","yrenrol_2013","yrenrol_2014","yrenrol_2015","yrenrol_2016",
         "TravelTime_below30min","TravelTime_over2h","TravelTime_1to2h","TravelTime_30to60min",
         "whostageinit_2","whostageinit_3","whostageinit_4","whostageinit_NA",
         "marriedinit_1","marriedinit_NA",
         "heightinit_RX","heightinit_R")

theX = c("new_cupd_1","log_cd4plus1") # LHS

theX_cov = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R") # RHS

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

#all(c(theA,theX,theY,theV) %in% colnames(newdat))

newdat = lagdat(newdat, c(theA, theX_cov, theY_ind)) 

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

#theY_lag1 = paste0(theY,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")


######################################################################################################################
### (3)
#[S1 | V, A0]
#[S2 | V, A1, S1]
#[S3 | V, A2, S2]

### Model building
maxt= 13 # for identifiability of parameters. No xfer or death in the observed data at t=15
mod =  vector("list", maxt) 

for(t in 1:maxt){
  if(t==1){
    fml = paste0(theY," ~ ",paste(c(theV,theA_lag1),collapse = " + "))
    mod[[t]] =  multinom(fml , data=newdat[time==200*t,]) 
    print(t)
  } else {
    fml = paste0(theY," ~ ",paste(c(theV,theA_lag1,theY_ind_lag1),collapse = " + "))
    mod[[t]] =  multinom(fml , data=newdat[time==200*t,]) 
    print(t)
  }
}

### Predictions
Vdist = newdat[time==0,c("myID", theV),with=F]
A = 1

Smat = matrix(NA, ncol=maxt, nrow=nrow(Vdist))


for(t in 1:maxt){
  if(t==1){

    dat.pred = data.frame(Vdist, A)
    colnames(dat.pred) = c("myID",theV,theA_lag1)
    cp = predict(mod[[t]], newdata = dat.pred, "probs")  
    Smat[,t] = rMultinom(cp, 1)   
   
    print(t)
  } else {
    
    # need to define theY_ind_lag1 from simulated values; indicator of engagement at t-1 (either eng or diseng)
    dat.pred = data.frame(Vdist, A, (Smat[,t-1] == 2)*1 )
    colnames(dat.pred) = c("myID",theV,theA_lag1,theY_ind_lag1)
    cp = predict(mod[[t]], newdata = dat.pred, "probs")  
    Smat[,t] = rMultinom(cp, 1)  
    
    done_ind = which(Smat[,t-1] %in% c(3,4)) # people who xferred or died from t-1
    Smat[done_ind, t] = Smat[done_ind, t-1]

    print(t)
  }

}


tmp = Smat
tmp = cbind(1,Smat)
tmp = as.factor(tmp)
levels(tmp) = c("diengaged","engaged","transferred out","died")
tmp = relevel(tmp, "engaged")

a = data.frame(time = rep(200*(0:maxt),each = nrow(Smat)), outcome = tmp)
counts <- table(a$outcome,a$time)

png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\SimEngA",A,".png"),width = 600, height = 400)
barplot(counts, main=paste0("Simulated Engagement Status at a=",A,"\n No Time-varying Covariates"),
        xlab="Days since baseline", ylab="Number of observations", col=1:length(levels(a$outcome)))
legend("bottomleft", 
       legend = rownames(counts), 
       fill = 1:4)
dev.off()

######################################################################################################################
### (4) 
#[L1 | V, A0, X0, L0]
#[X1 | V, A0, X0, L0, L1 = 1] # simulated Xt remains NA until an Lt == 1
#[S1 | V, A0, X1, L1]
#[L2 | V, A1, X1, L1, S1]
#[X2 | V, A1, X1, L1, S1, L2 = 1]
#[S2 | V, A1, S1, X2, L2]



fitmod = function(y, fml, fitdat){
  
  a = length(unique(y))
  
  if(is.factor(y) & a > 2){
    m =  multinom(fml , data=fitdat) 
  } else{
    if(a == 2){
      m = glm(fml, family=binomial(link='logit'), data=fitdat) 
    } else {
      m = lm(fml, data=fitdat)
    }
  }
  return(m)
}


predmod = function(y, mod, fitdat){
  
  a = length(unique(y))
  
  if(is.factor(y) & a > 2){
    cp = predict(mod, newdata = fitdat, "probs")  
    yhat = rMultinom(cp, 1) 
  } else{
    if(a == 2){
      cp = predict(mod, newdata = fitdat, type="response")  
      yhat =  rbinom(n=nrow(fitdat), size=1, prob=cp)
    } else {
      yhat = predict(mod, newdata = fitdat) 
      sighat = summary(mod)$sigma
      epsilon = rnorm(nrow(fitdat),mean = 0, sd = sighat)
      yhat = yhat+epsilon
    }
  }
  return(yhat)
}


### Model building
maxt= 13 # for identifiability of parameters. No xfer or death in the observed data at t=15
Smod = vector("list", maxt) 
Xmod = vector("list", maxt)

for(t in 1:maxt){
  Xmod[[t]] = vector("list", 2) 
  
  if(t==1){
    
    fml = paste0(theX[1]," ~ ",paste(c(theV,theA_lag1, theX_cov_lag1),collapse = " + "))
    a = newdat[,get(theX[1])]
    Xmod[[t]][[1]] = fitmod(a, fml, newdat[time==200*t,]) 
    
    fml = paste0(theX[2]," ~ ",paste(c(theV,theA_lag1, theX_cov_lag1),collapse = " + "))
    a = newdat[,get(theX[2])]
    Xmod[[t]][[2]] = fitmod(a, fml, newdat[time==200*t & get(theX[1])==1,])  # only fit for those current upd==1
    
    fml = paste0(theY," ~ ",paste(c(theV,theA_lag1,theX_cov),collapse = " + "))
    Smod[[t]] =  multinom(fml , data=newdat[time==200*t,]) 
    print(t)
  } else {

    fml = paste0(theX[1]," ~ ",paste(c(theV,theA_lag1, theX_cov_lag1, theY_ind_lag1),collapse = " + "))
    a = newdat[,get(theX[1])]
    Xmod[[t]][[1]] = fitmod(a, fml, newdat[time==200*t,]) 
    
    fml = paste0(theX[2]," ~ ",paste(c(theV,theA_lag1, theX_cov_lag1, theY_ind_lag1),collapse = " + "))
    a = newdat[,get(theX[2])]
    Xmod[[t]][[2]] = fitmod(a, fml, newdat[time==200*t & get(theX[1])==1,]) # only fit for those current upd==1
    
    fml = paste0(theY," ~ ",paste(c(theV,theA_lag1,theX_cov,theY_ind_lag1),collapse = " + "))
    Smod[[t]] =  multinom(fml , data=newdat[time==200*t,]) 
    print(t)
  }
}

### Predictions
Vdist = newdat[time==0,c("myID", theV),with=F]
A = 1 # Aprev, is a function of Xprev

Xm = vector("list", 2) 

for(i in 1:2){
  Xm[[i]] = matrix(NA, ncol=maxt, nrow=nrow(Vdist))
}

Sm = matrix(NA, ncol=maxt, nrow=nrow(Vdist))


for(t in 1:maxt){
  if(t==1){
    
    Xprev = newdat[time==0,theX_cov,with=F]
    
    # simulate L
    dat.pred = data.frame(Vdist, A, Xprev)
    colnames(dat.pred) = c("myID",theV,theA_lag1, theX_cov_lag1)
    a = newdat[,get(theX[1])]
    Xm[[1]][,t] = predmod(a, Xmod[[t]][[1]], dat.pred)
    
    # simulate X
    dat.pred = data.frame(Vdist, A, Xprev)
    colnames(dat.pred) = c("myID",theV,theA_lag1, theX_cov_lag1)
    a = newdat[,get(theX[2])]
    Xm[[2]][,t] = predmod(a, Xmod[[t]][[2]], dat.pred)  
    # only simulate X1 when L1 == 1
    ind = which(Xm[[1]][,t]==0)
    Xm[[2]][ind,t] = NA
    
    # define simulated theX_cov
    RX = Xm[[2]][,t]
    R = 1*is.na(RX)
    RX[is.na(RX)] = 0
    Xsim = data.frame(Xm[[1]][,t],  RX, R)
    
    # simulate S
    dat.pred = data.frame(Vdist, A, Xsim)
    colnames(dat.pred) = c("myID",theV,theA_lag1,theX_cov)
    a = newdat[,get(theY)]
    Sm[,t] = predmod(a, Smod[[t]], dat.pred)
    
    
    
    print(t)
  } else {
    
    # define simulated (Xprev,Sprev)
    RX = Xm[[2]][,t-1]
    R = 1*is.na(RX)
    RX[is.na(RX)] = 0
    
    XSprev = data.frame(Xm[[1]][,t-1], RX, R,(Sm[,t-1] == 2)*1)
    
    # simulate L
    dat.pred = data.frame(Vdist, A, XSprev)
    colnames(dat.pred) = c("myID",theV,theA_lag1, theX_cov_lag1, theY_ind_lag1)
    a = newdat[,get(theX[1])]
    Xm[[1]][,t] = predmod(a, Xmod[[t]][[1]], dat.pred)
    
    # simulate X
    dat.pred = data.frame(Vdist, A, XSprev)
    colnames(dat.pred) = c("myID",theV,theA_lag1, theX_cov_lag1, theY_ind_lag1)
    a = newdat[,get(theX[2])]
    Xm[[2]][,t] = predmod(a, Xmod[[t]][[2]], dat.pred)  
    # only simulate Xt when Rt == 1
    ind = which(Xm[[1]][,t]==0)
    Xm[[2]][ind,t] = Xm[[2]][ind,t-1]
    
    # define simulated theX_cov
    RX = Xm[[2]][,t]
    R = 1*is.na(RX)
    RX[is.na(RX)] = 0
    XSsim = data.frame(Xm[[1]][,t],  RX, R, (Sm[,t-1] == 2)*1)
    
    # simulate S
    dat.pred = data.frame(Vdist, A, XSsim)
    colnames(dat.pred) = c("myID",theV,theA_lag1,theX_cov,theY_ind_lag1)
    a = newdat[,get(theY)]
    Sm[,t] = predmod(a, Smod[[t]], dat.pred)
    
    done_ind = which(Sm[,t-1] %in% c(3,4)) # people who xferred or died from t-1
    Sm[done_ind, t] = Sm[done_ind, t-1]

    print(t)
  }
  
}


tmp = Sm
tmp = cbind(2,Sm)
tmp = as.factor(tmp)
levels(tmp) = c("diengaged","engaged","transferred out","died")
tmp = relevel(tmp, "engaged")

a = data.frame(time = rep(200*(0:maxt),each = nrow(Sm)), outcome = tmp)
counts <- table(a$outcome,a$time)

counts = simplify2array(lapply(1:ncol(counts), function(x) counts[,x]/sum(counts[,x])))
colnames(counts) = 200*(0:maxt)

png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\SimEngA",A,"_wCD4_maxt11.png"),width = 600, height = 400)
barplot(counts, main=paste0("Simulated Engagement Status at a=",A,"\n Time-varying Covariates: CD4 info"),
        xlab="Days since baseline", ylab="Number of observations", col=1:4)
legend("bottomleft", 
       legend = rownames(counts), 
       fill = 1:4)
dev.off()




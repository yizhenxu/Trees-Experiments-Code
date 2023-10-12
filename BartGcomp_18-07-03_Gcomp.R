#setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
#setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")

#screennum = 77313

# A = 1 143620
# A = 0 147610
# A = I{CD4 updates &  CD4<=350} 149262

#setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)
#load("I:\\Topic2\\B_PB_MNP_screen54559_ct33_nsub50000.RData")

load("B_PB_MNP_screen54559_ct33_nsub50000.RData")

maxt= 11
nd = 2000
nb = 100

nburn = nb
nthin = 0 # number of discarded samples between two valid samples
basey = "2" # engagement as reference level
nchain = 1 # number of MCMC chains
npost = nd 
ndraw = ceiling(nd/nchain) *(nthin+1) + nb

parell = F
######################################################################################################################
### ZEUS script::

library(GcompBART)

library(lubridate)# date processing
library(data.table)
library(zoo)# carry backward forward
library(foreign)#stata read in
library(nnet)# multinomial log-linear models
library(Hmisc)#rmultinom function

library("parallel")

source("Functions_new1.R")
######################################################################################################################
######################################################################################################################
theA = "onarv_cf"

theY = "outcome"

theV = c("ageinit",
         "male_1",
         "yrenrol",
         "TravelTime_RX","TravelTime_R",
         "whostageinit_RX","whostageinit_R",
         "marriedinit_1","marriedinit_NA",
         "heightinit_RX","heightinit_R",
         "log_weight_init_RX","log_weight_init_R",
         "log_vlplus1_init_RX","log_vlplus1_init_R","vl_0_init_1")


theX = c("new_cupd_1","log_cd4plus1") # LHS

theX_cov = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R") # RHS (from createind and createmis)

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

theA = "onarv_cf_1"

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

load("newdat.RData") # defined in 18-02-26 script
#load("./makeGridData/newdat.RData")

#OUTCOME: 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
newdat$outcome = as.numeric(newdat$outcome) #IMPORTANT! de-factorize outcome

######################################################################################################################

sink(paste0("output_",screennum,".txt")) #write all output to txt -- source('script.R',echo=T) to write input too

### optional: Gcomp on nsub people
seed = 1
nsub = ?
nn = max(newdat$myID)
set.seed(seed)
id = sample(1:nn, nsub, replace = F)

###regular
newdat = newdat[myID %in% id,]
### 

######################################################################################################################

### Outcome simulation 

A = 1 #theA_lag1
nlatent = length(unique(newdat$outcome)) - 1
nsub = sum(newdat$time==0)
ytab = table(newdat$outcome)
ylevels = as.numeric(names(ytab)[ytab > 0])

Smat = vector("list", maxt) 
Xmat = vector("list", maxt)

for(i in 1:maxt){
  Xmat[[i]] = array(NA, dim = c(nsub, nd, length(theX)))
}

Vdist = cbind(1,newdat[time==0,theV,with=F])
Vdist = Vdist[rep(1:nsub,nd),]

ptm <- proc.time()

ct = 1

for(t in 1:maxt){
  
  if(t==1){
    
    print(ct)
    #[L1 | V, A0, X0, L0]
    policy_lag = rep(A, nsub)
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, c(1,2,4), t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)

    Xmat[[t]][,,1] = predict_bart(mod[[ct]],dat.pred)$samp_y

    ct = ct+1
    
    print(ct)
    #[X1 | V, A0, X0, L0] 
    #dat.pred = makeX(Xname = theX, Vname = theV, time = t, modhere = mod[[ct]], simX = Xmat, policy_lag = A, dathere = newdat)
    # same RHS as in modeling L1
    postY = predict_bart(mod[[ct]],dat.pred)$samp_y
    
    L = Xmat[[t]][,,1]
    mis_ind = which(L == 0) # people who have no CD4 update
    postY[mis_ind] = NA
    Xmat[[t]][,,2] = postY
    remove(postY)
    
    ct = ct+1
    
    print(ct)
    #[S1 | V, A0, X1, L1]
    policy_lag = rep(A, nrow(Vdist))
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, time = t, modhere = modfml[[ct]], simX = Xmat, policy_lag = policy_lag, dathere = newdat)
    Smat[[t]] = MNPpredict(m = mod[[ct]], Xmat = dat.pred, nsub, basey, c(1,2,4), doparallel = parell)
    
    ct = ct+1
    
    print(paste0("t=",t))
    
  } else {
    
    print(ct)
    #[Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    #policy_lag = rep(A, nrow(Vdist))
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)
    Xmat[[t]][,,1] = predict_bart(mod[[ct]],dat.pred, Gcomp = TRUE)$samp_y
    
    ct = ct+1
    
    print(ct)
    #[Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    #policy_lag = rep(A, nrow(Vdist))
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)
    postY = predict_bart(mod[[ct]],dat.pred, Gcomp = TRUE)$samp_y
    
    L = Xmat[[t]][,,1]
    X_prev = Xmat[[t-1]][,,2]
    mis_ind = which(L == 0) # people who have no CD4 update keeps the old value
    postY[mis_ind] = X_prev[mis_ind]
    Xmat[[t]][,,2] = postY
    
    remove(postY)
    
    ct = ct+1
    
    print(ct)
    #[St | V, A_{t-1}, S_{t-1}, Xt, Lt]
    #policy_lag = rep(A, nrow(Vdist))
    dat.pred = makeX(Xname = theX, Xcovname = theX_cov, Vname = theV, repV = Vdist, yl = ylevels, time = t, modhere = modfml[[ct]], simX = Xmat, simY = Smat, policy_lag = policy_lag, dathere = newdat)
    Smat[[t]] = MNPpredict(m = mod[[ct]], Xmat = dat.pred, nsub, basey, ylevels, doparallel = parell)# basey is disengagement
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
    
    ct = ct+1
    
    print(paste0("t=",t))
    
    
  }
  
}

proc.time() - ptm
#
save(Xmat,Smat,file = paste0("Gcomp_screen",screennum,"_ct",ct,"_nsub",nsub,".RData"))

##################################################################################
library(mailR)
a =  proc.time()-ptm
status = as.character(a[3])

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

##################################################################################
##################################################################################

### BARPLOT

setwd("./Topic2")

load("Gcomp_screen143620_ct34_nsub30000.RData")#A=1
load("Gcomp_screen147610_ct34_nsub30000.RData")#A=0
load("Gcomp_screen149262_ct34_nsub30000.RData")#A=I{CD4<=350}
A = "1"
A = "0"
A = "CD4_350"

library(data.table)
maxt = 11

tmp = unlist(Smat)
#tmp = Smat
tmp1 = c(rep(2,length(Smat[[1]])),tmp)
#a = data.frame(time = rep(200*(0:maxt),each = length(Smat[[1]])), outcome = tmp1)
#counts <- table(a$outcome,a$time)
#longdat = setDT(a)
longdat = data.table(time = rep(200*(0:maxt),each = length(Smat[[1]])), outcome = tmp1)
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

save(counts,file = paste0("Counts_A",A,".RData"))

#load("I:\\Topic2\\tmp1.RData")


path = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"
load(paste0(path,"\\Counts_A",A,".RData"))

png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\BARTCD4SimEngA",A,".png"),width = 600, height = 400)
png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft9\\BPBMNPCD4SimEngA",A,".png"),width = 600, height = 400)

barplot(counts, main=paste0("Bayesian Simulated Engagement Status at a=",A,"\n Time-varying Covariates: CD4 info"),
        xlab="Days since baseline", ylab="Number of observations", col=1:4)
legend("bottomleft", 
       legend = rownames(counts), 
       fill = 1:4)
dev.off()

##################################################################################
##################################################################################

### MARGINAL PROBABILITY OF STATE MEMBERSHIP (WITH 95% CI)

CIvec = function(vec){
  v = sort(vec)
  n = length(vec)
  l = v[floor(n*0.025)]
  u = v[ceiling(n*0.975)]
  l = round(l,3)
  u = round(u,3)
  return(paste0(" (",l,", ",u,")"))
}

#maketab = function(mat){
#  longdat = data.table(draw = rep(1:ncol(mat),each = nrow(mat)), outcome = unlist(mat))
#  ldc = longdat[, .N, by = c("draw","outcome")]
#  return(ldc)
#}

sumt = function(i, Smat){
  print(i)
  res = apply(Smat[[i]],2, table)
  res1 = res/30000
  m = apply(res1,1,mean)
  m = round(m,3)
  ci = apply(res1,1,CIvec)
  return(paste0(m, ci))
}


tmp = lapply(1:maxt, function(i) sumt(i, Smat))
pt = matrix(unlist(tmp[2:maxt]), nrow = 4)
pt = cbind(rep(NA,4),pt)
pt[c(1,2,4),1] = tmp[[1]]
colnames(pt) = 200*(1:maxt)
rownames(pt) = c("Diseng","Eng","Xfer","Death")
print.noquote(pt)


#A = "1"
save(pt,file = paste0("Smat_Summ_A",A,".RData"))


##################################################################################
##################################################################################

### TRANSITION PROBABILITY FROM ENG/DISENG

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

A = "1"
A = "0"
A = "CD4_350"

maxt = 11

tmp = unlist(Smat)
#tmp = Smat
m = length(Smat[[1]])
tmp1 = c(rep(2,m),tmp[1:((maxt-1)*m)])

#a = data.frame(time = rep(200*(0:maxt),each = length(Smat[[1]])), outcome = tmp1)
#counts <- table(a$outcome,a$time)
#longdat = setDT(a)
longdat = data.table(time = rep(200*(1:maxt),each = length(Smat[[1]])), outcome_lag = tmp1, outcome = tmp)
ldc = longdat[, .N, by = c("time", "outcome_lag", "outcome")]
ldc_denom = longdat[, .N, by = c("time","outcome_lag")]

save(ldc,ldc_denom, file = "tmp.RData")
#load("tmp.RData")

ldc = ldc[-which(ldc$outcome_lag %in% c(3,4)),]
ldc = ldc[order(time, outcome_lag, outcome)]
ldc_denom = ldc_denom[-which(ldc_denom$outcome_lag %in% c(3,4)),]
ldc_denom = ldc_denom[order(time, outcome_lag)]

eldc = ldc[outcome_lag == 2,]
eldc = rbind(eldc,list(time = 200, outcome_lag = 2, outcome = 3, N = 0 ))
eldc = eldc[order(time, outcome, outcome_lag)]

dldc = ldc[outcome_lag == 1,]

eldc_denom = ldc_denom[outcome_lag == 2,]
dldc_denom = ldc_denom[outcome_lag == 1,]

engtp = matrix(eldc$N/rep(eldc_denom$N, each = 4), nrow = 4)
rownames(engtp) = c("diengaged","engaged","transferred out","died")
colnames(engtp) = unique(eldc$time)
dengtp = matrix(dldc$N/rep(dldc_denom$N, each = 4), nrow = 4)
rownames(dengtp) = c("diengaged","engaged","transferred out","died")
colnames(dengtp) = unique(dldc$time)

dd = t(engtp)
plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                     group = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
p1 = ggplot(plotdat, aes(time, transition_prob))+ geom_line(aes(group = group,col=group),size=1.2)+
  labs(title = "Transition from engagement")+ labs(y = "Transition Probabilities")+
  theme(plot.title = element_text(size=22),text = element_text(size=20),
        axis.text.x = element_text(size = 20))


dd = t(dengtp)
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

png(filename = paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft9\\TPtime_A",A,".png"), width =800, height = 300)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=1,widths=c(10, 2))
dev.off()

####################################################################################################

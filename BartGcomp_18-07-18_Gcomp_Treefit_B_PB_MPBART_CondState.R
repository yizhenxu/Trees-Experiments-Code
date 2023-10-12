#setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
#setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")

#screennum = 77313

# BART-PBART-MNP(ConS)
# A = 1 84286
# A = 0 88437
# A = I{CD4 updates &  CD4<=350} 88423 



#setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)
#load("I:\\Topic2\\B_PB_MNP_screen54559_ct33_nsub50000.RData")

load("CondSModel_screen122548_ct31_ntr40000.RData")

maxt= 10
basey1 = "2"
basey = c(rep("2", 2*(maxt-1)-1), "1") 
nd = 2000
nb = 100

ndraw = nd+nb 
nburn = nb
nthin = 0 # number of discarded samples between two valid samples
nchain = 1 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain

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

source("Functions_new.R")
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

### optional: pick nsub people
seed = 1
nsub = 90000
nn = max(newdat$myID)
set.seed(seed)
id = sample(1:nn, nsub, replace = F)

###regular
newdat = newdat[myID %in% id,]
### 

fitdat = newdat[time == 200,]
table(fitdat$outcome)
fitdat = newdat[time == 400,]
table(fitdat$outcome)


### Separate data into training set and test set

ntr = 40000 # floor(nsub * 0.7) # sample size of training set
nte = 10000
ngc = 40000

set.seed(1)
samp = sample(id,nsub)

trid = samp[1:ntr]  
trid = sort(trid)
trd = copy(newdat[myID %in% trid,])

teid = samp[(ntr+1):(ntr+nte)]  
teid = sort(teid)
ted = copy(newdat[myID %in% teid,])

gcid = samp[(ntr+nte+1):nsub]  
gcid = sort(gcid)
gcd = copy(newdat[myID %in% gcid,])


fitdat = trd[time == 200,]
table(fitdat$outcome)
fitdat = trd[time == 400,]
table(fitdat$outcome)

fitdat = ted[time == 200,]
table(fitdat$outcome)
fitdat = ted[time == 400,]
table(fitdat$outcome)
###

newdat = gcd

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
sct =1

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
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, c(1,2,4), t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)
    
    Smat[[t]] = predict_bart(mod[[ct]],dat.pred, Gcomp = TRUE)$samp_treefit
    
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
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, modfml[[ct]], Xmat, Smat, policy_lag, newdat)
    
    fitdat = trd[time==200*t,]
    fitdat1 =  fitdat[get(theY_ind_lag1)==0,]
    ytab = table(fitdat1$outcome)
    ylevels = as.numeric(names(ytab)[ytab > 0])
    
    Smat[[t]] = predict_bart(Smod[[sct]],dat.pred, Gcomp = TRUE)$samp_treefit
    sct = sct + 1
    
    fitdat1 =  fitdat[get(theY_ind_lag1)==1,]
    ytab = table(fitdat1$outcome)
    ylevels = as.numeric(names(ytab)[ytab > 0])
    tmp1 = predict_bart(Smod[[sct]],dat.pred, Gcomp = TRUE)$samp_treefit
    sct = sct + 1
    
    eng_ind = which(Smat[[t-1]] == 2)
    Smat[[t]][eng_ind] = tmp1[eng_ind]
    remove(tmp1)
    
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
    
    ct = ct+1
    
    print(paste0("t=",t))
    
    
  }
  
}

proc.time() - ptm
#
save(Xmat,Smat,file = paste0("CondSGcomp_screen",screennum,"_ct",ct,"_nsub",nsub,".RData"))

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

load("Gcomp_screen60786_ct34_ngc40000.RData")#A=1
load("Gcomp_screen60754_ct34_ngc40000.RData")#A=0
load("Gcomp_screen60772_ct34_ngc40000.RData")#A=I{CD4<=350}
type = "MNP"


load("Gcomp_screen60815_ct34_ngc40000.RData")#A=1
load("Gcomp_screen60464_ct34_ngc40000.RData")#A=0
load("Gcomp_screen60800_ct34_ngc40000.RData")#A=I{CD4<=350}
type = "MPBART"

load("CondSGcomp_screen84286_ct31_nsub40000.RData")#A=1
load("CondSGcomp_screen88437_ct31_nsub40000.RData")#A=0
load("CondSGcomp_screen88423_ct31_nsub40000.RData")#A=I{CD4<=350}
type = "CondS"
A = "1"
A = "0"
A = "CD4_350"

library(data.table)
maxt = 10

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

save(counts,file = paste0("Counts_",type,"_A",A,".RData"))

#load("I:\\Topic2\\tmp1.RData")


path = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"
load(paste0(path,"\\Counts_",type,"_A",A,".RData"))

png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft9\\Accuracy Comparison - Diag Plots\\07-13\\BARTCondSSimEngA",A,".png"),width = 600, height = 400)
#png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft9\\Accuracy Comparison - Diag Plots\\07-13\\BPBMNPSimEngA",A,".png"),width = 600, height = 400)
#png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft9\\Accuracy Comparison - Diag Plots\\07-13\\BARTSimEngA",A,".png"),width = 600, height = 400)

barplot(counts,
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

sumt = function(i, Smat, ngc){
  print(i)
  res = apply(Smat[[i]],2, table)
  res1 = res/ngc
  m = apply(res1,1,mean)
  m = round(m,3)
  ci = apply(res1,1,CIvec)
  return(paste0(m, ci))
}


tmp = lapply(1:maxt, function(i) sumt(i, Smat, ngc))
pt = matrix(unlist(tmp[2:maxt]), nrow = 4)
pt = cbind(rep(NA,4),pt)
pt[c(1,2,4),1] = tmp[[1]]
colnames(pt) = 200*(1:maxt)
rownames(pt) = c("Diseng","Eng","Xfer","Death")
print.noquote(pt)


#A = "1"
save(pt,file = paste0("Smat_Summ_A",A,".RData"))

engci = function(i, Smat, ngc){
  res = apply(Smat[[i]],2, table)
  res1 = res/ngc
  return(res1)
}
tmp = lapply(1:maxt, function(i) engci(i, Smat, ngc)) #in A = 1 and A = 1{CD4 < 350}  

###
k=4

tmp1 = tmp#A = 1
for(i in 1:maxt){
  tmp1[[i]] = tmp[[i]][k,]
}

###
#k=...
tmp2 = tmp#A = 1{CD4 < 350}
for(i in 1:maxt){
  tmp2[[i]] = tmp[[i]][k,]
}
save(tmp2, file = "tmp2.RData")
###

load("tmp2.RData")

nd = length(tmp1[[1]])
tmp3 = tmp1
tab = matrix(NA, 3, maxt)
for(i in 1:maxt){
  tmp3[[i]] = tmp1[[i]] - tmp2[[i]]
  tab[1,i] = mean(tmp3[[i]])
  a = sort(tmp3[[i]])
  tab[2,i] = a[floor(nd*0.025)]
  tab[3,i] = a[ceiling(nd*0.975)]
  
}

save(tab, file = "tab.RData")

#png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft9\\Accuracy Comparison - Diag Plots\\07-13\\BART_EngPercCI",A,".png"),width = 600, height = 400)
#png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft9\\Accuracy Comparison - Diag Plots\\07-13\\BPBMNP_EngPercCI",A,".png"),width = 600, height = 400)
png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft9\\Accuracy Comparison - Diag Plots\\07-13\\CondS_EngPercCI",A,"_k=",k,".png"),width = 600, height = 400)
#when k = 3, transffered, no first time point
#tab = tab[,-1]
#x = 200*(2:maxt)
x = 200*(1:maxt)
cname = c("Disengaged", "Engaged", "Transferred", "Died")
plot(x, tab[1,],
     ylim=range(c(tab)),
     pch=19, xlab="Time", ylab=paste0("Differece in Percent ",cname[k]),
     main="Treat Immediately v.s Treat when CD4 < 350"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(x, tab[2,], x, tab[3,], length=0.05, angle=90, code=3)
abline(h = 0, lty = 2)
dev.off()
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

maxt = 10

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
rownames(engtp) = c("disengaged","engaged","transferred out","died")
colnames(engtp) = unique(eldc$time)

dengtp = matrix(dldc$N/rep(dldc_denom$N, each = 3), nrow = 3)
rownames(dengtp) = c("disengaged","engaged","died")
colnames(dengtp) = unique(dldc$time)

library(car)
engtp1 = logit(engtp)
dengtp1 = logit(dengtp)

range(c(engtp1,dengtp1))

rg = c( -8.949624,  4.669631)

dd = t(engtp1)
plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                     group = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd) )
p1 = ggplot(plotdat, aes(time, transition_prob))+ geom_line(aes(group = group,col=group),size=1.2)+
  labs(title = "Transition from engagement")+ labs(y = "Transition Probabilities")+
  theme(plot.title = element_text(size=22),text = element_text(size=20),
        axis.text.x = element_text(size = 20))+ 
  scale_color_manual(values=c(4,2,1,3)) + ylim( -8.949624,  4.669631)
#plot(p1)

dd = t(dengtp1)
plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                     group = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
p2 = ggplot(plotdat, aes(time, transition_prob))+ geom_line(aes(group = group,col=group),size=1.2)+
  labs(title = "Transition from disengagement")+ labs(y = "Transition Probabilities")+
  theme(plot.title = element_text(size=22),text = element_text(size=20),
        axis.text.x = element_text(size = 20))+ 
  scale_color_manual(values=c(4,2,1))+ ylim( -8.949624,  4.669631)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1)

png(filename = paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft9\\Accuracy Comparison - Diag Plots\\07-13\\TPtime_A",A,".png"), width = 900, height = 300)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=1,widths=c(10, 2))
dev.off()

####################################################################################################

load("CondSModel_screen84204_ct30_ntr40000.RData")

#Outcomes model at t = 3 (day 600) sct = 3,4
cname = c("Disengaged", "Engaged","Transfer","Death")

m = Smod[[3]]
nlatent = 2 
yl = c(1,4)

m = Smod[[4]]
nlatent = 3 
yl = c(1,3,4)

s = apply(m, 2, sd)
mn = apply(m, 2, mean)

nhalfSig = nlatent*(nlatent - 1)/2 + nlatent
sigind = (ncol(m) - nhalfSig + 1): ncol(m)
parind = 1:(ncol(m) - nhalfSig)
ind = rep(1:nlatent,length(parind)/nlatent)

mntab = c()
for(j in 1:nlatent){
  mntab = cbind(mntab, mn[which(ind == j)])
}
colnames(mntab) = cname[yl]
rownames(mntab) = gsub(":1", "", rownames(mntab))

stab = c() 
for(j in 1:nlatent){
  stab = cbind(stab, s[which(ind == j)])
}

sigtab = 1*(abs(mntab/stab)>=2)

restab = round(mntab, 3)
for(i in 1:length(restab)){
  if(sigtab[i]==1) restab[i] = paste0(restab[i], "*")
}
noquote(restab)


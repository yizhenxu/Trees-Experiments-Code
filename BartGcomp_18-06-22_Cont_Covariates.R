#devtools::install_github("yizhenxu/GcompBART")
#remove.packages("GcompBART")

#setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")

screennum =150594
#BART ref death
#BART 20 trees 197009 (120283) ; Base engagement 197009 (120283)
#BART 50 trees 197598 (120325)
#BART 100 trees 1997643 (125505)
#LM ref disengagement 197719 (120347)
#LM ref death 36620
setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)
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
#setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")
#dat = fread("./makeGridData/my_analysis_data.csv")
dat = fread("my_analysis_data.csv")

dat = dat[,c("patient_id","time","apptdate","onarv_cf","trt_switch","ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit","new_cupd","new_wupd","new_vupd","outcome","log_weight","log_cd4plus1","vl_0","log_vlplus1"),with=F]
length(unique(dat$patient_id))

theV = c("ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit")
theX = c("new_wupd","log_weight",
         "new_cupd","log_cd4plus1",
         "new_vupd","vl_0","log_vlplus1")


# outcome: 0 disengage, 1 engage, 2 transfer_outcome, 3 death_outcome
# RELEVEL THE OUTCOME TO 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
dat$outcome = dat$outcome+1

### need to check data processing
### delete the 16 ppl who were xferred out at time==200
# check.... something wrong with xfer definition
todelete = dat[time == 200 & outcome == 3, "patient_id"]
delloc = which(dat$patient_id %in% todelete$patient_id)
dat = dat[-delloc,]

### SETTING UP VARIABLES

theA = "onarv_cf"

theY = "outcome"

# create new ID
num = dat[, .N,by=patient_id]
dat$myID = rep(1:length(num[[2]]), times = num[[2]])

dat[, log_weight_init := log_weight[1],by = patient_id]
dat[, log_cd4plus1_init := log_cd4plus1[1],by = patient_id]
dat[, log_vlplus1_init := log_vlplus1[1],by = patient_id]
dat[, vl_0_init := vl_0[1], by = patient_id]

if(0){
  for(j in 1:ncol(dat)){
    vname = colnames(dat)[j]
    tmp = length(unique(dat[,get(vname)]))
    if(tmp <= c){
      print(vname)
      print(table(dat[,get(vname)],useNA="always"))
    }
  }
}

#0 1: onarv_cf male new_cupd new_wupd new_vupd outcome vl_0
#0 1 NA: marriedinit
#-1 0 1: trt_switch
#cont: yrenrol 
#cont NA: TravelTime, whostageinit
### FACTORIZE CATEGORICAL VARIABLES WITH <= c 
newdat = factorizedata(dat,c = 4) # the 4 is not arbitrary...

a = rep(NA, nrow(newdat))
a[newdat$TravelTime == "< 30 minutes"] = 1
a[newdat$TravelTime == "30-60 minutes"] = 2
a[newdat$TravelTime == "1-2 hours"] = 3
a[newdat$TravelTime == "> 2 hours"] = 4
newdat$TravelTime = a

### create missing category (for factors)/ (RX,R) (for continuous variables)
newdat = createmis(newdat)

### create indicator variables for factors
newdat = createind(newdat, keep = c(theA,theY,theX)) # variables NOT in 'keep' discards the original factor variable when its indicators are generated


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

newdat = lagdat(newdat, c(theA, theX_cov, theY_ind)) 

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

newdat = newdat[,unique(c("myID","time",theY,theV,theX,theX_cov,theY_ind,theA,theA_lag1,theX_cov_lag1,theY_ind_lag1)),with=F]

save(newdat, file = "newdat.RData")

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

fitdat = newdat[time == 200,]
table(fitdat$outcome)
#1     2     4
#10903 73631  3664
fitdat = newdat[time == 400,]
table(fitdat$outcome)
#1     2     3     4
#26799 52493   639   871

#####################################
### optional: pick nsub people
seed = 1
nsub = 30000
nn = max(newdat$myID)
set.seed(seed)
id = sample(1:nn, nsub, replace = F)

###regular
newdat = newdat[myID %in% id,]
### OR
### Manually increase xfer (t==400) & death (both t) in trd and ted
id4t1 = newdat[time==200 & outcome == 4, "myID"]$myID
id34t2 = newdat[time==400 & outcome > 2, "myID"]$myID
id = unique(c(id,id4t1,id34t2))
newdat = newdat[myID %in% id,]
nsub = length(id)
###

fitdat = newdat[time == 200,]
table(fitdat$outcome)
fitdat = newdat[time == 400,]
table(fitdat$outcome)


### Separate data into training set and test set

ntr = floor(nsub*2/3) # floor(nsub * 0.7) # sample size of training set

set.seed(1)
samp = sample(id,nsub)

trid = samp[1:ntr]  
trid = sort(trid)
trd = copy(newdat[myID %in% trid,])

teid = samp[(ntr+1):nsub]  
teid = sort(teid)
ted = copy(newdat[myID %in% teid,])

fitdat = trd[time == 200,]
table(fitdat$outcome)
fitdat = trd[time == 400,]
table(fitdat$outcome)

fitdat = ted[time == 200,]
table(fitdat$outcome)
fitdat = ted[time == 400,]
table(fitdat$outcome)
#####################################

##############################################################
### Build models based on trd (training data)
#parell = TRUE # TRUE if on server, FALSE otherwise
#ncore = 5 # NULL if parall = FALSE

#maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15

### (1) Bayesian with CD4 information


tmp = (length(theX)+1)*maxt
mod =  vector("list", tmp) 
modfml =  vector("list", tmp) 


# initialize index count
ct = 1


ptm <- proc.time()
for(t in 1:maxt){
  fitdat = trd[time==200*t,]
  if(t==1){
    
    # [L1 | V, A0, X0, L0]
    
    RHS = c(theV, theA_lag1, theX_cov_lag1)
    # exclude lag X_R from models at t=1 because that is the same as lag L at baseline
    RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    tmp = model_bart(as.formula(fml), data = fitdat, type = "binary",
                     Prior = Prior_binary,
                     Mcmc = Mcmc_binary, make.prediction = FALSE)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod[[ct]] = tmp[-del]
    remove(tmp)
    ct = ct+1
    
    # [X1 | V, A0, X0, L0, L1 = 1] # simulated Xt remains NA until a first Lt == 1
    
    RHS = c(theV, theA_lag1, theX_cov_lag1)
    # exclude lag X_R from models at t=1 because that is the same as lag L at baseline
    RHS = RHS[-which(RHS=="log_cd4plus1_R_lag1")]
    RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
    fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    tmp = model_bart(as.formula(fml), data = fitdat[get(theX[1])==1,], type = "continuous",
                     Prior = Prior_cont,
                     Mcmc = Mcmc_cont, make.prediction = FALSE)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod[[ct]] = tmp[-del]
    remove(tmp)
    ct = ct+1
    
    # [S1 | V, A0, X1, L1]
    
    ytab = table(fitdat$outcome)
    ylevels = as.numeric(names(ytab)[ytab > 0])
    p = length(ylevels) # number of outcome categories
    
    RHS = c(theV,theA_lag1,theX_cov)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    tmp = model_bart(as.formula(fml), data = fitdat, type = "multinomial",
                     Prior = Prior_mult(p),
                     Mcmc = Mcmc_mult(p), make.prediction = FALSE)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod[[ct]] = tmp[-del]
    remove(tmp)
    ct = ct+1
    
    print(t)
  } else {
    
    # [Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    
    RHS = c(theV, theA_lag1, theX_cov_lag1,theY_ind_lag1)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theX[1]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    tmp = model_bart(as.formula(fml), data = fitdat, type = "binary",
                     Prior = Prior_binary,
                     Mcmc = Mcmc_binary, make.prediction = FALSE)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod[[ct]] = tmp[-del]
    remove(tmp)
    ct = ct+1
    
    # [Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}, Lt = 1]
    
    RHS = c(theV, theA_lag1, theX_cov_lag1, theY_ind_lag1)
    RHS = fmlfilter(fitdat[get(theX[1])==1,],RHS)
    fml = paste0(theX[2]," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    tmp = model_bart(as.formula(fml), data = fitdat[get(theX[1])==1,], type = "continuous",
                     Prior = Prior_cont,
                     Mcmc = Mcmc_cont, make.prediction = FALSE)
    del = which(names(tmp) %in% c("treefit","samp_y"))
    mod[[ct]] = tmp[-del]
    remove(tmp)
    ct = ct+1
    
    # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]
    
    ytab = table(fitdat$outcome)
    ylevels = as.numeric(names(ytab)[ytab > 0])
    p = length(ylevels) # number of outcome categories
    
    RHS = c(theV, theA_lag1, theX_cov, theY_ind_lag1)
    RHS = fmlfilter(fitdat,RHS)
    fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
    modfml[[ct]] = fml
    
    if(p>2){
      tmp = model_bart(as.formula(fml), data = fitdat, type = "multinomial",
                       Prior = Prior_mult(p),
                       Mcmc = Mcmc_mult(p), make.prediction = FALSE)
      del = which(names(tmp) %in% c("treefit","samp_y"))
      mod[[ct]] = tmp[-del]
      remove(tmp)
    }else if(p==2){
      
      Y = rep(0,nrow(fitdat))
      Y[fitdat$outcome == ylevels[2]] = 1
      fitdat$outcome = Y
      tmp = model_bart(as.formula(fml), data = fitdat, type = "binary",
                       Prior = Prior_binary,
                       Mcmc = Mcmc_binary, make.prediction = FALSE)
      del = which(names(tmp) %in% c("treefit","samp_y"))
      mod[[ct]] = tmp[-del]
      remove(tmp)
      #postY = mod[[ct]]$samp_y
      #postY[which(postY==0)] = ylevels[1]
      #postY[which(postY==1)] = ylevels[2]
      #mod[[ct]]$samp_y = postY
    }
    
    
    ct = ct+1
    print(t)
  }
}
proc.time() - ptm



### Evaluation --- per model on test set
# types and intypes give the evalulation output format
Eval_permod_bart = function(mod, modfml, ted, times, types, intypes, maxt){
  
  npost = mod[[1]]$ndraws
  lenmod = length(mod)
  nummod = length(types)
  
  E = vector("list", nummod)
  
  for(i in 1:nummod){
    if(types[i] == "mult"){
      E[[i]] = matrix(NA,5,maxt)
      rownames(E[[i]]) = c("y=1","y=2","y=3","y=4","overall")
      colnames(E[[i]]) = paste0("t=",1:maxt)
    } else {
      E[[i]] = rep(NA, maxt)
    }
  }
  
  #ylist = vector("list",lenmod) # list of observed y's
  #PM = vector("list",lenmod) # list of predictions' matrices
  
  for(ct in 1:lenmod){
    
    t = times[ct]
    
    fv = fmlvars(modfml[[ct]], itc = 0)
    
    dat.pred = ted[time==200*t,fv$RHS,with=F]
    y =  unlist(ted[time==200*t,fv$LHS,with = F])
    
    if(mod[[ct]]$type == "continuous"){
      upd = unlist(ted[time==200*t, theX[1],with = F])
      ind = which(upd==1)
      y = y[ind]
      dat.pred = dat.pred[ind,]
    }
    
    PM = predict_bart(mod[[ct]], dat.pred)$samp_y
    yM = matrix(rep(y, npost),ncol = npost)
    #ylist[[ct]] = y
    
    if(mod[[ct]]$type == "binary"){
      
      E[[intypes[ct]]][t] = mean(PM == yM)
      
    } else if(mod[[ct]]$type == "multinomial"){
      
      for(k in 1:4){
        ind = which(yM==k)
        E[[intypes[ct]]][k,t] = mean(PM[ind] == yM[ind])
      }
      E[[intypes[ct]]][5,t] = mean(PM == yM)
      
    } else {
      
      E[[intypes[ct]]][t] = mean((PM - yM)^2)
      
    }
    
    
    print(ct)
  }
  
  #res = list(Pred = PM, Y = ylist, Accuracy = E)
  res = E
  return(res)
}

# Bayesian G Comp
times = rep(1:maxt,each = 3) # time interval index for each ct model 
types = c("binary","cont","mult") # how many and which type should each E[[i]] be
intypes = rep(1:3, maxt) # which E[[i]] to go into

bart_allres = Eval_permod_bart(mod, modfml, ted, times, types, intypes, maxt)

##################################################################################
##################################################################################
library(mailR)
a =  proc.time()-ptm
status = paste0(screennum,",",as.character(a[3]))

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
save(mod, modfml, bart_allres, file = paste0("Bart_screen",screennum,"_ct",ct,"_nsub",nsub,".RData"))

##################################################################################

##################################################################################
######################################################################################################################
##################################################################################
######################################################################################################################

load("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData\\Bart_screen197009_ct7_nsub30000.RData")
load("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData\\Bart_screen197598_ct7_nsub30000.RData")
load("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData\\BLM_screen197719_ct7_nsub30000.RData")
DiagPlot(mod[[3]],1)
plot(mod[[6]]$sigmasample[3+9*(1:1000)])









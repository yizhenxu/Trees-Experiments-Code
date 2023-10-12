setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot")

load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS CODE/CondSModel_screen132394_ct30ntr50000.RData")

for(t in 1:10){
  ct = (t-1)*3+1
  png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot\\Binary_",t,".png"),width = 600, height = 400)
  DiagPlot(mod[[ct]],1)
  dev.off()
  
}

for(t in 1:10){
  ct = (t-1)*3+2
  png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot\\Continuous_",t,".png"),width = 600, height = 400)
  DiagPlot(mod[[ct]],1)
  dev.off()
}

for(t in 1:10){
  ct = (t-1)*3+2
  png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot\\Continuous_",t,"Sig.png"),width = 600, height = 400)
  plot(mod[[ct]]$sigmasample)
  dev.off()
}


for(t in 2:10){
  sct = (t-1)*2-1
  png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot\\MultinomialS0_",t,".png"),width = 600, height = 400)
  DiagPlot(Smod[[sct]],1)
  dev.off()
}

for(t in 2:10){
  sct = (t-1)*2-1
  png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot\\MultinomialS0_",t,"Sig.png"),width = 600, height = 400)
  plot(Smod[[sct]]$sigmasample)
  dev.off()
}

for(t in 1:10){
  sct = (t-1)*2
  if(t == 1){
    png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot\\MultinomialS1_",t,".png"),width = 600, height = 400)
    
    DiagPlot(mod[[3]],1)
  } else {
    png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot\\MultinomialS1_",t,".png"),width = 600, height = 600)
    
    DiagPlot(Smod[[sct]],1)
  }
  
  dev.off()
}


for(t in 1:10){
  sct = (t-1)*2
    if(t == 1){
    png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot\\MultinomialS1_",t,"Sig.png"),width = 600, height = 400)
      
    plot(mod[[3]]$sigmasample)
  } else {
    png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\18-09-09_convergence plot\\MultinomialS1_",t,"Sig.png"),width = 600, height = 600)
    
    plot(Smod[[sct]]$sigmasample)
  }
  
  dev.off()
}


######################################################################################################################
theA = "onarv_cf"

theY = "outcome"


theV = c("ageinit",
         "male_1",
         "yrenrol_2009","yrenrol_2010","yrenrol_2011","yrenrol_2012","yrenrol_2013","yrenrol_2014", "yrenrol_2015", "yrenrol_2016",
         "TravelTime_below30min","TravelTime_over2h","TravelTime_1to2h","TravelTime_30to60min",
         "whostageinit_2","whostageinit_3","whostageinit_4","whostageinit_NA",
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

#OUTCOME: 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
newdat$outcome = as.numeric(newdat$outcome) #IMPORTANT! de-factorize outcome

######################################################################################################################
######################################################################################################################
### optional: pick nsub people
seed = 1
nsub = 90000
nn = max(newdat$myID)
set.seed(seed)
id = sample(1:nn, nsub, replace = F)

###regular
newdat = newdat[myID %in% id,]
### 

### Separate data into training set and test set

ntr = 50000 # floor(nsub * 0.7) # sample size of training set
nte = 10000

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
ngc = length(unique(gcid))
##############################################################
# initialize index count
ct = 1
summ = vector("list", 2*9+1)

for(t in 1:10){
  fitdat = trd[time==200*t,]
  if(t==1){

    # [S1 | V, A0, X1, L1]
    
    ## trd
    summ[[ct]] = table(fitdat$outcome)
    names(summ)[ct] = "t=1"
    ##########################
    ct = ct+1
    
  } else {
    
    # [Lt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}]
    
    # [Xt | V, A_{t-1}, X_{t-1}, L_{t-1}, S_{t-1}, Lt = 1]
 
    # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]
 
    ### S_{t-1} = disengaged
    ## trd
    summ[[ct]] = table(fitdat[get(theY_ind_lag1)==0,]$outcome)
    names(summ)[ct] = paste0("t=",t,"_S0")
    ct = ct+1
    ### S_{t-1} = engaged
    summ[[ct]] = table(fitdat[get(theY_ind_lag1)==1,]$outcome)
    names(summ)[ct] = paste0("t=",t,"_S1")
    ct = ct+1
  }
}

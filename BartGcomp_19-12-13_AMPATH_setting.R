library(readstata13)
library(lubridate)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(data.table)
library(zoo)
library(nnet)

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
         "weightinit_RX","weightinit_R",
         "log_vlplus1_init_RX","log_vlplus1_init_R","vl_0_1")


theX = c("new_cupd_1","log_cd4plus1") # LHS

theX_cov = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R") # RHS (from createind and createmis)
# for X1, use (X[0],X[-1])
theX_cov_1 = c("new_cupd_1","log_cd4plus1_RX","log_cd4plus1_R","log_cd4preBplus1_RX","log_cd4preBplus1_R") # RHS (from createind and createmis)

theY_ind = "outcome_2" # RHS (2 is engagement. outcome 3(xfer) and 4(death) are absorbing state, hence no previous states in 3 or 4); theY for LHS

theA = "onarv_cf_1"

theA_lag1 = paste0(theA,"_lag1") # lagdat uses "_lag1" as name suffix for creating lagged variables

theX_cov_lag1 = paste0(theX_cov,"_lag1")
theX_cov_1_lag1 = c(theX_cov_lag1,"log_cd4preBplus1_RX","log_cd4preBplus1_R")

theY_ind_lag1 = paste0(theY_ind,"_lag1")

load("20-01-12_newdat.RData") # defined in BartGcomp_19-12-13_Reclean_AMPATH_data script

#OUTCOME: 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
newdat$outcome = as.numeric(newdat$outcome) #IMPORTANT! de-factorize outcome

######################################################################################################################

#fitdat = newdat[time == 200,]
#table(fitdat$outcome)
#1     2     4 old data
#10903 73631  3664
#1     2     4 current data
#14070 70426  3652 
#1     2     4 20-01-12
#6497 67462  2697
#fitdat = newdat[time == 400,]
#table(fitdat$outcome)
#1     2     3     4 old data
#26799 52493   639   871
#1     2     3     4 current data
#29423 49744   758   865 
#1     2     3     4 20-01-12
#21421 47982   766   812 
#####################################
seed = 1
nsub =   max(newdat$myID)

### Separate data into training set and test set

ntr = 50000 # floor(nsub * 0.7) # sample size of training set

set.seed(seed)
samp = sample(1:nsub, nsub, replace = F)

trid = samp[1:ntr]  
trid = sort(trid)
trd = copy(newdat[myID %in% trid,])

#teid = samp[(ntr+1):(ntr+nte)]  
teid = samp[(ntr+1):nsub]  
teid = sort(teid)
ted = copy(newdat[myID %in% teid,])

### Restrict to t = 1
t = 1

trd = trd[time==200*t,]

ted = ted[time==200*t,]

# rename outcome to y
yind = which(colnames(trd)==theY)
colnames(trd)[yind] = "y"
yind = which(colnames(ted)==theY)
colnames(ted)[yind] = "y"
theY = "y"

print(dim(trd))
print(dim(ted))

# [S1 | V, A0, X1, L1]
fmlfilter = function(fitdat, RHS){
  todelete  = c()
  for(i in 1:length(RHS)){
    a = fitdat[,get(RHS[i])]
    if(length(unique(a))==1) todelete = c(todelete,i)
  }
  if(length(todelete)>0) RHS = RHS[-todelete]
  return(RHS)
}

RHS = c(theV,theA_lag1,theX_cov)
RHS = fmlfilter(trd,RHS)
fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))

#####################################
### Build models based on trd (training data)
#parell = TRUE # TRUE if on server, FALSE otherwise
#ncore = 5 # NULL if parall = FALSE

#maxt= 11 # for identifiability of parameters. No xfer or death in the observed data at t=15

#-------------------------------

#maxt= 10 # cannot fit for t = 11, sample size not fittable
nd = 30000 # posterior draws
nb = 50000 # burn in
nt = 100 # number of trees
p=length(unique(trd$y)) # number of outcome levels

CatSeed = 888888

if(sim ==  1){KD = F; nuC = 1; Vm = diag(2); by = "1"; da = TRUE}
if(sim ==  2){KD = F; nuC = 1; Vm = diag(2); by = "1"; da = FALSE}
if(sim ==  3){KD = T; nuC = 1; Vm = diag(2); by = "1"}

if(sim ==  4){KD = F; nuC = 1; Vm = diag(2); by = "2"; da = TRUE}
if(sim ==  5){KD = F; nuC = 1; Vm = diag(2); by = "2"; da = FALSE}
if(sim ==  6){KD = T; nuC = 1; Vm = diag(2); by = "2"}

if(sim ==  7){KD = F; nuC = 1; Vm = diag(2); by = "4"; da = TRUE}
if(sim ==  8){KD = F; nuC = 1; Vm = diag(2); by = "4"; da = FALSE}
if(sim ==  9){KD = T; nuC = 1; Vm = diag(2); by = "4"}

##############################################################
##############################################################
### summary plot of newdat outcome
if(0){
  a = newdat$outcome
  levels(a) = c("Disengaged","Engaged","Transferred out","Died")
  a = relevel(a, "Engaged")
  counts <- table(a,newdat$time)
  barplot(counts, main="Observed Outcome Status",
          xlab="Days since baseline", ylab="Number of observations", col=1:length(levels(a)),
          legend = rownames(counts))
  
  library(ggplot2)
  cPalette <- c("#0072B2","#E69F00","#FF3333","#339966")
  pathto = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft3\\Accuracy Comparison - Diag Plots\\19-03-05\\Plots"
  
  tab = data.frame(Status = factor( rep(rownames(counts), each = ncol(counts)), levels = c("Died","Transferred out","Disengaged","Engaged") ),
                   X = rep(as.numeric(colnames(counts)), nrow(counts)),
                   Y = c(t(counts)))
  ggplot(data = tab, aes(x=X, y=Y, fill=Status)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values = c("Died"=cPalette[1],"Transferred out"=cPalette[2],"Disengaged"=cPalette[3],"Engaged"=cPalette[4])) +
    ggtitle("Observed Outcome Status")+theme(plot.title = element_text(hjust = 0.5))+
    xlab("Days since baseline") + ylab("Number of observations")
  ggsave(paste0(pathto,"\\ObaEng.png"),width = 5.5, height = 3.5)
}
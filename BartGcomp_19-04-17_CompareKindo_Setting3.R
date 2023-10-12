### AMPATH

#datpath = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\MPBART Computing\\UCI Contraceptive Choice"
#d = read.table(paste0(datpath,"\\cmc.data"),header = FALSE, sep=",")
source("BartGcomp_19-03-10_Athena0_prepare.R")

t = 1

trd = trd[time==200*t,]

ted = ted[time==200*t,]

yind = which(colnames(trd)==theY)
colnames(trd)[yind] = "y"
yind = which(colnames(ted)==theY)
colnames(ted)[yind] = "y"
theY = "y"

print(dim(trd))
print(dim(ted))

# [S1 | V, A0, X1, L1]
RHS = c(theV,theA_lag1,theX_cov)
RHS = fmlfilter(trd,RHS)
fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))


### Multinomial Probit BART Settings

set.seed(99)
nd = 2000 # posterior draws
nb = 500 # burn in
nt = 100 # number of trees
p=length(unique(trd$y)) # number of outcome levels

# seed to start with
CatSeed = 888888



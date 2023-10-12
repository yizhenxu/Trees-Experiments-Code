#setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
#setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
#setwd("./Topic2")

screennum = 76684
### ZEUS script::

sink(paste0("output_",screennum,".txt")) #write all output to txt -- source('script.R',echo=T) to write input too

#devtools::install_github("yizhenxu/GcompBART")
library(GcompBART)
library(MNP)

source("Functions_new.R")

n = 5000 #number of observations
ntr = 4000 #number of train
p=3
nd = 2000
nb = 1000

ptm <- proc.time()

f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)

set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y=Ey+sigma*rnorm(n)
z = (y-min(y))/(max(y)-min(y))
quantile(z)
y = 1*(z<0.4)+ 2*(z>=0.4 & z<0.6) + 3*(z>=0.6)
table(y)


newdat = data.frame(x,y)
trd = newdat[1:ntr,]
ted = newdat[(ntr+1):n,]
table(trd$y)
table(ted$y)

fml = "y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10"

set.seed(99)
bmpfit = model_bart(as.formula(fml), data = trd, type = "multinomial",
                    Prior = list(nu = p-1+3, V = diag(p-1),
                                 ntrees = 200,
                                 kfac = 2,
                                 pswap = 0.1, pbd = 0.5, pb = 0.25,
                                 alpha = 0.95, beta = 2.0,
                                 nc = 100, minobsnode = 10),
                    Mcmc = list(sigma0 = diag(p-1), burn = nb, ndraws = nd, nSigDr = 1))

DiagPlot(bmpfit,1)
# train prediction
TM = bmpfit$samp_y
yM = matrix(rep(trd$y, nd),ncol = nd)

E = rep(NA, p+1)
for(k in 1:p){
  ind = which(yM==k)
  E[k] = mean(TM[ind] == yM[ind])
}
E[p+1] = mean(TM == yM)
names(E) = c(paste0("y=",1:p),"overall")
TE = E
TE

### prediction as the most frequent yhat
yhat = function(vec){
  names(which.max(table(vec))) 
}
TMY = apply(TM, 1, yhat)
mean(TMY == trd$y)

# test prediction
PM = predict_bart(bmpfit, ted)$samp_y
yM = matrix(rep(ted$y, nd),ncol = nd)
 
E = rep(NA, p+1)
for(k in 1:p){
  ind = which(yM==k)
  E[k] = mean(PM[ind] == yM[ind])
}
E[p+1] = mean(PM == yM)
names(E) = c(paste0("y=",1:p),"overall")
E

PMY = apply(PM, 1, yhat)
mean(PMY == ted$y)
#-------------------------------

ndraw = nd+nb 
nburn = nb
nthin = 0 # number of discarded samples between two valid samples
basey = "3" # disengagement as reference level
nchain = 1 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain

mnpfit = mnp_chains(as.formula(fml), trd, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE)

d = ted
m = mnpfit
fv = fmlvars(fml)
dat.pred = d[, fv$RHS[-1]]
dat.pred = cbind(1, dat.pred)
colnames(dat.pred)[1] = fv$RHS[1]

ylevels = as.numeric(unique(trd$y))
latenty = sort(ylevels[-which(ylevels==basey)])
nlatent = length(latenty)

# simulate from multivariate normal (mu,Sigma)

nhalfSig = nlatent*(nlatent - 1)/2 + nlatent
sigmas = m[, (ncol(m) - nhalfSig + 1): ncol(m)]
params = m[, 1:(ncol(m) - nhalfSig)]

BbyY = getParamByY(params,nlatent)

mu = lapply(1:nlatent, function(yl) tcrossprod(as.matrix(dat.pred), BbyY[[yl]]) )
mu_array = simplify2array(mu) # 3 dims: nsub x npost x nlatent

postY = lapply(1:npost, function(j) getYhat(j, nlatent, basey, latenty, mu_array, sigmas))
PM = simplify2array(postY)

yM = matrix(rep(d$y, nd),ncol = nd)

E1 = rep(NA, p+1)
for(k in 1:p){
  ind = which(yM==k)
  E1[k] = mean(PM[ind] == yM[ind])
}
E1[p+1] = mean(PM == yM)
names(E1) = c(paste0("y=",1:p),"overall")
E1

### prediction as the most frequent yhat
yhat = function(vec){
  names(which.max(table(vec))) 
}
PMY = apply(PM, 1, yhat)
mean(PMY == d$y)
#-------------------------------

save(trd, ted, bmpfit, mnpfit, file = paste0("Simulation_screen",screennum,".RData"))

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

#> E
#y=1       y=2       y=3   overall
#0.6489614 0.5091692 0.4736429 0.5445145
#> E1
#y=1       y=2       y=3   overall
#0.6918858 0.5407462 0.6085393 0.6086975

load("I:\\Topic2\\Simulation_screen76684_nsub50000.RData")
DiagPlot(bmpfit,1)
DiagPlot(bmpfit,0)

nlatent = 2
nsig = (nlatent^2 -nlatent)/2 + nlatent
q = (ncol(mnpfit) - nsig)/2

library(coda)

autocorr.plot(as.mcmc(mnpfit[, (1:q)*2-1 ]))
autocorr.plot(as.mcmc(mnpfit[, (1:q)*2 ]))
autocorr.plot(as.mcmc(mnpfit[, (2*q+1):ncol(mnpfit) ]))
par(mfrow=c(3,4))
traceplot(as.mcmc(mnpfit[, (1:q)*2-1 ]))
par(mfrow=c(3,4))
traceplot(as.mcmc(mnpfit[, (1:q)*2 ]))
par(mfrow=c(3,4))
traceplot(as.mcmc(mnpfit[, (2*q+1):ncol(mnpfit) ]))


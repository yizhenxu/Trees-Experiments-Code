#Different simulation setting, similar to the MPBART paper
#setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
#setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART")
#setwd("./Topic2");source(paste0("script_",screennum,".R"),echo=T)

### ZEUS script::

sink(paste0("output_",screennum,".txt")) #write all output to txt -- source('script.R',echo=T) to write input too

#devtools::install_github("yizhenxu/GcompBART")
library(GcompBART)
library(MNP)

source("Functions_new.R")

n = 1000 #number of observations
ntr = 500 #number of train
p=3
nd = 3000
nb = 100

ptm <- proc.time()

f1 = function(x){
  15*sin(pi*x[,1]*x[,2]) + 2*(x[,3]-.5)^2 - 10*x[,4] + 5*x[,5]
}

f2 = function(x,v){
  (x[,3]-.5)^3 - 5*x[,4]*x[,5] + 4*v
}

set.seed(99)
u=matrix(runif(n*5),n,5) #10 variables, only first 5 matter
v = runif(n)*2

Sig = matrix(c(1, 0.5, 0.5, 1), 2, 2)

z1 = f1(u); z2 = f2(u,v)
z = mvrnorm(n, c(0,0), Sig)
z = z + cbind(z1, z2)
y = 1*(z[,2] >= z[,1]) +1
y[z[,1] < 0 & z[,2] < 0] = 3
table(y)



newdat = data.frame(u,v,y)
colnames(newdat)[6] = "X6"
trd = newdat[1:ntr,]
ted = newdat[(ntr+1):n,]
table(trd$y)
table(ted$y)

fml = "y ~ X1 + X2 + X3 + X4 + X5 + X6"


yhat = function(vec){
  names(which.max(table(vec))) 
}

mnppred = function(d, m, fml){
  
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
  return(PM)
}

set.seed(99)
bmpfit = model_bart(as.formula(fml), data = trd, type = "multinomial",
                    Prior = list(nu = p, V = diag(p-1),
                                 ntrees = 200,
                                 kfac = 2,
                                 pswap = 0.1, pbd = 0.5, pb = 0.25,
                                 alpha = 0.95, beta = 2.0,
                                 nc = 100, minobsnode = 10),
                    Mcmc = list(sigma0 = diag(p-1), burn = nb, ndraws = nd, nSigDr = 1),
                    correction = 0)

#DiagPlot(bmpfit,1)

# train prediction
BTM = bmpfit$samp_y
yM = matrix(rep(trd$y, nd),ncol = nd)

BTE = rep(NA, p+1)
for(k in 1:p){
  ind = which(yM==k)
  BTE[k] = mean(BTM[ind] == yM[ind])
}
BTE[p+1] = mean(BTM == yM)
names(BTE) = c(paste0("y=",1:p),"overall")
BTE


# test prediction
BPM = predict_bart(bmpfit, ted)$samp_y
yM = matrix(rep(ted$y, nd),ncol = nd)

BPE = rep(NA, p+1)
for(k in 1:p){
  ind = which(yM==k)
  BPE[k] = mean(BPM[ind] == yM[ind])
}
BPE[p+1] = mean(BPM == yM)
names(BPE) = c(paste0("y=",1:p),"overall")
BPE

### prediction as the most frequent yhat
BTY = apply(BTM, 1, yhat)
BPY = apply(BPM, 1, yhat)
a = c(BTE[p+1], mean(BTY == trd$y), BPE[p+1], mean(BPY == ted$y))
#-------------------------------

ndraw = nd+nb 
nburn = nb
nthin = 0 # number of discarded samples between two valid samples
basey = "3" # disengagement as reference level
nchain = 1 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain

mnpfit = mnp_chains(as.formula(fml), trd, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE)


# train prediction
d = trd
MTM = mnppred(d, mnpfit, fml)
yM = matrix(rep(d$y, nd),ncol = nd)

MTE = rep(NA, p+1)
for(k in 1:p){
  ind = which(yM==k)
  MTE[k] = mean(MTM[ind] == yM[ind])
}
MTE[p+1] = mean(MTM == yM)
names(MTE) = c(paste0("y=",1:p),"overall")
MTE


# test prediction
d = ted
MPM = mnppred(d, mnpfit, fml)
yM = matrix(rep(d$y, nd),ncol = nd)

MPE = rep(NA, p+1)
for(k in 1:p){
  ind = which(yM==k)
  MPE[k] = mean(MPM[ind] == yM[ind])
}
MPE[p+1] = mean(MPM == yM)
names(MPE) = c(paste0("y=",1:p),"overall")
MPE

### prediction as the most frequent yhat
MTY = apply(MTM, 1, yhat)
MPY = apply(MPM, 1, yhat)
b = c(MTE[p+1], mean(MTY == trd$y), MPE[p+1], mean(MPY == ted$y))
#-------------------------------

tab = rbind(a,b)
rownames(tab) = c("MPBART", "MNP")
colnames(tab) = c("Train Sim", "Train Argmax", "Test Sim", "Test Argmax")
tab
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


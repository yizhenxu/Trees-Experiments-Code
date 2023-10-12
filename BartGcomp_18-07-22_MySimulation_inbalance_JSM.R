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
redo = 20

priorhere = list(nu = p, V = diag(p-1),
                 ntrees = 100,
                 kfac = 2,
                 pswap = 0.1, pbd = 0.5, pb = 0.25,
                 alpha = 0.95, beta = 2.0,
                 nc = 100, minobsnode = 10)

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

Btab = matrix(NA, redo+2, 4)
aa = rep(NA, redo)

set.seed(1)

for(i in 1:redo){
  ptm <- proc.time()
  bmpfit = model_bart(as.formula(fml), data = trd, type = "multinomial",
                      Prior = priorhere,
                      Mcmc = list(sigma0 = diag(p-1), burn = nb, ndraws = nd, nSigDr = 1),
                      correction = 0)
  aa[i] =  (proc.time()-ptm)[3]
  #DiagPlot(bmpfit,1)
  
  # train prediction
  #BTM = bmpfit$samp_y
  BTM = predict_bart(bmpfit, trd)$samp_treefit
  yM = matrix(rep(trd$y, nd),ncol = nd)
  
  Btab[i,1] = mean(BTM == yM)
  
  BTY = apply(BTM, 1, yhat)
  Btab[i,2] = mean(BTY == trd$y)
  
  # test prediction
  #BPM = predict_bart(bmpfit, ted)$samp_y
  BPM = predict_bart(bmpfit, ted)$samp_treefit
  yM = matrix(rep(ted$y, nd),ncol = nd)
  
  Btab[i,3] = mean(BPM == yM)
  
  BPY = apply(BPM, 1, yhat)
  Btab[i,4] = mean(BPY == ted$y)
  
}

Btab[redo+1,] = apply(Btab[1:redo,], 2, mean)
Btab[redo+2,] = apply(Btab[1:redo,], 2, sd)
#-------------------------------

ndraw = nd+nb 
nburn = nb
nthin = 0 # number of discarded samples between two valid samples
basey = "1" # disengagement as reference level
nchain = 1 # number of MCMC chains
npost = floor((ndraw - nburn) / (nthin+1)) * nchain


Mtab = matrix(NA, redo+2, 4)
i = 1
len = 100
bb = rep(0, len)
for (ct in 1:len) {
  set.seed(ct)
  #ERROR HANDLING
  mnpfit <- tryCatch(
    mnp_chains(as.formula(fml), trd, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE),
    error=function(e) e
  )
  
  if(inherits(mnpfit, "error")) next
  
  # train prediction
  d = trd
  MTM = mnppred(d, mnpfit, fml)
  yM = matrix(rep(d$y, nd),ncol = nd)
  
  Mtab[i,1] = mean(MTM == yM)
  
  MTY = apply(MTM, 1, yhat)
  Mtab[i,2] = mean(MTY == trd$y)
  
  # test prediction
  d = ted
  MPM = mnppred(d, mnpfit, fml)
  yM = matrix(rep(d$y, nd),ncol = nd)
  
  Mtab[i,3] = mean(MPM == yM)
  
  MPY = apply(MPM, 1, yhat)
  Mtab[i,4] = mean(MPY == ted$y)
  
  bb[i] = 1
  if(sum(bb) == redo) break
  
  i = i+1
} 


Mtab[redo+1,] = apply(Mtab[1:redo,], 2, mean)
Mtab[redo+2,] = apply(Mtab[1:redo,], 2, sd)

#-------------------------------

tab = rbind(Btab[redo+1,], Mtab[redo+1,])
rownames(tab) = c("MPBART", "MNP")
colnames(tab) = c("Train Sim", "Train Argmax", "Test Sim", "Test Argmax")
tab
#-------------------------------

save(trd, ted, bmpfit, mnpfit, Btab, Mtab, file = paste0("Simulation_screen",screennum,".RData"))
##################################################################################
DiagPlot(bmpfit,1)
##################################################################################
library(mailR)

#status = as.character(aa[3])
status = as.character(Sys.time())

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



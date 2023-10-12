
### Data 

n = 10000 #number of observations, ntr+nte
ntr = 5000 #number of train

#ptm <- proc.time()

f1 = function(x){
  15*sin(pi*x[,1]*x[,2]) + (x[,3]-.5)^2 - 10*x[,4] - 5*x[,5]
}

f2 = function(x,v){
  (x[,3]-.5)^3 - 20*x[,4]*x[,5] + 4*v
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



### Multinomial Probit BART Settings

set.seed(99)
nd = 2000 # posterior draws
nb = 500 # burn in
nt = 100 # number of trees
p=length(unique(trd$y)) # number of outcome levels

# seed to start with
CatSeed = 888888



################################################################
if(0){
  
  source("BartGcomp_19-03-10_Athena0_prepare.R")
  
  t = 6
  
  td = ted[time==200*t,]
  
  ty =  unlist(ted[time==200*t, theY, with = F])
  
  res = matrix(NA, ncol=4, nrow=4)
  
  ###
  
  fitdat = trd[time==200*t,]
  print(dim(fitdat))
  
  Bmod2 = vector("list", 2)
  baseyl = rep(basey,2)
  ntL = vector("list", 2)
  
  # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]
  
  RHS = c(theV, theA_lag1, theX_cov)
  RHS = fmlfilter(fitdat,RHS)
  
  fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
  
  seedvec = rep(CatSeed, 2) # seeds for (prevS = 0, prevS = 1)
  
}
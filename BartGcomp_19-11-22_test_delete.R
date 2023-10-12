remove.packages("GcompBART")
devtools::install_github("yizhenxu/GcompBART")
library(GcompBART)

n=2000;ntr=1000
### Data 

n = 10000 #number of observations, ntr+nte
ntr = 5000 #number of train

#ptm <- proc.time()

f1 = function(x){
  15*sin(pi*x[,1]*x[,2]) + (x[,3]-.5)^2 - 10*x[,4] - 5*x[,5]
}

f2 = function(x,v){
  (x[,3]-.5)^2 - x[,4]*x[,5] + 4*v
}

set.seed(99)
u=matrix(runif(n*5),n,5) #10 variables, only first 5 matter
v = runif(n)*2

Sig = matrix(c(1, 0.5, 0.5, 1), 2, 2)

z1 = f1(u); z2 = f2(u,v)
library(MASS)
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
nuC = 1

# seed to start with
CatSeed = 888888

source("BartGcomp_19-11-22_Functions.R")

ylev = 1:3
bmpfit0 = vector("list",length(ylev))
bmpfit1 = vector("list",length(ylev))

nd = 200;nb = 200

j=1
bmpfit1[[j]] = model_bart(as.formula(fml), data = trd, type = "multinomial",
                          base = ylev[j],
                          Prior = Prior_mult(p = p, nu_C = nuC, ntree = nt),
                          Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                          correction = FALSE, Kindo = TRUE)

bmpfit0[[j]] = model_bart(as.formula(fml), data = trd, type = "multinomial",
                          base = ylev[j],
                          Prior = Prior_mult(p = p, nu_C = nuC, ntree = nt),
                          Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                          correction = FALSE, Kindo = FALSE)
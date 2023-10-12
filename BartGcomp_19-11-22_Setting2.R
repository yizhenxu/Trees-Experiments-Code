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
#nb 500 nd 2000
#nb 2000 nd 3000
set.seed(99)
nd = 3000 # posterior draws
nb = 2000 # burn in
nt = 100 # number of trees
p=length(unique(trd$y)) # number of outcome levels
nuC = 1
Vm = diag(p-1)

# seed to start with
CatSeed = 888888

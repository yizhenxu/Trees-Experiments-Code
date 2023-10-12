### Data 

n = 10000 #number of observations, ntr+nte
ntr = 5000 #number of train

#ptm <- proc.time()
#setting1 (more balanced); setting2 (extreme unbalanced)

setting = NULL; Vm = NULL; KD = NULL; nuC = NULL; by = NULL; 
if(sim ==  1){setting = 1; KD = F; nuC = 1; Vm = diag(2); by = 1}
if(sim ==  2){setting = 1; KD = T; nuC = 1; Vm = diag(2); by = 1}
if(sim ==  3){setting = 2; KD = F; nuC = 1; Vm = diag(2); by = 1}
if(sim ==  4){setting = 2; KD = T; nuC = 1; Vm = diag(2); by = 1}

if(sim ==  5){setting = 1; KD = F; nuC = 1; Vm = diag(2); by = 2}
if(sim ==  6){setting = 1; KD = T; nuC = 1; Vm = diag(2); by = 2}
if(sim ==  7){setting = 2; KD = F; nuC = 1; Vm = diag(2); by = 2}
if(sim ==  8){setting = 2; KD = T; nuC = 1; Vm = diag(2); by = 2}

if(sim ==  9){setting = 1; KD = F; nuC = 1; Vm = diag(2); by = 3}
if(sim == 10){setting = 1; KD = T; nuC = 1; Vm = diag(2); by = 3}
if(sim == 11){setting = 2; KD = F; nuC = 1; Vm = diag(2); by = 3}
if(sim == 12){setting = 2; KD = T; nuC = 1; Vm = diag(2); by = 3}

if(sim == 13){setting = 1; KD = F; nuC = 3; Vm = matrix(c(1,-0.5,-0.5,1),2,2); by = 3}
if(sim == 14){setting = 1; KD = T; nuC = 3; Vm = matrix(c(1,-0.5,-0.5,1),2,2); by = 3}
if(sim == 15){setting = 2; KD = F; nuC = 3; Vm = matrix(c(1,-0.5,-0.5,1),2,2); by = 3}
if(sim == 16){setting = 2; KD = T; nuC = 3; Vm = matrix(c(1,-0.5,-0.5,1),2,2); by = 3}

if(sim == 17){setting = 1; KD = F; nuC = 3; Vm = matrix(c(1,0.5,0.5,1),2,2); by = 3}
if(sim == 18){setting = 1; KD = T; nuC = 3; Vm = matrix(c(1,0.5,0.5,1),2,2); by = 3}
if(sim == 19){setting = 2; KD = F; nuC = 3; Vm = matrix(c(1,0.5,0.5,1),2,2); by = 3}
if(sim == 20){setting = 2; KD = T; nuC = 3; Vm = matrix(c(1,0.5,0.5,1),2,2); by = 3}


print(paste0("Simulation Setting ", setting))
print(paste0("Kindo ", KD))
print(Vm)
print(paste0("Reference Level ", by))

if(setting == 1){
  f1 = function(x){
    15*sin(pi*x[,1]*x[,2]) + (x[,3]-.5)^2 - 10*x[,4] - 5*x[,5]
  }
  
  f2 = function(x,v){
    (x[,3]-.5)^3 - 20*x[,4]*x[,5] + 4*v
  }
}

if(setting == 2){
  f1 = function(x){
    15*sin(pi*x[,1]*x[,2]) + (x[,3]-.5)^2 - 10*x[,4] - 5*x[,5]
  }
  
  f2 = function(x,v){
    (x[,3]-.5)^2 - x[,4]*x[,5] + 4*v
  }
}

set.seed(datseed)
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

p=length(unique(trd$y)) # number of outcome levels

fml = "y ~ X1 + X2 + X3 + X4 + X5 + X6"

if(0){ # multinomial logistic
  require("nnet")
  test = multinom(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X1*X2+X1*X3+X1*X4+X1*X5+X1*X6+X2*X3+X2*X4+X2*X5+X2*X6+X3*X4+X3*X5+X3*X6+X4*X5+X4*X6+X5*X6+X1*X1+X2*X2+X3*X3+X4*X4+X5*X5+X6*X6+X1*X1*X1+X2*X2*X2+X3*X3*X3+X4*X4*X4+X5*X5*X5+X6*X6*X6+X1*X2*X3+X1*X2*X4+X1*X2*X5+X1*X2*X6+X1*X3*X4+X1*X3*X5+X1*X3*X6+X1*X4*X5+X1*X4*X6+X1*X5*X6+X2*X3*X4+X2*X3*X5+X2*X3*X6+X2*X4*X5+X2*X4*X6+X2*X5*X6+X3*X4*X5+X3*X4*X6+X3*X5*X6+X4*X5*X6, data = trd)
  pp = predict(test, newdata = ted, "probs")

  ymode = function(vec){
    colnames(pp)[which.max(vec)]
  }
  yp = apply(pp, 1, ymode)
  mean(ted$y == as.numeric(yp))
}
### Multinomial Probit BART Settings
#nb 500 nd 2000
#nb 2000 nd 3000
#nb 15000 nd 5000
#nb 20000 nd 10000
#nb 30000 nd 20000
#nb 50000 nd 30000
nd = 30000 # posterior draws
nb = 50000 # burn in
nt = 100 # number of trees

# seed to start with
CatSeed = 888888

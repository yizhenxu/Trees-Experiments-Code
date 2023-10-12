devtools::install_github("yizhenxu/allBART")


library(allBART)


##simulate data (example from Friedman MARS paper)
f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
n = 5000 #number of observations
set.seed(8)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y=Ey+sigma*rnorm(n)

dd = data.frame(x,y)
dat = dd[1:floor(0.7*n),]
tdat = dd[(floor(0.7*n)+1):n,]

lmFit = lm(y ~ X1+X2, data = dd)
summary(lmFit)$sigma# 4.63

ndraw = 1000
nburn = 100
ntree = 100
ncut = 100
nsd = 20 # max repeated draw of Sigma for Mpbart
seed = 1

set.seed(seed)
fml = "y ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10"
bfit = BART_call(as.formula(fml), data = dat, test.data = tdat,
          Prior = list(nu = 3, sigq = 0.9,
                       ntrees = ntree,
                       kfac = 2,
                       pswap = 0.25, pbd = 0.5, pb = 0.25,
                       alpha = 0.95, beta = 2.0,
                       nc = ncut, minobsnode = 10),
          Mcmc = list(burn = nburn, ndraws = ndraw))


dyhat = bfit$samp_train
tyhat = bfit$samp_test
 
dmse = (dyhat - matrix(rep(dat$y,ndraw),ncol=ndraw))^2
dmse = apply(dmse, 2, mean)

tmse = (tyhat - matrix(rep(tdat$y,ndraw),ncol=ndraw))^2
tmse = apply(tmse, 2, mean)

#plot(dmse, ylim= range(c(tmse, dmse)), )
#points(tmse, col= "red")


bfit1 = BART_call(as.formula(fml), data = dat, test.data = tdat,
                 Prior = list(nu = 3, sigq = 0.9,
                              ntrees = ntree,
                              kfac = 2,
                              pswap = 0, pbd = 1, pb = 0.5,
                              alpha = 0.95, beta = 2.0,
                              nc = ncut, minobsnode = 10),
                 Mcmc = list(burn = nburn, ndraws = ndraw))


dyhat1 = bfit1$samp_train
tyhat1 = bfit1$samp_test

dmse1 = (dyhat1 - matrix(rep(dat$y,ndraw),ncol=ndraw))^2
dmse1 = apply(dmse1, 2, mean)

tmse1 = (tyhat1 - matrix(rep(tdat$y,ndraw),ncol=ndraw))^2
tmse1 = apply(tmse1, 2, mean)

#points(dmse1, col = "green")
#points(tmse1, col = "blue")

rng = range(c(tmse, dmse))
rng[2] = rng[2] + 0.1

png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\MSE_compare.png"),width = 600, height = 400)
png(paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-12-13_newData_multinomial_regression\\MSE_compare.png"),width = 600, height = 400)

plot(1, type="n", xlab="", ylab="Smoothed MSE", xlim=c(0, ndraw), ylim=rng)

lines(lowess(1:ndraw, dmse), col="red")
lines(lowess(1:ndraw, tmse), col="red", lty = "dashed")
lines(lowess(1:ndraw, dmse1), col="blue")
lines(lowess(1:ndraw, tmse1), col="blue", lty = "dashed")
legend(450, 3.3, legend=c("Grow/Prune/Swap/Change Train MSE", "Grow/Prune/Swap/Change Test MSE",
                         "Grow/Prune Train MSE", "Grow/Prune Test MSE"),
       col=c("red","red", "blue","blue"), lty=c(1,2,1,2), cex=1)

dev.off()
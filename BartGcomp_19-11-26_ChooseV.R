#install.packages("LaplacesDemon")
library("LaplacesDemon")

simy = function(nu_C){
  S = matrix(c(1,-0.5,-0.5,1),2,2)
  nu = 2+nu_C
  
  N=5000
  y = rep(NA, N)
  
  for(i in 1:N){
    
    x = rinvwishart(nu, S)
    alpha = (x[1,1]+x[2,2])/2
    y[i] = x[1,2]/alpha
  }
  return(y)
  
}

y1 = simy(1)
y2 = simy(3)
par(mfrow = c(1,2))
hist(y1, breaks = 20)
abline(v = mean(y1),col = "red")
hist(y2, breaks = 20)
abline(v = mean(y2),col = "red")


install.packages("TruncatedNormal")
library("TruncatedNormal")

d = 2;rho=0.5
u = c(0,0); l = rep(-Inf,d)
Sig=matrix(c(1,rho,rho,1),d,d)
mvNcdf(l,u,Sig,1e4) #prob # prob

1/4 + asin(rho)/(2*pi)

y = rep(NA, length(seq(-1,1,.05)))
i=1
for(rho in seq(-0.95,0.95,.05)){
  u = c(3,0); l = rep(-Inf,d)
  Sig=matrix(c(1,rho,rho,1),d,d)
  y[i]=mvNcdf(l,u,Sig,1e4)$prob #prob # prob
  i = i+1
}
plot(y);abline(h = 0)

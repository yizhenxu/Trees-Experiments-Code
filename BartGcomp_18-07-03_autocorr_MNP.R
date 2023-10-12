##################################################################################
library(coda)

m = mod[[3]]
nlatent = 2
nsig = (nlatent^2 -nlatent)/2 + nlatent
q = (ncol(m) - nsig)/2



autocorr.plot(as.mcmc(m[, (1:q)*2-1 ]))
autocorr.plot(as.mcmc(m[, (1:q)*2 ]))
autocorr.plot(as.mcmc(m[, (2*q+1):ncol(m) ]))

par(mfrow=c(3,4))
traceplot(as.mcmc(m[, (1:q)*2-1 ]))
par(mfrow=c(3,4))
traceplot(as.mcmc(m[, (1:q)*2 ]))
par(mfrow=c(3,4))
traceplot(as.mcmc(m[, (2*q+1):ncol(m) ]))


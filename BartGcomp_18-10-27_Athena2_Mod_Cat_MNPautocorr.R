library(coda)

m = mod
nlatent = 3
nsig = (nlatent^2 -nlatent)/2 + nlatent
q = (ncol(m) - nsig)/nlatent

autocorr.plot(as.mcmc(m[, (1:q)*nlatent-1 ]))
autocorr.plot(as.mcmc(m[, (1:q)*nlatent ]))
autocorr.plot(as.mcmc(m[, (nlatent*q+1):ncol(m) ]))

ln = length((1:q)*nlatent-1)
par(mfrow=c(3,floor(ln/3)))
traceplot(as.mcmc(m[, (1:q)*nlatent-2 ]))
par(mfrow=c(3,floor(ln/3)))
traceplot(as.mcmc(m[, (1:q)*nlatent-1 ]))
par(mfrow=c(3,floor(ln/3)))
traceplot(as.mcmc(m[, (1:q)*nlatent ]))

ln = length((nlatent*q+1):ncol(m))
par(mfrow=c(3,floor(ln/3)))
traceplot(as.mcmc(m[, (nlatent*q+1):ncol(m) ]))

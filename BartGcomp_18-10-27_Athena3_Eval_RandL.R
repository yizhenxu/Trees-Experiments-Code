#setwd("./Topic2")
##############################################################
### Third Part
### Evaluation plot  of simulated time-varying confounders R and L
### solid line as posterior mean of posterior estimations at each time point
### points as observed mean 
### shade as  95% posterior predictive intervals of the actual measurements
### on the test set
sink("check.txt") #write all output to txt -- source('script.R',echo=T) to write input too

source("BartGcomp_18-10-27_Athena0_prepare.R")

load("Mod18-10-27.RData")

summRL = function(M){# M: testn rows, nd cols
  vecp = apply(M, 2, function(x) mean(x))
  res = mean(vecp)
  
  svecp = sort(vecp)
  lp = length(svecp)
  ci = svecp[c(ceiling(lp*0.05), floor(lp*0.95))]
  res = c(res,ci)
  
  return(res)
}

res = array(NA, dim = c(4,maxt,2)) 
# the first dim: y, mean post, 5% post, 95% post
# the third dim: R, L

## R
for(t in 1:maxt){
  
  print(t)
  
  d = ted[time==200*t,]
  
  fv = fmlvars(modfml[[3*(t-1)+1]], itc = 0)
  y =  unlist(ted[time==200*t,fv$LHS,with = F])
  
  dat.pred = ted[time==200*t,fv$RHS,with=F]
  
  PM = predict_bart(mod1[[t]], dat.pred)$samp_y
  
  res[2:4, t, 1] = summRL(PM)
  
  res[1, t, 1] = mean(y)
  
}

## L
for(t in 1:maxt){
  
  print(t)
  
  d = ted[time==200*t,]
  
  fv = fmlvars(modfml[[3*(t-1)+2]], itc = 0)
  y =  unlist(ted[time==200*t,fv$LHS,with = F])
  
  dat.pred = ted[time==200*t,fv$RHS,with=F]
  
  upd = unlist(ted[time==200*t, theX[1],with = F])
  ind = which(upd==1)
  y = y[ind]
  dat.pred = dat.pred[ind,]
  PM = predict_bart(mod2[[t]], dat.pred)$samp_y
  
  res[2:4, t, 2] = summRL(PM)
  
  res[1, t, 2] = mean(y)
  
}
save(res, file = "resRL.RData")


sink()
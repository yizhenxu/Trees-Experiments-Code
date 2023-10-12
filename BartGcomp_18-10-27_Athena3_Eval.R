#sink(paste0("output_",screennum,".txt")) #write all output to txt -- source('script.R',echo=T) to write input too

source("BartGcomp_18-10-27_Athena0_prepare.R")

load("Mod18-10-27.RData")


##############################################################
### First Part
###BART-PBART-MPBART
yhat = function(vec){
  names(which.max(table(vec))) 
}

newAcc_bart = function(mod1, mod2, mod3_D, mod3_E, modfml, ted, times, types, intypes, maxt, basey1, basey){
  npost = mod[[1]]$ndraws
  lenmod = length(times)
  nummod = length(types)
  
  E = vector("list", nummod)
  
  for(i in 1:nummod){
    if(types[i] == "mult"){
      E[[i]] = matrix(NA,8,maxt)
      rownames(E[[i]]) = c("y=1","y=2","y=3","y=4","average","mode","tree_average","tree_mode")
      colnames(E[[i]]) = paste0("t=",1:maxt)
    } else if(types[i] == "binary"){
      E[[i]] = matrix(NA,2,maxt)
      rownames(E[[i]]) = c("average","mode")
    } else {
      E[[i]] = rep(NA, maxt)
    }
  }
  
  for(ct in 1:lenmod){
    
    t = times[ct]
    
    fv = fmlvars(modfml[[ct]], itc = 0)
    y =  unlist(ted[time==200*t,fv$LHS,with = F])
    
    if(ct %% 3 == 0){
      
      fitdat = trd[time==200*t,]
      d = ted[time==200*t,]
      
      if(t == 1){
        PB = predict_bart(mod3_E[[t]], d)
        PM = PB$samp_y
        PMT = PB$samp_treefit
      } else{
        PB = predict_bart(mod3_D[[t-1]], d)
        PM = PB$samp_y
        PMT = PB$samp_treefit
        
        PB1 = predict_bart(mod3_E[[t]], d)
        PM1 = PB1$samp_y
        PMT1 = PB1$samp_treefit
        
        ind = which(d[,theY_ind_lag1,with = F]==1)#check
        PM[ind,] = PM1[ind,]
        PMT[ind,] = PMT1[ind,]
      }
      
      
    } else {#cont or binary
      dat.pred = ted[time==200*t,fv$RHS,with=F]
      
      if(types[intypes[ct]] == "continuous"){ 
        upd = unlist(ted[time==200*t, theX[1],with = F])
        ind = which(upd==1)
        y = y[ind]
        dat.pred = dat.pred[ind,]
        PM = predict_bart(mod2[[t]], dat.pred)$samp_y
      } else {
        PM = predict_bart(mod1[[t]], dat.pred)$samp_y
      }
      
      
    }
    
    yM = matrix(rep(y, npost),ncol = npost)
    
    if( types[intypes[ct]] == "binary"){
      PY = apply(PM, 1, yhat)
      E[[intypes[ct]]][1,t] = mean(PM == yM)
      E[[intypes[ct]]][2,t] = mean(PY == y)
      
    } else if(types[intypes[ct]] == "mult"){
      
      for(k in 1:4){
        
        E[[intypes[ct]]][k,t] = mean(((PM != k) & (yM != k)) | (PM==yM))
        
      }
      E[[intypes[ct]]][5,t] = mean(PM == yM)
      
      PY = apply(PM, 1, yhat)
      E[[intypes[ct]]][6,t] = mean(PY == y)
      
      E[[intypes[ct]]][7,t] = mean(PMT == yM)
      
      PY = apply(PMT, 1, yhat)
      E[[intypes[ct]]][8,t] = mean(PY == y)
      
    } else {
      
      E[[intypes[ct]]][t] = mean((PM - yM)^2)
      
    }
    
    
    print(ct)
  }
  #res = list(Pred = PM, Y = ylist, Accuracy = E)
  return(E)
}



# Bayesian G Comp
times = rep(1:maxt,each = 3) # time interval index for each ct model 
types = c("binary","cont","mult") # how many and which type should each E[[i]] be
intypes = rep(1:3, maxt) # which E[[i]] to go into

allres = newAcc_bart(mod1, mod2, mod3_D, mod3_E, modfml, ted, times, types, intypes, maxt, basey1, basey)

##############################################################
### Second Part
### Evaluation plot stratified by outcome category
### solid line as posterior mean of proportions at each time point
### points as observed proportions 
### shade as  95% posterior predictive intervals of the proportions
### on the test set

#M = matrix(sample(1:4,36,replace=T),4,9)
#matrix(unlist(lapply(1:4, function(s) summ(M,s,1))), ncol=4)
#op = 1: return with 95% CI, 0 otherwise

source("BartGcomp_18-10-27_Athena0_prepare.R")

load("Mod18-10-27.RData")

summ = function(M, s, op){# M: testn rows, nd cols
  vecp = apply(M, 2, function(x) mean(x==s))
  res = mean(vecp)
  
  if(op){
    svecp = sort(vecp)
    lp = length(svecp)
    ci = svecp[c(ceiling(lp*0.05), floor(lp*0.95))]
    res = c(res,ci)
  }
  
  return(res)
}

taby = function(y){
  return(unlist(lapply(1:4, function(s) mean(y==s))))
}

res = array(NA, dim = c(4,maxt,4)) 
# the first dim: y, mean post, 5% post, 95% post
# the third dim: eng, diseng, xfer, death
for(t in 1:maxt){
  
  print(t)
  
  d = ted[time==200*t,]
  
  fv = fmlvars(modfml[[3*t]], itc = 0)
  y =  unlist(ted[time==200*t,fv$LHS,with = F])
  
  if(t == 1){
    PB = predict_bart(mod3_E[[t]], d)
    PM = PB$samp_y
    #PMT = PB$samp_treefit
  } else{
    PB = predict_bart(mod3_D[[t-1]], d)
    PM = PB$samp_y
    #PMT = PB$samp_treefit
    
    PB1 = predict_bart(mod3_E[[t]], d)
    PM1 = PB1$samp_y
    #PMT1 = PB1$samp_treefit
    
    ind = which(d[,theY_ind_lag1,with = F]==1)#check
    PM[ind,] = PM1[ind,]
    #PMT[ind,] = PMT1[ind,]
  }
  #yM = matrix(rep(y, npost),ncol = npost)
  res[2:4, t, ] = matrix(unlist(lapply(1:4, function(s) summ(PM,s,1))), ncol=4)
  
  res[1, t, ] = taby(y)

}
save(res, file = "res.RData")

#maxt = 10
pd = as.data.frame(matrix(NA, nrow = 4*maxt, ncol = 6))

colnames(pd) = c("Outcome", "Time", "obs", "Posteriors", "p5", "p95")
pd[,1] = rep(c("Engaged","Disengaged","Transferred","Died"), each = maxt)
pd[,2] = rep(1:maxt, 4)
pd[,3:6] = rbind(t(res[,,1]),t(res[,,2]),t(res[,,3]),t(res[,,4])) 
pd$Outcome = factor(pd$Outcome, levels = c("Engaged","Disengaged","Transferred","Died"))

library("ggplot2")

p=ggplot(pd, aes(Time)) + 
  geom_line(aes(y=Posteriors),colour="red") + 
  geom_point(aes(y=obs), colour="blue") +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.5)+  facet_wrap(~ Outcome, scales="free_y")

plot(p)

#facet_grid?


##############################################################
### Third Part
### Evaluation plot  of simulated time-varying confounders R and L
### solid line as posterior mean of posterior estimations at each time point
### points as observed mean 
### shade as  95% posterior predictive intervals of the actual measurements
### on the test set

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



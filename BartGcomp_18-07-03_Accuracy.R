#screennum = 76185
#load("I:\\Topic2\\B_PB_MNP_screen54559_ct33_nsub50000.RData")

### Evaluation --- per model on test set
# types and intypes give the evalulation output format
Eval_permod_bpbmnp = function(mod, modfml, ted, times, types, intypes, maxt){
  npost = mod[[1]]$ndraws
  lenmod = length(mod)
  nummod = length(types)
  
  E = vector("list", nummod)
  
  for(i in 1:nummod){
    if(types[i] == "mult"){
      E[[i]] = matrix(NA,5,maxt)
      rownames(E[[i]]) = c("y=1","y=2","y=3","y=4","overall")
      colnames(E[[i]]) = paste0("t=",1:maxt)
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
      ytab = table(fitdat$outcome)
      ylevels = as.numeric(names(ytab)[ytab > 0])
      p = length(ylevels) # number of outcome categories
      m = mod[[ct]]
      
      if(p>2){#do MNP
        fv = fmlvars(modfml[[ct]])
        dat.pred = ted[time==200*t,fv$RHS[-1],with=F]
        dat.pred = cbind(1, dat.pred)
        colnames(dat.pred)[1] = fv$RHS[1]
        
        latenty = sort(ylevels[-which(ylevels==basey)])
        nlatent = length(latenty)
        
        # simulate from multivariate normal (mu,Sigma)
        
        nhalfSig = nlatent*(nlatent - 1)/2 + nlatent
        sigmas = m[, (ncol(m) - nhalfSig + 1): ncol(m)]
        params = m[, 1:(ncol(m) - nhalfSig)]
        
        BbyY = getParamByY(params,nlatent)
        
        mu = lapply(1:nlatent, function(yl) tcrossprod(as.matrix(dat.pred), BbyY[[yl]]) )
        mu_array = simplify2array(mu) # 3 dims: nsub x npost x nlatent
        
        postY = lapply(1:npost, function(j) getYhat(j, nlatent, basey, latenty, mu_array, sigmas))
        PM = simplify2array(postY)
        remove(postY)
        
      } else {#multinomial reduces to binary
        PM = predict_bart(mod[[ct]], dat.pred)$samp_y
        
        PM[which(PM==0)] = ylevels[1]
        PM[which(PM==1)] = ylevels[2]
      }
    } else {#cont or binary
      dat.pred = ted[time==200*t,fv$RHS,with=F]
      
      if(mod[[ct]]$type == "continuous"){ 
        upd = unlist(ted[time==200*t, theX[1],with = F])
        ind = which(upd==1)
        y = y[ind]
        dat.pred = dat.pred[ind,]
      }
      
      PM = predict_bart(mod[[ct]], dat.pred)$samp_y
      
    }
    
    
    yM = matrix(rep(y, npost),ncol = npost)
    
    if( types[intypes[ct]] == "binary"){
      
      E[[intypes[ct]]][t] = mean(PM == yM)
      
    } else if(types[intypes[ct]] == "mult"){
      
      for(k in 1:4){
        ind = which(yM==k)
        E[[intypes[ct]]][k,t] = mean(PM[ind] == yM[ind])
      }
      E[[intypes[ct]]][5,t] = mean(PM == yM)
      
    } else {
      
      E[[intypes[ct]]][t] = mean((PM - yM)^2)
      
    }
    
    
    print(ct)
  }
  #res = list(Pred = PM, Y = ylist, Accuracy = E)
  res = E
  return(res)
}



# Bayesian G Comp
times = rep(1:maxt,each = 3) # time interval index for each ct model 
types = c("binary","cont","mult") # how many and which type should each E[[i]] be
intypes = rep(1:3, maxt) # which E[[i]] to go into

allres = Eval_permod_bpbmnp(mod, modfml, ted, times, types, intypes, maxt)

##################################################################################
library(coda)

m = Smod[[1]]
m = mod[[30]]
nlatent =2
nsig = (nlatent^2 -nlatent)/2 + nlatent
q = (ncol(m) - nsig)/nlatent



autocorr.plot(as.mcmc(m[, (1:q)*2-1 ]))
autocorr.plot(as.mcmc(m[, (1:q)*2 ]))
autocorr.plot(as.mcmc(m[, (2*q+1):ncol(m) ]))
par(mfrow=c(3,7))

traceplot(as.mcmc(m[, (1:q)*nlatent-1 ]))
par(mfrow=c(3,4))
traceplot(as.mcmc(m[, (1:q)*nlatent ]))
par(mfrow=c(3,4))
traceplot(as.mcmc(m[, (nlatent*q+1):ncol(m) ]))


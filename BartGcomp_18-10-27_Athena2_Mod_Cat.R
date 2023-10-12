# NOTE!!: basey="2" was used for all models, regardless of fitting MNP using "1" sometimes
# Hence, initial w and sig were not the ones that were planned to be inserted
# This was a code typo, it somehow worked... but the convergence cannot be compared
sink(paste0("output_",screennum,".txt")) #write all output to txt -- source('script.R',echo=T) to write input too

source("BartGcomp_18-10-27_Athena0_prepare.R")

baseyl = rep(0,2)

if(t==1){
  
  t=1
  fitdat = trd[time==200*t,]
  # [S1 | V, A0, X1, L1]
  
  ytab = table(fitdat$outcome)
  ylevels = as.numeric(names(ytab)[ytab > 0])
  p = length(ylevels) # number of outcome categories
  
  RHS = c(theV,theA_lag1,theX_cov)
  RHS = fmlfilter(fitdat,RHS)
  fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
  
  mod = mnp_chains(fml, fitdat, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE)
  ###
  fv = fmlvars(fml); d = fitdat
  dat.pred = d[, fv$RHS[-1], with=F]
  dat.pred = cbind(1, dat.pred)
  colnames(dat.pred)[1] = fv$RHS[1]
  
  latenty = sort(ylevels[-which(ylevels==basey)])
  nlatent = length(latenty)
  
  # simulate from multivariate normal (mu,Sigma)
  
  nhalfSig = nlatent*(nlatent - 1)/2 + nlatent
  sigmas = mod[, (ncol(mod) - nhalfSig + 1): ncol(mod)]
  params = mod[, 1:(ncol(mod) - nhalfSig)]
  
  BbyY = getParamByY(params,nlatent)
  
  mu = lapply(1:nlatent, function(yl) tcrossprod(as.matrix(dat.pred), BbyY[[yl]]) )
  mu_array = simplify2array(mu) # 3 dims: nsub x npost x nlatent
  
  m = mu_array[,npost,]
  sig = getSymat(sigmas[npost,], nlatent)
  w = mvrnorm(nrow(m), rep(0, nlatent), sig) + m
  
  #save(mod,w,m,sig, file = "tmp.RData")
  
  ###################################################
  ## trd
  
  ## diagnostic for all 
  Bmod2 = model_bart(as.formula(fml), data = fitdat, type = "multinomial2",
                     base = basey,
                     Prior = Prior_mult(p),
                     Mcmc = Mcmc_mult(p),
                     correction = FALSE,
                     w = w, sig = sig, fitMNP = nf)
  ###################################################
} else {
  ##############################################################
  
  
  fitdat = trd[time==200*t,]
  Bmod2 = vector("list",2)
  
  # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]
  
  #RHS = c(theV, theA_lag1, theX_cov, theY_ind_lag1)
  RHS = c(theV, theA_lag1, theX_cov)
  RHS = fmlfilter(fitdat,RHS)
  fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
  
  for(prevS in c(0,1)){ #### S_{t-1} = 0 - disengaged, 1 - engaged
    
    ind = prevS+1
    
    ## trd
    d = fitdat[get(theY_ind_lag1) == prevS,]
    ytab = table(d$outcome)
    ylevels = as.numeric(names(ytab)[ytab > 0])
    p = length(ylevels) # number of outcome categories
    
    #mod = mnp_chains(fml, d, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE)
    
    mod =  tryCatch(mnp_chains(fml, d, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE),
                    error = function(e) NULL )
    if(length(mod)==0){
      mod =  tryCatch(mnp_chains(fml, d, "1", ndraw, nburn, nthin, nchain, doparallel = FALSE),
                      error = function(e) NULL )
      if(length(mod)==0){
        baseyl[ind] = NA
        doit = 0
      } else {
        baseyl[ind] = "1"
        doit = 1
      }
    } else {
      baseyl[ind] = basey
      doit = 1
    }
    
    if(doit){
      ###
      fv = fmlvars(fml);
      dat.pred = d[, fv$RHS[-1], with=F]
      dat.pred = cbind(1, dat.pred)
      colnames(dat.pred)[1] = fv$RHS[1]
      
      latenty = sort(ylevels[-which(ylevels==basey)])
      nlatent = length(latenty)
      
      # simulate from multivariate normal (mu,Sigma)
      
      nhalfSig = nlatent*(nlatent - 1)/2 + nlatent
      sigmas = mod[, (ncol(mod) - nhalfSig + 1): ncol(mod)]
      params = mod[, 1:(ncol(mod) - nhalfSig)]
      
      BbyY = getParamByY(params,nlatent)
      
      mu = lapply(1:nlatent, function(yl) tcrossprod(as.matrix(dat.pred), BbyY[[yl]]) )
      mu_array = simplify2array(mu) # 3 dims: nsub x npost x nlatent
      
      m = mu_array[,npost,]
      sig = getSymat(sigmas[npost,], nlatent)
      w = mvrnorm(nrow(m), rep(0, nlatent), sig) + m
      
      ## diagnostic for all 
      Bmod2[[ind]] = model_bart(as.formula(fml), data = d, type = "multinomial2",
                                base = basey,
                                Prior = Prior_mult(p),
                                Mcmc = Mcmc_mult(p),
                                correction = FALSE,
                                w = w, sig = sig, fitMNP = nf)
      
    } #doit
  } # prevS
  
}

save(baseyl,fml, Bmod2, file = paste0("Cat_",t,".RData"))

sink()
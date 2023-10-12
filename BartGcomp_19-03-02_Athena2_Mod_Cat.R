args <- commandArgs(trailingOnly=TRUE)
t = eval(parse(text=args[1]))

# t - time interval index

##############################################################
source("BartGcomp_19-03-02_Athena0_prepare.R")

if(t==1){

  fitdat = trd[time==200*t,]
  print(dim(fitdat))
  
  # [S1 | V, A0, X1, L1]
  RHS = c(theV,theA_lag1,theX_cov)
  RHS = fmlfilter(fitdat,RHS)
  fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
  
  ytab = table(fitdat$outcome)
  print(ytab)
  
  ylevels = as.numeric(names(ytab)[ytab > 0])
  p = length(ylevels) # number of outcome categories
  
  yprop = ytab/sum(ytab)
  print(rbind(ytab,yprop))
  
  Yname = theY
  Xname = vector("list", p-1)
  for(i in 1:(p-1)){
    Xname[[i]] = RHS
  }
  
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
  
  ###################################################
  ## trd
  
  ## diagnostic for all 
  Bmod2 = model_bart(Yname = Yname, Xname = Xname, data = fitdat, type = "multinomial",
                            base = basey,
                            Prior = Prior_mult(p = p, ntree = rep(nt, p-1) ),
                            Mcmc = Mcmc_mult(p = p, w0 = w, sig0 = sig),
                            correction = FALSE,
                            fitMNP = nf, bylatent = blt)
  
  ##############################################################
} else {
  ##############################################################

  fitdat = trd[time==200*t,]
  print(dim(fitdat))
  
  Bmod2 = vector("list", 2)
  baseyl = rep(0,2)
  ntL = vector("list", 2)

  # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]
  
  RHS = c(theV, theA_lag1, theX_cov)
  RHS = fmlfilter(fitdat,RHS)
  
  fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
  
  for(prevS in c(0,1)){ #### S_{t-1} = 0 - disengaged, 1 - engaged
    
    ind = prevS+1
    
    ## trd
    d = fitdat[get(theY_ind_lag1) == prevS,]
    ytab = table(d$outcome)
    print(ytab)
    
    ylevels = as.numeric(names(ytab)[ytab > 0])
    p = length(ylevels) # number of outcome categories
    
    yprop = ytab/sum(ytab)
    print(rbind(ytab,yprop))
    
    Yname = theY
    Xname = vector("list", p-1)
    for(i in 1:(p-1)){
      Xname[[i]] = RHS
    }
 
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
      ntL[[ind]] = rep(nt, p-1)
      
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
      Bmod2[[ind]] = model_bart(Yname = Yname, Xname = Xname, data = d, type = "multinomial",
                                base = baseyl[ind],
                                Prior = Prior_mult(p = p, ntree = ntL[[ind]]),
                                Mcmc = Mcmc_mult(p = p, w0 = w, sig0 = sig),
                                correction = FALSE,
                                fitMNP = nf, bylatent = blt)
    } #doit
  } # prevS
  
}

if(t==1){
  save(fml, Bmod2, file = paste0("Cat19-03-02_",t,".RData"))
} else {
  save(fml, baseyl, Bmod2, file = paste0("Cat19-03-02_",t,".RData"))
}

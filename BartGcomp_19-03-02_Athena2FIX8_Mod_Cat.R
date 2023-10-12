args <- commandArgs(trailingOnly=TRUE)
s = eval(parse(text=args[1]))

print(s)

t = 8
##############################################################
source("BartGcomp_19-03-02_Athena0_prepare.R")

load("Cat19-03-02_8.RData")

fitdat = trd[time==200*t,]
print(dim(fitdat))

ntL = vector("list", 2)
ntL[[1]] = rep(nt, 2)

# [St | V, A_{t-1}, Xt, Lt, S_{t-1}]

RHS = c(theV, theA_lag1, theX_cov)
RHS = fmlfilter(fitdat,RHS)

fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))

#for(prevS in c(0,1)){ #### S_{t-1} = 0 - disengaged, 1 - engaged
  prevS = 1 
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
  
  set.seed(s)
  mod =  tryCatch(mnp_chains(fml, d, basey, ndraw, nburn, nthin, nchain, doparallel = FALSE),
                  error = function(e) NULL )
  if(length(mod)==0){
    set.seed(s)
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
    
    set.seed(s)
    Bmod2[[ind]] = model_bart(Yname = Yname, Xname = Xname, data = d, type = "multinomial",
                              base = baseyl[ind],
                              Prior = Prior_mult(p = p, ntree = ntL[[ind]]),
                              Mcmc = Mcmc_mult(p = p, w0 = w, sig0 = sig),
                              correction = FALSE,
                              fitMNP = nf, bylatent = blt)
  } #doit
#} # prevS


res = matrix(NA, ncol=4, nrow=4)
# the first dim: y, mean post, 5% post, 95% post
# the third dim: eng, diseng, xfer, death
d = ted[time==200*t,]

y =  unlist(ted[time==200*t, theY, with = F])

PB = predict_bart(Bmod2[[1]], d)
PM = PB$samp_y

PB1 = predict_bart(Bmod2[[2]], d)
PM1 = PB1$samp_y

ind = which(d[,theY_ind_lag1,with = F]==1)
PM[ind,] = PM1[ind,]

res[2:4, ] = matrix(unlist(lapply(1:4, function(s) summ(PM,s,1))), ncol=4)

res[1,] = taby(y)

print(res)

save(fml, res, baseyl, ntL, Bmod2, file = paste0("CatFIX19-03-02_",t,"_seed_",s,".RData"))

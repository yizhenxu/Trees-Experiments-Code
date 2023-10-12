args <- commandArgs(trailingOnly=TRUE)
t = eval(parse(text=args[1]))

# t - time interval index

##############################################################
source("BartGcomp_19-03-10_Athena0_prepare.R")

td = ted[time==200*t,]

ty =  unlist(ted[time==200*t, theY, with = F])

res = matrix(NA, ncol=4, nrow=4)

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
  
  ## trd
  ntL = rep(nt, p-1)
  baseyl = basey
  seedvec = CatSeed
  set.seed(seedvec)
  Bmod2 = model_bart(Yname = Yname, Xname = Xname, data = fitdat, type = "multinomial",
                     base = baseyl,
                     Prior = Prior_mult(p = p, ntree = ntL ),
                     Mcmc = Mcmc_mult(p = p),
                     correction = FALSE,
                     fitMNP = nf, bylatent = blt)
  
  ## ted
  PB = predict_bart(Bmod2, td)
  PM = PB$samp_y
  res[1,] = taby(ty)
  res[2:4, ] = matrix(unlist(lapply(1:4, function(s) summ(PM,s,1))), ncol=4) 
  ##############################################################
} else {
  ##############################################################
  
  fitdat = trd[time==200*t,]
  print(dim(fitdat))
  
  Bmod2 = vector("list", 2)
  baseyl = rep(basey,2)
  ntL = vector("list", 2)
  
  # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]
  
  RHS = c(theV, theA_lag1, theX_cov)
  RHS = fmlfilter(fitdat,RHS)
  
  fml = paste0(theY," ~ ",paste(RHS,collapse = " + "))
  
  seedvec = rep(CatSeed, 2) # seeds for (prevS = 0, prevS = 1)
  
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
    
    ## trd
    ntL[[ind]] = rep(nt, p-1)
    
    baseyl[ind] = names(which.max(ytab))[1]
    
    set.seed(seedvec[ind])
    Bmod2[[ind]] = model_bart(Yname = Yname, Xname = Xname, data = d, type = "multinomial",
                              base = baseyl[ind],
                              Prior = Prior_mult(p = p, ntree = ntL[[ind]]),
                              Mcmc = Mcmc_mult(p = p),
                              correction = FALSE,
                              fitMNP = nf, bylatent = blt)
    
    ## ted
    if(prevS == 0){
      PB0 = predict_bart(Bmod2[[1]], td)
      PM0 = PB0$samp_y
    } else { #prevS == 1
      
      res[1,] = taby(ty)
      tind = which(td[,theY_ind_lag1,with = F]==1)
      #Condition 1: convergence / Tree Depths after burnin < 3.5
      TD = Bmod2[[2]]$Tree_Depth
      TD = TD[, (ncol(TD)-nd+1):ncol(TD)] #tree depths of after-burnin samples
      diffTD = max(TD) - min(TD)
      cat(paste0("\n Tree Depth Max - Min:",diffTD,"\n\n\n"))
      
      #Condition 2: observed prop of xfer and death within 95% posterior CI
      PB1 = predict_bart(Bmod2[[2]], td)
      PM1 = PB1$samp_y
      
      PM = PM0; PM[tind,] = PM1[tind,]
      
      # the first dim: y, mean post, 5% post, 95% post
      # the second dim: eng, diseng, xfer, death
      res[2:4, ] = matrix(unlist(lapply(1:4, function(s) summ(PM,s,1))), ncol=4) 
      a = res[1,]-res[3,]; b = res[4,]-res[1,]
      cat(paste0("\n Accuracy status ",all(c(a[3:4],b[3:4])>=0),"\n"))
    }## ted end
  } # prevS
  
}



print(seedvec)
print(res)

save(seedvec, fml, res, baseyl, ntL, Bmod2, file = paste0("Cat19-03-10_",t,".RData"))

### Functions

#source("Functions_new.R")

source("G:/My Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/Functions_new.R")

Prior_mult = function(p, ntree){
  return(list(nu = p-1+3, V = diag(p-1),
              ntrees = ntree,
              kfac = 2,
              pswap = 0.1, pbd = 0.5, pb = 0.25,
              alpha = 0.95, beta = 2.0,
              nc = 100, minobsnode = 10))
}

Mcmc_mult = function(p,w0 = NULL,sig0 = NULL){
  res = list(sigma0 = diag(p-1), burn = nb, ndraws = nd, nSigDr = 50)  
  if(!is.null(w0)){
    res = append(res, list("w0" = w0), length(res))
  } 
  if(!is.null(sig0)){
    res$sigma0 = sig0
  }
  return(res)
}

#get the outcome level of the largest frequency 
ymode = function(vec){
  names(which.max(table(vec))) 
}

ylev = function(vec){
  tab = table(vec); all_lev = names(tab)
  lev = all_lev[-which.max(tab)]
  return(lev)
}

TDPlot <- function(bfit, trdy){
  nd = bfit$burn + bfit$ndraws
  lev = bfit$releveled
  
  for(i in 1:pm1){
    a = bfit$Tree_Depth
    scatter.smooth(1:nd, a[i,], lpars = list(col = "red", lwd = 3, lty = 3),
                   xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
                   main = paste0("Latent for Y=",lev[i]) )
  }
  #par(mfrow=c(1,1))
}


accuracyfun = function(bmpfit,testd,trdy,tedy){
  # posterior prediction by the entire distribution
  # mode accuracy, mean accuracy
  res = rep(NA, 4); names(res) = c("Train Mean", "Train Mode", "Test Mean", "Test Mode")
  nd = bmpfit$ndraws
  # train prediction
  BTr = bmpfit$samp_y
  yTr = matrix(rep(trdy, nd),ncol = nd)
  PmodeTr = apply(BTr, 1, ymode)
  res[1] = mean(BTr == yTr)
  res[2] = mean(PmodeTr == trdy)
  # test prediction
  BTe = predict_bart(obj = bmpfit, newdata = testd)$samp_y
  yTe = matrix(rep(tedy, nd),ncol = nd)
  PmodeTe = apply(BTe, 1, ymode)
  res[3] = mean(BTe == yTe)
  res[4] = mean(PmodeTe == tedy)
  return(res)
}

### Generate summary for the estimated Sigma
SigSumm <- function(bfit){
  nd = bfit$ndraws
  pm1 = bfit$ndim
  sig = bfit$sigmasample
  ss = matrix(sig, nrow = pm1^2)
  sm = matrix(apply(ss,1,mean), nrow = pm1)
  
  tmp = apply(ss,1,function(x) sort(x)[c(ceiling(nd*0.05), floor(nd*0.95))])
  s05 = matrix(tmp[1,], nrow = pm1)
  s95 = matrix(tmp[2,], nrow = pm1)
  
  sm = round(sm,3);s05 = round(s05,3);s95 = round(s95,3)
  res = matrix(NA, pm1, pm1)
  for(i in 1:pm1){
    for(j in 1:pm1){
      res[i,j] = paste0(sm[i,j], " (",s05[i,j],", ",s95[i,j],")")
    }
  }
  return(res)
}



### FACTORIZE CATEGORICAL VARIABLES WITH <= C LEVELS
factorizedata = function(dd, fl){
  d = copy(dd)
  for(j in 1:length(fl)){
    vname = fl[j]
    tmp = unique(d[,get(vname)])
    if(length(tmp)==1){
      d[,(vname) := NULL]  
    } else if(NA %in% tmp){
      d[,(vname) := addNA(get(vname))]  
    } else {
      d[,(vname) := as.factor(get(vname))]
    }
    
    print(vname)
  }
  return(d)
}


### CREATE INDICATOR VARIABLES FOR ALL FACTORS
# discard indicators with <=c number of 1's
# returns the data
createind = function(d,c){
  library(nnet)
  newd = copy(d)
  cname = colnames(d)
  for(j in 1:ncol(d)){
  
    if(is.factor(d[,get(cname[j])])){
      print(cname[j])
      tmp = as.data.frame(class.ind(d[,get(cname[j])]))# class.ind is function in package nnet
      tmpname = colnames(tmp)
      del = as.numeric(which(apply(tmp,2,function(x) sum(x))<=c))
      
      if(length(del)==0) del = 1 #if not deleting, delete 1st col as ref
      
      tmpname = tmpname[-del];  dtmp = tmp[,-del]
      if(ncol(tmp)>0){
        cnametmp = gsub(" ", "", tmpname, fixed = TRUE)# remove space from string for var names
        newd[, (paste0(cname[j], "_",cnametmp)) := dtmp]  
      }
      newd[,(cname[j]) := NULL]
    }#if
    
  }#j
  return(newd)
}

createmis = function(d,nl){
  newd = copy(d)
  cname = colnames(d)
  for(j in 1:length(nl)){
    cname = nl[j]
    a = d[,get(cname)]
    
    if(sum(is.na(a))>0){
      
      RX = a
      R = 1*(is.na(a)) # indicator of missing
      RX[is.na(RX)] = 0 # replace missing by 0
      tmp = data.frame(RX, R)
      newd[, (paste0(cname, "_",c("RX","R"))) := tmp]
      newd[,(cname) := NULL]
    }#if missing exists
    
  }#for all covariates
  
  return(newd)
  
}
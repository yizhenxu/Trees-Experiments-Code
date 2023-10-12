### Functions

source("Functions_new.R")
#source("/Volumes/GoogleDrive/My Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/Functions_new.R")

#source("G:/My Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/Functions_new.R")
#source("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/Functions_new.R")

# p is the number of outcome levels, d = p-1 is the number of latent utilities
# nu essentially acts as the number of observations we had observed prior to collecting the data, or, alternatively,
# the number of observations on which our prior sums-of-squares matrix (sum ZZ^T) is based.

# nu_C should be positive, nu = d + nu_C -- should try nu_C = 1, 2, 3
# diagonal V and nu_C = 1 results in marginal uniform distribution for all correlations, 
# still a proper prior, BUT nu=d+1 does not give closed analytic form of prior mean.
# when nu_C > 1, the prior mean has closed analytic form, V/(nu_C-1)
Prior_mult = function(p, nu_C, Vmat, ntree){
  return(list(nu = p-1+nu_C, V = Vmat,
              ntrees = ntree,
              kfac = 2,
              pswap = 0.1, pbd = 0.5, pb = 0.25,
              alpha = 0.95, beta = 2.0,
              nc = 100, minobsnode = 10))
}

Mcmc_mult = function(p,w0 = NULL,sig0 = NULL,nb, nd){
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

scplot = function(x,y, sp = 0.7, xlab = NULL, ylab = NULL, main = NULL){
  plot(x, y, type = "l", xlab = xlab, ylab = ylab, main = main)
  lines(loess.smooth(x, y, span = sp), col = "red", lwd = 3, lty = 3)
}

TDPlot <- function(bfit, j, all=T){ # j is the reference level
  lev = bfit$releveled
  a = bfit$Tree_Depth
  nd = bfit$ndraws; nb = bfit$burn
  tnd = nb + nd
  
  if(all){
    
    for(i in 1:pm1){
      scplot(1:tnd, a[i,], 
                     xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
                     main = paste0("L",lev[i]," - L",j) )
    }
    
  } else {
    
    for(i in 1:pm1){
      scplot(1:nd, a[i,(nb+1):tnd], 
                     xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
                     main = paste0("L",lev[i]," - L",j) )
    }
    
  }
  
}


ALLPlot <- function(bfit, j, all=T){ # j is the reference level
  lev = bfit$releveled
  nd = bfit$ndraws; nb = bfit$burn
  tnd = nb + nd
  
  if(all){
    
    for(i in 1:pm1){
      
      scplot(1:tnd, bfit$Percent_Acceptance[i,], 
             xlab = "MCMC Iteration", ylab = "Tree Acceptance Rate", 
             main = paste0("L",lev[i]," - L",j) )      
      scplot(1:tnd, bfit$Tree_Num_Nodes[i,], 
             xlab = "MCMC Iteration", ylab = "Average Num of Tree Nodes", 
             main = paste0("L",lev[i]," - L",j) )    
      scplot(1:tnd, bfit$Tree_Num_Leaves[i,], 
             xlab = "MCMC Iteration", ylab = "Average Num of Tree Leaves", 
             main = paste0("L",lev[i]," - L",j) )  
      scplot(1:tnd, bfit$Tree_Depth[i,], 
             xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
             main = paste0("L",lev[i]," - L",j) )
    }
    
  } else {
    
    for(i in 1:pm1){
      scplot(1:nd, bfit$Percent_Acceptance[i,(nb+1):tnd], 
             xlab = "MCMC Iteration", ylab = "Tree Acceptance Rate", 
             main = paste0("L",lev[i]," - L",j) )      
      scplot(1:nd, bfit$Tree_Num_Nodes[i,(nb+1):tnd], 
             xlab = "MCMC Iteration", ylab = "Average Num of Tree Nodes", 
             main = paste0("L",lev[i]," - L",j) )    
      scplot(1:nd, bfit$Tree_Num_Leaves[i,(nb+1):tnd], 
             xlab = "MCMC Iteration", ylab = "Average Num of Tree Leaves", 
             main = paste0("L",lev[i]," - L",j) )  
      scplot(1:nd, bfit$Tree_Depth[i,(nb+1):tnd], 
             xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
             main = paste0("L",lev[i]," - L",j) )
      
    }
    
  }
  
}

scplot2 = function(x,y1,y2, sp = 0.7, xlab = NULL, ylab = NULL, main = NULL){
  rg = range(c(y1,y2))
  plot(x, y1, type = "l",col='#00000088', ylim = rg, xlab = xlab, ylab = ylab, main = main)# black (88 for opacity, '#FF000088' red)
  lines(loess.smooth(x, y1, span = sp), col="green", lwd = 3, lty = 3)
  lines(x, y2, type="l", col='#FF000088')# red (blue '#0080FF88')
  lines(loess.smooth(x, y2, span = sp), col="blue", lwd = 3, lty = 3)
}

# KD=T, KD=F in one plot, same scale
TDPlot2 <- function(bfit0, bfit1, j, all=T){ # j is the reference level
  lev = bfit0$releveled
  a0 = bfit0$Tree_Depth; a1 = bfit1$Tree_Depth
  nd = bfit0$ndraws; nb = bfit0$burn
  tnd = nb + nd
  
  if(all){
    
    for(i in 1:pm1){
      scplot2(1:tnd, a0[i,], a1[i,], 
             xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
             main = paste0("L",lev[i]," - L",j) )
    }
    
  } else {
    
    for(i in 1:pm1){
      scplot2(1:nd, a0[i,(nb+1):tnd], a1[i,(nb+1):tnd], 
              xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
              main = paste0("L",lev[i]," - L",j) )
    }
    
  }
  
}



ALLPlot2 <- function(bfit0, bfit1, j, all=T){ # j is the reference level
  lev = bfit0$releveled
  nd = bfit0$ndraws; nb = bfit0$burn
  tnd = nb + nd
  
  if(all){
    
    for(i in 1:pm1){
      
      scplot2(1:tnd, bfit0$Percent_Acceptance[i,], bfit1$Percent_Acceptance[i,], 
             xlab = "MCMC Iteration", ylab = "Tree Acceptance Rate", 
             main = paste0("L",lev[i]," - L",j) )      
      scplot2(1:tnd, bfit0$Tree_Num_Nodes[i,], bfit1$Tree_Num_Nodes[i,],
             xlab = "MCMC Iteration", ylab = "Average Num of Tree Nodes", 
             main = paste0("L",lev[i]," - L",j) )    
      scplot2(1:tnd, bfit0$Tree_Num_Leaves[i,], bfit1$Tree_Num_Leaves[i,],
             xlab = "MCMC Iteration", ylab = "Average Num of Tree Leaves", 
             main = paste0("L",lev[i]," - L",j) )  
      scplot2(1:tnd, bfit0$Tree_Depth[i,], bfit1$Tree_Depth[i,], 
             xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
             main = paste0("L",lev[i]," - L",j) )
    }
    
  } else {
    
    for(i in 1:pm1){
      scplot2(1:nd, bfit0$Percent_Acceptance[i,(nb+1):tnd], bfit1$Percent_Acceptance[i,(nb+1):tnd], 
             xlab = "MCMC Iteration", ylab = "Tree Acceptance Rate", 
             main = paste0("L",lev[i]," - L",j) )      
      scplot2(1:nd, bfit0$Tree_Num_Nodes[i,(nb+1):tnd], bfit1$Tree_Num_Nodes[i,(nb+1):tnd], 
             xlab = "MCMC Iteration", ylab = "Average Num of Tree Nodes", 
             main = paste0("L",lev[i]," - L",j) )    
      scplot2(1:nd, bfit0$Tree_Num_Leaves[i,(nb+1):tnd], bfit1$Tree_Num_Leaves[i,(nb+1):tnd], 
             xlab = "MCMC Iteration", ylab = "Average Num of Tree Leaves", 
             main = paste0("L",lev[i]," - L",j) )  
      scplot2(1:nd, bfit0$Tree_Depth[i,(nb+1):tnd], bfit1$Tree_Depth[i,(nb+1):tnd], 
             xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
             main = paste0("L",lev[i]," - L",j) )
      
    }
    
  }
  
}

##############################
scplot3 = function(x,y1,y2,y3, sp = 0.65, xlab = NULL, ylab = NULL, main = NULL){
  rg = range(c(y1,y2))
  plot(x, y1, type = "l",col='#00000088', ylim = rg, xlab = xlab, ylab = ylab, main = main)# black (88 for opacity, '#FF000088' red)
  lines(x, y2, type="l", col='#FF000088')# red 
  lines(x, y3, type="l", col='#0080FF88')# blue 
  
  lines(loess.smooth(x, y1, span = sp), col="yellow", lwd = 3, lty = 3)
  lines(loess.smooth(x, y2, span = sp), col="blue", lwd = 3, lty = 3)
  lines(loess.smooth(x, y3, span = sp), col="red", lwd = 3, lty = 3)
}

# KD=T, KD=F mg=T, KD=F mg=F in one plot, same scale
TDPlot3 <- function(bfit0, bfit1, bfit2, j, all=T){ # j is the reference level
  lev = bfit0$releveled
  a0 = bfit0$Tree_Depth; a1 = bfit1$Tree_Depth; a2 = bfit2$Tree_Depth
  nd = bfit0$ndraws; nb = bfit0$burn
  tnd = nb + nd
  
  if(all){
    
    for(i in 1:pm1){
      scplot3(1:tnd, a0[i,], a1[i,], a2[i,], 
              xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
              main = paste0("L",lev[i]," - L",j) )
    }
    
  } else {
    
    for(i in 1:pm1){
      scplot3(1:nd, a0[i,(nb+1):tnd], a1[i,(nb+1):tnd], a2[i,(nb+1):tnd], 
              xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
              main = paste0("L",lev[i]," - L",j) )
    }
    
  }
  
}



ALLPlot3 <- function(bfit0, bfit1, bfit2, j, all=T){ # j is the reference level
  lev = bfit0$releveled
  nd = bfit0$ndraws; nb = bfit0$burn
  tnd = nb + nd
  
  if(all){
    
    for(i in 1:pm1){
      
      scplot3(1:tnd, bfit0$Percent_Acceptance[i,], bfit1$Percent_Acceptance[i,], bfit2$Percent_Acceptance[i,], 
              xlab = "MCMC Iteration", ylab = "Tree Acceptance Rate", 
              main = paste0("L",lev[i]," - L",j) )      
      scplot3(1:tnd, bfit0$Tree_Num_Nodes[i,], bfit1$Tree_Num_Nodes[i,], bfit2$Tree_Num_Nodes[i,],
              xlab = "MCMC Iteration", ylab = "Average Num of Tree Nodes", 
              main = paste0("L",lev[i]," - L",j) )    
      scplot3(1:tnd, bfit0$Tree_Num_Leaves[i,], bfit1$Tree_Num_Leaves[i,], bfit2$Tree_Num_Leaves[i,],
              xlab = "MCMC Iteration", ylab = "Average Num of Tree Leaves", 
              main = paste0("L",lev[i]," - L",j) )  
      scplot3(1:tnd, bfit0$Tree_Depth[i,], bfit1$Tree_Depth[i,], bfit2$Tree_Depth[i,], 
              xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
              main = paste0("L",lev[i]," - L",j) )
    }
    
  } else {
    
    for(i in 1:pm1){
      scplot3(1:nd, bfit0$Percent_Acceptance[i,(nb+1):tnd], bfit1$Percent_Acceptance[i,(nb+1):tnd], bfit2$Percent_Acceptance[i,(nb+1):tnd], 
              xlab = "MCMC Iteration", ylab = "Tree Acceptance Rate", 
              main = paste0("L",lev[i]," - L",j) )      
      scplot3(1:nd, bfit0$Tree_Num_Nodes[i,(nb+1):tnd], bfit1$Tree_Num_Nodes[i,(nb+1):tnd], bfit2$Tree_Num_Nodes[i,(nb+1):tnd], 
              xlab = "MCMC Iteration", ylab = "Average Num of Tree Nodes", 
              main = paste0("L",lev[i]," - L",j) )    
      scplot3(1:nd, bfit0$Tree_Num_Leaves[i,(nb+1):tnd], bfit1$Tree_Num_Leaves[i,(nb+1):tnd], bfit2$Tree_Num_Leaves[i,(nb+1):tnd], 
              xlab = "MCMC Iteration", ylab = "Average Num of Tree Leaves", 
              main = paste0("L",lev[i]," - L",j) )  
      scplot3(1:nd, bfit0$Tree_Depth[i,(nb+1):tnd], bfit1$Tree_Depth[i,(nb+1):tnd], bfit2$Tree_Depth[i,(nb+1):tnd], 
              xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
              main = paste0("L",lev[i]," - L",j) )
      
    }
    
  }
  
}
##############################

if(0){
  TDPlot <- function(bfit, j, all=T){ # j is the reference level
    lev = bfit$releveled
    a = bfit$Tree_Depth
    nd = bfit$ndraws; nb = bfit$burn
    tnd = nb + nd
    
    if(all){
      
      for(i in 1:pm1){
        scatter.smooth(1:tnd, a[i,], lpars = list(col = "red", lwd = 3, lty = 3),
                       xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
                       main = paste0("L",lev[i]," - L",j) )
      }
      
    } else {
      
      for(i in 1:pm1){
        scatter.smooth(1:nd, a[i,(nb+1):tnd], lpars = list(col = "red", lwd = 3, lty = 3),
                       xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
                       main = paste0("L",lev[i]," - L",j) )
      }
      
    }
    
  }
  
  
  ALLPlot <- function(bfit, j, all=T){ # j is the reference level
    lev = bfit$releveled
    nd = bfit$ndraws; nb = bfit$burn
    tnd = nb + nd
    
    if(all){
      
      for(i in 1:pm1){
        
        scatter.smooth(1:tnd, bfit$Percent_Acceptance[i,], lpars = list(col = "red", lwd = 3, lty = 3),
                       xlab = "MCMC Iteration", ylab = "Tree Acceptance Rate", 
                       main = paste0("L",lev[i]," - L",j) )      
        scatter.smooth(1:tnd, bfit$Tree_Num_Nodes[i,], lpars = list(col = "red", lwd = 3, lty = 3),
                       xlab = "MCMC Iteration", ylab = "Average Num of Tree Nodes", 
                       main = paste0("L",lev[i]," - L",j) )    
        scatter.smooth(1:tnd, bfit$Tree_Num_Leaves[i,], lpars = list(col = "red", lwd = 3, lty = 3),
                       xlab = "MCMC Iteration", ylab = "Average Num of Tree Leaves", 
                       main = paste0("L",lev[i]," - L",j) )  
        scatter.smooth(1:tnd, bfit$Tree_Depth[i,], lpars = list(col = "red", lwd = 3, lty = 3),
                       xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
                       main = paste0("L",lev[i]," - L",j) )
      }
      
    } else {
      
      for(i in 1:pm1){
        scatter.smooth(1:nd, bfit$Percent_Acceptance[i,(nb+1):tnd], lpars = list(col = "red", lwd = 3, lty = 3),
                       xlab = "MCMC Iteration", ylab = "Tree Acceptance Rate", 
                       main = paste0("L",lev[i]," - L",j) )      
        scatter.smooth(1:nd, bfit$Tree_Num_Nodes[i,(nb+1):tnd], lpars = list(col = "red", lwd = 3, lty = 3),
                       xlab = "MCMC Iteration", ylab = "Average Num of Tree Nodes", 
                       main = paste0("L",lev[i]," - L",j) )    
        scatter.smooth(1:nd, bfit$Tree_Num_Leaves[i,(nb+1):tnd], lpars = list(col = "red", lwd = 3, lty = 3),
                       xlab = "MCMC Iteration", ylab = "Average Num of Tree Leaves", 
                       main = paste0("L",lev[i]," - L",j) )  
        scatter.smooth(1:nd, bfit$Tree_Depth[i,(nb+1):tnd], lpars = list(col = "red", lwd = 3, lty = 3),
                       xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
                       main = paste0("L",lev[i]," - L",j) )
        
      }
      
    }
    
  }
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
  
  sm = round(sm,2);s05 = round(s05,2);s95 = round(s95,2)
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

drawprocess = function(x){#x = post$prob.test[1,]
  tmp = matrix(x, nrow = 3 )
  r = apply(tmp, 2, which.max)
  return(r)
}

drawagreement = function(pmat, y){#pmat = post$prob.test; y = newdat$y
  r = pmat[,(0:(length(y)-1))*3+y]
  tmp = apply(r,2,mean)
  return(tmp)
}



accuracyfun_mbart = function(post, y, ntr, n){
  
  # posterior prediction by the entire distribution
  # mode accuracy, mean accuracy
  res = rep(NA, 4); names(res) = c("Train Mean", "Test Mean", "Train Mode", "Test Mode")
  
  # agreement
  ag = drawagreement(post$prob.test, y)
  res[1] = mean(ag[1:ntr]) # train prediction
  res[2] = mean(ag[(ntr+1):n]) # test prediction
  
  # mode
  modep = drawprocess(post$prob.test.mean)
  macc = 1*(modep == y)
  res[3] = mean(macc[1:ntr]) # train prediction
  res[4] = mean(macc[(ntr+1):n]) # test prediction
  
  return(res)
}
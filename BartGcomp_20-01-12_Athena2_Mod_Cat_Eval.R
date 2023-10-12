#!/bin/bash
#SBATCH -J CAT
#SBATCH --time 72:00:00
#SBATCH --array=1-10
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.err
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

#Rscript BartGcomp_20-01-12_Athena2_Mod_Cat_Eval.R

rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer

t = as.numeric(slurm_arrayid)#1:10
setting = 4 

#sim = as.numeric(slurm_arrayid)#1:40
#setting = ceiling(sim/10)
#t = sim - 10*(setting-1)

# t - time interval index

##############################################################
source("BartGcomp_20-01-12_Athena0_prepare.R")

td = ted[time==200*t,]

ty =  unlist(ted[time==200*t, theY, with = F])

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
  
  ## trd

  baseyl = basey
  set.seed(CatSeed)
  Bmod2 = model_bart(as.formula(fml), data = fitdat, type = "multinomial",
                     base = baseyl,
                     Prior = Prior_mult(p = p, ntree = nt),
                     Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                     correction = FALSE,
                     Kindo = FALSE, do_alpha2_prior = FALSE)

  ##############################################################
} else {
  ##############################################################
  
  fitdat = trd[time==200*t,]
  print(dim(fitdat))
  
  Bmod2 = vector("list", 2)
  baseyl = rep(basey,2)
  
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
    
    ## trd
    
    baseyl[ind] = names(which.max(ytab))[1]
    
    set.seed(CatSeed)
    if(p==2){
      Bmod2[[ind]] = model_bart(as.formula(fml), data = d, type = "binary",
                 Prior = Prior_binary,
                 Mcmc = Mcmc_binary)
    } else {
      Bmod2[[ind]] = model_bart(as.formula(fml), data = d, type = "multinomial",
                                base = baseyl[ind],
                                Prior = Prior_mult(p = p, ntree = nt),
                                Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                                correction = FALSE,
                                Kindo = FALSE, do_alpha2_prior = FALSE)
    }
    
    
    #Condition 1: convergence / Tree Depths after burnin < 3.5
    #TD = Bmod2[[ind]]$Tree_Depth
    #TD = TD[, (ncol(TD)-nd+1):ncol(TD)] #tree depths of after-burnin samples
    #diffTD = max(TD) - min(TD)
    #cat(paste0("\n Tree Depth Max - Min:",diffTD,"\n\n\n"))    
      
  } # prevS
  
}


save( fml, baseyl,  Bmod2, file = paste0("Cat20-01-12_t",t,"_setting",setting,".RData"))
#save( fml, res, baseyl,  Bmod2, file = paste0("Cat20-01-12_t",t,"_setting",setting,".RData"))

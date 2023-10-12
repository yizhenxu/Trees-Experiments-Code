# run BartGcomp_19-11-26_mbart50_Athena.R by mbash2.sh for 50 replicates
# run BartGcomp_19-11-26_mbart_Athena.R by mbash1.sh for 1 run
#setwd("./MPBART")

nd = 30000 # posterior draws
nb = 50000 # burn in
nt = 100 # number of trees

sim = 1

book = matrix(NA, ncol = 2, nrow = 100)
colnames(book) = c("sim", "datseed")
book[,1] = rep(c(1,3) ,each = 50); book[,2] = rep(0:49, 2);

resl = vector("list", 2)

for(simj in 1:2){
  sim = c(1,3)[simj]
  datseed = 1; source("BartGcomp_19-11-26_Setting.R") #uniform prior for correlation
  for(datseed in 0:49){
    sfn = paste0("./res/CK11-26_b",nb,"d",nd,"s",sim,"_ct",datseed,"_mbart.RData")
    if(file.exists(sfn)){
      load(sfn)
      resl[[simj]] = rbind(resl[[simj]], tmp)
      rm(tmp)
    }
  }
}


summ = function(x){
  r = paste0(round(mean(x),3),' (',round(sd(x),3),')')
  return(r)
}

noquote(rbind(apply(resl[[1]], 2, summ),
        apply(resl[[2]], 2, summ)) )

Train Mean    Test Mean     Train Mode    Test Mode    
0.869 (0.003) 0.855 (0.003) 0.937 (0.004) 0.916 (0.004)
0.882 (0.003) 0.87 (0.003)  0.939 (0.003) 0.92 (0.004) 

load("./res/CK11-26_b50000d30000s1_mbart.RData");tmp
Train Mean  Test Mean Train Mode  Test Mode 
0.8683298  0.8520842  0.9334000  0.9142000


load("./res/CK11-26_b50000d30000s3_mbart.RData");tmp
Train Mean  Test Mean Train Mode  Test Mode 
0.8815273  0.8661162  0.9368000  0.9130000 


load("./res/CK_AMPATHacc_b10000d10000_AMPATH.RData");tmp
Train Mean  Test Mean Train Mode  Test Mode 
0.8143889  0.8131551  0.8885307  0.8883731 
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################




#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

if(0){
  
  accuracyfun_mbart1 = function(BTr,BTe,trdy,tedy){
    # posterior prediction by the entire distribution
    # mode accuracy, mean accuracy
    res = rep(NA, 4); names(res) = c("Train Mean", "Train Mode", "Test Mean", "Test Mode")
    
    nd = ncol(BTr)
    # train prediction
    yTr = matrix(rep(trdy, nd),ncol = nd)
    PmodeTr = apply(BTr, 1, ymode)
    res[1] = mean(BTr == yTr)
    res[2] = mean(PmodeTr == trdy)
    # test prediction
    yTe = matrix(rep(tedy, nd),ncol = nd)
    PmodeTe = apply(BTe, 1, ymode)
    res[3] = mean(BTe == yTe)
    res[4] = mean(PmodeTe == tedy)
    return(res)
  }
  
  py = apply(post$prob.test, 1, drawprocess) # Ntest x npost
  pytrn = py[1:ntr, (nb+1):(nb+nd)] # train prediction
  pyted = py[(ntr+1):n, (nb+1):(nb+nd)] # test prediction
  
  accuracyfun_mbart1(pytrn, pyted, trd$y, ted$y)
}

if(0){
  library(BART)
  
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
  
  
  
  accuracyfun_mbart = function(post, y, nd, nb, ntr, n){
    
    # posterior prediction by the entire distribution
    # mode accuracy, mean accuracy
    res = rep(NA, 4); names(res) = c("Train Mean", "Test Mean", "Train Mode", "Test Mode")
    
    # agreement
    ag = drawagreement(post$prob.test[(nb+1):(nb+nd),], y)
    res[1] = mean(ag[1:ntr]) # train prediction
    res[2] = mean(ag[(ntr+1):n]) # test prediction
    
    # mode
    modep = drawprocess(post$prob.test.mean)
    macc = 1*(modep == y)
    res[3] = mean(macc[1:ntr]) # train prediction
    res[4] = mean(macc[(ntr+1):n]) # test prediction
    
    return(res)
  }
  
  set.seed(99)
  post <- mbart(x.train = trd[,1:6], y.train = trd[,7], x.test = newdat[,1:6], ndpost = nd+nb, nskip = nb, keepevery = 1, ntree=nt)
  
  
  tmp = accuracyfun_mbart(post, newdat$y, nd, nb, ntr, n)
  
}
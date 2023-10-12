## run  BB.sh for parallelizing k posterior draws: BartGcomp_20-01-12_Athena4_BayesianBootstrap_batch.R 
# cd ./MPBART
# qsub BB.sh

# run joinBB.sh to join posterior draws: BartGcomp_20-01-12_Athena4_BayesianBootstrap_joinbatch.R 
# set round number rdnum to 1 with array 1:120 and then rdnum to 2 with array 1:3
# generate res_A.RData

rm(list=ls()) 


# coerce the value to an integer
taskid = as.numeric(Sys.getenv('SGE_TASK_ID'))

skip = 4; ndthin = 1000#floor(nd/(skip+1))

book = expand.grid(c(0,1,2), 1:ndthin)
# coerce the value to an integer
A = book[taskid, 1]; k = book[taskid, 2]
# theA_lag1 A = 0,1,2; A==2 means A = 1{CD4<=350}
print(paste0("A is ",A))

#setwd("./MPBART")
source("Functions_new.R")

setting = 4; maxt = 10; 
source("BartGcomp_20-01-12_Athena0_prepare.R")

##############################################################
### Outcome simulation 

#nlatent = length(unique(gcd$outcome)) - 1


Smat = vector("list", maxt) 
Xmat = vector("list", maxt)

#Xmatk = array(NA, dim = c(maxt, ngc, length(theX)))
#Smatk = matrix(NA, nrow = maxt, ncol = ngc)

for(t in 1:maxt){
  Xmat[[t]] = array(NA, dim = c(ngc, ndthin, length(theX)))
  Smat[[t]] = matrix(NA, nrow = ngc, ncol = ndthin)
}

set.seed(1)
all = lapply(1:ndthin, function(i){
  gcid = sample(1:nsub, ngc, replace = F)
  gcid = sort(gcid)
  return(which(newdat$myID %in% gcid))
})


for(t in 1:maxt){
  
  if(t==1){
    
    #[A0 | X(-1)]
    if(!(A %in% c(0,1))){
      #theA_lag1 as a function of theX_lag1 
      I0 = gcd[time==0, log_cd4preBplus1_R] #missing indicator
      C0 = gcd[time==0, log_cd4preBplus1_RX]
      policy_lag = 1*(I0==0 & (C0<= log(351)))    #rep(A, ngc)
    } else {
      policy_lag = rep(A, ngc)
    }
    
    #[R1 | V, A0, X0, X(-1)]
    load( paste0("ContBin20-01-12_t",t,"_setting",setting,".RData") )  
    #[S1 | V, A0, R1, L1]
    load(paste0("Cat20-01-12_t",t,"_setting",setting,".RData")) 
    
    #for(k in 1:ndthin){
      gcd = newdat[all[[k]],]#copy(newdat[myID %in% gcid,])
      
      ytab = table(gcd$outcome)
      ylevels = as.numeric(names(ytab)[ytab > 0])
      
      X0 = gcd[time ==0 ,log_cd4plus1_RX]
      X0[gcd[time ==0,log_cd4plus1_R] == 1] = NA
      
      Vdist = cbind(1,gcd[time==0,theV,with=F])
      
      #[R1 | V, A0, X0, X(-1)]
      theV0 = c(theV,"log_cd4preBplus1_RX","log_cd4preBplus1_R")
      dat.pred = makeX_bart_boot(theX, theX_cov, theV0, Vdist, ylevels, t,fml1, Xmat, Smat, policy_lag, gcd,k)
      
      tmp = mod1; tmp$ndraws = k+skip*(k-1)
      Xmat[[t]][,k,1] = predict_bart(tmp,skip, dat.pred)$samp_y[,k]
      
      # [L1 | V, A0, X0, X(-1), R1 = 1] 
      # same RHS as in modeling L1
      tmp = mod2; tmp$ndraws = k+skip*(k-1)
      postY = predict_bart(tmp,skip, dat.pred)$samp_y[,k]
      
      R = Xmat[[t]][,k,1]
      # X(1) = X(0) for people who have no CD4 update at time==0 (days (0,12])
      mis_ind = which(R == 0) 
      postY[mis_ind] = X0[mis_ind] 
      Xmat[[t]][,k,2] = postY
      
      remove(postY)
      
      #[S1 | V, A0, R1, L1]
      dat.pred = makeX_bart_boot(theX, theX_cov, theV, Vdist, ylevels, t, fml, Xmat, Smat, policy_lag, gcd,k)
      if(k==1){
        tmp = Bmod2; tmp$ndraws = 2#k+skip*(k-1)
        tmp$sigmasample = tmp$sigmasample[tmp$ndraw * (tmp$ndim)^2]
        Smat[[t]][,k] = predict_bart(tmp, 0,dat.pred, Gcomp = FALSE)$samp_y[,1]
      } else {
        tmp = Bmod2; tmp$ndraws = k+skip*(k-1)
        tmp$sigmasample = tmp$sigmasample[tmp$ndraw * (tmp$ndim)^2]
        Smat[[t]][,k] = predict_bart(tmp, skip,dat.pred, Gcomp = FALSE)$samp_y[,1]
      }
      
    #}
    
    print(paste0("t=",t))
    remove(fml1, fml2, mod1, mod2);
    remove(fml, baseyl,  Bmod2); gc()
    
  } else {# t>1
    load( paste0("ContBin20-01-12_t",t,"_setting",setting,".RData") ) 
    load(paste0("Cat20-01-12_t",t,"_setting",setting,".RData")) 
    #for(k in 1:ndthin){
      gcd = newdat[all[[k]],]
      Vdist = cbind(1,gcd[time==0,theV,with=F])
      
      ytab = table(gcd$outcome)
      ylevels = as.numeric(names(ytab)[ytab > 0])
      
      #[A(t-1) | A(t-2), L(t-2), S(t-1)]
      if(!(A %in% c(0,1))){
        A_prev = policy_lag
        I0 = Smat[[t-1]][,k]
        if(t==2){
          X0 = gcd[time ==0 ,log_cd4plus1_RX]
          X0[gcd[time ==0,log_cd4plus1_R] == 1] = NA
          C0 = X0
        } else {
          C0 = Xmat[[t-2]][,k,2]
        }
        dim(I0) = NULL
        dim(C0) = NULL
        policy_lag = 1*(A_prev== 1) + 1*(A_prev==0 & I0==2 & !is.na(C0) & (C0<= log(351)))    #rep(A, nrow(Vdist))
      }
      
      #[Rt | V, A_{t-1}, R_{t-1}, L_{t-1}, S_{t-1}]
      dat.pred = makeX_bart_boot(theX, theX_cov, theV, Vdist, ylevels, t, fml1, Xmat, Smat, policy_lag, gcd,k)
      tmp = mod1; tmp$ndraws = k+skip*(k-1)
      Xmat[[t]][,k,1] = predict_bart(tmp,skip, dat.pred)$samp_y[,k]
      
      #[Lt | V, A_{t-1}, R_{t-1}, L_{t-1}, S_{t-1}]
      tmp = mod2; tmp$ndraws = k+skip*(k-1)
      postY = predict_bart(tmp, skip, dat.pred)$samp_y[,k]
      
      R = Xmat[[t]][,k,1]
      X_prev = Xmat[[t-1]][,k,2]
      mis_ind = which(R == 0) 
      postY[mis_ind] = X_prev[mis_ind]
      Xmat[[t]][,k,2] = postY
      
      
      #[St | V, A_{t-1}, S_{t-1}, Rt, Lt]
      dat.pred = makeX_bart_boot(theX, theX_cov, theV, Vdist, ylevels, t, fml, Xmat, Smat, policy_lag, gcd,k)
      if(k==1){
        tmp = Bmod2[[1]]; tmp$ndraws = 2
        tmp$sigmasample = tmp$sigmasample[tmp$ndraw * (tmp$ndim)^2]
        Smat[[t]][,k] = predict_bart(tmp, 0, dat.pred)$samp_y[,k]
        
        remove(tmp); tmp = Bmod2[[2]]; tmp$ndraws = 2
        tmp$sigmasample = tmp$sigmasample[tmp$ndraw * (tmp$ndim)^2]
        tmp1 = predict_bart(tmp, 0, dat.pred)$samp_y[,k]
      } else {
        tmp = Bmod2[[1]]; tmp$ndraws = k+skip*(k-1)
        tmp$sigmasample = tmp$sigmasample[tmp$ndraw * (tmp$ndim)^2]
        Smat[[t]][,k] = predict_bart(tmp, skip, dat.pred)$samp_y[,k]
        
        remove(tmp); tmp = Bmod2[[2]]; tmp$ndraws = k+skip*(k-1)
        tmp$sigmasample = tmp$sigmasample[tmp$ndraw * (tmp$ndim)^2]
        tmp1 = predict_bart(tmp, skip, dat.pred)$samp_y[,k]
      }
      
      eng_ind = which(Smat[[t-1]][,k] == 2)
      Smat[[t]][eng_ind,k] = tmp1[eng_ind]
      
      done_ind = which(Smat[[t-1]][,k] %in% c(3,4)) # people who xferred or died from t-1
      Smat[[t]][done_ind,k] = Smat[[t-1]][done_ind,k]
      
    #}
    print(paste0("t=",t))
    remove(fml1, fml2, mod1, mod2);
    remove(fml, baseyl,  Bmod2); gc()
    gc()
  }#if t
  
}#for t

save(Xmat,Smat,file = paste0("Gcomp_A",A,"_k",k,"_22-05-04.RData"))

###############################################################################


if(0){
  
  for(A in 0:2){
    load(paste0("res_A",A,".RData"))
    
    ### TP plot
    tmp = unlist(Smat)
    
    m = length(Smat[[1]])
    tmp1 = c(rep(2,m),tmp[1:((maxt-1)*m)])
    
    longdat = data.table(time = rep(200*(1:maxt),each = length(Smat[[1]])), outcome_lag = tmp1, outcome = tmp)
    ldc = longdat[, .N, by = c("time", "outcome_lag", "outcome")]
    ldc_denom = longdat[, .N, by = c("time","outcome_lag")]
    
    ldc = ldc[-which(ldc$outcome_lag %in% c(3,4)),]
    ldc = ldc[order(time, outcome_lag, outcome)]
    ldc_denom = ldc_denom[-which(ldc_denom$outcome_lag %in% c(3,4)),]
    ldc_denom = ldc_denom[order(time, outcome_lag)]
    
    eldc = ldc[outcome_lag == 2,]
    eldc = rbind(eldc,list(time = 200, outcome_lag = 2, outcome = 3, N = 0 ))
    eldc = eldc[order(time, outcome, outcome_lag)]
    
    dldc = ldc[outcome_lag == 1,]
    
    eldc_denom = ldc_denom[outcome_lag == 2,]
    dldc_denom = ldc_denom[outcome_lag == 1,]
    
    engtp = matrix(eldc$N/rep(eldc_denom$N, each = 4), nrow = 4)
    rownames(engtp) = c("diengaged","engaged","transferred out","died")
    colnames(engtp) = unique(eldc$time)
    dengtp = matrix(dldc$N/rep(dldc_denom$N, each = 3), nrow = 3)
    rownames(dengtp) = c("diengaged","engaged","died")
    colnames(dengtp) = unique(dldc$time)
    
    ##
    #tmp = unlist(Smat)
    tmp1 = c(rep(2, m),tmp)
    
    longdat = data.table(time = rep(200*(0:maxt),each = m), outcome = tmp1)
    ldc = longdat[, .N, by = c("time","outcome")]
    
    tt = sum(ldc$time==200)
    if(tt == 3){
      ldc = rbind(ldc,list(time = c(0,0,0,200),outcome = c(1,3,4,3),N=c(0,0,0,0) ))
    } else {
      ldc = rbind(ldc,list(time = c(0,0,0),outcome = c(1,3,4),N=c(0,0,0) ))
    }
    
    ldc = ldc[order(time,outcome)]
    counts = matrix(ldc$N,nrow = 4)
    
    ctmp = counts[1,]
    counts[1,] = counts[2,]
    counts[2,] = ctmp
    
    #rownames(counts) = c("diengaged","engaged","transferred out","died")
    counts = simplify2array(lapply(1:ncol(counts), function(x) counts[,x]/sum(counts[,x])))
    colnames(counts) = 200*(0:maxt)
    save(counts,file = paste0("Counts_A",A,".RData"))
    
  }
  ### start here
  
  library(RColorBrewer)
  library(ggplot2)
  p = vector("list",3)
  Tt = list("Never Treat","Treat Immediately","Treat when CD4 Drops below 350 cells/mm3")
  for(A in 0:2){
    load( paste0("Counts_A",A,".RData"))
    
    rownames(counts) = c("Engaged","Disengaged","Transferred out","Died")
    
    tab = data.frame(Status = factor( rep(rownames(counts), each = ncol(counts)), levels = c("Died","Transferred out","Disengaged","Engaged") ),
                     Time = rep(as.numeric(colnames(counts)), nrow(counts)),
                     Proportion = c(t(counts)))
    
    colourCount = 4
    getPalette = colorRampPalette(brewer.pal(4, "Set2"))
    p[[A+1]] = ggplot(data = tab, aes(x=Time, y=Proportion, fill=Status)) +
      geom_bar(stat="identity")+
      scale_fill_manual(values = rev(getPalette(colourCount)))+
      theme(axis.title.y = element_text(size = 20, angle = 90),
            axis.title.x = element_text(size = 20),
            legend.text=element_text(size=20),
            legend.title=element_text(size=20),
            plot.title = element_text(size=20))+
      ggtitle(Tt[[A+1]])
    pdf(paste0("CondSEngA",A,".pdf"), width =9, height = 5)
    plot(p[[A+1]])
    dev.off()
  }
  
  mylegend<-g_legend(p[[1]])
  
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(reshape2)
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  mylegend<-g_legend(p[[A+1]])
  pdf("CondSEngA.pdf", width =9, height = 15)
  
  grid.arrange(arrangeGrob(p[[1]] + theme(legend.position="none"),
                           p[[2]] + theme(legend.position="none"),
                           p[[3]] + theme(legend.position="none"),
                           ncol=1),
               mylegend, nrow=1,widths=c(5.5, 2))
  dev.off()
  
  
}




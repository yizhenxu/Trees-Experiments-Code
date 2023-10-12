#!/bin/bash
#SBATCH -J gcomp
#SBATCH --time 72:00:00
#SBATCH --array=0-2
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.err
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

#Rscript BartGcomp_20-01-12_Athena4_Gcomp_thin1000.R

rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
A = as.numeric(slurm_arrayid)
# theA_lag1 A = 0,1,2; A==2 means A = 1{CD4<=350}
print(paste0("A is ",A))

# TP plot
# Gcomp barplot
# diffCI

setting = 4; maxt = 10; 
source("BartGcomp_20-01-12_Athena0_prepare.R")
skip = 4; ndthin = floor(nd/(skip+1))

X0 = gcd[time ==0 ,log_cd4plus1_RX]
X0[gcd[time ==0,log_cd4plus1_R] == 1] = NA

X0Mat = X0[rep(1:ngc,ndthin)]
##############################################################
### Outcome simulation 

nlatent = length(unique(gcd$outcome)) - 1
ytab = table(gcd$outcome)
ylevels = as.numeric(names(ytab)[ytab > 0])

Smat = vector("list", maxt) 
Xmat = vector("list", maxt)

for(i in 1:maxt){
  Xmat[[i]] = array(NA, dim = c(ngc, ndthin, length(theX)))
}

Vdist = cbind(1,gcd[time==0,theV,with=F])
Vdist = Vdist[rep(1:ngc,ndthin),]

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
    
    load( paste0("ContBin20-01-12_t",t,"_setting",setting,".RData") )  
    
    #[R1 | V, A0, X0, X(-1)]
    
    theV0 = c(theV,"log_cd4preBplus1_RX","log_cd4preBplus1_R")
    dat.pred = makeX_bart(theX, theX_cov, theV0, Vdist, c(1,2,4), t,fml1, Xmat, Smat, policy_lag, gcd)
    
    Xmat[[t]][,,1] = predict_bart(mod1,skip, dat.pred)$samp_y
    
    # [L1 | V, A0, X0, X(-1), R1 = 1] 
    #dat.pred = makeX(Xname = theX, Vname = theV, time = t, modhere = mod[[ct]], simX = Xmat, policy_lag = A, dathere = gcd)
    # same RHS as in modeling L1
    postY = predict_bart(mod2,skip, dat.pred)$samp_y
    
    remove(fml1, fml2, mod1, mod2); gc()
    
    R = Xmat[[t]][,,1]
    # X(1) = X(0) for people who have no CD4 update at time==0 (days (0,12])
    mis_ind = which(R == 0) 
    postY[mis_ind] = X0Mat[mis_ind] 
    Xmat[[t]][,,2] = postY
    
    remove(postY)
    
    #[S1 | V, A0, R1, L1]
    
    load(paste0("Cat20-01-12_t",t,"_setting",setting,".RData")) 
    
    if(!(A %in% c(0,1))){
      policy_lag = policy_lag[rep(1:ngc,ndthin)]
    } else {
      policy_lag = rep(A, nrow(Vdist))
    }
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, c(1,2,4), t, fml, Xmat, Smat, policy_lag, gcd)
    
    Smat[[t]] = predict_bart(Bmod2, skip,dat.pred, Gcomp = TRUE)$samp_y
    
    remove(fml, res, baseyl,  Bmod2); gc()
    
    print(paste0("t=",t))
    
  } else {# t>1
    
    #[A(t-1) | A(t-2), L(t-2), S(t-1)]
    if(!(A %in% c(0,1))){
      
      A_prev = policy_lag
      I0 = Smat[[t-1]]
      if(t==2){
        C0 = X0Mat
        remove(X0Mat); gc()
      } else {
        C0 = Xmat[[t-2]][,,2]
      }
      dim(I0) = NULL
      dim(C0) = NULL
      policy_lag = 1*(A_prev== 1) + 1*(A_prev==0 & I0==2 & !is.na(C0) & (C0<= log(351)))    #rep(A, nrow(Vdist))
      
    }
    
    load( paste0("ContBin20-01-12_t",t,"_setting",setting,".RData") ) 
    
    #[Rt | V, A_{t-1}, R_{t-1}, L_{t-1}, S_{t-1}]
    
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, fml1, Xmat, Smat, policy_lag, gcd)
    Xmat[[t]][,,1] = predict_bart(mod1, skip,dat.pred, Gcomp = TRUE)$samp_y
    
    #[Lt | V, A_{t-1}, R_{t-1}, L_{t-1}, S_{t-1}]
    
    #dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, fml2, Xmat, Smat, policy_lag, gcd)
    postY = predict_bart(mod2, skip,dat.pred, Gcomp = TRUE)$samp_y
    
    remove(fml1, fml2, mod1, mod2); gc()
    
    R = Xmat[[t]][,,1]
    X_prev = Xmat[[t-1]][,,2]
    mis_ind = which(R == 0) 
    postY[mis_ind] = X_prev[mis_ind]
    Xmat[[t]][,,2] = postY
    
    remove(postY)
    
    #[St | V, A_{t-1}, S_{t-1}, Rt, Lt]
    
    load(paste0("Cat20-01-12_t",t,"_setting",setting,".RData")) 
    
    dat.pred = makeX_bart(theX, theX_cov, theV, Vdist, ylevels, t, fml, Xmat, Smat, policy_lag, gcd)
    
    Smat[[t]] = predict_bart(Bmod2[[1]], skip, dat.pred, Gcomp = TRUE)$samp_y
    
    tmp1 = predict_bart(Bmod2[[2]], skip, dat.pred, Gcomp = TRUE)$samp_y
    
    remove(fml, res, baseyl,  Bmod2); gc()
    
    eng_ind = which(Smat[[t-1]] == 2)
    Smat[[t]][eng_ind] = tmp1[eng_ind]
    remove(tmp1)
    
    done_ind = which(Smat[[t-1]] %in% c(3,4)) # people who xferred or died from t-1
    Smat[[t]][done_ind] = Smat[[t-1]][done_ind]
    
    print(paste0("t=",t))
    
    
  }#if t
  
  
  save(Xmat,Smat,file = paste0("Gcomp_A",A,"_t",t,"_setting",setting,"_20-01-12.RData"))
  
}#for t

save(Xmat,Smat,file = paste0("Gcomp_A",A,"_setting",setting,"_20-01-12.RData"))

####################################################################

if(0){
  for(A in 0:2){
    load(paste0("Gcomp_A",A,"_setting",setting,"_20-01-12.RData"))
    
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
    
    save(engtp,dengtp,file = paste0("TP_A",A,".RData"))
    
    
   
    ####################################################################
    ###Gcomp bartplot
    
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
    ####################################################################
    ### Prepare for diffCI
    summ2 = function(vec){
      return(unlist(lapply(1:4, function(x) mean(vec==x))))
    }
    engci = function(i, Smat){
      res = apply(Smat[[i]],2, summ2)
      return(res)
    }
    
    mat = lapply(1:maxt, function(i) engci(i, Smat)) #in A = 1 and A = 1{CD4 < 350}  
    
    save(mat,file = paste0("CImat_A",A,".RData"))
  }
  
  ####################################################################
   ### TPplot
 
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(reshape2)
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  p1 = vector("list",3)
  p2 = vector("list",3)
  for(A in 0:2){
    load(paste0("TP_A",A,".RData"))
    dd = t(engtp)
    colnames(dd) = c("Disengaged","Engaged","Transferred Out","Died")
    plotdat = data.frame(Time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                         Status = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
    p1[[A+1]] = ggplot(plotdat, aes(Time, transition_prob))+ geom_line(aes(group = Status,col=Status),size=1.2)+
      labs(title = "Transition from Engagement")+ labs(y = "Transition Probabilities")+
      theme(axis.title.y = element_text(size = 20, angle = 90),
            axis.title.x = element_text(size = 20),
            legend.text=element_text(size=20),
            legend.title=element_text(size=20),
            plot.title = element_text(size=20))
    
    
    dd = t(dengtp)
    colnames(dd) = c("Disengaged","Engaged","Died")
    plotdat = data.frame(Time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                         Status = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
    p2[[A+1]] = ggplot(plotdat, aes(Time, transition_prob))+ geom_line(aes(group = Status,col=Status),size=1.2)+
      labs(title = "Transition from Disengagement")+ labs(y = "Transition Probabilities")+
      theme(axis.title.y = element_text(size = 20, angle = 90),
            axis.title.x = element_text(size = 20),
            legend.text=element_text(size=20),
            legend.title=element_text(size=20),
            plot.title = element_text(size=20))
  
    mylegend<-g_legend(p1[[A+1]])
    pdf(paste0("TPtime_A",A,".pdf"), width =15, height = 5)
    
    grid.arrange(arrangeGrob(p1[[A+1]] + theme(legend.position="none"),
                             p2[[A+1]] + theme(legend.position="none"),
                             nrow=1),
                 mylegend, nrow=1,widths=c(8, 2))
    dev.off()
    
  }
  
  mylegend<-g_legend(p1[[1]])
  pdf("TPtime.pdf", width =15, height = 15)
  
  grid.arrange(arrangeGrob(p1[[1]] + theme(legend.position="none"),
                           p2[[1]] + theme(legend.position="none"),
                           p1[[2]] + theme(legend.position="none"),
                           p2[[2]] + theme(legend.position="none"),
                           p1[[3]] + theme(legend.position="none"),
                           p2[[3]] + theme(legend.position="none"),
                           nrow=3),
               mylegend, nrow=1,widths=c(10, 2))
  dev.off()
  
  #####################################################################
  ### Barplot
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
  pdf("CondSEngA.pdf", width =9, height = 15)
  
  grid.arrange(arrangeGrob(p[[1]] + theme(legend.position="none"),
                           p[[2]] + theme(legend.position="none"),
                           p[[3]] + theme(legend.position="none"),
                           ncol=1),
               mylegend, nrow=1,widths=c(5.5, 2))
  dev.off()
  
  mylegend<-g_legend(p[[1]])
  pdf("CondSEngAThesis.pdf", width =18, height = 8)
  
  grid.arrange(arrangeGrob(p[[2]] + theme(legend.position="none"),
                           p[[3]] + theme(legend.position="none"),
                           ncol=2),
               mylegend, nrow=1,widths=c(5.5, 1))
  dev.off()
    #####################################################################
    ### Diff Confidence interval plots
    CImat = vector("list",3)
    for(A in c(0,1,2)){
      load(paste0("CImat_A",A,".RData"))
      CImat[[A+1]] = mat
    }
    
    
    cp = matrix(c(2,3,2,1,3,1),ncol=2,byrow = T)
    Type = list(expression(plain("Treat Immediately v.s. Treat when CD4<350 cells/mm")^3),
                "Treat Immediately vs No Treatment",
                expression(paste(plain("Treat when CD4<350 cells/mm")^3, plain(" v.s. No Treatment"))) )
    Title = list("A1-A2","A1-A0","A2-A0") 
    
    #cp = c(2,3) #Treat immediately vs Treat when CD4 >= 350
    #cp = c(2,1) #Treat immediately vs No treatment
    #cp = c(3,1) #Treat when CD4 >= 350 vs No treatment
    
    
    l = 1
    nd = ncol(CImat[[1]][[1]])
    tab = vector("list",4) #one for each outcome category
    
    for(k in 1:4){
      tab[[k]] = matrix(NA, 3, maxt)
      
      for(i in 1:maxt){
        dif = CImat[[ cp[l,1] ]][[i]][k,] - CImat[[ cp[l,2] ]][[i]][k,]
        tab[[k]][1,i] = mean(dif)
        a = sort(dif)
        tab[[k]][2,i] = a[floor(nd*0.025)]
        tab[[k]][3,i] = a[ceiling(nd*0.975)]
        
      }
    }
    
    pdf(paste0("DiffCI_",Title[[l]],".pdf"), width =9, height = 9)
    
    par(mfrow =c(2,2))
    for(k in 1:4){
      pt = tab[[k]]
      x = 200*(1:maxt)
      cname = c("Disengaged", "Engaged", "Transferred", "Died")
      if(k==3){#when k = 3, transffered, no first time point
        pt = pt[,-1]
        x = 200*(2:maxt)
      }
      plot(x, pt[1,],
           ylim=range(c(pt,0)),
           pch=19, xlab="Time", ylab=paste0("Differece in Percent ",cname[k]) ,cex.lab=1.4)
      # hack: we draw arrows but with very special "arrowheads"
      arrows(x, pt[2,], x, pt[3,], length=0.05, angle=90, code=3)
      abline(h = 0, lty = 2)
    }
    
   # mtext(Type[[l]], side = 3, line = -2.7, outer = TRUE,cex = 2)
    
    dev.off()
    
    #####################################################################
    ### Diff Confidence interval plots Separated
    
    
    l = 1
    nd = ncol(CImat[[1]][[1]])
    tab = vector("list",4) #one for each outcome category
    
    for(k in 1:4){
      tab[[k]] = matrix(NA, 3, maxt)
      
      for(i in 1:maxt){
        dif = CImat[[ cp[l,1] ]][[i]][k,] - CImat[[ cp[l,2] ]][[i]][k,]
        tab[[k]][1,i] = mean(dif)
        a = sort(dif)
        tab[[k]][2,i] = a[floor(nd*0.025)]
        tab[[k]][3,i] = a[ceiling(nd*0.975)]
        
      }
    }
    
    
    cname = c("Disengaged", "Engaged", "Transferred", "Died")
    
    for(k in 1:4){
      pdf(paste0("DiffCI_",Title[[l]],"_",cname[k],".pdf"), width =9, height = 7)
      
      pt = tab[[k]]
      x = 200*(1:maxt)
      
      if(k==3){#when k = 3, transffered, no first time point
        pt = pt[,-1]
        x = 200*(2:maxt)
      }
      plot(x, pt[1,],
           ylim=range(c(pt,0)),
           pch=19, xlab="Time", ylab=paste0("Differece in Percent ",cname[k]) ,cex.lab=1.5)
      # hack: we draw arrows but with very special "arrowheads"
      arrows(x, pt[2,], x, pt[3,], length=0.05, angle=90, code=3)
      abline(h = 0, lty = 2)
      dev.off()
      
    }
    
  #############################################################
  ### CHECK engagement at time=200
  a = newdat[time==200,]
  tab = matrix(0,2,3); rownames(tab) = c("A0=0","A0=1"); colnames(tab) = c("CD4<=350","CD4>350 or missing","all")
  for(A0 in 0:1){
    for(j in 1:3){
      if(j == 1) a1 = a[onarv_cf_1_lag1 == A0 & log_cd4plus1_R == 0 & log_cd4plus1_RX<=log(351),]
      if(j == 2) a1 = a[onarv_cf_1_lag1 == A0 & (log_cd4plus1_R == 1 | (log_cd4plus1_R == 0 & log_cd4plus1_RX>log(351))),]
      if(j == 3) a1 = a[onarv_cf_1_lag1 == A0, ]
      tab[A0+1, j] = paste0(round(mean(a1$outcome==2),3), " (",nrow(a1),")")
    }
  }
  print(noquote(tab))
  
    a1 = a[onarv_cf_1_lag1 == 0 & log_cd4plus1_R == 0 & log_cd4plus1_RX<=log(351),]
  mean(a1$outcome==2);dim(a1)
  
  
}
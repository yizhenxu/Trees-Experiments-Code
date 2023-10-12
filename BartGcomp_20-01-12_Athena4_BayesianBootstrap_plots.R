setwd("./Estimation RData")

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(data.table)
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
  
  #mylegend<-g_legend(p1[[A+1]])
  #pdf(paste0("TPtime_A",A,".pdf"), width =15, height = 5)
  
  #grid.arrange(arrangeGrob(p1[[A+1]] + theme(legend.position="none"),
  #                         p2[[A+1]] + theme(legend.position="none"),
  #                         nrow=1),
  #             mylegend, nrow=1,widths=c(8, 2))
  #dev.off()
  
}

mylegend<-g_legend(p1[[1]])
pdf("TPtime.pdf", width =15, height = 15)

grid.arrange(arrangeGrob(p1[[1]] + theme(legend.position="none"),
                         p2[[1]] + theme(legend.position="none"),
                         p1[[3]] + theme(legend.position="none"),
                         p2[[3]] + theme(legend.position="none"),
                         p1[[2]] + theme(legend.position="none"),
                         p2[[2]] + theme(legend.position="none"),
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
  #pdf(paste0("CondSEngA",A,".pdf"), width =9, height = 5)
  #plot(p[[A+1]])
  #dev.off()
}

mylegend<-g_legend(p[[1]])
pdf("CondSEngA.pdf", width =9, height = 15)

grid.arrange(arrangeGrob(p[[1]] + theme(legend.position="none"),
                         p[[3]] + theme(legend.position="none"),
                         p[[2]] + theme(legend.position="none"),
                         ncol=1),
             mylegend, nrow=1,widths=c(5.5, 2))
dev.off()

#mylegend<-g_legend(p[[1]])
#pdf("CondSEngAThesis.pdf", width =18, height = 8)

#grid.arrange(arrangeGrob(p[[2]] + theme(legend.position="none"),
#                         p[[3]] + theme(legend.position="none"),
#                         ncol=2),
#             mylegend, nrow=1,widths=c(5.5, 1))
#dev.off()
#####################################################################
### Diff Confidence interval plots
CImat = vector("list",3)
for(A in c(0,1,2)){
  load(paste0("CImat_A",A,".RData"))
  CImat[[A+1]] = mat
}


cp = matrix(c(2,3,2,1,3,1),ncol=2,byrow = T)
Type = list(expression(plain("Treat Immediately v.s. Treat when CD4<350 cells/mm")^3),
            expression("Treat Immediately vs No Treatment"),
            expression(paste(plain("Treat when CD4<350 cells/mm")^3, plain(" v.s. No Treatment"))) )
Title = list("A1-A2","A1-A0","A2-A0") 

#cp = c(2,3) #Treat immediately vs Treat when CD4 >= 350
#cp = c(2,1) #Treat immediately vs No treatment
#cp = c(3,1) #Treat when CD4 >= 350 vs No treatment


for(l in 1:3){
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
    if(l!=1) rg = range(c(pt,0))
    if(l==1) rg = c(-0.065, 0.065)
    plot(x, pt[1,],
         ylim = rg,
         pch=19, xlab="Time", ylab=paste0("Differece in P( ",cname[k]," | Data)") ,cex.lab=1.4)
    # hack: we draw arrows but with very special "arrowheads"
    arrows(x, pt[2,], x, pt[3,], length=0.05, angle=90, code=3)
    abline(h = 0, lty = 2)
  }
  
  mtext(Type[[l]], side = 3, line = -2.7, outer = TRUE,cex = 2)
  
  dev.off()
  
}

pdf("DiffCI.pdf")#, width =9, height = 25)
par(mfrow=c(3,1))
gap=25
for(l in 1:3){
  
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
  
  for(k in 1:4){
    pt = tab[[k]]
    x = 200*(1:maxt)
    cname = c("Disengaged", "Engaged", "Transferred", "Died")
    if(k==3){#when k = 3, transffered, no first time point
      pt = pt[,-1]
      x = 200*(2:maxt)
    }
    rg = c(-0.27, 0.27)
    if(k==1){
      plot(x, pt[1,],
           ylim = rg,cex=0.5,
           pch=19, xlab="Time", ylab=expression(Delta*Probability),
           main = Type[[l]], xlim = c(190, 2000+gap*4),cex.lab=1.4)
      # hack: we draw arrows but with very special "arrowheads"
      arrows(x, pt[2,], x, pt[3,], length=0.05, angle=90, code=3)
    } else {
      points(x+gap*k, pt[1,], pch=19, col=k,cex=0.5)
      # hack: we draw arrows but with very special "arrowheads"
      arrows(x+gap*k, pt[2,], x+gap*k, pt[3,], length=0.03, angle=90, code=3,col=k)
    }
    
  }
  abline(h = 0, lty = 2)
}
legend(1, 95, legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8,
       box.lty=0)
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


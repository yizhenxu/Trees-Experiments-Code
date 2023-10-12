# theA_lag1 A = 0,1,2; A==2 means A = 1{CD4<=350}

source("BartGcomp_18-10-27_Athena0_prepare.R")

#pathfrom = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"
pathto = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft3\\Accuracy Comparison - Diag Plots\\19-03-05\\Plots"
pathfrom = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"

#####################################################################
### Gcomp counts barplots
#0,1,2; 2 for 1{CD4<=350}

library(ggplot2)

#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cPalette <- c("#0072B2","#E69F00","#FF3333","#339966")

for(A in c(0,1,2)){
  load(paste0(pathfrom,"\\Counts_A",A,".RData"))
  rownames(counts) = c("Engaged","Disengaged","Transferred out","Died")
  
  tab = data.frame(Status = factor( rep(rownames(counts), each = ncol(counts)), levels = c("Died","Transferred out","Disengaged","Engaged") ),
                   Time = rep(as.numeric(colnames(counts)), nrow(counts)),
                   Proportion = c(t(counts)))
  
  #png(paste0(pathto,"\\CondSEngA",A,".png"),width = 600, height = 400)
  ggplot(data = tab, aes(x=Time, y=Proportion, fill=Status)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values = c("Died"=cPalette[1],"Transferred out"=cPalette[2],"Disengaged"=cPalette[3],"Engaged"=cPalette[4]))
  ggsave(paste0(pathto,"\\CondSEngA",A,".png"),width = 6, height = 3.5)
  #dev.off()  
}

#####################################################################
### Transition probability plots

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

for(A in c(0,1,2)){
  load(paste0(pathfrom,"\\TP_A",A,".RData"))
  
  rownames(engtp) = c("Disengaged","Engaged","Transferred out","Died")
  rownames(dengtp) = c("Disengaged","Engaged","Died")
  
  dd = t(engtp)
  plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                       group = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
  p1 = ggplot(plotdat, aes(time, transition_prob))+ geom_line(aes(group = group,col=group),size=1.2)+
    labs(title = "Transition from engagement")+ labs(y = "Transition Probabilities")+
    theme(plot.title = element_text(size=22),text = element_text(size=20),
          axis.text.x = element_text(size = 20))
  
  
  dd = t(dengtp)
  plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                       group = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
  p2 = ggplot(plotdat, aes(time, transition_prob))+ geom_line(aes(group = group,col=group),size=1.2)+
    labs(title = "Transition from disengagement")+ labs(y = "Transition Probabilities")+
    theme(plot.title = element_text(size=22),text = element_text(size=20),
          axis.text.x = element_text(size = 20))
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  mylegend<-g_legend(p1)
  
  png(filename = paste0(pathto,"\\TPtime_A",A,".png"), width =800, height = 300)
  
  grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                           p2 + theme(legend.position="none"),
                           nrow=1),
               mylegend, nrow=1,widths=c(10, 2))
  dev.off()
  
}

#####################################################################
### Diff Confidence interval plots
CImat = vector("list",3)
for(A in c(0,1,2)){
  load(paste0(pathfrom,"\\CImat_A",A,".RData"))
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

png(filename = paste0(pathto,"\\DiffCI_",Title[[l]],".png"), width =800, height = 650)

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
       pch=19, xlab="Time", ylab=paste0("Differece in Percent ",cname[k]) ,cex.lab=1.5)
  # hack: we draw arrows but with very special "arrowheads"
  arrows(x, pt[2,], x, pt[3,], length=0.05, angle=90, code=3)
  abline(h = 0, lty = 2)
}

mtext(Type[[l]], side = 3, line = -2.7, outer = TRUE,cex = 2)

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
  png(filename = paste0(pathto,"\\DiffCI_",Title[[l]],"_",cname[k],".png"), width =600, height = 400)
  
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

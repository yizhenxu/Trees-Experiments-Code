### upload models for categorical BART
t=10
sc = c(211025,
       210968,
       210996,
       210953,
       210921,
       211039,
       210939,
       211067,
       211053,
       211011)

screennum = sc[t]

setwd("./Topic2");source("BartGcomp_18-10-27_Athena2_Mod_Cat.R",echo=T)

### upload models for continuous and binary BART
screennum = 210982
setwd("./Topic2");source("BartGcomp_18-10-27_Athena1_Mod_ContBin.R",echo=T)

##################################################################################

### combine all models into "Mod18-10-27.RData"
# at screen 211011
#rm(list = ls())

load("ContBin.RData")

mod3_D =  vector("list", maxt-1) # categorical models for those previously disengaged
mod3_E =  vector("list", maxt) # categorical models for those previously engaged
baseylist = rep(NA, 2*(maxt-1)+1)

for(t in 1:9){
  load(paste0("Cat_",t,".RData"))
  
  modfml[[3*t]] = fml
  
  if(t==1){
    baseylist[1] = 2
    mod3_E[[1]] = Bmod2
  } else {
    baseylist[(2*(t-1)):(2*(t-1)+1)] = baseyl
    mod3_D[[t-1]] = Bmod2[[1]]
    mod3_E[[t]] = Bmod2[[2]]
  }
  remove(Bmod2)
  remove(fml)
  remove(baseyl)
}

save(mod1, mod2, mod3_D, mod3_E, baseylist, modfml, file = "Mod18-10-27.RData")
##################################################################################

### convergence plot

## Run BartGcomp_18-10-27_Athena6_ConvergencePlot.R in ZEUS

##################################################################################

### accuracy plot for S using second part of BartGcomp_18-10-27_Athena3_Eval and loadtest.sh
### plot from accuracy

path = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"
load(paste0(path,"\\res.RData"))

#maxt = 10
pd = as.data.frame(matrix(NA, nrow = 4*maxt, ncol = 6))

colnames(pd) = c("Outcome", "Time", "obs", "Posteriors", "p5", "p95")
pd[,1] = rep(c("Engaged","Disengaged","Transferred","Died"), each = maxt)
pd[,2] = rep(1:maxt, 4)
pd[,3:6] = rbind(t(res[,,1]),t(res[,,2]),t(res[,,3]),t(res[,,4])) 
pd$Outcome = factor(pd$Outcome, levels = c("Engaged","Disengaged","Transferred","Died"))

library("ggplot2")

p=ggplot(pd, aes(Time)) + 
  geom_line(aes(y=Posteriors),colour="red") + 
  geom_point(aes(y=obs), colour="blue") +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.5)+  facet_wrap(~ Outcome, scales="free_y")

png("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft2\\Accuracy Comparison - Diag Plots\\11-05\\Accuracy.png",width = 600, height = 400)
plot(p)
dev.off()

### R eps does not support half transparency -- turn plot from png to eps manually online
#setEPS()
#postscript("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft2\\Accuracy Comparison - Diag Plots\\11-05\\Accuracy.eps")
#plot(p)
#dev.off()

### accuracy plot for S using second part of BartGcomp_18-10-27_Athena3_Eval and loadtest.sh
### plot from accuracy
path = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"
load(paste0(path,"\\resRL.RData"))

#maxt = 10
pd = as.data.frame(matrix(NA, nrow = 2*maxt, ncol = 6))

colnames(pd) = c("Outcome", "Time", "obs", "Posteriors", "p5", "p95")
pd[,1] = rep(c("Rt","Lt"), each = maxt)
pd[,2] = rep(1:maxt, 2)
pd[,3:6] = rbind(t(res[,,1]),t(res[,,2])) 
pd$Outcome = factor(pd$Outcome, levels = c("Engaged","Disengaged","Transferred","Died"))

library("ggplot2")

p=ggplot(pd, aes(Time)) + 
  geom_line(aes(y=Posteriors),colour="red") + 
  geom_point(aes(y=obs), colour="blue") +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.5)+  facet_wrap(~ Outcome, scales="free_y")

png("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft3\\Accuracy Comparison - Diag Plots\\11-05\\AccuracyRL.png",width = 600, height = 400)
png("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft3\\Accuracy Comparison - Diag Plots\\11-05\\AccuracyRL.png",width = 600, height = 250)
plot(p)
dev.off()
##################################################################################

### Gcomp simulation
### Gcomp.sh

#!/bin/bash
#SBATCH -J Name
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --time 72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

i=1

for A in `seq 0 1 2`; do

echo tmp${i}.out

R CMD BATCH --vanilla "--args $A"  BartGcomp_18-10-27_Athena4_Gcomp.R  ./tmp${i}.out

i=$((i+1))

done

### Gcomp barplot
### Gcomp.sh
  
#!/bin/bash
#SBATCH -J Name
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --time 72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

i=1

for A in `seq 0 1 2`; do

echo tmp${i}.out

R CMD BATCH --vanilla "--args $A"  BartGcomp_18-10-27_Athena5_GcompPlot.R  ./tmp${i}.out

i=$((i+1))

done


### plot from Gcomp counts

A = #0,1,2; 2 for 1{CD4<=350}
path = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"
load(paste0(path,"\\Counts_A",A,".RData"))

png(paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft2\\Accuracy Comparison - Diag Plots\\11-05\\CondSEngA",A,".png"),width = 600, height = 400)
barplot(counts,
        xlab="Days since baseline", ylab="Number of observations", col=1:4)
legend("bottomleft", 
       legend = rownames(counts), 
       fill = 1:4)
dev.off()

####################################################################################################
### TRANSITION PROBABILITY FROM ENG/DISENG

### Gcomp.sh

#!/bin/bash
#SBATCH -J Name
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --time 72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

i=1

for A in `seq 0 1 2`; do

echo tmp${i}.out

R CMD BATCH --vanilla "--args $A"  BartGcomp_18-10-27_Athena5_TPPlot.R  ./tmp${i}.out

i=$((i+1))

done

### Plot in R
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

A = #0,1,2; 2 for 1{CD4<=350}
path = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"
load(paste0(path,"\\TP_A",A,".RData"))

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

png(filename = paste0("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft2\\Accuracy Comparison - Diag Plots\\11-05\\TPtime_A",A,".png"), width =800, height = 300)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=1,widths=c(10, 2))
dev.off()

####################################################################################################
#!/bin/bash
#SBATCH -J Name
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --time 72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

R CMD BATCH --vanilla BartGcomp_18-10-27_Athena5_DiffCI.R DiffCI.out


### plot using tab from DiffCI.R
### MARGINAL PROBABILITY OF STATE MEMBERSHIP (WITH 95% CI)
### "Treat Immediately v.s Treat when CD4 < 350"

path = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"
load(paste0(path,"\\tab.RData"))

png(filename = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft2\\Accuracy Comparison - Diag Plots\\11-05\\DiffCIwTitle.png", width =800, height = 650)

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
       ylim=range(c(pt)),
       pch=19, xlab="Time", ylab=paste0("Differece in Percent ",cname[k]) ,cex.lab=1.5)
  # hack: we draw arrows but with very special "arrowheads"
  arrows(x, pt[2,], x, pt[3,], length=0.05, angle=90, code=3)
  abline(h = 0, lty = 2)
}

mtext(expression(plain("Treat Immediately v.s. Treat when CD4<350 cells/mm")^3), side = 3, line = -2.7, outer = TRUE,cex = 2)

dev.off()

##################################################################################

### freq table of outcomes in the training set for each model  

taby = function(x){
  return(unlist(lapply(1:4,function(s) sum(x==s))))
}

tab = matrix(NA,nrow = 19, ncol = 6)
colnames(tab) = c("t","sum","1","2","3","4")
tab[,1] = c(1,rep(2:10,each = 2))

for(t in 1:10){
  if(t==1){
    
    fitdat = trd[time==200*t,]
    # [S1 | V, A0, X1, L1]
    
    tab[t,3:6] = taby(fitdat$outcome)
    tab[t,2] = sum(tab[t,3:6])
    
  } else {
    
    fitdat = trd[time==200*t,]
    
    # [St | V, A_{t-1}, Xt, Lt, S_{t-1}]
    
    for(prevS in c(0,1)){ #### S_{t-1} = 0 - disengaged, 1 - engaged
      
      ## trd
      d = fitdat[get(theY_ind_lag1) == prevS,]
      loc = 2*(t-1)+prevS
      tab[loc,3:6] = taby(d$outcome)
      tab[loc,2] = sum(tab[loc,3:6])
    } # prevS
    
  }
}


save(tab,file = "tab.RData")

##################################################################################

# making report for 02-23_DiffTreesDiffP

## create tables
LL = c(3,1,2) # Table 1
LL = c(7,1,8,9,4,5,6) # Table 2

for(k in LL){
  load(paste0("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS RData/res",k,".RData"))
  assign(paste0("res",k), res)
}

vec = matrix(NA,length(LL)+1,4)
for(k in 1:length(LL)){
  
  res = get(paste0("res",LL[k]))
  res = round(res*100,2)
  
  for(j in 1:4){
    vec[k+1,j] = paste0(res[2,j],"% (", res[3,j],"%, ",res[4,j],"%)")
  }
  
}

vec[1,] = paste0(round(res1[1,]*100,2),"%")
print(vec,quote = F)

## create convergence plots
library(GcompBART)

LL = c(3,1,2) # Table 1
LL = c(7,1,8,9,4,5,6) # Table 2
for(s in LL){
  load(paste0("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS RData/Mod3Run",s,".RData"))
  png(paste0("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/reports/Exploratory/19-02-23_DiffTreesDiffP/pngs/Diag_",s,"_D.png"),width = 600, height = 400)
  DiagPlot(Bmod2[[1]],1)
  dev.off()
  png(paste0("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/reports/Exploratory/19-02-23_DiffTreesDiffP/pngs/Diag_",s,"_E.png"),width = 600, height = 600)
  DiagPlot(Bmod2[[2]],1)
  dev.off()
}

## pull the "All 20 Trees" model and do variable selection for Eng model
### Plot
load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS RData/Mod3Run7.RData")
Bmod2[[2]]$Inclusion_Proportions
bfit = Bmod2[[2]]
pm1 = bfit$ndim

png("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/reports/Exploratory/19-02-23_DiffTreesDiffP/pngs/IncProp.png",width = 800, height = 500)

par(mfrow=c(1,pm1),mar=c(14,4,4,4))

for(i in 1:pm1){
  tmp = bfit$Inclusion_Proportions[[i]]
  barplot(sort(tmp,decreasing = T), ylab = "Inclusion Proportion", las=2, 
          cex.axis=1.5, cex.names=1.5, main=paste0("latent ",i))
  if(i > 1)  abline(h=c(0.03,0.05), col="red")
}
dev.off()
par(mfrow=c(1,1))

### create convergence plots
library(GcompBART)

LL = c(10,11) # Table 3

for(s in LL){
  load(paste0("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS RData/Mod3Run",s,".RData"))

  png(paste0("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/reports/Exploratory/19-02-23_DiffTreesDiffP/pngs/Diag_",s,"_E.png"),width = 600, height = 600)
  DiagPlot(Bmod2[[2]],1)
  dev.off()
}

### make accuracy table for engagement model
for(k in LL){
  load(paste0("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS RData/res",k,".RData"))
  assign(paste0("res",k), res)
}

vec = matrix(NA,length(LL)+1,4)
for(k in 1:length(LL)){
  
  res = get(paste0("res",LL[k]))
  res = round(res*100,2)
  
  for(j in 1:4){
    vec[k+1,j] = paste0(res[2,j]," (", res[3,j],", ",res[4,j],")")
  }
  
}

vec[1,] = res[1,]
print(vec,quote = F)

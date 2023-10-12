### Plots

### Functions
scplot = function(x,y, sp = 0.7, xlab = NULL, ylab = NULL, main = NULL){
  plot(x, y, type = "l", xlab = xlab, ylab = ylab, main = main)
  lines(loess.smooth(x, y, span = sp), col = "red", lwd = 3, lty = 3)
}


TDPlot <- function(bfit0, j, all=T){ # j is the reference level
  lev = bfit0$releveled
  a0 = bfit0$Tree_Depth
  nd = bfit0$ndraws; nb = bfit0$burn
  tnd = nb + nd
  
  if(all){
    
    for(i in 1:pm1){
      scplot(1:tnd, a0[i,],
             xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
             main = paste0("L",lev[i]," - L",j) )
    }
    
  } else {
    
    for(i in 1:pm1){
      scplot(1:nd, a0[i,(nb+1):tnd], 
             xlab = "MCMC Iteration", ylab = "Average Tree Depth", 
             main = paste0("L",lev[i]," - L",j) )
    }
    
  }
  
}


### 1. Tree Depth MCMC npost
### 2. Bar Histogram of Observed Outcomes
### 3. Observed Transtion Probability Table
### 4. Observed Transtion Probability Over Time
### 5. Observed Treatment Status

source("BartGcomp_20-01-12_Athena0_prepare.R")
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
#install.packages("RColorBrewer")
library(RColorBrewer)
#options(scipen = 10)

### 1. Tree depth MCMC npost 

setting = 4# 1:4 in prepare.R

for(t in 1:10){
  load(paste0("Cat20-01-12_t",t,"_setting",setting,".RData"))
  
  if(t == 1){
    
    mod = Bmod2
    pm1 = mod$ndim
    pdf(paste0("TDplot_t",t,".pdf"),width = pm1*3.5, height = 3.5)
    par(mfrow=c(1,pm1))
    TDPlot(mod, baseyl, all=F)
    dev.off()
    
  } else {
    
    for(prevS in 1:2){
      
      mod = Bmod2[[prevS]]
      pm1 = mod$ndim
      pdf(paste0("TDplot_t",t,"prevS",prevS,".pdf"),width = pm1*3.5, height = 3.5)
      par(mfrow=c(1,pm1))
      TDPlot(mod, baseyl[[prevS]], all=F)
      dev.off()

    }
  }
  
  remove(baseyl, Bmod2)
}

### 2. Bar Histogram of Observed Outcomes

a = newdat$outcome
a = as.factor(a)
a = relevel(a, "2")
levels(a) = c("Engaged","Disengaged","Transferred Out","Died")
counts = table(a,newdat$time)
pdf("ObsEng.pdf",width =8, height = 6)
barplot(counts, main="Observed Outcome Status",
        xlab="Days since Baseline", ylab="Number of Observations", col=brewer.pal(n = 4, name = "Set2"),
        legend = rownames(counts))#col=1:length(levels(a))
#display.brewer.pal(n = 4, name = 'Set2')
dev.off()

### 3. Observed Transtion Probability Table

tab = table(newdat$outcome_2_lag1,newdat$outcome)
tab[1,] = tab[1,]/sum(tab[1,])
tab[2,] = tab[2,]/sum(tab[2,])

prop = table(newdat$outcome)
prop = prop/sum(prop)
tp = matrix(0,5,4)
rownames(tp) = c("Overall","Disengaged","Engaged","Transferred Out","Died")
colnames(tp) = c("Disengaged","Engaged","Transferred Out","Died")
tp[1,] = prop
tp[2:3,] = tab
tp[4,3] = tp[5,4] = 1
round(tp,3)

### new table 3

tab = table(newdat$outcome_2_lag1,newdat$outcome)
tab[1,] = tab[1,]/sum(tab[1,])
tab[2,] = tab[2,]/sum(tab[2,])

tp = matrix(0,4,5)
rownames(tp) = c("Disengaged","Engaged","Transferred Out","Died")
colnames(tp) = c("Disengaged","Engaged","Transferred Out","Died","Count")
tp[1:2,1:4] = tab
tp[3,3] = tp[4,4] = 1
tp[1:2,5] = table(newdat$outcome_2_lag1)
tp[3:4,5] = table(newdat$outcome)[3:4]
round(tp,3)

nrow(newdat); length(unique(newdat$myID))


### 4. Observed Transtion Probability Over Time

tmp = unique(newdat$time)[-1]

engtp = matrix(NA,length(tmp),4)
colnames(engtp) = c("Disengaged","Engaged","Transferred Out","Died")
rownames(engtp) = tmp
dengtp = matrix(NA,length(tmp)-1,4)
colnames(dengtp) = c("Disengaged","Engaged","Transferred Out","Died")
rownames(dengtp) = tmp[-1]


for(t in 1:length(tmp)){

  td = newdat[time == tmp[t] & outcome_2_lag1==1,]
  tab = table(td$outcome_2_lag1,td$outcome)
  engtp[t,as.numeric(colnames(tab))]= tab/sum(tab)
 
  if(t > 1){
    td = newdat[time == tmp[t] & outcome_2_lag1==0,]
    tab = table(td$outcome_2_lag1,td$outcome)
    dengtp[t-1,as.numeric(colnames(tab))] = tab/sum(tab)
  }
  
}

engtp[is.na(engtp)] = 0; dengtp[is.na(dengtp)] = 0;

dd = engtp
plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                     group = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
p1 = ggplot(plotdat, aes(time, transition_prob))+ geom_line(aes(group = group,col=group),size=1.2)+
  labs(title = "Transition from Engagement")+ labs(y = "Transition Probabilities")+
  labs(x = "Days since Baseline")+
  theme(plot.title = element_text(size=22),text = element_text(size=20),
        axis.text.x = element_text(size = 20))


dd = dengtp
plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                     group = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
p2 = ggplot(plotdat, aes(time, transition_prob))+ geom_line(aes(group = group,col=group),size=1.2)+
  labs(title = "Transition from Disengagement")+ labs(y = "Transition Probabilities")+
  labs(x = "Days since Baseline")+
  theme(plot.title = element_text(size=22),text = element_text(size=20),
        axis.text.x = element_text(size = 20))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1)

pdf("obsTPtime.pdf",width =12, height = 4)
grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=1,widths=c(10, 2))
dev.off()

### 5. Observed Treated Proportion

tmp = unique(newdat$time)[-1]

tab = matrix(NA,3, length(tmp))
for(t in 1:length(tmp)){
  # overall on trt
  td = newdat[time == tmp[t] , .N , onarv_cf_1]
  a = td$N / sum(td$N)
  tab[1,t] = a[which(td[,1]==1)]   
  
  # on trt | prevS = Diseng
  if(t == 1){
    tab[2,t] = NA
  } else {
    td = newdat[time == tmp[t] & outcome_2_lag1==0, .N , onarv_cf_1]
    a = td$N / sum(td$N)
    tab[2,t] = a[which(td[,1]==1)]   
  }
  
  # on trt | prevS = Eng
  td = newdat[time == tmp[t] & outcome_2_lag1==1, .N , onarv_cf_1]
  a = td$N / sum(td$N)
  tab[3,t] = a[which(td[,1]==1)]   
  
}

dd = t(tab); rownames(dd) = tmp; colnames(dd) = c("Overall", "Disengaged at t-1", "Engaged at t-1")
plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                     group = rep(colnames(dd),each = nrow(dd)),onarv = c(dd))
p1 = ggplot(plotdat, aes(time, onarv))+ geom_line(aes(group = group,col=group),size=1.2)+
  labs(title = "Proportion of On Treatment")+ labs(y = "Probabilities")+
  labs(x = "Days since Baseline")+
  theme(plot.title = element_text(size=22),text = element_text(size=20),
        axis.text.x = element_text(size = 20))+
  scale_color_brewer(palette="Set2")

# barplot
#p1 = ggplot(plotdat, aes(x=time, y=onarv,fill=group)) +
#  geom_bar(stat="identity", position=position_dodge())+
#  scale_fill_brewer(palette="Set2")+
#  theme_minimal()

pdf("OnARV.pdf",width =8, height = 4)
plot(p1)
dev.off()
####################################################################################################

### 6. Summary table of covariates 
#path = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART"
path = "/Volumes/GoogleDrive/My Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART"

source(paste0(path,"/BartGcomp_19-12-13_Reclean_AMPATH_data_Functions.R"))

load("/Volumes/GoogleDrive/My Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/makeGridData/my working datasets/AMPATH_12-13-2019_dt4.RData")

varlist =c("gender","ageinit","marriedinit","yrenrol","TravelTime","whostageinit","weightinit","heightinit","viralloadinit",
           "cd4_preB","cd4_0","onarv_0")
tmp = dt4[time == 0, ]
tab = TabCovSumm2(varlist, tmp)

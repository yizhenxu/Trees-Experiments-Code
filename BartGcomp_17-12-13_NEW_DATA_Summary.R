### make transition matrix averaged over time
### outcome bar histogram over time
### transition prob over time
### from observed data

setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")
setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")

### PACKAGES

library(lubridate)# date processing
library(data.table)
library(zoo)# carry backward forward
library(foreign)#stata read in
library(nnet)# multinomial log-linear models
library(Hmisc)#rmultinom function
# plotting
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)


options(scipen = 10)

#setwd("./Topic2")
dat = fread("my_analysis_data.csv")
length(unique(dat$patient_id))
dim(dat)
colnames(dat)

# outcome: 0 disengage, 1 engage, 2 transfer_outcome, 3 death_outcome
######################################################################################################################
### DEFINE LAG OUTCOME
#dat = dt6
theY = "outcome"

### DEFINE LAGGED THEX AND THEY
# USE DATA.TABLE
lagfun = function(data,vname){
  vname1 = paste0(vname,"_1")
  data[, (vname1) := c(NA,get(vname)[-.N]), by = patient_id]
  return(data[,get(vname1)])
}

lagX = lagfun(dat,theY)  
vname1 = paste0(theY,"_1")
dat[,(vname1) := lagX]


dat$outcome = as.factor(dat$outcome)
levels(dat$outcome) =c("disengaged","engaged","transferred out","died")

dat$outcome_1 = as.factor(dat$outcome_1)
levels(dat$outcome_1) =c("disengaged","engaged","transferred out","died")

######################################################################################################################
### BAR HIST OF OBSERVED OUTCOME
#png(filename = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft8\\ObsEng.png", width =500, height = 400)
#png("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft7\\ObsEng.png", width =500, height = 400)
png(filename = "C:\\Users\\Yizhen Xu\\Desktop\\Draft8\\ObsEng.png", width =550, height = 400)

a = dat$outcome
a = relevel(a, "engaged")
counts <- table(a,dat$time)
barplot(counts, main="Observed Outcome Status",
        xlab="Days since baseline", ylab="Number of observations", col=1:length(levels(a)),
        legend = rownames(counts))
dev.off()
######################################################################################################################
# transition prob

tab = table(dat$outcome_1,dat$outcome)
tab[1,] = tab[1,]/sum(tab[1,])
tab[2,] = tab[2,]/sum(tab[2,])

prop = table(dat$outcome)
prop = prop/sum(prop)
tp = matrix(0,5,4)
rownames(tp) = c("overall","diengaged","engaged","transferred out","died")
colnames(tp) = c("diengaged","engaged","transferred out","died")
tp[1,] = prop
tp[2:3,] = tab
tp[4,3] = tp[5,4] = 1
round(tp,3)


# transition prob over time

tmp = unique(dat$time)[-1]

engtp = matrix(NA,length(tmp),4)
colnames(engtp) = c("diengaged","engaged","transferred out","died")
rownames(engtp) = tmp
dengtp = matrix(NA,length(tmp)-1,4)
colnames(dengtp) = c("diengaged","engaged","transferred out","died")
rownames(dengtp) = tmp[-1]

erestab = vector("list", length = length(tmp))
drestab = vector("list", length = length(tmp))


for(t in 1:length(tmp)){
  
  td = dat[time == tmp[t] & outcome_1=="engaged",]
  tab = table(td$outcome_1,td$outcome)
  
  ptab = lapply(1:4, function(k) tab[k,]/sum(tab[k,]))
  erestab[[t]]  = t(simplify2array(ptab))
  
  td = dat[time == tmp[t] & outcome_1=="disengaged",]
  tab = table(td$outcome_1,td$outcome)
  
  ptab = lapply(1:4, function(k) tab[k,]/sum(tab[k,]))
  drestab[[t]]  = t(simplify2array(ptab))
  
  engtp[t,] = erestab[[t]][2,]
  
  if(t>1) dengtp[t-1,] = drestab[[t]][1,]
  
}

dd = engtp
plotdat = data.frame(time = rep(as.numeric(rownames(dd)),ncol(dd)), 
                     group = rep(colnames(dd),each = nrow(dd)),transition_prob = c(dd))
p1 = ggplot(plotdat, aes(time, transition_prob))+ geom_line(aes(group = group,col=group),size=1.2)+
    labs(title = "Transition from engagement")+ labs(y = "Transition Probabilities")+
    theme(plot.title = element_text(size=22),text = element_text(size=20),
          axis.text.x = element_text(size = 20))


dd = dengtp
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

png(filename = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft8\\TPtime.png", width =800, height = 300)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=1,widths=c(10, 2))
dev.off()

####################################################################################################

### FIGURE 2.4: gap between RVD and Record dates
dat[,Rj_Aj := returnvisitdate - apptdate]

png(filename = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft8\\Rj_Aj.png", width =500, height = 300)
hist(as.numeric(dt1$Rj_Aj),xlim = c(0,200),xlab = "Days",breaks = 1000,main = "Scheduled Next Return Time - Current Visit Time")
dev.off()



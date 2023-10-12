setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\Stata Process")
setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\Stata Process")

library(readstata13)
library(lubridate)
library(data.table)

# merged data of demog and visits, done in STATA, with abnormal death eliminated
dat = read.dta13("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\denom_all.dta")
#dat = read.dta13("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\denom_all.dta")

colnames(dat)
var_enrol = c("state_province", "initclinic", "male","yrenrol","ageatenroll")
var_time  = c("arvstart","apptdate","returnvisitdate","death","transfer", "numadmeas","arvperfectad", "clinic","weight","height",  "cd4", 
              "whostage", "onarv","married", "age", "viralload", "visittype" )

dat = dat[,c("patient_id","enroldate",var_enrol,var_time)] # rearrange the order of the variables


#################################### # process state province
# everyone starts trt at/after enrollment visit
all(dat$enroldate <= dat$arvstart,na.rm = T) # TRUE


###########################################################
### distribution of visit number , length of follow-up, their ratio, and gaps between appts

# define last vist indicator, visit or record number, length of follow up 
dt = data.table(dat)
dt[, num := 1:.N, by=patient_id] # record number
dt[, lvisit := c(rep(0,.N-1),1), by=patient_id] # last visit indicator
visitnum =  dt[,num[.N],by=patient_id]
latime = dt[,apptdate[.N],by=patient_id] # last return visit date
fatime = dt[,apptdate[1],by=patient_id] # first appointment date
followup = as.numeric(latime$V1 - fatime$V1) # length of followup
#sum(!is.na(followup))
#sum(followup < 0 ,na.rm=T)#0 problematic follow-up periods

dt[,gapA := c(NA,apptdate[-.1]-apptdate[-.N]),by=patient_id] # gap between appointment dates

# histogram of record number, followup period, average gap between records
png('C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-09-22_Data Processing and Description\\Raw1.png',
    width = 900,height=750)
par(mfrow=c(2,2))
hist(visitnum$V1,breaks=80,xlab = "Number of Visits",
     main="Number of Visits within Patient",cex.main=1.8,cex.lab=1.5)
text(60,10000,"Range = 1:115\n 96,045 patients in the plot",cex = 1.5)
hist(followup,breaks=80,xlab = "Days",
     main="Length of Follow-up",cex.main=1.8,cex.lab=1.5)
#range(followup) #0 3000
text(1800, 14000, "Last Appt Date - First Appt Date\n Range = (0, 3000)\n 96,045 patients in the plot",cex = 1.5)
hist(followup/visitnum$V1,breaks=100,xlab = "Days",
     main="Average Gap Time between Visits",cex.main=1.8,cex.lab=1.5)
#range(followup/visitnum$V1,na.rm = T)# 0 1474
abline(v=c(50,200))
text(900, 20000, "Length of Follow-up / Number of Visits \n Range = (0, 1474)\n Vertical lines at 50 & 200\n 96,045 patients in the plot",cex = 1.5)
gapAhist = hist(dt$gapA[!is.na(dt$gapA) & dt$gapA<300],breaks=seq(0,300,1), xlab="Days",
     main="Gap Time between Adjacent Appointment Dates",cex.main=1.8,cex.lab=1.5)
#sum(dt$gapA>300,na.rm = T) #12399
#gapAhist$breaks[tail(order(gapAhist$counts), 5)]
#range(dt$gapA,na.rm=T)
text(200, 200000, paste0("Appt Date j+1 - Appt Data j \n Range = (1, 2948)\n",sum(dt$gapA>300,na.rm = T)," records have gaps > 300 days\n 5 modes (days): 41 13 83 55 27" ),cex = 1.5)
dev.off()
par(mfrow=c(1,1))

table(visitnum$V1)

numtab = matrix(0,nrow=4,ncol=10)
numtab[1,] = 1:10
numtab[3,] = 11:20
numtab[2,] = table(visitnum$V1)[1:10]
numtab[4,] = table(visitnum$V1)[11:20]
rownames(numtab) = c("Total Visits","Patient #","Total Visits","Patient #")
###########################################################
### Problematic Ri's

sum(is.na(dt$returnvisitdate))#125338
mean(is.na(dt$returnvisitdate)) #0.07427811

#dt = data.table(dat)

#check Rj-Aj
dt[, Rj_Aj := returnvisitdate - apptdate]
sum(dt$Rj_Aj < 0,na.rm=T) #88
sum(dt$Rj_Aj == 0,na.rm=T) #0
bad1 = dt[dt$patient %in% dt$patient_id[which(dt$Rj_Aj < 0 & !is.na(dt$Rj_Aj))],]
dt[, flag:= 1*(Rj_Aj < 0)]

dt[, Alag1 := c(NA,apptdate[-.N]),by=patient_id]
dt[, Rlag1 := c(NA,returnvisitdate[-.N]),by=patient_id]
dt[, gapR:= c(NA,returnvisitdate[-.1]-returnvisitdate[-.N]), by=patient_id]
dt[, gapRA:= c(NA,as.numeric(apptdate[-.1])-as.numeric(returnvisitdate[-.N])), by=patient_id]

sum(dt$gapR==0,na.rm = T)
sum(dt$gapR < 0,na.rm = T)
sum(dt$gapRA == 0,na.rm = T)
sum(dt$gapRA < 0,na.rm = T)

# how many patients have R_j missing for non-enrollment visit
temp = dt[num>1,]
tempID = unique(temp$patient_id[is.na(temp$returnvisitdate)]) 
length(tempID)#46975
46975/length(unique(dt$patient_id))#0.4890937
#show = dat[dat$patient_id==496337,] # Figure 2
#show = dat[dat$patient_id==4233,] # Figure 3

###########################################################
### replace those Rj smaller than Aj with NA

dt$returnvisitdate[dt$Rj_Aj < 0] = NA
#check
#Rj_Aj_del = as.numeric(dt$returnvisitdate - dt$apptdate)
#sum(Rj_Aj_del < 0,na.rm=T) #0
###########################################################
### how many R1 are missing?
R1 = dt$returnvisitdate[dt$num==1] # length 96045
A1 = dt$apptdate[dt$num==1]
sum(is.na(R1))# 41655/96045 missing R1
sum(R1 < A1, na.rm=T) # 5 patients have R1<A1, find out their records
tempID = dt$patient_id[dt$num==1]
showbad = dt[dt$patient_id %in% tempID[which(R1<A1)],] # view showbad
sum(R1 > A1, na.rm=T) #54385/96045 R1>A1
sum(R1 == A1, na.rm=T) #0

###########################################################
### distribution of R1-A1 A2-A1 A2-R1 

R1_A1all = dt[,returnvisitdate[1] - apptdate[1],by=patient_id]
A2_A1all = dt[,apptdate[2] - apptdate[1],by=patient_id]
A2_R1all = dt[,apptdate[2] - returnvisitdate[1],by=patient_id]
R1_A1 = as.numeric(R1_A1all$V1)
A2_A1 = as.numeric(A2_A1all$V1)
A2_R1 = as.numeric(A2_R1all$V1)


png('C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-09-22_Data Processing and Description\\A1R1.png',
    width = 950,height=300)
par(mfrow=c(1,3))
hist(A2_A1[!is.na(A2_A1) & A2_A1>=0 & A2_A1<100 ],breaks=100,xlab = "Days",
     main="A2-A1 in [0,100)",cex.main=2,cex.lab=1.5)
# range(A2_A1,na.rm = T)#1 2948
#mean(A2_A1 >=100, na.rm=T)#0.07267048
#mean(A2_A1 <3, na.rm=T)#0.1386124
#quantile(A2_A1[!is.na(A2_A1)])
text(60, 13000, "Quantiles: \n 0%  25%  50%  75% 100%\n 1     5     14      28     2948\n 7.3% over 100\n 13.9% over 50",cex = 1.8)

hist(R1_A1[!is.na(R1_A1) & R1_A1>=0 & R1_A1<100 ],breaks=100,xlab = "Days",
     main="R1-A1 in [0,100)",cex.main=2,cex.lab=1.5)
# range(R1_A1,na.rm = T)#-345  728
#sum(R1_A1<0,na.rm =T)
#mean(R1_A1 >=100, na.rm=T)#0.01104983
text(70, 20000, "Range = (-345, 728)\n 54,390 records\n 5 negative values\n 1.1% over 100",cex = 1.8)

hist(A2_R1[!is.na(A2_R1) & A2_R1>=-50 & A2_R1<100],breaks=100,xlab = "Days",
     main="A2-R1 in [0,100)",cex.main=2,cex.lab=1.5)
#sum(A2_R1<0,na.rm =T) # 13,906
#mean(A2_R1<0,na.rm =T) # 0.2912373
#mean(A2_R1==0,na.rm =T) #0.474156
#range(A2_R1[!is.na(R2_R1)]) # -622 2941
#mean(A2_R1 >= 100,na.rm=T)#0.03013486
text(60, 20000, "Range = (-622, 2941)\n 47,752 records\n 13,906 negative values\n 3% over 100",cex = 1.8)
dev.off()

###########################################################
### distribution of Rj-Aj 

png('C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-09-22_Data Processing and Description\\R_A.png',
    width = 1150,height=1200)
par(mfrow=c(5,4))
for(j in 1:20){
  R = dt$returnvisitdate[dt$num==j] # length 96045
  A = dt$apptdate[dt$num==j]
  R_A = as.numeric(R)-as.numeric(A)
  hist(R_A[!is.na(R_A) & R_A>=0 & R_A <120],breaks=seq(1,120,1),xlab = "Days",
       main=paste0("R",j,"-A",j," in [0,120)"),cex.main=2.5,cex.lab=1.5,freq = F,
       ylim=c(0,0.5))
  print(range(R_A[!is.na(R_A)]))
}
dev.off()

R_Ahist20 = hist(R_A[!is.na(R_A) & R_A>=0 & R_A <120],breaks=seq(1,120,1),xlab = "Days",
                 main=paste0("R",j,"-A",j," in [0,120)"),cex.main=1.5,cex.lab=1.5,freq = F,
                 ylim=c(0,0.5))
R_Ahist20$breaks[tail(order(R_Ahist20$counts), 5)]
###########################################################
### clean state province
a = dat$state_province
a1 = names(table(a))
RV = c(grep("^R",a1, perl=TRUE, value=TRUE)[-c(8,19)] , grep("^r",a1, perl=TRUE, value=TRUE)[-17],
       " Rift Valley"," RIFTVALLEY","`R/VALLEY")
Bun = c("bungoma","Bungoma", "BUNGOMA","Bungoma North")
Busia = c("Busia","BUSIA","Buusia")
Cent = c("central","Central","CENTRAL","Central Province")
Coast = c("coast","Coast","COAST","COAST PROVINCE")
Eastern = c("Eastern","EASTERN","N/EASTEN","ESTREN")
Eld = c("ELDORET EAST")
Nair = c("nairobi","Nairobi","NAIROBI")
Nandi = c("nandi" ,"Nandi" ,"NANDI","Nandi North")
Nak = c( "Nakuru","NAKURU")
Nyan = c("NYANVZA","NYANXA","nyanza","Nyanza","NYANZA")
Sia = c( "Siaya","SIAYA")
Kisu = grep("^KISU|Kisu|kisu",a1, perl=TRUE, value=TRUE)
Trans = grep("Tr|TRAN|tran",a1, perl=TRUE, value=TRUE)
Gishu = grep("GISHU|Gishu|gishu|UASIN|Uasin|UG",a1, perl=TRUE, value=TRUE)[-11]
Ugan = grep("Ugan|UGAN",a1, perl=TRUE, value=TRUE)
Vihi = grep("Vihi|VIHI",a1, perl=TRUE, value=TRUE)
Poko = grep("Poko|POKO",a1, perl=TRUE, value=TRUE)
West = grep("WEAT|^Wes|^WES|^wes|^WET|^WST",a1, perl=TRUE, value=TRUE)
Kit = grep("Kit|KIT|kit",a1, perl=TRUE, value=TRUE)
Kaka = grep("kaka|Kaka|KAKA",a1, perl=TRUE, value=TRUE)
a[which(a %in% RV)] = "Rift Valley"
a[which(a %in% Bun)] = "Bungoma"
a[which(a %in% Busia)] = "Busia"
a[which(a %in% Cent)] = "Cent"
a[which(a %in% Coast)] = "Coast"
a[which(a %in% Eastern)] = "Eastern"
a[which(a %in% Eld)] = "Eldoret"
a[which(a %in% Nair)] = "Nairobi"
a[which(a %in% Nandi)] = "Nandi"
a[which(a %in% Nak)] = "Nakuru"
a[which(a %in% Nyan)] = "Nyanza"
a[which(a %in% Sia)] = "Siaya"
a[which(a %in% Kisu)] = "Kisumu"
a[which(a %in% Trans)] = "Trans Nzoia"
a[which(a %in% Gishu)] = "Uasin Gishu"
a[which(a %in% Ugan)] = "Uganda"
a[which(a %in% Vihi)] = "Vihiga"
a[which(a %in% Poko)] = "Pokot"
a[which(a %in% West)] = "Western"
a[which(a %in% Kit)] = "Kitale"
a[which(a %in% Kaka)] = "Kakamega"
a[which(a =="LODWAR")] = "Lodwar"
a[which(a =="TURKANA")] = "Turkana"
#table(a)
dat$state_province = a


###########################################################
### TWO TABLES: VAR AVAILABLE AT BASELINE / TIME-VARYING VAR

# data of var at enrollment / time-varying
dat_enrol = dat[,c("patient_id",var_enrol)]
dat_time  = dat[,c("patient_id",var_time)]

# summary function
summ = function(id,x){
  temp = tryCatch(all(x==round(x),na.rm = TRUE),  error = function(e) print(2)) #FALSE(0) if float, TRUE(1) if integer, 2 if categorical
  res = rep(NA,9)
  if(!is.Date(x)){
    if(temp==0){
      res[1] = length(unique(id[which(is.na(x))]))
      res[2] = res[1]/length(unique(id))*100
      res[3] = sum(is.na(x))
      res[4] = mean(is.na(x))*100
      res[5] = NA
      res[6] = mean(x,na.rm=TRUE)
      res[7:8] = quantile(x, probs = c(0.25,0.75), na.rm = TRUE)
      res[9] = sd(x,na.rm=TRUE)
    } else {
      if(temp==2){
        res[1] = length(unique(id[which(x=="")]))
        res[2] = res[1]/length(unique(id))*100
        res[3] = sum(x=="")
        res[4] = mean(x=="")*100
        res[5] = length(unique(x[x!=""]))
      } else{
        if(temp==1){
          res[1] = length(unique(id[which(is.na(x))]))
          res[2] = res[1]/length(unique(id))*100
          res[3] = sum(is.na(x))
          res[4] = mean(is.na(x))*100
          res[5] = length(unique(x[!is.na(x)]))
        } 
      }
    }
  } else {
    res[1] = length(unique(id[which(is.na(x))]))
    res[2] = res[1]/length(unique(id))*100
    res[3] = sum(is.na(x))
    res[4] = mean(is.na(x))*100
  }
  
  return(res)
}

# make summary table 
# data at enrollment
tab_enrol = matrix(NA,ncol=9,nrow = ncol(dat_enrol)-1)
colnames(tab_enrol) = c("missingbyperson#","missingbyperson%","missing#","missing%","Levels","Mean","25%","75%","SD")
rownames(tab_enrol) = colnames(dat_enrol)[-1]
for(i in 2:ncol(dat_enrol)){
  tab_enrol[i-1,] = summ(dat_enrol[,1],dat_enrol[,i])
}


# data of time-varying vars
tab_time = matrix(NA,ncol=9,nrow = ncol(dat_time)-1)
colnames(tab_time) = c("missingbyperson#","missingbyperson%","missing#","missing%","Levels","Mean","25%","75%","SD")
rownames(tab_time) = colnames(dat_time)[-1]
for(i in 2:ncol(dat_time)){
  tab_time[i-1,] = summ(dat_enrol[,1],dat_time[,i])
}

tab_enrol
tab_time

# frequency tables for death transfer male arvperfectad  married, BY ON ARV |||  whostage , BY ON ARV
twoltab = matrix(NA, nrow=5,ncol=4)
colnames(twoltab) = c("0,arv0","1,arv0","0,arv1","1,arv1")
rownames(twoltab) = c("death", "transfer", "male", "arvperfectad", "married")

temp = setDF(dt)[,rownames(twoltab)]
for(k in 1:nrow(twoltab)){
  if(k<4){
    # time fixed
    twoltab[k,] = c( table(temp[dt$lvisit==1,k],dt$onarv[dt$lvisit==1]))
  } else {
    # time varying
    twoltab[k,] = c( table(temp[,k],dt$onarv))
  }
  
}

twoltab

whoarv = table(dt$whostage,dt$onarv)
whoarv[,1] = whoarv[,1]/sum(whoarv[,1])
whoarv[,2] = whoarv[,2]/sum(whoarv[,2])
whoarv

# visittype by missing Rj
vtarv = table(dt$visittype,dt$missR,useNA = "always")
vtarv[,1] = vtarv[,1]/sum(vtarv[,1])
vtarv[,2] = vtarv[,2]/sum(vtarv[,2])
round(vtarv*100,4)

# CD4 by arv
plot(dt$cd4[dt$onarv==1])
a = dt$cd4[dt$onarv==1]
d <- density(a[!is.na(a)]) # returns the density data 
plot(d) # plots the results
polygon(d, col="red", border="blue")

library(ggplot2)

#log CD4 by arv
a = dt$cd4[dt$onarv==1]
a = log(a[!is.na(a)])
b = dt$cd4[dt$onarv==0]
b = log(b[!is.na(b)])
temp <- data.frame(logCD4 = c(a,b)
                  , onARV =c(rep("1",length(a)), rep("0",length(b))))
p1 = ggplot(temp, aes(x = logCD4, fill = onARV)) + geom_density(alpha = 0.5) + theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.text=element_text(size=12),
                                                                                     axis.title=element_text(size=18,face="bold"))


#VL by arv
a = dt$viralload[dt$onarv==1]
a = log(a[!is.na(a)])
b = dt$viralload[dt$onarv==0]
b = log(b[!is.na(b)])
temp <- data.frame(logVL = c(a,b)
                   , onARV =c(rep("1",length(a)), rep("0",length(b))))
p2 = ggplot(temp, aes(x = logVL, fill = onARV)) + geom_density(alpha = 0.5)+ theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.text=element_text(size=12),
                                                                                   axis.title=element_text(size=18,face="bold"))
require(gridExtra)
png('C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-09-22_Data Processing and Description\\logCD4VLbyARV.png',
    width = 700,height=300)
grid.arrange(p1,p2, ncol=2)
dev.off()

# 
table(dt$onarv[dt$lvisit==1],dt$num[dt$lvisit==1])

###########################################################
### visits plot

pid = 496337


timeofcare = as.numeric(dt$apptdate - dt$enroldate)
dat = cbind(dat,timeofcare)

onePlot <- function(dat.visit, dat.cd4, dat.state, pat, legend = T) {
  # if (legend == T) {
  #   par(mar=c(5.1, 4.1, 2.1, 8.1)) # , xpd=F)
  # } else {
  #   par(mar=c(5.1, 4.1, 2.1, 2.1)) #, xpd=F)
  # }
  
  if (dim(dat.cd4)[1]==0) {
    stop('No cd4 info!')
  }
  colors <- c("#78c679", # engaged
              "#fcae91", "#fb6a4a", "#de2d26", # disengaged
              "#E69F00", # transferred-out
              "#2c7fb8") # deseased
  # x and y limits
  # xlim <- c(-50, max(dat.state$gridtime)*200)
  xlim <- c(-50, 3000)
  ylim <- c(0, 1050)
  # cd4 plot
  par(xpd = F)
  plot(dat.cd4, xlim = xlim, xlab = 'Days since enrollment', ylab = 'CD4 count', 
       ylim = ylim, cex.axis=0.8, pch = 19, xaxt='n') #, 
  # main = paste0('pid: ', pat))
  axis(1, at = seq(0, xlim[2], by = 200), cex.axis = 0.8) #, las=2)
  # visit
  abline(v = dat.visit$date[dat.visit$onarv==1], col = 'grey')
  abline(v = dat.visit$date[dat.visit$onarv==0], lty = 2, col = 'grey')
  points(dat.cd4, pch = 19, xaxt='n')
  # time in care
  timeincare <- dat.visit$timeincare[1]
  abline(v = timeincare)
  # state
  for (i in 1:max(dat.state$gridtime)) {
    col <- colors[which(dat.state$state[i] == c(1, 21, 22, 23, 3, 5))]
    rect(200*(i-1), -50, 200*i, 0, col = col, border = NA)
    rect(200*(i-1), -50, 200*(i-1), 0)
    rect(200*(i), -50, 200*(i), 0)
  }
  # legend
  if (legend == T) {
    # par(xpd = T)
    legend(2000, 900, fill = colors, title="State",
           legend=c("Engaged","Disengaged 1", 'Disengaged 2', 
                    'Disengaged 3+', 'Transferred-out', 'Deceased'), 
           cex = 0.7, bty = "n")
  }
}


###########################################################
### Plot CD4 count value at treatment initiation by enrolment year

library(zoo)
dt[,arvdiff := c(NA,onarv[-1] - onarv[-.N]),by=patient_id]
dt[,cd4cb := na.locf(cd4,na.rm=F,fromLast=T),by=patient_id]
temp = dt[dt$arvdiff==1,c("yrenrol","cd4cb")]
setDF(temp)
temp$yrenrol = as.factor(temp$yrenrol)
temp1 = temp[temp$cd4cb<2000,]
library(ggplot2)

p <- ggplot(temp1, aes(x = cd4cb, colour = yrenrol))
p + geom_density()


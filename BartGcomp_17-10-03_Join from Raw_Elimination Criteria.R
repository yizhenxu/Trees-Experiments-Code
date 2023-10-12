setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")
setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")

library(readstata13)
library(lubridate)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

# merged data of demog and visits, done in STATA, with abnormal death eliminated
dat = read.dta13(".\\denom_all.dta")

colnames(dat)
var_enrol = c("state_province", "initclinic", "male","yrenrol","ageatenroll")
var_time  = c("arvstart","apptdate","returnvisitdate","death","transfer", "numadmeas","arvperfectad", "clinic","weight","height",  "cd4", 
              "whostage", "onarv","married", "age", "viralload", "visittype","lost2fup" )

dat = dat[,c("patient_id","enroldate","hivstatus",var_enrol,var_time)] # rearrange the order of the variables

dt = setDT(dat)

# number of cd4/VL measurement, by patient
dt[,cd4orvl := sum(!is.na(cd4)) + sum(!is.na(viralload)),by=patient_id]
# is there a cd4/VL measurement? 
dt[,cd4vlind := 1*(cd4orvl>=1),by=patient_id] 
# adjacent arv status change. 1 switch on, -1 switch off, 0 unchanged, NA missing
#dt[,arvdiff := c(NA,onarv[-1] - onarv[-.N]),by=patient_id] 
# is there arv initiation?
#dt[,arvind := 1*(max(arvdiff,na.rm=T)==1),by=patient_id] 
# is there cd4/VL measurement OR arv initiation?
#dt[,cd4vlarvind := 1*(cd4orvl>=1 || max(arvdiff,na.rm=T)==1),by=patient_id] 
# define visit num
dt[,num := 1:.N, by=patient_id]

# processing hivstatus -- convert both "" and unknown to NA
#sum(dt$hivstatus== "")
dt$hivstatus[dt$hivstatus == ""] = NA
#sum(dt$hivstatus== "unknown",na.rm=T)
dt$hivstatus[dt$hivstatus == "unknown"] = NA

# is hivstatus time varying?
# hivstatuschange 0=no known hiv status, 1=one known hiv status, larger than 1=hiv status is time varying
#dt[,hivstatuschange:= length(unique(hivstatus[!is.na(hivstatus)])),by=patient_id] 
#dt[,.N,by=hivstatuschange]
# hiv status is NOT time-varying

# tables for report
#tmp = table(dt$hivstatus[dt$num==1],dt$cd4vlind[dt$num==1],useNA="always")
#rbind(tmp,apply(tmp,2,sum))
#tmp = table(dt$hivstatus[dt$num==1],dt$cd4vlarvind[dt$num==1],useNA="always")
#rbind(tmp,apply(tmp,2,sum))

# number of cd4 measurement, by patient
#dt[,cd4allind := (sum(!is.na(cd4))>=1)*1,by=patient_id]
#tmp = table(dt$hivstatus[dt$num==1],dt$cd4allind[dt$num==1],useNA="always")
#rbind(tmp,apply(tmp,2,sum))

# delete patients with no CD4/VL measurements AND hivstate negative or NA
#length(unique(dt$patient_id))#96045
dt$hivnegmis = 1*(((is.na(dt$hivstatus)) + (!is.na(dt$hivstatus) & dt$hivstatus=="NEGATIVE"))>=1)
dt[,todelete := 1*(cd4vlind==0 & hivnegmis==1)]
#length(unique(dt$patient_id[dt$todelete==1]))#11013 = 10731+282

dt = dt[todelete==0]
#dim(dt)#1648259      33
#length(unique(dt$patient_id)) # 85032


dt[, Rj_Aj := returnvisitdate - apptdate]
mean(dt$Rj_Aj>120,na.rm = T)*100
tmp=as.numeric(dt$Rj_Aj[!is.na(dt$Rj_Aj)])

hist(tmp[tmp>0],breaks=seq(1,730,1),xlab = "Days", main="Rj-Aj",cex.main=1.5,cex.lab=1.5,freq = F)
text(400, 0.2, "Last obvious mode at 24 weeks (168 days)",cex = 1.2)

A2_A1all = dt[,apptdate[2] - apptdate[1],by=patient_id]
A2_A1 = as.numeric(A2_A1all$V1)
R1 = dt$returnvisitdate[dt$num==1]
R1_A1all = dt[,returnvisitdate[1] - apptdate[1],by=patient_id]
R1_A1 = as.numeric(R1_A1all$V1)
R2_A2all = dt[,returnvisitdate[2] - apptdate[2],by=patient_id]
R2_A2 = as.numeric(R2_A2all$V1)


R1_missing = 1*is.na(R1)
a = A2_A1[R1_missing == 1]
a = a[!is.na(a) & a<=100]
b = A2_A1[R1_missing == 0]
b = b[!is.na(b) & b<=100]
temp <- data.frame(A2_A1 = c(a,b), R1_missing =c(rep("1",length(a)), rep("0",length(b))))
p1 = ggplot(temp, aes(x = A2_A1, fill = R1_missing)) + geom_density(alpha = 0.5) + theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.text=element_text(size=12),
                                                                                     axis.title=element_text(size=18,face="bold")) + geom_vline(xintercept = c(5,10))+ggtitle("Verticle Lines at 5 & 10")
plot(p1)


a = R2_A2[R1_missing == 1]
a = a[!is.na(a) & a<=100 & a>=0]
b = R1_A1[R1_missing == 0]
b = b[!is.na(b) & b<=100 & b>=0]
c = R2_A2[R1_missing == 0]
c = b[!is.na(b) & b<=100 & b>=0]
temp <- data.frame(Days = c(a,b), Type =c(rep("R2-A2|R1=NA",length(a)), rep("R1-A1|R1!=NA",length(b))))
p2 = ggplot(temp, aes(x = Days, fill = Type)) + geom_density(alpha = 0.5) + theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.text=element_text(size=12),
                                                                                         axis.title=element_text(size=18,face="bold")) 

temp1 <- data.frame(Days = c(b,c), Type =c( rep("R1-A1|R1!=NA",length(b)), rep("R2-A2|R1!=NA",length(c))))
p3 = ggplot(temp1, aes(x = Days, fill = Type)) + geom_density(alpha = 0.5) + theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.text=element_text(size=12),
                                                                                  axis.title=element_text(size=18,face="bold")) 


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend2<-g_legend(p2)
mylegend3<-g_legend(p3)
#blank<-grid.rect(gp=gpar(col="white"))

png(filename = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-03_Log_Join from Raw\\R_AcondR1mis.png", width =700, height = 700,res=83)

grid.arrange(p2+ theme(legend.position="none"),
             arrangeGrob(blank,mylegend2,blank,nrow=1,widths=c(0.01,0.68,0.01)), 
             p3 + theme(legend.position="none"),
             arrangeGrob(blank,mylegend3,blank,nrow=1,widths=c(0.01,0.68,0.01)), 
             ncol=2,widths=c(1,0.3))
dev.off()
#############################################################

#mean(A2_A1 %in% c(1,2,3,4,5), na.rm=T)
#dt$nummis = apply(dt,1, function(x) sum(is.na(x)))
#hist(dt$nummis[dt$num==1])
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")

library(readstata13)
library(lubridate)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(data.table)

# merged data of demog and visits, done in STATA, with abnormal death eliminated
dat = read.dta13(".\\denom_all.dta")

colnames(dat)
var_enrol = c("state_province", "initclinic", "male","yrenrol","ageatenroll")
var_time  = c("arvstart","apptdate","returnvisitdate","death","transfer", "numadmeas","arvperfectad", "clinic","weight","height",  "cd4", 
              "whostage", "onarv","married", "age", "viralload", "visittype","lost2fup" )

dat = dat[,c("patient_id","enroldate","hivstatus",var_enrol,var_time)] # rearrange the order of the variables

dt = setDT(dat)

# load 2 csv files
invisit = fread(".\\raw SAS data\\hivvisit2_ET.csv")
indemog = fread(".\\raw SAS data\\hivdemog2_travel&tribe.csv")

# merge travel time and tribe
setkey(indemog, patient_id)
setkey(dt, patient_id)
dt = merge(dt, indemog, all.x = TRUE)

# merge encounter_type
invisit$apptdate = as.Date(invisit$apptdate, "%d%b%Y")
dt[,tmpid:=paste0(patient_id,"_",apptdate)]
invisit[,tmpid:=paste0(patient_id,"_",apptdate)]
invisit[, patient_id:= NULL]
invisit[, apptdate:=NULL]
setkey(dt,tmpid)
setkey(invisit, tmpid)
dt = merge(dt, invisit, all.x = TRUE)

# define visit number 
dt[, num := 1:.N, by=patient_id]


################################################################################
# more Rj-Aj summaries

dt[, Rj_Aj := returnvisitdate - apptdate]
mean(dt$Rj_Aj>120,na.rm = T)*100
tmp=as.numeric(dt$Rj_Aj[!is.na(dt$Rj_Aj)])

par(mfrow=c(1,2))
hist(tmp[tmp>0],breaks=seq(1,730,1),xlab = "Days", main="Rj-Aj",cex.main=1,cex.lab=1,freq = F)
text(400, 0.2, "Last obvious mode\n at 168 days (24 weeks)",cex = 1)

hist(tmp[tmp>200],breaks=100,main="Rj-Aj in [200,730]",xlab="Days",cex.main=1,cex.lab=1)
text(500,450,"Peak at 252 days\n (9 months)")
################################################################################
### ET Tables

# percent of Rj = NA in each category of ET
percRjNA = dt[, mean(is.na(returnvisitdate)), by = encounter_type]
# percent of HIV positive
dt$hivstatus[dt$hivstatus == "unknown"] = ""
tabhiv = table(dt$encounter_type,dt$hivstatus)
for(k in 1:nrow(tabhiv)){
  tabhiv[k,] = tabhiv[k,]/sum(tabhiv[k,])
}
tabhiv = round(as.data.frame.matrix(tabhiv)*100)

# make overall ET table
a = table(dt$encounter_type)
ETtab = data.frame(ET = names(a), Count = as.numeric(a), RjNA = round(percRjNA$V1[match(names(a),percRjNA$encounter_type)]*100))
ETtab = cbind(ETtab, tabhiv)
print(ETtab,row.names = F)

# table of hiv status
tabhivstatus = matrix(NA,nrow=2,ncol = 5)
tabhivstatus[1,1:4] = table(dt$hivstatus)
tabhivstatus[1,5] = sum(table(dt$hivstatus))
tabhivstatus[2,1:4] = table(dt$hivstatus[dt$num==1])
tabhivstatus[2,5] = sum(table(dt$hivstatus[dt$num==1]))
tabhivstatus = as.data.frame(tabhivstatus)
colnames(tabhivstatus)[1:4] = names(table(dt$hivstatus))
colnames(tabhivstatus)[5] = "Sum"
rownames(tabhivstatus) = c("# Records", "# Patients")

# how many patients have existence of missing ET?
length(unique(dt$patient_id[dt$encounter_type==""]))# 57822


### ET table for the first visit
# percent of HIV positive for the first records by ET
percRjNA1 = dt[num==1, mean(is.na(returnvisitdate)), by = encounter_type]
# percent of HIV positive for the first records by ET
tabhiv1 = table(dt$encounter_type[dt$num==1],dt$hivstatus[dt$num==1])
for(k in 1:nrow(tabhiv1)){
  tabhiv1[k,] = tabhiv1[k,]/sum(tabhiv1[k,])
}
tabhiv1 = round(as.data.frame.matrix(tabhiv1)*100)
# percent of clinic NA  for the first records by ET
percCLNA1 = dt[num==1, sum(is.na(clinic)), by = encounter_type]
# percent of CD4 measured  for the first records by ET
percCD4 = dt[num==1, sum(!is.na(cd4)), by = encounter_type]
# make ET table for first records
a = table(dt$encounter_type[dt$num==1])
ETtab1 = data.frame(ET = names(a), Count = as.numeric(a), percRjNA = round(percRjNA1$V1[match(names(a),percRjNA1$encounter_type)]*100),
                    CLNA = round(percCLNA1$V1[match(names(a),percCLNA1$encounter_type)]),
                    CD4yes = round(percCD4$V1[match(names(a),percCD4$encounter_type)]))
print(ETtab1,row.names = F)

# enrollment year of those first visits without Return visit date

table(dt$yrenrol[dt$num==1 & dt$encounter_type==""])
table(dt$yrenrol[dt$encounter_type==""])



# Pick some people to display

library(gridExtra)
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.5)),
  colhead = list(fg_params=list(cex = 0.5)),
  rowhead = list(fg_params=list(cex = 0.5)))


tmp = dt[dt$patient_id %in% dt$patient_id[dt$encounter_type=="PHCT INITIAL"],]
tmp = dt[dt$patient_id %in% dt$patient_id[dt$encounter_type=="PEP INITIAL"],]
tmp = dt[dt$patient_id %in% dt$patient_id[dt$encounter_type=="PEDS RETURN UNDER FIVE"],]
tmp = dt[dt$patient_id %in% dt$patient_id[dt$encounter_type=="116"],]

look = dt[dt$patient_id == 100573,c("patient_id","apptdate","returnvisitdate","encounter_type","lost2fup","clinic","onarv","whostage","weight","height","cd4","viralload","death")]
myt <- gridExtra::tableGrob(look, theme = mytheme)
pdf(file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient100573.pdf")
grid.draw(myt)
dev.off()

look = dt[dt$patient_id == 104421,c("patient_id","apptdate","returnvisitdate","encounter_type","clinic","onarv","whostage","weight","height","cd4","viralload","death")]
look = look[1:31,]
myt <- gridExtra::tableGrob(look, theme = mytheme)
pdf(file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient104421.pdf")
grid.draw(myt)
dev.off()


look = dt[dt$patient_id %in% c(671377,690246,693660,696564),c("patient_id","apptdate","returnvisitdate","encounter_type","clinic","onarv","whostage","weight","height","cd4","viralload","death")]
myt <- gridExtra::tableGrob(look, theme = mytheme)
pdf(file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient671377.pdf")
grid.draw(myt)
dev.off()


look = dt[dt$patient_id == 107298,c("patient_id","apptdate","returnvisitdate","encounter_type","clinic","onarv","whostage","weight","height","cd4","viralload","death")]
look = look[1:31,]
myt <- gridExtra::tableGrob(look, theme = mytheme)
pdf(file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient107298.pdf")
grid.draw(myt)
dev.off()


# 699972 IMPORTANT

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.45)),
  colhead = list(fg_params=list(cex = 0.45)),
  rowhead = list(fg_params=list(cex = 0.45)))


look = dt[dt$patient_id == 699972,c("patient_id","apptdate","returnvisitdate","encounter_type","clinic","onarv","whostage","weight","height","cd4","viralload","transfer")]
myt <- gridExtra::tableGrob(look, theme = mytheme)
pdf(file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient699972.pdf")
grid.draw(myt)
dev.off()


# why are all 116 forms the last viist?
dt[,lvisit := c(rep(0,.N-1),1),by=patient_id]
dt[, mean(lvisit==1), by = encounter_type]
mean(dt$lvisit[dt$encounter_type=="AMPATH transfer form"]==1 )
table(dt$yrenrol[dt$encounter_type=="AMPATH transfer form"])
table(dt$encounter_type[dt$lvisit==1], dt$transfer[dt$lvisit==1])
table(dt$transfer[dt$encounter_type=="AMPATH transfer form"],useNA="always")

# triage
tmp = dt[dt$encounter_type =="Triage - Vital Signs Form",]
look = dt[dt$patient_id==207628,]
table(dt$lvisit[dt$encounter_type =="Triage - Vital Signs Form"])
mean((dt$lvisit[dt$encounter_type =="Triage - Vital Signs Form"]))
#########################################################################

### section adult initial

# does each person only have one adult init? NO
tmp = dt[, sum(encounter_type == "adult initial"), by=patient_id]
table(tmp$V1)
#0     1     2     3 
#2623 92880   537     5                                            
tmp1 = dt[dt$patient_id %in% tmp$patient_id[tmp$V1 == 3],]

### Pick some people to display

look1 = dt[dt$patient_id == 108251,c("patient_id","apptdate","returnvisitdate","encounter_type","lost2fup","clinic","onarv","whostage","weight","height","cd4","viralload","death")]
look1 = look1[1:31,]
#write.csv(look1, file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient108251.csv",row.names=FALSE, na="")

myt1 <- gridExtra::tableGrob(look1, theme = mytheme)

pdf(file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient108251.pdf")
grid.draw(myt1)
dev.off()

look2 = dt[dt$patient_id == 417703,c("patient_id","apptdate","returnvisitdate","encounter_type","clinic","onarv","whostage","weight","height","cd4","viralload","death")]
look2 = look2[1:30,]
myt2 <- gridExtra::tableGrob(look2, theme = mytheme)
pdf(file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient417703.pdf")
grid.draw(myt2)
dev.off()

look3 = dt[dt$patient_id %in% c(665993,665995),c("patient_id","apptdate","returnvisitdate","encounter_type","clinic","onarv","whostage","weight","height","cd4","viralload","death")]
#look3 = look3[1:30,]
myt3 <- gridExtra::tableGrob(look3, theme = mytheme)
pdf(file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient665993&665995.pdf")
grid.draw(myt3)
dev.off()

mean(is.na(dt$height[dt$encounter_type == "Express Care Stable Patient Form"]))#0.998956
mean(is.na(dt$weight[dt$encounter_type == "Express Care Stable Patient Form"]))#0.9885479
mean(is.na(dt$cd4[dt$encounter_type == "Express Care Stable Patient Form"]))
mean(!is.na(dt$cd4[dt$encounter_type == ""]))
mean(!is.na(dt$viralload[dt$encounter_type == ""]))


# distribution of gaps between consecutive adult initial encounters
dt$ind = 1*(dt$encounter_type == "adult initial")
dt[, twotog := c(NA,ind[-1] - ind[-.N]) , by=patient_id]
table(dt$twotog) # ind == 1 & twotog == 0 : the second or the third of consecutive adult initial encounters
dt[, gap := c(NA, as.numeric(apptdate[-1] - apptdate[-.N])), by=patient_id]
tmp = dt$gap[dt$ind==1 & (!is.na(dt$ind)) & dt$twotog == 0 & (!is.na(dt$twotog))]
range(tmp,na.rm=T)
length(tmp)
png("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\gapsAI.png",width = 450,height=350)
hist(tmp,breaks = seq(1, 1758, by=1),freq = T,main="Gaps between Consecutive Adult Initial Encounters")
abline(v=c(30,100),col="red")
text(1000, 15, "#{second or third consecutive adult initial encounters} = 193 \n Range (1, 1758) \n Red vertical lines at 30 and 100 days",cex = 0.8)
dev.off()
# how many second/third adult initial encounters are from switching clinics?
dt[, clsw := c(NA,clinic[-1] - clinic[-.N]) , by=patient_id] # clinic switch, 0 = no switch, all others mean switching
sum(dt$clsw !=0 & !is.na(dt$clsw) & dt$ind == 1 & dt$twotog == 0,na.rm=T)
sum( dt$ind == 1 & dt$twotog == 0,na.rm=T)


# for those with ONE adult init, which visit does the adult init happen?
tmp = dt[, sum(encounter_type == "adult initial"), by=patient_id]
table(tmp$V1)
dt1ai = dt[dt$patient_id %in% tmp$patient_id[tmp$V1==1],] 
# dim 1653119      34, 92880 patients
tmp = dt1ai$num[dt1ai$encounter_type == "adult initial"]
table(tmp)

###################################################################
### PICK EXAMPLE TO DISPLAY

# what happened to those with ai encounter at the 5th - 10th visit ?
tmp = dt[dt$patient_id %in% dt$patient_id[dt$num==10 & dt$encounter_type=="adult initial"],]

look = dt[dt$patient_id == 528561,c("patient_id","apptdate","returnvisitdate","encounter_type","clinic","onarv","whostage","weight","height","cd4","viralload","death")]
look = look[1:30,]
myt <- gridExtra::tableGrob(look, theme = mytheme)
pdf(file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient528561.pdf")
grid.draw(myt)
dev.off()

# what about PHCT initial?
tmp = dt[dt$patient_id %in% dt$patient_id[dt$encounter_type=="PHCT INITIAL"],]


look = dt[dt$patient_id == 412509,c("patient_id","apptdate","returnvisitdate","encounter_type","clinic","onarv","whostage","weight","height","cd4","viralload","death")]
myt <- gridExtra::tableGrob(look, theme = mytheme)
pdf(file = "C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\reports\\Exploratory\\17-10-05_Encounter Type\\Patient412509.pdf")
grid.draw(myt)
dev.off()
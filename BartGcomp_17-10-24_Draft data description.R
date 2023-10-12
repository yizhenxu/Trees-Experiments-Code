setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData")

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

# define last visit
dt[,lvisit := c(rep(0,.N-1),1),by=patient_id]

### !!!!!!!!!!!!!!!!!!!!!!!!!!!
### update the transfer variable
dt$transfer[dt$encounter_type =="AMPATH transfer form"] = 2

length(unique(dt$patient_id))#96045

### SELECTION CRITERIA
### patients with only one adult initial
tmp = dt[, sum(encounter_type == "adult initial"), by=patient_id]
keepid_1init = tmp$patient_id[tmp$V1 == 1] # those IDs with only 1 adult init
dt1 = dt[patient_id %in% keepid_1init,]
### patients with only one adult initial in the first four records
tmp = dt1[, num[encounter_type == "adult initial"], by = patient_id]
keepid_in4 = tmp$patient_id[tmp$V1 < 5] 
dt1 = dt1[patient_id %in% keepid_in4,]

length(unique(dt1$patient_id))#92628

###########################################################
### clean state province
a = dt1$state_province
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
dt1$state_province = a


###########################################################
### !!!!!!!!!!!!!!!!!!!!!!!!
### SKIP THIS TRANSFER SECTION
### Redefine transfer

table(dt1$transfer)
# 1 for transfer within ampath
# 2 for transfer out of ampath

# one example of multiple tranfer outs
tmp = dt1[,sum(!is.na(transfer) & transfer == 2), by = patient_id]
table(tmp$V1)
which(tmp$V1==3)[1]
tmp[25589,]
dt1[patient_id == 349176,]

# look at transfer and death from the lvisit
tmp = dt1[lvisit==1,]
table(tmp$transfer)

# redefine transfer with encounter_type == "AMPATH transfer form"
table(dt1$transfer)
#1    2 
#2193 3665 
dt1$transfer[dt1$encounter_type =="AMPATH transfer form"] = 2
table(dt1$transfer)
#1    2 
#2193 3966 

tmp = dt[dt$transfer==2,]
length(unique(tmp$patient_id)) # 3670 people previously
tmp = dt1[dt1$transfer==2,]
length(unique(tmp$patient_id)) #  3864 people now

###########################################################
### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### CLEAN AND FLAG WEIGHT, HEIGH, CD4
dt1$cd4[(dt1$cd4 > 2000) & !is.na(dt1$cd4)] = -99

dt1[, wmed := median(weight,na.rm=T), by = patient_id]
dt1[, hmed := median(height,na.rm=T), by = patient_id]

dt1$weight[xor((dt1$weight < dt1$wmed-25) , (dt1$weight > dt1$wmed+25)) & !is.na(dt1$weight)] = -99
dt1$height[xor((dt1$height < dt1$hmed-20) , (dt1$height > dt1$hmed+20)) & !is.na(dt1$height)] = -99

summary(dt1$height[dt1$height >-99])
summary(dt1$weight[dt1$weight >-99])

sum(dt1$height[dt1$height >-99 & !is.na(dt1$height)] < 100, na.rm=T)

sum(dt1$weight[dt1$weight >-99 & !is.na(dt1$weight)] < 25, na.rm=T)

###########################################################
### Remake table 2.3

colnames(dt1)
varlist = c("gender","ageatenroll","married","yrenrol","initclinic","TravelTime","state_province",
            "whostage","onarv",
            "weight","height","cd4","viralload")

dt1_bl = dt1[num==1, ]


# summary function
summ = function(x){
  temp = tryCatch(all(x==round(x),na.rm = TRUE),  error = function(e) print(2)) #FALSE(0) if float, TRUE(1) if integer, 2 if categorical
  
  # float
  if(temp==0){
    res = matrix(NA,nrow = 1, ncol = 9)
    res[1] = ""
    res[2] = sum(!is.na(x))
    res[3] = round(mean(!is.na(x))*100,2)
    res[4:9] = round(as.numeric(summary(x)[1:6]),2)
    res = as.data.frame(res)
  }   
  # categorical or integer
  if(temp %in% 1:2){
    x[x == ""] = NA
    tt = table(x,useNA="always")
    if(length(tt)<=20){
      res = matrix(NA,ncol = 9, nrow = length(tt))
      res[,1] = names(tt)
      res[,2] = as.numeric(tt)
      res[,3] = round(as.numeric(res[,2])/length(x)*100,2)
      res = as.data.frame(res)
    } else {
      res = matrix(NA,nrow = 1, ncol = 9)
      res[1] = "Levels"
      res[2] = length(tt)
      res = as.data.frame(res)
    }
    
  }
  return(res)
}

# make summary table 
# data at enrollment
for(i in 1:length(varlist)){
  if(varlist[i] == "gender"){
    x = factor(dt1_bl$male,labels = c("Female","Male"))
  }
  if(varlist[i] %in% c("weight","height","cd4")){
    x = dt1_bl[[which(colnames(dt1_bl)==varlist[i])]]
    x[x == -99] = NA
  }
  if(varlist[i] %in% c("ageatenroll","married","yrenrol","initclinic","TravelTime","state_province",
                       "whostage","onarv","viralload")){
    x = dt1_bl[[which(colnames(dt1_bl)==varlist[i])]]
  }
  if(i == 1){
    tab = summ(x)
    tab = cbind(c(varlist[i],rep("",nrow(tab)-1)), tab) # add the var name col
    colnames(tab) = c("","","n","percent","min","25%","50%","mean","75%","max")
  } else {
    res = summ(x)
    res = cbind(c(varlist[i],rep("",nrow(res)-1)), res) # add the var name col
    colnames(res) = colnames(tab)
    tab = rbind.data.frame(tab,res)
  }

}

tab

####################################################################
### continuous variables time varying

library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(reshape2)

tmp = dt1[!is.na(cd4) & (cd4 > -99),c("cd4","onarv")]
a = tmp$cd4
b = tmp$cd4[tmp$onarv==1]
c = tmp$cd4[tmp$onarv==0]
data <- data.frame(CD4=c(a, b, c), 
                   Type = as.factor(c(rep("Overall",length(a)),rep("On ARV",length(b)),rep("Not On ARV",length(c)))))
p1 = ggplot(data,aes(x=Type, y=CD4)) + geom_boxplot() +
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold"))+
      scale_y_continuous(limits = c(-100, max(a))) + annotate("text", x = 2:1, y = -100, label = paste0("n =  ", c(length(b),length(c))))
plot(p1)

tmp = dt1[!is.na(viralload),c("viralload","onarv")]
a = tmp$viralload
b = tmp$viralload[tmp$onarv==1]
c = tmp$viralload[tmp$onarv==0]
data <- data.frame(ViralLoad=c(a, b, c), 
                   Type = as.factor(c(rep("Overall",length(a)),rep("On ARV",length(b)),rep("Not On ARV",length(c)))))
p2 = ggplot(data,aes(x=Type, y=ViralLoad)) + geom_boxplot() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  scale_y_continuous(limits = c(-500000, max(a))) + annotate("text", x = 2:1, y = -500000, label = paste0("n =  ", c(length(b),length(c))))
plot(p2)

tmp = dt1[!is.na(weight) & (weight > -99),c("weight","onarv")]
a = tmp$weight
b = tmp$weight[tmp$onarv==1]
c = tmp$weight[tmp$onarv==0]
data <- data.frame(weight=c(a, b, c), 
                   Type = as.factor(c(rep("Overall",length(a)),rep("On ARV",length(b)),rep("Not On ARV",length(c)))))
p3 = ggplot(data,aes(x=Type, y=weight)) + geom_boxplot() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  scale_y_continuous(limits = c(-15, max(a))) + annotate("text", x = 2:1, y = -15, label = paste0("n =  ", c(length(b),length(c))))
plot(p3)


tmp = dt1[!is.na(height) & (height > -99),c("height","onarv")]
a = tmp$height
b = tmp$height[tmp$onarv==1]
c = tmp$height[tmp$onarv==0]
data <- data.frame(height=c(a, b, c), 
                   Type = as.factor(c(rep("Overall",length(a)),rep("On ARV",length(b)),rep("Not On ARV",length(c)))))
p4 = ggplot(data,aes(x=Type, y=height)) + geom_boxplot() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  scale_y_continuous(limits = c(-15, max(a))) + annotate("text", x = 2:1, y = -15, label = paste0("n =  ", c(length(b),length(c))))
plot(p4)


blank<-grid.rect(gp=gpar(col="white"))

#png(filename = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft6\\cont_var.png", width =700, height = 700)

png(filename = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft66\\cont_var.png", width =700, height = 650)

grid.arrange(p1,p2,p3,p4, 
             ncol=2,nrow = 2)
dev.off()

### visit number
visitnum =  dt1[,num[.N],by=patient_id]
table(visitnum$V1)

png(filename = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft66\\visitnum.png", width =500, height = 300)
hist(visitnum$V1, main = "Total Number of Records Within Patient",xlab = "Number of Records")
dev.off()

### gap between RVD and Record dates
dt1[,Rj_Aj := returnvisitdate - apptdate]

png(filename = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft66\\Rj_Aj.png", width =500, height = 300)
hist(as.numeric(dt1$Rj_Aj),xlim = c(0,200),xlab = "Days",breaks = 1000,main = "Scheduled Next Return Time - Current Visit Time")
dev.off()

mean(as.numeric(dt1$Rj_Aj) < 200, na.rm=TRUE) # 99.8% smaller than 200
mean(!is.na(dt1$Rj_Aj)) # 92.7% non-missing
max(dt1$Rj_Aj,na.rm=T) # max 730

####################################################################
### follow-up period by patient

fup = dt1[, apptdate[.N] - apptdate[1], by = patient_id]

range(fup$V1)

### baseline height and weight

tmp = dt1[num==4,]
mean(is.na(tmp$weight))
mean(is.na(tmp$height))

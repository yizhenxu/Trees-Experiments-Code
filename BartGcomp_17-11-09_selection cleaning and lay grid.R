setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\my working datasets")

setwd("G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData\\my working datasets")


library(readstata13)
library(lubridate)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(data.table)
library(zoo)

#d1 = read.dta13(".\\Stata Process\\hivdemoggic2008_2016Stata.dta")
#d2 = read.dta13(".\\Stata Process\\hivvisitgic2008_2016Stata.dta")
#setDT(d1)
#setDT(d2)
#setkey(d1, patient_id)
#setkey(d2, patient_id)
#d3 = merge(d1,d2, all.y=TRUE)
#tmp = d3[dod <= enroldate,]
#length(unique(tmp$patient_id))

# merged data of demog and visits, done in STATA, with abnormal death eliminated
dat = read.dta13(".\\denom_all.dta")

colnames(dat)
var_enrol = c("state_province", "initclinic", 
              "male","yrenrol","ageatenroll")
var_time  = c("arvstart","apptdate","returnvisitdate","death","transfer", 
              "numadmeas","arvperfectad", "clinic","weight","height",  "cd4", 
              "whostage", "onarv","married", "age", "viralload", 
              "visittype","lost2fup" )

dat = dat[,c("patient_id","dod","enroldate","hivstatus",var_enrol,var_time)] # rearrange the order of the variables

dt = setDT(dat)

# load 2 csv files
invisit = fread("hivvisit2_ET.csv")
indemog = fread("hivdemog2_travel&tribe.csv")

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

# define record number 
dt[, num := 1:.N, by=patient_id]


### !!!!!!!!!!!!!!!!!!!!!!!!!!!
### update the transfer variable
# 1 for transfer within ampath
# 2 for transfer out of ampath
dt$transfer[dt$encounter_type =="AMPATH transfer form"] = 2

length(unique(dt$patient_id))#96045 the very raw

### SELECTION CRITERIA
# patients with only one adult initial
tmp = dt[, sum(encounter_type == "adult initial"), by=patient_id]
keepid_1init = tmp$patient_id[tmp$V1 == 1] # those IDs with only 1 adult init
dt1 = dt[patient_id %in% keepid_1init,]

dim(dt1)#1653119  
length(unique(dt1$patient_id))# 92880

### 1% of adult init record (578 ppl) have missing return visit date ---> DELETE THESE
tmp = dt1[encounter_type == "adult initial" & is.na(returnvisitdate)]
dt1 = dt1[!(patient_id %in% tmp$patient_id)]

dim(dt1)#1642075
length(unique(dt1$patient_id))#92299

# patients with only one adult initial in the first four records
tmp = dt1[, num[encounter_type == "adult initial"], by = 
            patient_id]
keepid_in4 = tmp$patient_id[tmp$V1 < 5] 
dt1 = dt1[patient_id %in% keepid_in4,]

dim(dt1)
length(unique(dt1$patient_id))#92050


###########################################################
### TABLE 2.1 -- RECORD INDEX OF ADULT INITIAL ENCOUNTER AFTER SELECTION CRITERIA
tmp = dt1$num[dt1$encounter_type == "adult initial"]
table(tmp)

###########################################################
### clean state province
a = dt1$state_province
a1 = names(table(a))
RV = c(grep("^R",a1, perl=TRUE, value=TRUE)[-c(8,19)] , grep
       ("^r",a1, perl=TRUE, value=TRUE)[-17],
       " Rift Valley"," RIFTVALLEY","`R/VALLEY")
Bun = c("bungoma","Bungoma", "BUNGOMA","Bungoma North")
Busia = c("Busia","BUSIA","Buusia","BUSSIA")
Cent = c("central","Central","CENTRAL","Central Province")
Coast = c("coast","Coast","COAST","COAST PROVINCE")
Eastern = c("Eastern","EASTERN","N/EASTEN","ESTREN")
Eld = c("ELDORET EAST")
EM = c("Elgeyo Marakwet", "Elgeiyo Marakwet", "elgeyo marakwet", "ELGEYO MARAKWET")
HB = c("Homa Bay", "HomaBay")
Nair = c("nairobi","Nairobi","NAIROBI")
Nandi = c("nandi" ,"Nandi" ,"NANDI","Nandi North")
Nak = c( "Nakuru","NAKURU")
Nyan = c("NYANVZA","NYANXA","nyanza","Nyanza","NYANZA")
Sia = c( "Siaya","SIAYA")
Kisu = grep("^KISU|Kisu|kisu",a1, perl=TRUE, value=TRUE)
Trans = grep("Tr|TRAN|tran",a1, perl=TRUE, value=TRUE)
Gishu = grep("GISHU|Gishu|gishu|UASIN|Uasin|UG",a1, perl=TRUE, 
             value=TRUE)[-11]
Ugan = grep("Ugan|UGAN",a1, perl=TRUE, value=TRUE)
Vihi = grep("Vihi|VIHI",a1, perl=TRUE, value=TRUE)
Poko = grep("Poko|POKO",a1, perl=TRUE, value=TRUE)
West = grep("WEAT|^Wes|^WES|^wes|^WET|^WST",a1, perl=TRUE, 
            value=TRUE)
Kit = grep("Kit|KIT|kit",a1, perl=TRUE, value=TRUE)
Kaka = grep("kaka|Kaka|KAKA",a1, perl=TRUE, value=TRUE)
a[which(a %in% RV)] = "Rift Valley"
a[which(a %in% Bun)] = "Bungoma"
a[which(a %in% Busia)] = "Busia"
a[which(a %in% Cent)] = "Cent"
a[which(a %in% Coast)] = "Coast"
a[which(a %in% Eastern)] = "Eastern"
a[which(a %in% Eld)] = "Eldoret"
a[which(a %in% EM)] = "Elgeyo Marakwet"
a[which(a %in% HB)] = "HomaBay"
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
# clean the global environment and keep only list / dataframe 
rm(list = ls.str(mode = 'character'))
rm(list = ls.str(mode = 'numeric'))

# save temporary data
save(dt1, file = "rawdat92050_1.RData")

###########################################################
### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### CLEAN AND FLAG WEIGHT, HEIGH, CD4 (BY WITHIN PATIENT MEDIANS)
dt1$cd4[(dt1$cd4 > 2000) & !is.na(dt1$cd4)] = -99

dt1[, wmed := median(weight,na.rm=T), by = patient_id]
dt1[, hmed := median(height,na.rm=T), by = patient_id]

dt1$weight[xor((dt1$weight < dt1$wmed-25) , (dt1$weight > dt1$wmed+25)) & !is.na(dt1$weight)] = -99
dt1$height[xor((dt1$height < dt1$hmed-20) , (dt1$height > dt1$hmed+20)) & !is.na(dt1$height)] = -99

#summary(dt1$height[dt1$height >-99])
#summary(dt1$weight[dt1$weight >-99])
#sum(dt1$height[dt1$height >-99 & !is.na(dt1$height)] < 100, na.rm=T)
#sum(dt1$weight[dt1$weight >-99 & !is.na(dt1$weight)] < 25, na.rm=T)

###########################################################
# IS PATIENT ID IN ORDER?
#is.unsorted(dt1$patient_id, na.rm = TRUE, strictly = FALSE)
# SORT BY PATIENT ID AND NUM TO REDUCE POSSIBLE MISTAKES IN THE FOLLOWING CLEANING PROCESS
dt1 = dt1[order(patient_id,num)]

### DEFINE NEW INDEX NEW_NUM, ADULT INITIAL RECORD AS 0, AND RECORDS BEFORE THAT AS NEGATIVE
dt1[, ai := which(encounter_type=="adult initial"), by = patient_id]
dt1$new_num = dt1$num - dt1$ai
dt1[, ai:= NULL]

###########################################################
### DEFINE RECORD_TYPE, INDICATOR OF ENCOUNTER, 1 WHEN EXISTS A RETURN VISIT DATE

### DEFINE MEASUREMENT POINTS AND ENCOUNTER PROCESS
# RECORD_TYPE =  1 AS ENCOUNTER, 0 AS NON-ENCOUNTER MEASUREMENT, -1 AS INVALID MEASUREMENT (NA IN WEIGHT, CD4, VL), 0.5 AS FORCED INIT ENCOUNTER (THOSE WITHOUT RVD)
dt1$record_type = 1*(!is.na(dt1$returnvisitdate))

#table(dt1$record_type)
#0       1 
#119015 1516620 

###########################################################
# how many patients have returnvisitdate for every record?
look = dt1[, mean(record_type), by = patient_id]
mean(look$V1==1) #32.4%,  29809 ppl out of 92050
allenc = look$patient_id[look$V1==1]
allencdat = dt1[patient_id %in% allenc,  ]
dim(allencdat)# 269788 33

### frequency of vl measurements
numvl = allencdat[, sum(!is.na(viralload)), by=patient_id]
summary(numvl$V1)
hist(numvl$V1)

### patients' behavior
allencdat[, Rj := c(NA,returnvisitdate[-.N]),by=patient_id]
allencdat[, Aj_Rj:= apptdate - Rj] # gap between actual show up and the scheduled date
a = as.numeric(allencdat$Aj_Rj)
#hist(a[a<=100 & a>1],freq=F,breaks = 100)
#hist(as.numeric(a[a<0]),freq=F,breaks=30) # early arrivals
#hist(as.numeric(a[a>0]),freq=F,breaks=30) # late ones
a[a < 1] = 1 # do not look at early arrivals for now, too complicated
#look at actual cases
hist(a[a<=20],freq=F,breaks = 100)


look1 = allencdat[, .N, by=patient_id]
b= look1$patient_id[order(look1$N,decreasing = T)]

pn = 50
sn = 100
#png(file = paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\Aj-Rj for people with complete RVD\\plot",sn,".png"),height=800,width=1000)
plot(c(0,1000), c(0,pn+1), 
      yaxt='n', ann=F,
     type="n")
for(i in (sn+1):(sn+pn)){
  look2 = a[allencdat$patient_id == b[i]]
  look2 = look2[!is.na(look2)]
  points(cumsum(look2),rep(i-sn,length(look2)),pch=22,cex = 0.3)
}
#dev.off()

### clinicians' behavior
allencdat[, Rj_Aj:= returnvisitdate - apptdate] # gap between actual show up and the scheduled date
a = as.numeric(allencdat$Rj_Aj)
#hist(as.numeric(a),freq=F)
a[a < 0] = NA # backward scheduling is wrong input
#look at actual cases
look1 = allencdat[, .N, by=patient_id]
b= look1$patient_id[order(look1$N,decreasing = T)]

pn = 50
sn = 100
#png(file = paste0("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\Aj-Rj for people with complete RVD\\plot",sn,".png"),height=800,width=1000)
plot(c(0,3000), c(0,pn+1), 
     yaxt='n', ann=F,
     type="n")
for(i in (sn+1):(sn+pn)){
  look2 = a[allencdat$patient_id == b[i]]
  look2 = look2[!is.na(look2)]
  points(cumsum(look2),rep(i-sn,length(look2)),pch=22,cex = 0.3)
}
#dev.off()


###########################################################

### CARRY FORWARD MEASUREMENTS TO ADULT INITIAL RECORD FROM NEG INDEX RECORDS, 
### IFF MEASUREMENT WAS MADE WITHIN 31 DAYS, FOR HEIGHT WEIGHT CD4 AND VL

# 1. extract patient_id and their records with new_num<=0. keep only those having new_num<0.
# 2. for each HWCV, get the carried forward value and the correponding date at which the value was measured
# 3. for day gap (adult ini date - measurement date of the carried forward value) > 30, set to NA for the corresponding HWCV at new_num==0
# 4. update to dt1, check if the nonmissing percentage increased (should be)

cfdata = function(data, tmp_before, colnum, gap){
  var = colnames(data)[colnum]
  
  # cf
  # mt for measurement time of the non-missing values of var
  tmp_before$mt = tmp_before$apptdate
  tmp_before$mt[is.na(tmp_before[,get(var)])] = NA 
  # var1 : carried forward values of var
  tmp_before[, var1 := na.locf(get(var),na.rm = F), by = patient_id]
  # mt1 : carried forward values of mt
  tmp_before[, mt1 := na.locf(mt,na.rm=F), by = patient_id] # mt1 tells you which date this carried forward record was measured
  # we will void those carried forwarded number when timediff > gap days
  timediff = as.numeric(tmp_before$apptdate - tmp_before$mt1) 
  tmp_before[new_num==0 & timediff > gap, var1 := NA]
  # update these cf meaurements to new_num==0 in data
  tmp_before1 = tmp_before[tmp_before$new_num==0] 
  data[(data$patient_id %in% tmp_before1$patient_id) & new_num==0, (var):= tmp_before1$var1]
  return(data[,get(var)])  
  
}



cbdata = function(data, tmp_after, colnum, gap){
  var = colnames(data)[colnum]
  # mt for measurement time of the non-missing values of var
  tmp_after$mt = tmp_after$apptdate
  tmp_after$mt[is.na(tmp_after[,get(var)])] = NA 
  # var1 : carried backward values of var
  tmp_after[, var1 := na.locf(get(var),na.rm = F,fromLast = TRUE), by = patient_id]
  # mt1 : carried backward values of mt
  tmp_after[, mt1 := na.locf(mt,na.rm=F, fromLast = TRUE), by = patient_id] 
  # we will void those carried forwarded number when timediff > gap days
  timediff = as.numeric(tmp_after$mt1 - tmp_after$apptdate) 
  tmp_after[new_num==0 & timediff > gap, var1 := NA]
  # update these cb meaurements to new_num==0 in data
  tmp_after1 = tmp_after[tmp_after$new_num==0] 
  data[(data$patient_id %in% tmp_after1$patient_id) & new_num==0, (var):= tmp_after1$var1]
  return(data[,get(var)])  
}


data = dt1
varlist = c("height" ,"weight" , "cd4" ,"viralload","whostage", "married" )
cfcbcols = which(colnames(data) %in% varlist)

# prepare for carry forward
tmp_before = data[new_num <= 0 ]
tmp_before[tmp_before == -99] = NA
tmp_before[, prev_num := c(NA,new_num[-.N]), by = patient_id ]
tmp_before = tmp_before[-which(tmp_before$new_num==0 & is.na(tmp_before$prev_num)),] # exclude from target those with no negative new_num

# carry forward all V and X
cfdt = lapply(cfcbcols, function(k) cfdata(data, tmp_before, k,gap=31))
vname1 = paste0(colnames(data)[cfcbcols],1)
data[, (vname1) := cfdt]

# prepare for carry backward

tmp_after = data[new_num >= 0 ]
tmp_after[tmp_after == -99] = NA
tmp_after[, post_num := c(new_num[-1],NA), by = patient_id ]
tmp_after = tmp_after[-which(tmp_after$new_num==0 & is.na(tmp_after$post_num)),] # exclude from target those with no records after adult init record

# only carry backward V
vname2 = c("height1", "whostage1",  "married1")
cfcbcols1 = which(colnames(data) %in% vname2)
#colnames(data)[cfcbcols1]

cbdt = lapply(cfcbcols1, function(k) cbdata(data, tmp_after, k,gap=31))
vname22 = colnames(data)[cfcbcols1]
data[, (vname22) := cbdt]



rm(list = ls.str(mode = 'numeric'))

# DEFINE INIT HEIGHT, WHOSTAGE, MARRIED, AGE

dt1[, heightinit := rep(height1[new_num==0],.N), by=patient_id]

dt1[, whostageinit := rep(whostage1[new_num==0],.N), by=patient_id]

dt1[, marriedinit := rep(married1[new_num==0],.N), by=patient_id]

dt1[, ageinit := rep(age[new_num==0],.N), by=patient_id]

###########################################################


save(dt1, file = "rawdat92050_2.RData")
load("rawdat92050_2.RData")



#dd = setDT(dat)
#vn = c("patient_id","apptdate","dod")
#ddd = dd[, vn,with=F]
#ddd[,tmpid:=paste0(patient_id,"_",apptdate)]

#head(ddd)
#ddd[,patient_id := NULL]
#ddd[,apptdate:= NULL]

#setkey(dt1,tmpid)
#setkey(ddd, tmpid)
#dt1 = merge(dt1, ddd, all.x = TRUE)

#save(dt1, file = "rawdat92050_2.RData")

###########################################################
### FIGURE 2.3, TABLE 2.4: RECORD NUMBER HISTOGRAM
visitnum =  dt1[,num[.N],by=patient_id]
table(visitnum$V1[visitnum$V1<21])

png(filename = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft7\\visitnum.png", width =500, height = 300)
hist(visitnum$V1, main = "Total Number of Records Within Patient",xlab = "Number of Records")
dev.off()

### TABLE 2.5: ENCOUNTER NUMBER TABLE
encnum = dt1[, sum(!is.na(returnvisitdate)), by=patient_id]
table(encnum$V1[encnum$V1<21])

### FIGURE 2.4: gap between RVD and Record dates
dt1[,Rj_Aj := returnvisitdate - apptdate]

png(filename = "C:\\Users\\Yizhen Xu\\Desktop\\Draft8\\Rj_Aj.png", width =500, height = 300)
#png(filename = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft7\\Rj_Aj.png", width =500, height = 300)
hist(as.numeric(dt1$Rj_Aj),xlim = c(0,200),xlab = "Days",breaks = 1000,main = "Scheduled Next Return Time - Current Visit Time")
dev.off()

mean(as.numeric(dt1$Rj_Aj) < 200, na.rm=TRUE) # 99.8% smaller than 200
mean(!is.na(dt1$Rj_Aj)) # 92.7% non-missing
max(dt1$Rj_Aj,na.rm=T) # max 730
###########################################################
### TABLE 2.3 SUMMARY OF ALL VARIABLES AT BASELINE

colnames(dt1)
varlist = c("gender","ageinit","marriedinit","yrenrol","initclinic","TravelTime","state_province",
            "whostageinit","onarv",
            "weight1","heightinit","cd41","viralload1")

dt1_bl = dt1[new_num==0, ]


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
# data at baseline
for(i in 1:length(varlist)){
  if(varlist[i] == "gender"){
    x = factor(dt1_bl$male,labels = c("Female","Male"))
  }
  if(varlist[i] %in% c("weight1","heightinit","cd41")){
    x = dt1_bl[[which(colnames(dt1_bl)==varlist[i])]]
    x[x == -99] = NA
  }
  if(varlist[i] %in% c("ageinit","marriedinit","yrenrol","initclinic","TravelTime","state_province",
                       "whostageinit","onarv","viralload1")){
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
write.csv(tab,"G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft7\\table2point3_new.csv",row.names=F)

### Figure 2.2: BOX PLOT OF HEIGHT WEIGHT CD4 AND VL

library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(reshape2)

tmp = dt1_bl[!is.na(cd41) & (cd41 > -99)]
colnames(tmp)[which(colnames(tmp)=="cd41")] = "CD4"
p1 = ggplot(data = tmp, aes(x = "", y = CD4)) + geom_boxplot() +
  theme(axis.title.y=element_text(size=14,face="bold"),axis.title.x=element_blank())

tmp = dt1_bl[!is.na(viralload1) ]
sum(tmp$viralload1 <= 50,na.rm=T)
tmp$viralload1[tmp$viralload1==0] = NA
tmp = tmp[!is.na(viralload1) ]
tmp$logVL = log(tmp$viralload1)
p2 = ggplot(data = tmp, aes(x = "", y = logVL)) + geom_boxplot() +
  theme(axis.title.y=element_text(size=14,face="bold"),axis.title.x=element_blank())

tmp = dt1_bl[!is.na(weight1) & (weight1 > -99)]
colnames(tmp)[which(colnames(tmp)=="weight1")] = "Weight"
p3 = ggplot(data = tmp, aes(x = "", y = Weight)) + geom_boxplot() +
  theme(axis.title.y=element_text(size=14,face="bold"),axis.title.x=element_blank())

tmp = dt1_bl[!is.na(height1) & (height1 > -99)]
colnames(tmp)[which(colnames(tmp)=="height1")] = "Height"
p4 = ggplot(data = tmp, aes(x = "", y = Height)) + geom_boxplot() +
  theme(axis.title.y=element_text(size=14,face="bold"),axis.title.x=element_blank())


png(filename = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft7\\cont_var.png", width =550, height = 400)

#png(filename = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft66\\cont_var.png", width =700, height = 650)

grid.arrange(p1,p2,p3,p4, ncol=2,nrow = 2)
dev.off()

###########################################################
# delete records after dod (259 records)
dt1 = dt1[-which(!is.na(dod) & apptdate>dod),] #1635635 ->1635376   

# define LAST ENCOUNTER & LAST RECORD
# LENC SHOULD BE THE LAST ENCOUNTER (RECORD_TYPE==1) BUT NOT THE LAST RECORD!!!!
dt1[, whichlast := last(which(record_type==1)),by=patient_id] # index of the last encounter
dt1[, lenc := ifelse(!is.na(whichlast), c(rep(0,whichlast[.N]-1),1,rep(0,.N - whichlast[.N])), 0 ) ,by=patient_id]
dt1[, whichlast:= NULL]

dt1[, lrec := c(rep(0,.N-1),1),by = patient_id]
#tmp1 = dt1[patient_id==73548]
#tmp1 = dt1[patient_id==75099]
sum(dt1$lrec)#92050
sum(dt1$lenc)#92050
###########################################################
### CLEAN XFER RECORDS
# dd = dt1

# tranfer should be pushed to the last encounter and not conflict with death
#    (a)if both death and transfer happened, it's death
#    (b)set transfer to 2 at lenc == 1 if transfer happened before AND death!=1
#    (c)set transfer to 0 for non-last-encounter records (lenc!=1)  
#look = dt1[patient_id == 78349]
#look = dt1[patient_id == 74358]

dt1$transfer[is.na(dt1$transfer)] = 0

#(a)
dt1$transfer[dt1$death==1] = 0

dt1[,transfer1 := max(transfer,na.rm = T),by=patient_id] # 1: transfer within, 2: transfer out

#(b)
dt1$transfer[dt1$lenc == 1 & dt1$transfer1==2] = 2

#(c)
dt1$transfer[dt1$lenc != 1] = 0
###########################################################


save(dt1, file = "rawdat92050_3.RData")
load("rawdat92050_3.RData")


###########################################################
### KEEP ONLY THOSE NEW_NUM >= 0 ---> define and use dt2 for data analysis 

dt2 = dt1[new_num >=0,c("patient_id","onarv","arvstart","apptdate","returnvisitdate",
                        "death", "dod", "transfer","transfer1", "lenc", "lrec","lost2fup",
                        "ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit",
                        "weight1", "cd41", "viralload1","new_num","record_type")]
#dt2 = dt1[new_num>=0,]
length(unique(dt2$patient_id)) # 92050

dt2[ dt2 == -99]  = NA

#mean(!is.na(dt2$whostage[dt2$new_num==0])) # 0.7062792
#mean(!is.na(dt2$married[dt2$new_num==0])) # 0.9500054

###########################################################
### LAY GRID
int_days = 200

# DEFINE TIME
dt2$apptdate1 = as.numeric(dt2$apptdate)
dt2[, time := apptdate1 - apptdate1[1], by=patient_id]

# DEFINE DATA CLOSURE DATE
dt2$closure = as.numeric(max(dt2$apptdate)) # 2016-8-24

# HOW MANY GRID POINT FOR EACH PERSON?
length_in_care = dt2[, closure[new_num==0] - apptdate1[new_num==0],by=patient_id]
dt2[, apptdate1 := NULL]

num_grid = floor(length_in_care$V1/int_days) # range: 0-15

# DEFINE GRID DATA
grid_id = rep(length_in_care$patient_id, times = num_grid+1) # plus 1 for the zeros in num_grid (baselines)
grid_times = unlist(lapply(num_grid,function(x) (0:x)*int_days))
grid_ind = unlist(lapply(num_grid,function(x) (0:x)))
grid_data = as.data.frame(cbind(patient_id = grid_id, time = 
                                  grid_times, grid_ind = grid_ind))
grid_data = setDT(grid_data)

# SET KEY AND MERGE 
dt2[,tmpid:=paste0(patient_id,"_",time)]
grid_data[,tmpid:=paste0(patient_id,"_",time)]
setkey(dt2,tmpid)
setkey(grid_data, tmpid)
dt3 = merge(dt2, grid_data, all = TRUE)

# CLEAN SOME VARIABLES FROM MERGING, AND SORT
dt3$patient_id  = dt3$patient_id.x
dt3$patient_id[is.na(dt3$patient_id)] = dt3$patient_id.y[is.na(dt3$patient_id)]

dt3$time  = dt3$time.x
dt3$time[is.na(dt3$time)] = dt3$time.y[is.na(dt3$time)]

dt3 = dt3[order(patient_id,time)]

dt3[, tmpid:= NULL]
dt3[, time.x:= NULL]
dt3[, time.y:= NULL]
dt3[, patient_id.x:= NULL]
dt3[, patient_id.y:= NULL]

rm(list = ls.str(mode = 'numeric'))

###########################################################

save(dt3, file = "rawgriddat_dt3.RData")
load("rawgriddat_dt3.RData")

###########################################################
### BE CAREFUL DEFINING XFER AND DEATH ON GRID
# 1. KEEP the measurement records after the last encounter for people who xferred out
#    careful with X values in LGRIID
# 2. death is the only outcome category that depends on date of death INSTEAD of apptdate
#    careful with LGRID
#tmp1 = dt3[patient_id == 771827]
#tmp1 = dt3[patient_id == 73552]
#tmp1 = dt3[patient_id == 277813]
#tmp1 = dt3[patient_id == 75022]
#tmp1 = dt3[patient_id == 78349]

# DEFINE OUTCOME RELATED VARIABLES BASED ON RECORD_TYPE AND DOD
# RECORD_TYPE =  1 AS ENCOUNTER, 0 AS NON-ENCOUNTER records/grid
#table(dt3$record_type) # only 1 and 0
#sum(dt3$lenc,na.rm=TRUE) #92050

# DEFINE TMPID = PATIENT_ID X GRID_ID 
dt3[, grid_ind_cb := na.locf(grid_ind,na.rm=F,fromLast = TRUE), by = patient_id] # grid intervals (,]
dt3$tmpgridid = paste0(dt3$patient_id,"_",dt3$grid_ind_cb)
dt3[,totrec := sum(!is.na(apptdate)), by = tmpgridid]

# ENGAGE
dt3[,totenc := sum(record_type,na.rm=T), by = tmpgridid]
dt3$engage = 1*(dt3$totenc > 0)

### DEFINE LAST GRID INDICATORS
# engagement and xfer are based on encounter -- LGRID
# while death is based on dod -- LGRID_DEATH

# LAST GRID BASED ON DOD: LGRID_DEATH

#tmp1 = dt3[patient_id == 95122,]
#tmp1 = dt3[patient_id == 75022]
#tmp1 = dt3[patient_id == 73741]
dt3[, dod1 := dod[1],by = patient_id]
dt3[, appt1:=apptdate[1],by=patient_id]
dt3[,timetodeath := as.numeric(dod1) -as.numeric(appt1)]# time to death from baseline
dt3[, dod1 :=NULL];dt3[, appt1 :=NULL];
dt3[, lgrid_death:=ifelse(timetodeath!=0, 1*(!is.na(timetodeath) & (time/timetodeath)>=1), 1) ] #timetodeath at 0 are those died at the first clinic encounter (baseline)
# lgrid_death is 0 for non-last-grids, 1s for all lines at or after last time
dt3[, lgrid_death := lgrid_death - c(0,lgrid_death[-.N]),by = patient_id ]
dt3$lgrid_death[is.na(dt3$lgrid_death)] = 0
dt3[, lgrid_death := max(lgrid_death,na.rm = T), by = tmpgridid]
# lgrid_death is 0 for non-last-grids, 1s for only lines in the last grid (where dod is)

# LAST GRID BASED ON LENC: LGRID

dt3$lenc[is.na(dt3$lenc)] = 0
dt3[, lgrid := max(lenc,na.rm = T), by = tmpgridid]

################# REDEFINE LAST GRIDS FOR TRNASFER
# when XFER happens --- 
# grid t is the last record, then grid t is engaged, grid t+1 is xfer out 
# WE ALREADY PUSHED ALL XFERs TO THE LAST ENCOUNTER (LENC==1) see line 492

#look = dt3[patient_id == 78349]
#look = dt3[patient_id == 73747]
#look = dt3[patient_id == 73566]

# all patients have no or only 1 transfer==2 so far. nobody has multiple transfer==2 (i.e.1 at lenc==1 and 1 at measurement record after lenc==1)
dt3[, numxfer := sum(transfer==2,na.rm=T), by=patient_id]
table(dt3$numxfer)
dt3[,numxfer := NULL]


# update death and transfer1 to all lines (a lot of NAs from combining grid)
dt3[,death:= max(death,na.rm = T), by=patient_id]
dt3$transfer[is.na(dt3$transfer)] = 0 # transfer from inserted grids were NA, set to 0
dt3[,transfer1 := max(transfer,na.rm = T),by=patient_id] # update transfer1; 1: transfer within, 2: transfer out

# DEFINE LGRID_XFER: INDICATOR OF LAST GRID OF XFER
# grid t is NOT xfer out, grid t+1 is xfer out 
# from numxfer, we can directly push lgrid[transfer==2] to the next record (definitely a inserted grid, not in the obs data) if it's not the last record(lrec==0)

dt3[, lgrid_xfer_tmp := 1*(transfer1==2)*lgrid] # only 1s for last grid lines of XFER, 0 o.w.
dt3[, lgrid_xfer_tmp := lgrid_xfer_tmp - c(0,lgrid_xfer_tmp[-.N]), by=patient_id]
dt3[, lgrid_xfer := 1*(lgrid_xfer_tmp == (-1))] # some might be censored from data closure
dt3[, lgrid_xfer_tmp := NULL]

# DEATH
dt3$death_outcome = dt3$lgrid_death * dt3$death # 1 for those in the grid with death occurence (last encounter), 0 o.w.
dt3$transfer_outcome = dt3$lgrid_xfer * 1 * (dt3$transfer1==2)

save(dt3, file = "rawgriddat_dt3_tmp.RData")
load("rawgriddat_dt3_tmp.RData")
###########################################################
### REFINE LAST GRID & 
### (a) EXCLUDE THOSE BEYOND THE MAXIMUM GRID OF OBSERVATIONS (ALL lgrid==0 for these ppl)
### (b) DEFINE KEEPOBS. CONSIDER 4 CASES: DEATH, XFER, NO LGRID DUE TO CENSORING, NO DEATH OR XFER WITH LGRID==1
### KEEP ALL IF CENSORED, NO DEATH AND NO XFER BEFORE MAXGRID
### (c) DELETE THOSE WE DON'T WANT TO KEEP

# (a) grid_ind_cb == NA ---> going to discard these due to data closure censoring, and PULL LGRID TO THE PREVIOUS GRID
#tmp1 = dt3[patient_id ==  803978]
#tmp1 =  dt3[patient_id == 48647]
#tmp1 = dt3[patient_id == 672531]
dt3[, closure_censor := (sum(is.na(grid_ind_cb))>0 )*1, by = patient_id] # 1 for patients who have records beyond maximum grid
# FOR CLOSURE_CENSOR==1, SHOULD LET LGRID_CENSOR := !IS.NA(GRID_IND_CB)
dt3[, lgrid_censor := 1*(closure_censor==1)*(!is.na(grid_ind_cb)), by=patient_id]


# someone who had death in the closure censor part
#sum(dt3$lgrid_death[dt3$lgrid_censor==0 & dt3$closure_censor==1])
#which(dt3$lgrid_death==1 & dt3$lgrid_censor==0 & dt3$closure_censor==1)[3]
#tmp1 = dt3[patient_id==807674]

# DELETE THOSE WITH CLOSURE_CENSOR==1 & LGRID_CENSOR==0; we are going to lost some deaths and xfers here
dt4 = dt3[-which(closure_censor==1 & lgrid_censor==0),]

sum(dt4$closure_censor==1 & dt4$lgrid_death==1)#0 
sum(dt4$closure_censor==1 & dt4$lgrid_xfer==1)#166 patients have both xfer and maxgrid censoring

ddt4 =dt4
#dt4 = ddt4
# (b) define keepobs from the four cases
# have to define xfer and death after the closure_censor --> in case xfer or death happened before maxgrid for those having cencored obs
dt4[, keepobs := 1] # keep all if neither dead nor xfer
dt4[closure_censor == 1, keepobs := lgrid_censor]
dt4[death == 1, keepobs := lgrid_death]
dt4[transfer1 == 2, keepobs := lgrid_xfer]
# SPECIAL CASES OF (b)
# (b.1) xfer record happened in the maxgrid -- xfer is censored -- all lgrid_xfer == 0 -- keep all
#tmp1 = dt4[patient_id ==73892]
dt4[,xfercensored := (sum(lgrid_xfer)==0)*1, by=patient_id]
dt4[transfer1 ==2 & xfercensored == 1, keepobs := xfercensored]
# (b.2) death date censored after maxgrid -- death is censored -- all lgrid_death == 0 -- keep all
dt4[,dodcensored := (sum(lgrid_death)==0)*1, by=patient_id]
dt4[death == 1 & dodcensored == 1, keepobs := dodcensored]


dt4[,closure_censor :=NULL]; dt4[,lgrid_censor :=NULL];
dt4[, dodcensored:=NULL]; dt4[,xfercensored:=NULL];


# (c)
dt4$keepobs[dt4$keepobs==0] = NA # prepare for carry backward 
dt4[, keepobs := na.locf(keepobs,na.rm = F,fromLast = TRUE), by = patient_id]
#table(dt4$keepobs)

dt4 = dt4[!is.na(keepobs)] # rows 2321631 ->  2208145 

length(unique(dt4$patient_id)) # 92050 ...FINALLY!!! this is correct...
###############################################################################
save(dt4, file = "rawgriddat_dt4_tmp.RData")
load("rawgriddat_dt4_tmp.RData")
 # lenc IS THE LAST ENCOUNTER (RECORD_TYPE==1) BUT NOT THE LAST RECORD!!!!

# CHECK:
#tmp1 = dt4[patient_id==73500]   ---> death. okay.
#tmp1 = dt4[patient_id==73598]  ---> 1 xfer, push lgrid to the next grid
#tmp1 = dt4[patient_id == 48647] ---> 2 descard obs
#tmp1 = dt4[patient_id == 73548] ---> 3 redefine lenc, selection criteria record_type != -1 caused problem, maybe exclude these later
#tmp1 =dt4[patient_id == 117314] ---> only baseline record. okay.
#tmp1 =dt4[patient_id == 130]    ---> keep all if closure_censor==0 & death1==0 & transfer1==0

###########################################################
### CARRY FORWARD FEATURES -> EXTRACT GRIDS -> DEFINE OUTCOME -> DEFINE COVARIATES FOR ANALYSIS

# features to carry forward: ageinit, male, yrenrol, TravelTime, whostageinit, marriedinit, heightinit   
dt4[, ageinit := na.locf(ageinit,na.rm = F), by = patient_id]
dt4[, male := na.locf(male,na.rm = F), by = patient_id]
dt4[, yrenrol := na.locf(yrenrol,na.rm = F), by = patient_id]
dt4[, TravelTime := na.locf(TravelTime,na.rm = F), by = patient_id]
dt4[, whostageinit := na.locf(whostageinit,na.rm = F), by = patient_id]
dt4[, marriedinit := na.locf(marriedinit,na.rm = F), by = patient_id]
dt4[, heightinit := na.locf(heightinit,na.rm = F), by = patient_id]

# DEFINE TIME VARYING COVARIATES 

# CONTINUOUS -- LOG WEIGHT, LOG CD4, LOG VL (MOST RECENT) -- SET LOG VL TO NA WHEN VL == 0
# INDICATOR OF VL==0
dt4[, weight_cf := na.locf(weight1,na.rm = F), by = patient_id]
dt4[, cd4_cf := na.locf(cd41,na.rm = F), by = patient_id]
dt4[, vl_cf := na.locf(viralload1,na.rm = F), by = patient_id]

################################################################
##### DO WE WANT AN INDICATOR OF VL==0 OR OF UNDETECTABLE VL?
##### Weldhagen (2008): VL<50 RNA COPIES/ML IS CONSIDERED AS UNDETECTALE IN AMPATH
##### SET CARRIED FORWARD LOG VL TO NA WHEN UNDETECTABLE?

# INDICATORS -- IS THERE UPDATE ON WEIGHT, CD4, VL
dt4[, wupd := max((sum(!is.na(weight1)) > 0)*1), by = tmpgridid] # 1 if there exists at least one update in the grid, 0 o.w.
dt4[, cupd := max((sum(!is.na(cd41)) > 0)*1), by = tmpgridid] 
dt4[, vupd := max((sum(!is.na(viralload1)) > 0)*1), by = tmpgridid] 


# DEFINE FUNCTIONAL VARIABLES
# INDICATOR OF TRT INITIATION, TIME OF TRT INITIATION, TIME OF CD4 MEASUREMENT IN THAT SAME GRID
# TREATMENT SWITCH
dt4[, onarv_cf := na.locf(onarv,na.rm = F), by = patient_id]
dt4[, trt_switch := onarv_cf - c(0,onarv_cf[-.N]), by = patient_id] 
#table(dt4$trt_switch,useNA = "always") # 1: start trt, 0: trt status unchanged, -1: switch off trt
dt4[, trt_init_time := last(time[trt_switch==1]), by = tmpgridid]

dt4[, last_cd4_time := last(time[which(!is.na(cd41))]), by = tmpgridid ]
#dt4[, first_cd4_time := time[which(!is.na(cd41))][1], by = tmpgridid ]

################################################################
# time of the last cd4 measurement update before treatment initiation in the same grid
dt4[, new_cd4_time := last(time[which(!is.na(cd41) & (time <= trt_init_time) )]), by = tmpgridid ]
# the exact value of that last cd4 update before treatment initiation 
dt4[, new_cd4 := cd41[time == new_cd4_time[1]], by = tmpgridid]

# CASES IN A GRID:
# 1. last cd4 in grid was <= trt init -- okay
# 2. last cd4 in grid was > trt init, no cd4 <= trt init in grid -> current grid cd4 NA; if next grid has no CD4, this last cd4 goes to the next grid
# 3. last cd4 in grid was > trt init, another cd4 in grid <= trt init -> that cd4 immediately before trt init is the current grid cd4; if next grid has no CD4, this last current grid cd4 after trt init goes to the next grid

# how to update cd4_cf to new_cd4_cf ?
dt4[,new_cd4_cf := cd4_cf]
### CURRENT GRID CD4
# what grid records need to change cd4_cf? --- (2)(3) --- last_cd4_time > trt_init_time 
# (2) when there is no cd4 <= trt init in grid --- new_cd4 and new_cd4_time are NA --- let all new_cd4_cf be the cd4_cf[trt_init_time] (NA)
dt4[last_cd4_time > trt_init_time & is.na(new_cd4_time), new_cd4_cf := cd4_cf[time == trt_init_time], by=tmpgridid]
# (3) when there is cd4 <= trt init in grid --- exist new_cd4 and new_cd4_time for the whole grid (not NA) --- last_cd4_time!=new_cd4_time, update new_cd4 to the whole grid
# existence of new_cd4  and last_cd4_time!= new_cd4_time implies last_cd4_time > trt_init_time
dt4[last_cd4_time!=new_cd4_time, new_cd4_cf := new_cd4, by=tmpgridid]

# how to update cupd?
# what grid records need to change cupd? --- last_cd4_time[previous grid] > trt_init_time[previous grid], and cupd[current grid] == 0
# --- update new_cupd[curent grid] (from 0) to 1 AND, new_cupd[previous grid] from 1 to 0 if there was no cd4 <= trt init in previous grid

dt4[, cd4_after_trt:= (last_cd4_time > trt_init_time)*1, by=tmpgridid]
dt4$cd4_after_trt[is.na(dt4$cd4_after_trt)] = 0
dt4[, find := cd4_after_trt - c(0,cd4_after_trt[-.N]), by=patient_id]
dt4[, find := 1*(find==-1)]
dt4$find[dt4$find==0] = NA #prepare for carry forward by tmpgridid
dt4[, upd_next_grid := na.locf(find, na.rm=F), by=tmpgridid] # the grid lines with last_cd4_time[previous grid] > trt_init_time[previous grid]
dt4[, find := NULL];dt4[, cd4_after_trt:=NULL]

dt4[, new_cupd := cupd]
#last_cd4_time[previous grid] > trt_init_time[previous grid], and cupd[current grid] == 0
#cd4_cf[current grid] is automatically the last cd4 from previous grid if cupd[current grid] == 0
#update new_cupd[curent grid] (from 0) to 1
dt4[cupd == 0 & upd_next_grid == 1, new_cupd := 1, by = patient_id]
#new_cupd[previous grid] from 1 to 0 if there was no cd4 <= trt init in previous grid
dt4[last_cd4_time > trt_init_time & is.na(new_cd4_time), new_cupd := 0, by=patient_id]


##############

# update new_weight_cf and new_wupd
dt4[, last_weight_time := last(time[which(!is.na(weight1))]), by = tmpgridid ]

dt4[, new_weight_time := last(time[which(!is.na(weight1) & (time <= trt_init_time) )]), by = tmpgridid ]
dt4[, new_weight := weight1[time == new_weight_time[1]], by = tmpgridid]

dt4[,new_weight_cf := weight_cf]
dt4[last_weight_time > trt_init_time & is.na(new_weight_time), new_weight_cf := weight_cf[time == trt_init_time], by=tmpgridid]
dt4[last_weight_time!=new_weight_time, new_weight_cf := new_weight, by=tmpgridid]

dt4[, weight_after_trt:= (last_weight_time > trt_init_time)*1, by=tmpgridid]
dt4$weight_after_trt[is.na(dt4$weight_after_trt)] = 0
dt4[, find := weight_after_trt - c(0,weight_after_trt[-.N]), by=patient_id]
dt4[, find := 1*(find==-1)]
dt4$find[dt4$find==0] = NA 
dt4[, upd_next_grid := na.locf(find, na.rm=F), by=tmpgridid] 
dt4[, find := NULL]
dt4[, weight_after_trt:=NULL]

dt4[, new_wupd := wupd]
dt4[wupd == 0 & upd_next_grid == 1, new_wupd := 1, by = patient_id]
dt4[last_weight_time > trt_init_time & is.na(new_weight_time), new_wupd := 0, by=patient_id]

# update new_vl_cf and new_vupd
dt4[, last_vl_time := last(time[which(!is.na(viralload1))]), by = tmpgridid ]

dt4[, new_vl_time := last(time[which(!is.na(viralload1) & (time <= trt_init_time) )]), by = tmpgridid ]
dt4[, new_vl := viralload1[time == new_vl_time[1]], by = tmpgridid]

dt4[,new_vl_cf := vl_cf]
dt4[last_vl_time > trt_init_time & is.na(new_vl_time), new_vl_cf := vl_cf[time == trt_init_time], by=tmpgridid]
dt4[last_vl_time!=new_vl_time, new_vl_cf := new_vl, by=tmpgridid]

dt4[, vl_after_trt:= (last_vl_time > trt_init_time)*1, by=tmpgridid]
dt4$vl_after_trt[is.na(dt4$vl_after_trt)] = 0
dt4[, find := vl_after_trt - c(0,vl_after_trt[-.N]), by=patient_id]
dt4[, find := 1*(find==-1)]
#dt4$find[dt4$find==0] = NA 
dt4[find == 0, find := NA]
dt4[, upd_next_grid := na.locf(find, na.rm=F), by=tmpgridid] 
dt4[, find := NULL]
dt4[, vl_after_trt:=NULL]

dt4[, new_vupd := vupd]
dt4[vupd == 0 & upd_next_grid == 1, new_vupd := 1, by = patient_id]
dt4[last_vl_time > trt_init_time & is.na(new_vl_time), new_vupd := 0, by=patient_id]



# TEST
#tmp1 = dt4[patient_id==73641]
#tmp1 = dt4[patient_id==73548]
#tmp1 = dt4[patient_id==73744]

save(dt4, file = "rawdt4.RData")

########################################################
# DELETE LINES OF IS.NA(GRID_IND)
dt5 = dt4[!is.na(grid_ind)]

length(unique(dt5$patient_id))

dim(dt5)#772868     66
########################################################
# COMPARE TRT_INIT_TIME VERSUS LAST_CD4_TIME
tmp = dt5[!is.na(trt_init_time) & !is.na(last_cd4_time),]
dim(tmp) #29704 grids have both trt init & cd4 measurement 
sum(tmp$trt_init_time < tmp$last_cd4_time) # 10654 grids have cd4 measured AFTER trt init
mean(tmp$trt_init_time < tmp$last_cd4_time) # 35.9%

sum((tmp$trt_init_time+30) < tmp$last_cd4_time) # 9964 grids have cd4 measured AFTER one month post trt init

png(filename = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\Draft7\\Rplot.png", width =500, height = 300)
hist(tmp$last_cd4_time - tmp$trt_init_time, main="Last CD4 Time minus Trt Init Time", xlab = "Days")
dev.off()

sum(tmp$trt_init_time == tmp$last_cd4_time)#9257
mean(tmp$trt_init_time == tmp$last_cd4_time)#31.16%

sum(!is.na(dt5$trt_init_time)) # 72680 trt init grids in total
########################################################

dt6 = dt5[,c("patient_id","time","apptdate","onarv_cf","trt_switch",
             "ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit",
             "new_cd4_cf","new_cupd","new_weight_cf","new_wupd","new_vl_cf","new_vupd",
             "totrec","totenc","engage","transfer_outcome","death_outcome")]


#define OUTCOME: 0 disengage, 1 engage, 2 transfer_outcome, 3 death_outcome
dt6$outcome = dt6$engage
dt6$outcome[dt6$transfer_outcome==1] = 2
dt6$outcome[dt6$death_outcome==1] = 3
# some xfers might lost due to its pushed-back last grid and data closure

dt6$log_weight = log(dt6$new_weight_cf)
dt6$log_cd4plus1 = log(dt6$new_cd4_cf + 1)
dt6$vl_0 = (dt6$new_vl_cf<=50 & !is.na(dt6$new_vl_cf))*1 # indicator of vl exists and <=50
dt6$new_vl_cf[dt6$vl_0==1] = 0 # set undetectable but existing vl to 0
dt6$log_vlplus1 = log(dt6$new_vl_cf + 1)

dt6$log_weight_RX = dt6$log_weight
dt6$miss_weight = (is.na(dt6$log_weight))*1
dt6$log_weight_RX[dt6$miss_weight==1] = 0

dt6$log_cd4plus1_RX = dt6$log_cd4plus1
dt6$miss_cd4 = (is.na(dt6$log_cd4plus1))*1
dt6$log_cd4plus1_RX[dt6$miss_cd4==1] = 0

dt6$log_vlplus1_RX = dt6$log_vlplus1
dt6$miss_vl = (is.na(dt6$log_vlplus1))*1
dt6$log_vlplus1_RX[dt6$miss_vl==1] = 0

#A series of "ad-hoc" methods have been suggested to address missing data. One approach for continuous variables involves recoding missing values to some common value, creation of an indicator of missingness as a new variable, and including both these variables along with their interaction in the regression model. A similar approach for categorical variables involves the creation of an additional category for missing values. These ad-hoc approaches have the potential to induce bias and are not recommended (Jones 1996, Greenland & Finkle 1995).

# COVARIATES: log_weight, log_cd4plus1, vl_0, log_vlplus1, miss_weight, miss_cd4, miss_vl

# save dt6
wd = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\makeGridData"
fwrite(dt6, paste0(wd,"\\my_analysis_data.csv"))

# read data
#dt = fread("my_analysis_data1.csv")

###############################################################################




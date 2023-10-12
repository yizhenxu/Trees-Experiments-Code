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

summIQR = function(x){
  temp = tryCatch(all(x==round(x),na.rm = TRUE),  error = function(e) print(2)) #FALSE(0) if float, TRUE(1) if integer, 2 if categorical
  
  # float
  if(temp==0){
    res = matrix(NA,nrow = 1, ncol = 6)
    res[1] = ""
    res[2] = sum(!is.na(x))
    res[3] = round(mean(!is.na(x))*100,1)
    res[4:6] = round(as.numeric(quantile(x, probs = c(0.25,0.5,0.75), na.rm=T)),1)
    res = as.data.frame(res)
  }   
  # categorical or integer
  if(temp %in% 1:2){
    x[x == ""] = NA
    tt = table(x,useNA="always")
    if(length(tt)<=20){
      res = matrix(NA,ncol = 6, nrow = length(tt))
      res[,1] = names(tt)
      res[,2] = as.numeric(tt)
      res[,3] = round(as.numeric(res[,2])/length(x)*100,1)
      res = as.data.frame(res)
    } else {
      res = matrix(NA,nrow = 1, ncol = 6)
      res[1] = "Levels"
      res[2] = length(tt)
      res = as.data.frame(res)
    }
    
  }
  return(res)
}

# make summary table 
# data at enrollment
TabCovSumm = function(varlist, d){
  for(i in 1:length(varlist)){
    if(varlist[i] == "gender"){
      x = factor(d$male,labels = c("Female","Male"))
    }
    if(varlist[i] %in% c("weight","height","cd4")){
      x = d[[which(colnames(d)==varlist[i])]]
      x[x == -99] = NA
    }
    if(varlist[i] %in% c("ageatenroll","married","yrenrol","initclinic","TravelTime","state_province",
                         "whostage","onarv","viralload")){
      x = d[[which(colnames(d)==varlist[i])]]
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
  
  tab[,2] = as.character(tab[,2])
  tab[which(is.na(tab[,2])), 2] = "Missing"
  tab = tab[-which(tab[,4] == 0),]
  return(tab)
}


TabCovSumm1 = function(varlist, d){
  for(i in 1:length(varlist)){
    if(varlist[i] == "gender"){
      x = factor(d$male,labels = c("Female","Male"))
    } else {
      x = d[[which(colnames(d)==varlist[i])]]
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
  
  tab[,2] = as.character(tab[,2])
  tab[which(is.na(tab[,2])), 2] = "Missing"
  tab = tab[-which(tab[,4] == 0),]
  return(tab)
}


TabCovSumm2 = function(varlist, d){
  for(i in 1:length(varlist)){
    if(varlist[i] == "gender"){
      x = factor(d$male,labels = c("Female","Male"))
    } else {
      x = d[[which(colnames(d)==varlist[i])]]
    }
    
    if(i == 1){
      tab = summIQR(x)
      tab = cbind(c(varlist[i],rep("",nrow(tab)-1)), tab) # add the var name col
      colnames(tab) = c("","","n","percent","25%","50%","75%")
    } else {
      res = summIQR(x)
      res = cbind(c(varlist[i],rep("",nrow(res)-1)), res) # add the var name col
      colnames(res) = colnames(tab)
      tab = rbind.data.frame(tab,res)
    }
  }
  
  tab[,2] = as.character(tab[,2])
  tab[which(is.na(tab[,2])), 2] = "Missing"
  tab = tab[-which(tab[,4] == 0),]
  return(tab)
}

### 3. clean state province record with text parsing

text_state = function(a){
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
  
  return(a)
}


movefun = function(data, var, gap, anchor){
  
  if(gap[1]>0 || gap[2]<0) return("wrong")
  
  tmp = copy(data)
  
  tmp[, tag1 := 1*(time >= (anchor+gap[1]) & time <= anchor+gap[2])]
  tmp = tmp[tag1 == 1,]
  tmp$mt = tmp$time; tmp$mt[is.na(tmp[,get(var)])] = NA 
  if(gap[1] < 0){
    # var1 : carried forward values of var
    tmp[, var1 := na.locf(get(var),na.rm = F), by = patient_id]
    # mt1 : carried forward values of mt; tells which date this carried forward record was measured
    tmp[, mt1 := na.locf(mt,na.rm=F), by = patient_id] 
  }
  if(gap[2] > 0){
    # var2 : carried backward values of var
    tmp[, var2 := na.locf(get(var),na.rm = F,fromLast = TRUE), by = patient_id]
    # mt2 : carried backward values of mt
    tmp[, mt2 := na.locf(mt,na.rm=F, fromLast = TRUE), by = patient_id] 
  }
  
  if(gap[1] == 0 & gap[2] > 0){ #[0, gap[2]]
    tmp[time == anchor & is.na(mt) & !is.na(mt2), (var) := var2]
  } else if(gap[1] < 0 & gap[2] == 0){ #[gap[1], 0]
    tmp[time == anchor & is.na(mt) & !is.na(mt1), (var) := var1]
  } else {
    tmp[time == anchor & is.na(mt) & (!is.na(mt1)+!is.na(mt2))==2, (var) := ifelse(time - mt1 <= mt2 - time, var1, var2)]
    tmp[time == anchor & is.na(mt) & (!is.na(mt1)+!is.na(mt2))==1, (var) := ifelse(!is.na(mt1), var1, var2)]
  }
  
  print(paste0("Non-Missing Proportion of ",var,": ",sum(!is.na(data[time == anchor, get(var)]))," -> ", sum(!is.na(tmp[time == anchor,get(var)]))))
  #return(tmp[time == anchor, c("patient_id", var), with=F])
  return(tmp[time == anchor, get(var)])
  
}


if(0){# carry forward backward code
  
  cfdata = function(data, tmp_before, colnum, gap){
    var = colnames(data)[colnum]
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
    data[(patient_id %in% tmp_before1$patient_id) & new_num==0, (var):= tmp_before1$var1]
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
    data[(patient_id %in% tmp_after1$patient_id) & new_num==0, (var):= tmp_after1$var1]
    return(data[,get(var)])  
  }
  
  
  data = copy(dt1)
  varlist = c("height" ,"weight" , "cd4" ,"viralload","whostage", "onarv", "married","hivstatus" )
  cfcbcols = which(colnames(data) %in% varlist)
  varlist = colnames(data)[cfcbcols] # reorder varlist as it is in data
  
  # prepare for carry forward
  tmp_before = data[new_num <= 0 ]
  tmp_before[, prev_num := c(NA,new_num[-.N]), by = patient_id ]
  tmp_before = tmp_before[-which(tmp_before$new_num==0 & is.na(tmp_before$prev_num)),] # exclude those with no negative new_num
  
  # carry forward all before adult initial within 31 days: V, A, X
  cfdt = lapply(cfcbcols, function(k) cfdata(data, tmp_before, k, gap=31))
  data[, (varlist) := cfdt]
  
  # only carry backward V to adult initial within 31 days
  varlistb = c("height" ,"weight" ,"viralload","whostage", "married","hivstatus" )
  cfcbcolsb = which(colnames(data) %in% varlistb)
  varlistb = colnames(data)[cfcbcolsb] # reorder varlist as it is in data
  
  # prepare for carry backward
  tmp_after = data[new_num >= 0 ]
  tmp_after[, post_num := c(new_num[-1],NA), by = patient_id ]
  tmp_after = tmp_after[-which(tmp_after$new_num==0 & is.na(tmp_after$post_num)),] # exclude from target those with no records after adult init record
  
  
  # carry forward all before adult initial within 31 days: V, A, X
  cbdt = lapply(cfcbcolsb, function(k) cbdata(data, tmp_after, k, gap=31))
  data[, (varlistb) := cbdt]
  
  ##############################################
  
  ### 7. carry forward non-adult-return measurements of TIME-VARYING CONFOUNDERS to the nearest next adult return
  
  ## carry backward patient_id x adult return apptdate as temporary ID
  dt1[, tmp := apptdate]
  dt1[! encounter_type %in% c("adult initial","adult return"), tmp := NA]
  dt1[, tmp1 := na.locf(tmp,na.rm = F,fromLast = TRUE), by = patient_id]
  dt1[, tmp := paste0(patient_id,"_",tmp1)]; dt1[, tmp1:=NULL]
  
  ## carry forward all TIME-VARYING CONFOUNDERS forward by the temporary ID
  #tvar = c("cd4", "weight","height","viralload") # only "cd4" is sufficient, or just use tvar = varlist
  tvar = c("cd4")
  for(i in 1:length(tvar)){
    var = tvar[i]
    dt1[, (var) :=  na.locf(get(var),na.rm = F), by = tmp]
  }
  dt1[,tmp:=NULL]
  
  ##############################################
  
  ### DECRIPTION 2. Destribution of treatment initiation time
  png(filename = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft6\\DTrtInitTime.png", width =800, height = 300)
  
  hist(a[a<=100 & a >= (0)], main = "Treatment Initiation Time", breaks = 200, xlab="Day (Initial Encounter at Day 0)",col="red")
  abline(v = c(7,30), col = "blue")
  dev.off()
  
  range(a, na.rm=T)
  mean(a < 0, na.rm=T); sum(a < 0, na.rm=T); sum(a <= (-200), na.rm=T); mean(a==0, na.rm=T)
  
  ### DECRIPTION 3. Destribution of CD4 measurement time in [0, 200]
  
  png(filename = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft6\\DCD4Time.png", width =500, height = 300)
  
  hist(b[b<=100 & b > (0)],add=T,col="black",breaks = 200)
  hist(b[b<=100 & b > (0)], main = "CD4 Measurement Time", breaks = 200, xlab="Day (Initial Encounter at Day 0)",col="black")
  abline(v = c(-30,-7,7,30), col = "red")
  dev.off()
}

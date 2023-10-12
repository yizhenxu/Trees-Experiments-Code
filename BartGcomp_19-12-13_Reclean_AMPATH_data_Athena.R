
library(readstata13)
library(lubridate)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(data.table)
library(zoo)

#dt2 = dt3[,c("patient_id","onarv","arvstart","apptdate","returnvisitdate",
#             "death", "dod", "transfer",
#             "ageinit","male","yrenrol","TravelTime","whostageinit","marriedinit","heightinit",
#             "weightinit", "cd4", "viralload","new_num")]

### Lay grids of 200 days

### LAY GRID
int_days = 200

# DEFINE DATA CLOSURE DATE
dt2[,closure := max(apptdate)] # 2016-8-24

### 8. Include only people with baseline at least 200 days before the closure date
length_in_care = dt2[, as.numeric(closure - apptdate[1L]), by=patient_id]
include = 1*(length_in_care$V1<200)
dt2 = dt2[include == 0,]

dim(dt2)# 1510686
length(unique(dt2$patient_id))#  76748

### 9. Define grid data and merge

# HOW MANY GRID POINT FOR EACH PERSON?
length_in_care = dt2[, as.numeric(closure[1L] - apptdate[1L]), by=patient_id]
num_grid = floor(length_in_care$V1/int_days) # range(num_grid): 1-15

# DEFINE GRID DATA
grid_id = rep(length_in_care$patient_id, num_grid+2) # plus 2 for the zeros (baseline) and 12 (care initiation period end)
grid_times = unlist(lapply(num_grid,function(x) c(0,12,(1:x)*int_days)))

# day 12 as grid ind 0.5
grid_ind = unlist(lapply(num_grid,function(x) c(0,0.5,1:x))) 
grid_data = as.data.frame(cbind(patient_id = grid_id, time = 
                                  grid_times, grid_ind = grid_ind))
grid_data = setDT(grid_data)
remove(grid_id,grid_times,grid_ind,num_grid,length_in_care)

# SET KEY AND MERGE 
dt2[,tmpid:=paste0(patient_id,"_",time)]
grid_data[,tmpid:=paste0(patient_id,"_",time)]
setkey(dt2,tmpid)
setkey(grid_data, tmpid)

#memory.limit(size= 10000)
dt2 = merge(dt2, grid_data, all = TRUE)

dim(dt2)# 2235500      55

remove(grid_data);gc()

# CLEAN SOME VARIABLES FROM MERGING, AND SORT
dt2[,patient_id  := patient_id.x]
dt2$patient_id[is.na(dt2$patient_id)] = dt2$patient_id.y[is.na(dt2$patient_id)]
#dt2[is.na(patient_id), patient_id := patient_id.y]

dt2[,time := time.x]
dt2[is.na(time), time := time.y]

dt2 = dt2[order(patient_id,time)]

dt2[, tmpid:= NULL]
dt2[, time.x:= NULL]
dt2[, time.y:= NULL]
dt2[, patient_id.x:= NULL]
dt2[, patient_id.y:= NULL]

rm(list = ls.str(mode = 'numeric'))

#tmp2 = dt2[patient_id == 277813]
#print(tmp2)

###########################################################

# DEFINE TMPID = PATIENT_ID X GRID_ID 
dt2[, grid_ind_cb := na.locf(grid_ind,na.rm=F,fromLast = TRUE), by = patient_id] # grid intervals (,]
dt2[, tmpgridid := paste0(patient_id,"_",grid_ind_cb)]
dt2[, totrec := sum(!is.na(apptdate)), by = tmpgridid]
dt2[, totenc := sum(record_type,na.rm=T)+ sum(transfer==2,na.rm=T), by = tmpgridid]
#dt2$engage = 1*(dt2$totenc > 0)

### Mark all records in a grid where dod falls in, indicator lgrid_death

# There is already no record after death --> only joined grid after time of the death record
dt2[, timetodeath := as.numeric(dod[1L] - apptdate[1L]),by = patient_id]# time to death from baseline
dt2[, lgrid_death:=ifelse(timetodeath!=0, 1*(!is.na(timetodeath) & (time/timetodeath)>=1), 1) ] #timetodeath at 0 are those died at the first clinic encounter (baseline)
# lgrid_death is all NA's for non-death individuals 
# otherwise, 0 for lines with time < timetodeath, 1 for all lines at or after the death record time
dt2[is.na(lgrid_death), lgrid_death := 0]
dt2[, lgrid_death := lgrid_death - c(0,lgrid_death[-.N]),by = patient_id ]
dt2[, lgrid_death := max(lgrid_death,na.rm = T), by = tmpgridid]
# lgrid_death is 1 for only lines in the grid where dod is, otherwise 0
dt2[, .N, by = lgrid_death]

#a = dt2[patient_id==73500,]
#a[17:21,]

dt2[, timetodeath := NULL]

### Mark the grid record after the grid of xfer==2, as having lgrid_xfer == 1 

# grid t containing record with transfer == 2 is NOT xfer out, grid t+1 is xfer out (if grid t+1 happens before closure) 
# lgrid_xfer is indicator of last grid of xfer, definitely an inserted grid, not in the obs data
dt2[, xfer2 := 1*(transfer == 2)]
dt2[, lgrid_xfer_tmp := max(c(0, xfer2),na.rm = T), by = tmpgridid] 
# lgrid_xfer_tmp = indicator of records in the grid of transfer==2
dt2[, lgrid_xfer_tmp := lgrid_xfer_tmp - c(0,lgrid_xfer_tmp[-.N]), by=patient_id]
dt2[, lgrid_xfer := 1*(lgrid_xfer_tmp == (-1))] # some might be censored from data closure
dt2[, lgrid_xfer := max(lgrid_xfer,na.rm = T), by = tmpgridid]
dt2[, .N, by = lgrid_xfer]

#a = dt2[patient_id==73598,]
# a[,c("time","transfer","lgrid_xfer"),with=F]

dt2[, lgrid_xfer_tmp :=NULL]
dt2[, xfer2 := NULL]

### 10. Define outcome: 0 disengaged, 1 engaged, 2 xfer out, 3 death
dt2[, outcome := 1*(totenc >= 1), by = tmpgridid]
dt2[lgrid_xfer == 1, outcome :=2]
dt2[lgrid_death == 1, outcome:=3]

### IMPORTANT! the grid where xfer happened is engaged!
dt2[transfer==2,.N,outcome] # should be all ones!

### Mark the last grid based on lenc

dt2$lenc[is.na(dt2$lenc)] = 0
dt2[, lgrid := max(lenc,na.rm = T), by = tmpgridid]





save(dt2, file = "AMPATH_12-13-2019_dt2.RData")
load("AMPATH_12-13-2019_dt2.RData")

### 19. exclude lines after the last grid, lgrid_death == 1, or lgrid_xfer == 1
dim(dt2)# 2235500

dt2 = dt2[-which(is.na(grid_ind_cb)),]

dt2[lgrid_death==0, lgrid_death := NA]
dt2[, lgrid_death1 := na.locf(lgrid_death,na.rm = F), by = patient_id]
dt2[, ind := 1*(is.na(lgrid_death) & lgrid_death1 == 1)]
# ind: indicator of not in the dod grid but are lines after the dod grid

dt2[lgrid_xfer==0, lgrid_xfer := NA]
dt2[, lgrid_xfer1 := na.locf(lgrid_xfer,na.rm = F), by = patient_id]
dt2[is.na(lgrid_xfer) & lgrid_xfer1 == 1, ind := 1]

dt2 = dt2[ind == 0 | is.na(ind),]
dim(dt2)  # 2138641
  
dt2[, c("lgrid_death","lgrid_death1","lgrid_xfer","lgrid_xfer1","ind"):=NULL]

###########################################################
### CARRY FORWARD FEATURES -> EXTRACT GRIDS -> DEFINE OUTCOME -> DEFINE COVARIATES FOR ANALYSIS

# features to carry forward: ageinit, male, yrenrol, TravelTime, whostageinit, marriedinit, heightinit   
dt2[, ageinit := na.locf(ageinit,na.rm = F), by = patient_id]
dt2[, male := na.locf(male,na.rm = F), by = patient_id]
dt2[, yrenrol := na.locf(yrenrol,na.rm = F), by = patient_id]
dt2[, TravelTime := na.locf(TravelTime,na.rm = F), by = patient_id]
dt2[, whostageinit := na.locf(whostageinit,na.rm = F), by = patient_id]
dt2[, marriedinit := na.locf(marriedinit,na.rm = F), by = patient_id]
dt2[, heightinit := na.locf(heightinit,na.rm = F), by = patient_id]
dt2[, weightinit := na.locf(weightinit,na.rm = F), by = patient_id]
dt2[, viralloadinit := na.locf(viralloadinit,na.rm = F), by = patient_id]


# DEFINE TIME VARYING COVARIATES 

# CONTINUOUS -- LOG WEIGHT, LOG CD4, LOG VL (MOST RECENT) -- SET LOG VL TO NA WHEN VL == 0
# INDICATOR OF VL==0
colnames(dt2)
dt2[, original_cd4 := cd4]
dt2[time == 0, cd4:=cd4_preB]
dt2[, cd4_cf := na.locf(cd4,na.rm = F), by = patient_id]

# INDICATORS -- 
# IS THERE UPDATE ON CD4
dt2[, cupd := max((sum(!is.na(cd4)) > 0)*1), by = tmpgridid] 
# ARV
dt2[, onarv_cf := na.locf(onarv,na.rm = F), by = patient_id]

########################################################
########################################################
dt3 = dt2[!is.na(grid_ind)]

length(unique(dt3$patient_id))#76748

dim(dt3)# 755579     62

### Align onarv and cd4_cf on grid_ind 0 and 0.5
# CD4

dt3[, new_cd4_cf := cd4_cf]
dt3[, tmp := cd4_cf[2L], by = patient_id]
dt3[time==0, new_cd4_cf := tmp, by= patient_id]

dt3[, new_cupd := cupd]
dt3[, tmp := cupd[2L], by = patient_id]
dt3[time==0, new_cupd := tmp, by= patient_id]
dt3[, tmp:=NULL]

dt3[time==0, onarv_cf := onarv_0, by= patient_id]

dt4 = dt3[grid_ind != 0.5,]

dim(dt4)#678831     64
length(unique(dt4$patient_id))#76748

### Delete patients who transferred out during (0, 12]
a = dt4[time==200 & outcome == 2,]
dt4 = dt4[!(patient_id %in% a$patient_id),]

dim(dt4)#678815     64
length(unique(dt4$patient_id))#76740

########################################################
########################################################

### Process variables

##### DO WE WANT AN INDICATOR OF VL==0 OR OF UNDETECTABLE VL?
##### Weldhagen (2008): VL<50 RNA COPIES/ML IS CONSIDERED AS UNDETECTALE IN AMPATH
##### SET CARRIED FORWARD LOG VL TO NA WHEN UNDETECTABLE?
dt4[, vl_0 := 1*(!is.na(viralloadinit) & viralloadinit<=50) ] # indicator of vl exists and <=50
dt4[vl_0 == 1, viralloadinit := 0] # set undetectable but existing vl to 0
dt4[, log_vlplus1_init := log(viralloadinit +1)]

dt4[, log_cd4plus1 := log(new_cd4_cf + 1)]

dt4[, cd4_preB := na.locf(cd4_preB,na.rm = F), by = patient_id]
dt4[, log_cd4preBplus1 := log(cd4_preB+1)]

dt4[, log_weight_init := log(weightinit),by = patient_id]# not neccessary

dt4[, log_height_init := log(heightinit),by = patient_id]# not neccessary

### RELEVEL THE OUTCOME TO 1 (DISENG), 2 (ENG), 3 (XFER), 4 (DEATH)
### outcome: 0 disengage, 1 engage, 2 transfer_outcome, 3 death_outcome
dt4$outcome = dt4$outcome+1
dt4[,.N,outcome]

#   outcome      N
#1:       2 419662
#2:       1 250899
#3:       4   5450
#4:       3   2804


### Relevel TravelTime to numeric
a = dt4$TravelTime
a[a == "< 30 minutes"] = 1
a[a == "30-60 minutes"] = 2
a[a == "1-2 hours"] = 3
a[a == "> 2 hours"] = 4
dt4$TravelTime = as.numeric(a)


save(dt4, file = "AMPATH_12-13-2019_dt4.RData")
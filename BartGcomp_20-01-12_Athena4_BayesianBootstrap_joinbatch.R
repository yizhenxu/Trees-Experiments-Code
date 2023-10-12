# round 1: set rdnum = 1, bash array 1:120
# round 2: set rdnum = 2, bash array 1:3

rdnum = 2

taskid = as.numeric(Sys.getenv('SGE_TASK_ID')) # 1 to 36 for batch size 30 rd1, 0 to 2 for rd 2

setting = 4; maxt = 10; 
source("BartGcomp_20-01-12_Athena0_prepare.R")
skip = 4; ndthin = floor(nd/(skip+1))

###

L = list.files(getwd())

book = cbind(expand.grid(c(0,1,2), 1:ndthin), 0)
for(i in 1:nrow(book)){
  A = book[i,1]; k = book[i,2]
  if(file.exists( paste0("Gcomp_A",A,"_k",k,"_22-05-04.RData") )) book[i,3] = 1
}
npost_here = min(max(book[which(book[,3]==1 & book[,1] ==0) ,2]),
                 max(book[which(book[,3]==1 & book[,1] ==1),2]),
                 max(book[which(book[,3]==1 & book[,1] ==2),2]))

###

bs = 25 # combine batch size
rd = floor(npost_here/bs)
rbs = npost_here - bs*rd
totalrd = ifelse(npost_here/bs == rd, rd, rd+1)
idxlist = vector("list", totalrd)
for(l in 1:rd){
  idxlist[[l]] = 1:bs + (l-1)*bs
}
if(totalrd > rd) idxlist[[totalrd]] = 1:rbs + rd*bs

totalnum = rd # should be totalrd

if(rdnum == 1){ # taskid 1:3*totalrd, 36 now
  newbook = expand.grid(c(0,1,2), 1:totalnum)
  
  if(taskid <= nrow(newbook)){
    
    A = newbook[taskid, 1]
    rdid = newbook[taskid, 2]
    
    fn <- paste0("tmp_A",A,"_rd1_",rdid,".RData")
    if (!file.exists(fn)){
      Smatall = vector("list", maxt) 
      Xmatall = vector("list", maxt)
      
      for(t in 1:maxt){
        Xmatall[[t]] = array(NA, dim = c(ngc, ndthin, length(theX)))
        Smatall[[t]] = matrix(NA, nrow = ngc, ncol = ndthin)
      }
      for(k in idxlist[[rdid]]){
        load(paste0("Gcomp_A",A,"_k",k,"_22-05-04.RData") )
        for(t in 1:maxt){
          Xmatall[[t]][,k,] = Xmat[[t]][,k,]
          Smatall[[t]][,k] = Smat[[t]][,k]                      
        }
        rm(Xmat,Smat)
      }
      
      save(Xmatall,Smatall,file = fn)
    } 
    
  }
}

if(rdnum == 2){ # taskid 0:2
  if(taskid %in% 1:3){
    A =  taskid - 1
    
    Smat = vector("list", maxt) 
    Xmat = vector("list", maxt)
    
    for(t in 1:maxt){
      Xmat[[t]] = array(NA, dim = c(ngc, ndthin, length(theX)))
      Smat[[t]] = matrix(NA, nrow = ngc, ncol = ndthin)
    }
    for(rdid in 1:totalnum){
      load(paste0("tmp_A",A,"_rd1_",rdid,".RData"))
      idx = idxlist[[rdid]]
      for(t in 1:maxt){
        Xmat[[t]][,idx,] = Xmatall[[t]][,idx,]
        Smat[[t]][,idx] = Smatall[[t]][,idx]
      }
      rm(Smatall, Xmatall)
    }
    
    save(Xmat,Smat,file = paste0("res_A",A,".RData"))
  }
  
}


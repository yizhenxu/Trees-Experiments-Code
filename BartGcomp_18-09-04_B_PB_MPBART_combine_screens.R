setwd("./Topic2")
screens = c(132394, 133012, 174198)
screens = c(133134, 133152, 174216)
ct = c(10, 19, 31)

load(paste0("CondSModel_screen",screens[1],"_ct",ct[1],".RData"))
mod1 = mod
modfml1 = modfml
Smod1 = Smod

load(paste0("CondSModel_screen",screens[2],"_ct",ct[2],".RData"))
for(i in 10:18){
  mod1[i] = mod[i]
  modfml1[i] = modfml[i]
}
for(i in 5:10){
  Smod1[i] = Smod[i]
}

load(paste0("CondSModel_screen",screens[3],"_ct",ct[3],".RData"))
for(i in 19:(3*maxt)){
  mod1[i] = mod[i]
  modfml1[i] = modfml[i]
}
for(i in 11:(2*(maxt-1))){
  Smod1[i] = Smod[i]
}

mod = mod1
modfml = modfml1
Smod = Smod1


save(mod,Smod,modfml,file = paste0("CondSModel_screen",screens[1],"_ct",30,"ntr",ntr,".RData"))




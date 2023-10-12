
source("BartGcomp_19-03-02_Athena0_prepare.R")

load("ContBin.RData")

maxt = 10
mod3_D =  vector("list", maxt-1) # categorical models for those previously disengaged
mod3_E =  vector("list", maxt) # categorical models for those previously engaged
baseylist = vector("list", maxt)
ntlist = vector("list", maxt)
seedlist = vector("list", maxt)

resCat = array(NA, dim = c(4,maxt,4)) # Accuracy - Categorical Models 
# the first dim: y, mean post, 5% post, 95% post
# the second dim: time from 1 to maxt
# the third dim: eng, diseng, xfer, death

for(t in 1:maxt){
  load(paste0("Cat19-03-05_",t,".RData"))
  
  modfml[[3*t]] = fml
  
  resCat[, t, ] = res
  baseylist[[t]] = baseyl
  ntlist[[t]] = ntL
  seedlist[[t]] = seedvec
  
  if(t==1){
    mod3_E[[1]] = Bmod2
  } else {
    mod3_D[[t-1]] = Bmod2[[1]]
    mod3_E[[t]] = Bmod2[[2]]
    remove(baseyl)
  }
  
  remove(Bmod2)
  remove(fml)
  remove(res)
  remove(ntL)
  remove(seedvec)
}

print(seedlist)
print(baseylist)
resCat

save(mod1, mod2, mod3_D, mod3_E, baseylist, ntlist, seedlist, modfml, file = "Mod19-03-05.RData")
save(resCat, file = "resCat.RData")








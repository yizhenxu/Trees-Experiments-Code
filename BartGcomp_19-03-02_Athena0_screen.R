### Categorical models
### bash

#!/bin/bash
#SBATCH -J Cat1
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --time 72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu


for t in `seq 1 1 10`; do

R CMD BATCH --vanilla "--args $t"  BartGcomp_19-03-02_Athena2_Mod_Cat.R  ./cat${t}.out

done

##################################################################################

### combine all models into "Mod19-03-4.RData"
### mympbartmod3 for fitting categorical models in Athena
### continuous and binary models using the same ones from 18-10-27

#!/bin/bash
#SBATCH -J Eval
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time 72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

R CMD BATCH --vanilla BartGcomp_19-03-02_Athena3_JoinAndEval.R eval.out


### plot accuracy pd for categorical models
load("resCat.RData")
maxt = 10
res = resCat
pd = as.data.frame(matrix(NA, nrow = 4*maxt, ncol = 6))

colnames(pd) = c("Outcome", "Time", "obs", "Posteriors", "p5", "p95")
pd[,1] = rep(c("Engaged","Disengaged","Transferred","Died"), each = maxt)
pd[,2] = rep(1:maxt, 4)
pd[,3:6] = rbind(t(res[,,1]), t(res[,,2]), t(res[,,3]), t(res[,,4]) ) 
pd$Outcome = factor(pd$Outcome, levels = c("Engaged","Disengaged","Transferred","Died"))
pd = pd[-which(pd[,1]=="Transferred" & pd[,2]==1),]

library("ggplot2")

p=ggplot(pd, aes(Time)) + 
  geom_line(aes(y=Posteriors),colour="red") + 
  geom_point(aes(y=obs), colour="blue") +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.5)+  facet_wrap(~ Outcome, scales="free_y")
plot(p)

path = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft3\\Accuracy Comparison - Diag Plots\\19-03-05\\"

png(paste0(path,"AccuracyFix.png"),width = 600, height = 400)
plot(p)
dev.off()

if(0){
  library(GcompBART)
  load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART/ZEUS RData/Cat19-03-02_8.RData")
  png(paste0(path,"Diag_8_D.png"),width = 600, height = 400)
  DiagPlot(Bmod2[[1]],1)
  dev.off()
  png(paste0(path,"Diag_8_E.png"),width = 600, height = 600)
  DiagPlot(Bmod2[[2]],1)
  dev.off()
}
##################################################################################

### CatFix.sh

#!/bin/bash
#SBATCH -J FIX
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time 72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

for seed in `seq 0 10 20`; do

R CMD BATCH --vanilla "--args $seed"  BartGcomp_19-03-02_Athena2FIX8_Mod_Cat.R  ./FIX${seed}.out

done

##################################################################################

### Gcomp.sh

#!/bin/bash
#SBATCH -J Gcomp
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time 72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

for A in `seq 0 1 2`; do

R CMD BATCH --vanilla "--args $A"  BartGcomp_19-03-02_Athena4_Gcomp.R  ./Gcomp${A}.out

done

### plot from Gcomp counts

A = #0,1,2; 2 for 1{CD4<=350}
path = "G:\\My Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"
pathto = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\Draft\\paper2Draft3\\Accuracy Comparison - Diag Plots\\19-03-05\\"
path = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\CODE\\G COMPUTATION BART\\ZEUS RData"
load(paste0(path,"\\Counts_A",A,".RData"))

png(paste0(pathto,"CondSEngA",A,".png"),width = 600, height = 400)
barplot(counts,
        xlab="Days since baseline", ylab="Number of observations", col=1:4)
legend("bottomleft", 
       legend = rownames(counts), 
       fill = 1:4)
dev.off()


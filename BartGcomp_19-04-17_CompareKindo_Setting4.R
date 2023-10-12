### Contraceptive Method Choice

#datpath = "C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\HIV Treatment Cascade\\MPBART Computing\\UCI Contraceptive Choice"
#d = read.table(paste0(datpath,"\\cmc.data"),header = FALSE, sep=",")
d = read.table("cmc.data",header = FALSE, sep=",")

colnames(d) = c("age","edu","husbedu","children","islam","work","husboccp","liveind","media","y")
head(d)

table(d$Y)
dim(d)

yind = which(colnames(d)=="Y")
colnames(d)[yind] = "y"

set.seed(1)
rnd = sample(1:nrow(d)); ntrd = ceiling(nrow(d)*0.7)
trd = d[rnd[1:ntrd],]; ted = d[rnd[(ntrd+1):nrow(d)],]

fml = "y~age+edu+husbedu+children+islam+work+husboccp+liveind+media"



### Multinomial Probit BART Settings

set.seed(99)
nd = 2000 # posterior draws
nb = 500 # burn in
nt = 100 # number of trees
p=length(unique(trd$y)) # number of outcome levels

# seed to start with
CatSeed = 888888



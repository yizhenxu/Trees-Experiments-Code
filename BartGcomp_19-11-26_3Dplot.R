#install.packages("plotly")
#install.packages("TruncatedNormal")
library("TruncatedNormal")
library("plotly")
packageVersion('plotly')



### ### mu1 = 0 (plot in paper)

mu1 = 0 #mu2 = mu1 + x
y =  seq(-0.4,0.9,0.1)#rho
x = seq(-0.5,0.5,0.1)#mu2-mu1
z1 = matrix(NA, nrow = length(y), ncol = length(x))
z2 = matrix(NA, nrow = length(y), ncol = length(x))
z3 = matrix(NA, nrow = length(y), ncol = length(x))
for( i in 1:length(x)){
  for( j in 1:length(y)){
    z1[j,i] = mvNcdf(l = c(-mu1,x[i]), u = c(Inf,Inf), Sig = matrix(c(1,1-y[j],1-y[j],2),2))$prob
    z2[j,i] = mvNcdf(l = c(-(mu1+x[i]),-x[i]), u = c(Inf,Inf), Sig = matrix(c(1,1-y[j],1-y[j],2),2))$prob
    z3[j,i] = mvNcdf(l = c(-Inf,-Inf), u = c(-mu1, -(mu1+x[i])), Sig = matrix(c(1,y[j],y[j],1),2))$prob
  }
}


a = seq(min(x)+mu1, max(x)+mu1,length.out =11)
axx <- list(
  ticketmode = 'array',
  ticktext = a,
  tickvals = seq(0,length(x)-1,length.out =11),
  title = "mu2"
)
a = seq(min(y), max(y),length.out =11)
axy <- list(
  ticketmode = 'array',
  ticktext = a,
  tickvals = seq(0,length(y)-1,length.out =11),
  title = "rho"
)

hl = -0.4; loc = which(round(x+mu1,2)==hl);   

p <- plot_ly(showscale = FALSE) %>%
  add_surface(z = ~z1, colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)")))%>%
  add_surface(z = ~z2, opacity = 0.98, colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)"))) %>%
  add_surface(z = ~z3, opacity = 0.98, colorscale = list(c(0, 1), c("rgb(184,253,155)", "rgb(102,204,0)"))) %>%
  layout(scene = list(xaxis=axx,yaxis=axy, zaxis = list(title="Probability")))%>%
  add_paths(x = rep(loc-1,length(y)),
             y = 0:(length(y)-1),
             z = apply(cbind(z1[,loc],z2[,loc],z3[,loc]),1,max),
             line = list(color = "black", width = 5, dash = 'dash'))

p # 1 red 2 blue 3 green

### mu1 = other numbers

mu1 = 0.3 #mu2 = mu1 + x
y =  seq(-0.4,0.9,0.1)#rho # 2>(1-rho)^2  ->  rho \in (1-sqrt(2), 1+sqrt(2))
x = seq(-1,0.5,0.1)#mu2-mu1
z1 = matrix(NA, nrow = length(y), ncol = length(x))
z2 = matrix(NA, nrow = length(y), ncol = length(x))
z3 = matrix(NA, nrow = length(y), ncol = length(x))
for( i in 1:length(x)){
  for( j in 1:length(y)){
    z1[j,i] = mvNcdf(l = c(-mu1,x[i]), u = c(Inf,Inf), Sig = matrix(c(1,1-y[j],1-y[j],2),2))$prob
    z2[j,i] = mvNcdf(l = c(-(mu1+x[i]),-x[i]), u = c(Inf,Inf), Sig = matrix(c(1,1-y[j],1-y[j],2),2))$prob
    z3[j,i] = mvNcdf(l = c(-Inf,-Inf), u = c(-mu1, -(mu1+x[i])), Sig = matrix(c(1,y[j],y[j],1),2))$prob
  }
}


a = seq(min(x)+mu1, max(x)+mu1,length.out =11)
axx <- list(
  ticketmode = 'array',
  ticktext = a,
  tickvals = seq(0,length(x)-1,length.out =11),
  title = "mu2"
)
a = seq(min(y), max(y),length.out =11)
axy <- list(
  ticketmode = 'array',
  ticktext = a,
  tickvals = seq(0,length(y)-1,length.out =11),
  title = "rho"
)

hl = 0.4; loc = which(round(x+mu1,2)==hl);    # mu2 = hl line

p <- plot_ly(showscale = FALSE) %>%
  add_surface(z = ~z1, colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)")))%>%
  add_surface(z = ~z2, opacity = 0.98, colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)"))) %>%
  add_surface(z = ~z3, opacity = 0.98, colorscale = list(c(0, 1), c("rgb(184,253,155)", "rgb(102,204,0)"))) %>%
  layout(scene = list(xaxis=axx,yaxis=axy, zaxis = list(title="Probability")))%>%
  add_paths(x = rep(loc-1,length(y)),
            y = 0:(length(y)-1),
            z = apply(cbind(z1[,loc],z2[,loc],z3[,loc]),1,max),
            line = list(color = "black", width = 5, dash = 'dash'))

p # 1 red 2 blue 3 green

##############################################################################################
##############################################################################################
##############################################################################################

### 2D plot
#install.packages("latex2exp")
library("latex2exp")
library("TruncatedNormal")

doplot = function(mu1,mu2){
  y =  seq(-0.4,0.9,0.1)#rho
  
  z1 = rep(NA, length(y))
  z2 = rep(NA, length(y))
  z3 = rep(NA, length(y))
  
  for( j in 1:length(y)){
    z1[j] = mvNcdf(l = c(-mu1, mu2-mu1), u = c(Inf,Inf), Sig = matrix(c(1,1-y[j],1-y[j],2),2))$prob
    z2[j] = mvNcdf(l = c(-mu2, mu1-mu2), u = c(Inf,Inf), Sig = matrix(c(1,1-y[j],1-y[j],2),2))$prob
    z3[j] = mvNcdf(l = c(-Inf,-Inf), u = c(-mu1, -mu2), Sig = matrix(c(1,y[j],y[j],1),2))$prob
  }
  
  z = apply(cbind(z1,z2,z3),1,max)
  
  d1 = data.frame(y,z1)
  d2 = data.frame(y,z2)
  d3 = data.frame(y,z3)
  
  cols = c("sig", "P")
  colnames(d1) = colnames(d2) = colnames(d3) = cols
  
  d4 = data.frame(Probability = rep(c("P(Y=1)","P(Y=2)","P(Y=3)"),each = length(y)),rbind(d1,d2,d3))
  
  #p=ggplot(d4,aes(y = P,x = sig,color = Probability)) + 
  #  geom_line(size=2) + geom_vline(xintercept = c(-0.3,0.5,0.8))+
  #  xlab(TeX("$\\rho$"))+ ylab("Probability")+
  #  ggtitle(TeX(paste0("$\\mu_1 = ", mu1, ", \\mu_2 = ", mu2,"$")))
  p=ggplot(d4,aes(y = P,x = sig,color = Probability)) + 
    geom_line(size=2) +
    xlab(TeX("$\\rho$"))+ ylab("Probability")+
    ggtitle(TeX(paste0("$\\mu_1 = ", mu1, ", \\mu_2 = ", mu2,"$")))
  
  return(p)
}


library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p = vector("list",3)
p[[1]] = doplot(0,0)
p[[2]] = doplot(0,-0.2)
p[[3]] = doplot(-0.1,0.3)
mylegend<-g_legend(p[[1]])
grid.arrange(arrangeGrob(p[[1]] + theme(legend.position="none"),
                         p[[2]] + theme(legend.position="none"),
                         p[[3]] + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=1,widths=c(10, 1))

#############################################################

 try = function(a,r,x){
   return(pnorm(sqrt((1+r)/2)*(a+x/2))-abs(x)/(1-r))
 }

r = seq(-0.95,0.95,0.05)
x = seq(-1,1,0.1)
a = seq(-10,10,0.1)
res = array(NA, dim = c(length(a),length(r),length(x)) )
for(k in 1:length(a)){
  for(i in 1:length(r)){
    for(j in 1:length(x)){
      res[k,i,j] = try(a[k], r[i], x[j])   
      
    }
  }
}

compup = compdn = matrix(NA, length(a),length(x))
for(k in 1:length(a)){
  for(j in 1:length(x)){
    tt = res[k,,j]
    compup[k,j] = 1*(round(min(tt),2) == round(tt[1],2)) + 1*(round(max(tt),2) == round(tt[length(r)],2))
    compdn[k,j] = 1*(round(min(tt),2) == round(tt[length(r)],2)) + 1*(round(max(tt),2) == round(tt[1],2) || round(max(tt),2) == round(tt[2],2)|| round(max(tt),2) == round(tt[3],2))
  }
}
table(compup,compdn)
table(compup+compdn)

loc = which(compdn==1 & compup == 2)[1]
which(compdn == 1)
loc = 2135
cl = floor(loc/length(a))+1
rw = loc - (cl-1)*length(a)
plot(res[rw, , cl])

p = doplot(1.99, 2.01)
plot(p)

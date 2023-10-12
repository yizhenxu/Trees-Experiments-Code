setwd("C:/Users/yxu/Google Drive/Desktop/2015 summer/HIV Treatment Cascade/CODE/G COMPUTATION BART")

setwd("I:")
load("m0ndraws1000a1age.RData")

out[1:3,]
hist(post1)
hist(post2)

post = data.frame(time200 = post1, time400 = post2)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
p1 = ggplot(data=post, aes(time200)) + 
  geom_histogram(breaks=seq(min(post1), max(post1), by =0.0008), 
                 col="red", 
                 aes(y = (..count..)/sum(..count..),fill=..count..)) +
  scale_fill_gradient("Count", low = "green", high = "red")+ 
  scale_y_continuous(labels=percent)+ theme(legend.position="none") +
  xlab('Posterior Engagement Rate (1 ~ 200 Days)') +
  ylab('Percentage')+ annotate("text", x = 0.932, y = 0.115, label = c("a = 1"), colour = "blue", size = 6.5)#+geom_vline(xintercept = 0.8643095)



p2 = ggplot(data=post, aes(time400)) + 
  geom_histogram(breaks=seq(min(post2), max(post2), by =0.0014), 
                 col="red", 
                 aes(y = (..count..)/sum(..count..),fill=..count..)) +
  scale_fill_gradient("Count", low = "green", high = "red")+ 
  scale_y_continuous(labels=percent)+ theme(legend.position="none") +
  xlab('Posterior Engagement Rate (201 ~ 400 Days)') +
  ylab('Percentage')+ annotate('text', x = 0.77, y = 0.12, label ='a[1] == 1', colour = "blue", size = 6,parse=TRUE)+geom_vline(xintercept = 0.6649428)

plot(p2)

t1 = arrangeGrob(p1 + theme(legend.position="none"),ncol=1,nrow=1)
t2 = arrangeGrob(p2 + theme(legend.position="none"),ncol=1,nrow=1)

png(filename = "HISTPLOT.png",width = 550, height = 275)
g4 = grid.arrange(t1,t2, 
                  nrow=1,ncol = 2, widths = c(1,1))
dev.off()
out

mean(post1)
mean(post2)
save(post1,post2,post,p1,p2,out, file = "readytogo.RData")

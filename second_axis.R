d <- data.frame(name=c("zhao","qian","sun","li"),weight=c(62,58,79,60),height=c(178,169,180,173))
x <- d$name
y1 <- d$weight
y2 <- d$height


###########
###  plotrix
###########

#png("2y_1.png")
library(plotrix)
xpos <- 1:4
twoord.plot(xpos,y1,xpos,y2,xlim=c(0,5),lylim=c(0,100),rylim=c(100,190), lcol=4,rcol=2,xlab="name",ylab="Weight (kg)",rylab="Height (cm)",type=c("bar","b"),xticklab=x,halfwidth=0.2)
#dev.off()


###########
## TeachingDemos
###########

#png("2y_2.png")
library(TeachingDemos)
par(mar=c(5,5,4,5)+0.1)
bar <- barplot(y1,ylim=c(0,100),ylab="Weight (kg)",col="blue",col.axis="blue",col.lab="blue")
updateusr(1:2,range(0,100),1:2,range(100,190))
lines(bar,y2,type="b",col="red")
axis(4,col="red",col.ticks="red",col.axis="red")
mtext("Heigth (cm)",side=4,line=3,col="red")
mtext(x,side=1,line=1,at=bar,col="black")
#dev.off()


########
## no package
########
#png("2y_3.png")
par(mar=c(5,5,4,5)+0.1)
bar <- barplot(y1,xlim=c(0,5),ylim=c(0,100),ylab="Weight (kg)",col="blue",col.axis="blue",col.lab="blue")

mtext(x,side=1,line=1,at=bar,col="black")
mtext("Name",side=1,line=3,col="black")
par(new=T)
plot(bar,y2,axes=F,xlim=c(0,5),ylim=c(100,190),xlab="",ylab="",col="red",type="b")
axis(4,col="red",col.ticks="red",col.axis="red")
mtext("Heigth (cm)",side=4,line=3,col="red")
#dev.off()





#load the different results

load("partitionBHmax101to200.RData")
partitionmaxscores<-maxyscores
load("edgerevBHmax101to200.RData")
edgerevmaxscores<-maxyscores
load("partitionBHedgerevmax101to200.RData")
combinedmaxscores<-maxyscores

#png("BHmaxscoredensity.png", width=7.5, height=5, unit="in",res=100,bg="transparent")
pdf("BHmaxscoredensity.pdf", width=7.5, height=5, onefile=F, pointsize=10,  paper="special")

# plotting parameters

par(mar=c(2.25,2.75,0.5,0.5))
par(mgp=c(3,1,0))
par(cex.axis=1.25)
par(cex.lab=1.5)

# kernel density plots 

edged<-density(edgerevmaxscores,bw=1,adjust=1,from=-20,to=4)
partd<-density(partitionmaxscores,bw=1,adjust=1,from=-20,to=4)
combd<-density(combinedmaxscores,bw=1,adjust=1,from=-20,to=4)
plot(edged,col=rgb(0.8,0.4,0),lwd=2.5,lty=5,ylab="",xlab="",main="")
polygon(partd,col=rgb(0,0.4,0.8,0.35),border=rgb(0,0.4,0.8))
polygon(edged,col=rgb(0.8,0.4,0,0.35),border=rgb(0.8,0.4,0))
polygon(combd,col=rgb(0.5,0.2,0.5,0.35),border=rgb(0.5,0.2,0.5))
abline(h=0)
lines(partd,col=rgb(0,0.4,0.8),lwd=3)
lines(edged,col=rgb(0.8,0.4,0),lwd=2.5,lty=5)
lines(combd,col=rgb(0.5,0.2,0.5),lwd=4,lty=3)

# save results

dev.off()

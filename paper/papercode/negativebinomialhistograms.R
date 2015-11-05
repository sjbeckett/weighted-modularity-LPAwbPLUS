setEPS()
postscript("negativebinomial.eps",height=5.5,width=9)

par(mfrow=c(1,2))

hist(rnbinom(10000,size=0.5,mu=4),ylim=c(0,10000),main=NULL,100,xlab="size=0.5")
hist(rnbinom(10000,size=2.5,mu=4),ylim=c(0,10000),main=NULL,100,xlab="size=2.5")

dev.off()

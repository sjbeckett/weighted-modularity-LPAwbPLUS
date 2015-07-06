setEPS()
postscript("negativebinomial.eps",height=5.5,width=9)

par(mfrow=c(1,2))

hist(rnbinom(10000,size=0.2,mu=2),ylim=c(0,10000),main=NULL,100,xlab="size=0.2")
hist(rnbinom(10000,size=1,mu=2),ylim=c(0,10000),main=NULL,100,xlab="size=1")

dev.off()

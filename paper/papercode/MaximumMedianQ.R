
QQ = read.csv("output/summary/summaryQuaQBM.csv")
QL = read.csv("output/summary/summaryQuaLPAwb+.csv")
QE = read.csv("output/summary/summaryQuaEXLPAwb+.csv")

BQ = read.csv("output/summary/summaryBinQBM.csv")
BL = read.csv("output/summary/summaryBinLPAwb+.csv")
BE = read.csv("output/summary/summaryBinEXLPAwb+.csv")


MAX_QW = apply(cbind(QQ$MAX_QW_QQ, QL$MAX_QW_QL, QE$MAX_QW_QE), 1, max)
MAX_QB = apply(cbind(BQ$MAX_QW_BQ, BL$MAX_QW_BL, BE$MAX_QW_BE), 1, max)


setEPS()
postscript("MedianMaximum.eps",height=5.5,width=9)
par(mfrow = c(1,2),oma=c(2.5,0,0,0))

plot(0,0, type='n', ylab = expression("Median Q"[B]*" per algorithm"), xlab = expression("Consensus Maximum Q"[B]), xlim=c(0,0.7),ylim=c(0,0.7),bty="n")
axis(1,col="grey",labels=FALSE,lwd=2)
axis(2,col="grey",labels=FALSE,lwd=2)
lines(0:1,0:1,lty=3,lwd=4,col="grey80")
points(MAX_QB,BQ$MEDIAN_QW_BQ)
points(MAX_QB,BL$MEDIAN_QW_BL,pch=3,col="grey30")
points(MAX_QB,BE$MEDIAN_QW_BE,pch=4,col="grey65")



plot(0,0,type='n', ylab=expression("Median Q"[W]*" per algorithm"), xlab = expression("Consensus Maximum Q"[W]), xlim=c(0,0.7),ylim=c(0,0.7),bty="n")
axis(1,col="grey",labels=FALSE,lwd=2)
axis(2,col="grey",labels=FALSE,lwd=2)
lines(0:1,0:1,lty=3,lwd=4,col="grey80")
points(MAX_QW,QQ$MEDIAN_QW_QQ)
points(MAX_QW,QL$MEDIAN_QW_QL,pch=3,col="grey30")
points(MAX_QW,QE$MEDIAN_QW_QE,pch=4,col="grey65")



par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("bottom",c("QuanBiMo","LPAwb+","DIRTLPAwb+"),pch=c(1,3,4),col=c("black","grey30","grey65"),bty="n",horiz=TRUE)

#dev.copy2eps(file="Maximumvs.Median.eps")
dev.off()


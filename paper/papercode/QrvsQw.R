##Binary data

A=read.csv("output/summary/summaryBinEXLPAwb+.csv",row.names=1)
B=read.csv("output/summary/summaryBinLPAwb+.csv",row.names=1)
C=read.csv("output/summary/summaryBinQBM.csv",row.names=1)


BinMod = cbind(A$MAX_QW_BE,B$MAX_QW_BL,C$MAX_QW_BQ)

BQW = apply(BinMod,1,max)

BQR = c()


for(each in 1:23) {
	H=which(BinMod[each,] == BQW[each])[1]
	if(H==1) {
		BQR[each] = A$QR_BE[each]
	} else if(H==2) {
		BQR[each] = B$QR_BL[each]
	} else {
		BQR[each] = C$QR_BQ[each]
	}
}

## Quantitative data

A=read.csv("output/summary/summaryQuaEXLPAwb+.csv",row.names=1)
B=read.csv("output/summary/summaryQuaLPAwb+.csv",row.names=1)
C=read.csv("output/summary/summaryQuaQBM.csv",row.names=1)


QuaMod = cbind(A$MAX_QW_QE,B$MAX_QW_QL,C$MAX_QW_QQ)

QQW = apply(QuaMod,1,max)

QQR = c()


for(each in 1:23) {
	H=which(QuaMod[each,] == QQW[each])[1]
	if(H==1) {
		QQR[each] = A$QR_QE[each]
	} else if(H==2) {
		QQR[each] = Q$QR_BL[each]
	} else {
		QQR[each] = C$QR_QQ[each]
	}
}

#### Modularity vs. Realised Modularity

plot(0,0,type='n',xlim=c(0.2,1),ylim=c(-0.2,1),ylab=expression("Q'"[R]),xlab="Modularity",bty="n")
axis(1,col="grey",labels=FALSE,lwd=2)
axis(2,col="grey",labels=FALSE,lwd=2)
lines(0:1,rep(0,2),lty=3,lwd=4,col="grey80")
points(BQW,BQR)
for( aa in 1:23 ) {
lines(c(BQW[aa],QQW[aa]),c(BQR[aa],QQR[aa]),col="indianred1")
}
points(QQW,QQR,pch=4)
legend(0.2,0.95,c(expression("Binary Q"[B]),expression("Quantitative Q"[W])),pch=c(1,4),bty='n')

dev.copy2eps(file="ModRealisedModularity.eps")
dev.off()




#### Normalised Modularity vs. Realised Modularity

A=read.csv("output/summary/MAXMODoutput.csv",row.names=1)

MBM_B = A$MBM_B #maximum binary modularity
MBM_Q = A$MBM_Q # maximum weighted modularity

NQB = BQW/MBM_B # normalised binary modularity
NQW = QQW/MBM_Q # normalised weighted modularity


plot(0,0,type='n',xlim=c(0.2,1),ylim=c(-0.2,1),ylab=expression("Q'"[R]),xlab="Normalised Modularity",bty="n")
axis(1,col="grey",labels=FALSE,lwd=2)
axis(2,col="grey",labels=FALSE,lwd=2)
lines(0:1,rep(0,2),lty=3,lwd=4,col="grey80")

for( aa in 1:23 ) {
lines(c(NQW[aa],NQB[aa]),c(QQR[aa],BQR[aa]),col="indianred1")
}
points(NQB,BQR)
points(NQW,QQR,pch=4)


legend(0.25,0.95,c(expression("Binary Q"[B]^norm),expression("Quantitative Q"[W]^norm)),pch=c(1,4),bty='n')


dev.copy2eps(file="NormQvsRealisedModularity.eps")
dev.off()














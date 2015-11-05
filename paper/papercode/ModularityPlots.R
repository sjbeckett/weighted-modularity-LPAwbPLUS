
B_LPA = read.csv("output/summary/Bin_Modscore_LPA_wb_plus.csv",row.names=1)
B_EXL = read.csv("output/summary/Bin_Modscore_EXLPA_wb_plus.csv",row.names=1)
B_QBM = read.csv("output/summary/Bin_Modscore_QuanBiMo.csv",row.names=1)

Q_LPA = read.csv("output/summary/Qua_Modscore_LPA_wb_plus.csv",row.names=1)
Q_EXL = read.csv("output/summary/Qua_Modscore_EXLPA_wb_plus.csv",row.names=1)
Q_QBM = read.csv("output/summary/Qua_Modscore_QuanBiMo.csv",row.names=1)


BL = apply(B_LPA,1,max,na.rm=TRUE)
BE = apply(B_EXL,1,max,na.rm=TRUE)
BQ = apply(B_QBM,1,max,na.rm=TRUE)

QL = apply(Q_LPA,1,max,na.rm=TRUE)
QE = apply(Q_EXL,1,max,na.rm=TRUE)
QQ = apply(Q_QBM,1,max,na.rm=TRUE)





colw="grey"
setEPS()
postscript("MaxModularity.eps",height=5.5,width=9)
par(mfrow = c(1,2),oma=c(2.5,0,0,0))


plot(seq(0,1,0.1),seq(0,1,0.1),type='n',col=colw,ylim=c(0,0.7),xlim=c(0,0.7),ylab = expression("New algorithms maximum Q"[B]),xlab =expression("Maximum QuanBiMo Q"[B]),bty="n")
lines(0:1,0:1,lty=3,lwd=4,col="grey80")
points(BQ,BL,pch=3, col = "grey30")
points(BQ,BE,pch=4, col = "grey65")
axis(1,col=colw,labels=FALSE,lwd=2)
axis(2,col=colw,labels=FALSE,lwd=2)

plot(seq(0,1,0.1),seq(0,1,0.1),type="n",lwd=2,col=colw,ylim=c(0,0.7),xlim=c(0,0.7),ylab = expression("New algorithms maximum Q"[W]),xlab =expression("Maximum QuanBiMo Q"[W]),bty="n")
lines(0:1,0:1,lty=3,lwd=4,col="grey80")
points(QQ,QL,pch=3, col = "grey30")
points(QQ,QE,pch=4, col = "grey65")
axis(1,col=colw,labels=FALSE,lwd=2)
axis(2,col=colw,labels=FALSE,lwd=2)


par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("bottom",c("LPAwb+","DIRTLPAwb+"),pch=c(3,4),col=c("grey30","grey65"),bty="n",horiz=TRUE)


#dev.copy2eps(file="Modularity.eps")
dev.off()

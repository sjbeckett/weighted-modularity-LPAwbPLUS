 BQ=read.csv("output/summary/summaryBinQBM.csv",row.names=1)
 BL=read.csv("output/summary/summaryBinLPAwb+.csv",row.names=1)
 BE=read.csv("output/summary/summaryBinEXLPAwb+.csv",row.names=1)

 QQ=read.csv("output/summary/summaryQuaQBM.csv",row.names=1)
 QL=read.csv("output/summary/summaryQuaLPAwb+.csv",row.names=1)
 QE=read.csv("output/summary/summaryQuaEXLPAwb+.csv",row.names=1)


 t1 = BQ$AVG_time_BQ
 t2 = BL$AVG_time_BL
 t3 = BE$AVG_time_BE

 t4 = QQ$AVG_time_QQ
 t5 = QL$AVG_time_QL
 t6 = QE$AVG_time_QE


 T1 = log(t1)
 T2 = log(t2)
 T3 = log(t3)

 T4 = log(t4)
 T5 = log(t5)
 T6 = log(t6)

 MIN=min(T1,T2,T3,T4,T5,T6)
 MAX=max(T1,T2,T3,T4,T5,T6)

 NETWORKS = list("Safariland","barrett1987","bezerra2009","elberling1999","inouye1988","junker2013","kato1990","kevan1970","memmott1999","mosquin1967","motten1982","olesen2002aigrettes","olesen2002flores","ollerton2003","schemske1978","small1976","vazarr","vazcer","vazllao","vazmasc","vazmasnc","vazquec","vazquenc")

h = 1:23

 setEPS()
 postscript("Timings.eps",height=6,width=12)
 par(mfrow = c(1,2),oma=c(2.5,0,0,0))

 par(mar = c(7, 4, 4, 2) + 0.1)
 
 plot(0,0,type='n',ylab="log(computation time (seconds))", xaxt = "n", xlab="", ylim=c(MIN,MAX), xlim=c(1,23) ,bty='n')
 axis(1, at=1:23,labels = FALSE,col="grey",lwd=2)
 axis(2, labels = FALSE, col="grey", lwd=2)
 text(1:23, par("usr")[1] - 6.5, srt = 45, adj = c(1,0.5),
     labels = NETWORKS, xpd = TRUE)
 points(h,T1)
 points(h,T2,pch=3,col="grey30")
 points(h,T3,pch=4,col="grey65")

 plot(0,0,type='n',ylab="log(computation time (seconds))", xaxt= "n", xlab="", ylim=c(MIN,MAX), xlim=c(1,23),bty='n')
 axis(1, at=1:23,labels = FALSE, col="grey", lwd=2)
 axis(2, labels = FALSE, col="grey", lwd=2)
 text(1:23, par("usr")[1] - 6.5, srt = 45, adj = c(1,0.5),
     labels = NETWORKS, xpd = TRUE)
 points(h,T4)
 points(h,T5,pch=3,col="grey30")
 points(h,T6,pch=4,col="grey65")

 par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
 plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
 legend("bottom",c("QuanBiMo","LPAwb+","DIRTLPAwb+"),pch=c(1,3,4),col=c("black","grey30","grey65"),bty="n",horiz=TRUE)

 dev.off()


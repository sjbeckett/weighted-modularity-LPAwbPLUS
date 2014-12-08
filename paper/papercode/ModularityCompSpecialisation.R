
library(bipartite)

DATALIST = list(Safariland,barrett1987,bezerra2009,elberling1999,inouye1988,junker2013,kato1990,kevan1970,memmott1999,mosquin1967,motten1982,olesen2002aigrettes,olesen2002flores,ollerton2003,schemske1978,small1976,vazarr,vazcer,vazllao,vazmasc,vazmasnc,vazquec,vazquenc)


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

MODBIN = apply(cbind(BL,BE,BQ),1,max)
MODQUA = apply(cbind(QL,QE,QQ),1,max)

H2QUA = c()
num = 1

for(network in DATALIST){
	print(num)
	Qua=network[rowSums(network)>0,colSums(network)>0]
	Bin=1*(Qua>0)
	H2QUA[num] = H2fun(Qua)[1]
	
	num = num+1
}

#H2 needs quantitative data!

plot(H2QUA,MODQUA,xlab="Complementary specialisation H2'",ylab="Modularity",bty="n",xlim=c(0,1),ylim=c(0.2,0.8))
points(H2QUA,MODBIN,pch=3)
legend(0.05,0.7,c(expression("Q"[W]),expression("Q"[B])),pch=c(1,3),box.col="grey90",bg="grey90")


source("NormalisedMutualInformation.R")
source("MaximumBipartiteModularity.R")


DATALIST = list(Safariland,barrett1987,bezerra2009,elberling1999,inouye1988,junker2013,kato1990,kevan1970,memmott1999,mosquin1967,motten1982,olesen2002aigrettes,olesen2002flores,ollerton2003,schemske1978,small1976,vazarr,vazcer,vazllao,vazmasc,vazmasnc,vazquec,vazquenc)


Networks =c("Safariland","barrett1987","bezerra2009","elberling1999","inouye1988","junker2013","kato1990","kevan1970","memmott1999","mosquin1967","motten1982","olesen2002aigrettes","olesen2002flores","ollerton2003","schemske1978","small1976","vazarr","vazcer","vazllao","vazmasc","vazmasnc","vazquec","vazquenc")

FRONT = "output/summary/"
##
filestart = "Qua"
front1 = "Q"

#import data
A = as.matrix(read.csv(paste(FRONT,filestart,"_Modscore_QuanBiMo.csv",sep=""), row.names=1))
B = as.matrix(read.csv(paste(FRONT,filestart,"_Modscore_LPA_wb_plus.csv",sep=""), row.names=1))
C = as.matrix(read.csv(paste(FRONT,filestart,"_Modscore_EXLPA_wb_plus.csv",sep=""), row.names=1))	

#find maximum modularity for each row and position
AA = apply(A,1,max,na.rm =TRUE)
BB = apply(B,1,max,na.rm =TRUE)
CC = apply(C,1,max,na.rm =TRUE)

#list of run(s) where maximum Modularity achieved.
LSA_q=list()
LSB_q=list()
LSC_q=list()
for(data in 1:dim(A)[1]) {
	LSA_q = append(LSA_q,list(which(A[data,]==AA[data])))
	LSB_q = append(LSB_q,list(which(B[data,]==BB[data])))
	LSC_q = append(LSC_q,list(which(C[data,]==CC[data])))
}

maxesQua = cbind(AA,BB,CC)
QW = apply(maxesQua,1,max)   #Best weighted Modularity



##
filestart = "Bin"
front1 = "bin"

#import data
A = as.matrix(read.csv(paste(FRONT,filestart,"_Modscore_QuanBiMo.csv",sep=""), row.names=1))
B = as.matrix(read.csv(paste(FRONT,filestart,"_Modscore_LPA_wb_plus.csv",sep=""), row.names=1))
C = as.matrix(read.csv(paste(FRONT,filestart,"_Modscore_EXLPA_wb_plus.csv",sep=""), row.names=1))

#find maximum modularity for each row ignoring NA's
AA = apply(A,1,max,na.rm =TRUE)
BB = apply(B,1,max,na.rm =TRUE)
CC = apply(C,1,max,na.rm =TRUE)


#list of run(s) where maximum Modularity achieved.
LSA_b=list()
LSB_b=list()
LSC_b=list()
for(data in 1:dim(A)[1]) {
	LSA_b = append(LSA_b,list(which(A[data,]==AA[data])))
	LSB_b = append(LSB_b,list(which(B[data,]==BB[data])))
	LSC_b = append(LSC_b,list(which(C[data,]==CC[data])))
}

maxesBin = cbind(AA,BB,CC)
QB = apply(maxesBin,1,max) #Best binary Modularity



NMI=c() # to store NMI scores
MBM_B=c() # to store maximum binary bipartite modularity
MBM_Q=c() # to store maximum weighted bipartite modularity


for(data in 1:dim(A)[1]) {

	aa=as.matrix(DATALIST[[data]])
	Q <- as.matrix(aa) #network
	#convert to binary network
	A = 1*(Q>0)
	#Remove empty Rows/Columns from analysis
	if(length(which(rowSums(A)==0))>0) {
		Q = Q[-which(rowSums(A)==0),]
		A = A[-which(rowSums(A)==0),] }
	if(length(which(colSums(A)==0))>0) {
		Q = Q[,-which(colSums(A)==0)]
		A = A[,-which(colSums(A)==0)] }



	DIREFRONT = "output/"
	

	Algorithm_q = which(maxesQua[data,] == QW[data])[1]
	if(Algorithm_q == 1) {
		ROWS = as.matrix(read.csv(paste(DIREFRONT,"QrowsQBM",Networks[data],".csv",sep=""),row.names=1))
		Rchoice = ROWS[LSA_q[[data]][[1]],]
		COLS = as.matrix(read.csv(paste(DIREFRONT,"QcolsQBM",Networks[data],".csv",sep=""),row.names=1))
		Cchoice = COLS[LSA_q[[data]][[1]],]
	} else if(Algorithm_q==2) {
		ROWS = as.matrix(read.csv(paste(DIREFRONT,"QrowsLPAwb+",Networks[data],".csv",sep=""),row.names=1))
		Rchoice = ROWS[LSB_q[[data]][[1]],]
		COLS = as.matrix(read.csv(paste(DIREFRONT,"QcolsLPAwb+",Networks[data],".csv",sep=""),row.names=1))
		Cchoice = COLS[LSB_q[[data]][[1]],]
	} else {
		ROWS = as.matrix(read.csv(paste(DIREFRONT,"QrowsEXLPAwb+",Networks[data],".csv",sep=""),row.names=1))
		Rchoice = ROWS[LSC_q[[data]][[1]],]
		COLS = as.matrix(read.csv(paste(DIREFRONT,"QcolsEXLPAwb+",Networks[data],".csv",sep=""),row.names=1))
		Cchoice = COLS[LSC_q[[data]][[1]],]
	}

	QuanList = c(Rchoice,Cchoice)

	MBM_Q[data] = MAXIMUM_BIP_MOD( Q,Rchoice,Cchoice)
	
	Algorithm_b = which(maxesBin[data,] == QB[data])[1]
	if(Algorithm_b == 1) {
		ROWS = as.matrix(read.csv(paste(DIREFRONT,"binrowsQBM",Networks[data],".csv",sep=""),row.names=1))
		Rchoice = ROWS[LSA_b[[data]][[1]],]
		COLS = as.matrix(read.csv(paste(DIREFRONT,"bincolsQBM",Networks[data],".csv",sep=""),row.names=1))
		Cchoice = COLS[LSA_b[[data]][[1]],]
	} else if(Algorithm_b==2) {
		ROWS = as.matrix(read.csv(paste(DIREFRONT,"binrowsLPAwb+",Networks[data],".csv",sep=""),row.names=1))
		Rchoice = ROWS[LSB_b[[data]][[1]],]
		COLS = as.matrix(read.csv(paste(DIREFRONT,"bincolsLPAwb+",Networks[data],".csv",sep=""),row.names=1))
		Cchoice = COLS[LSB_b[[data]][[1]],]
	} else {
		ROWS = as.matrix(read.csv(paste(DIREFRONT,"binrowsEXLPAwb+",Networks[data],".csv",sep=""),row.names=1))
		Rchoice = ROWS[LSC_b[[data]][[1]],]
		COLS = as.matrix(read.csv(paste(DIREFRONT,"bincolsEXLPAwb+",Networks[data],".csv",sep=""),row.names=1))
		Cchoice = COLS[LSC_b[[data]][[1]],]
	}

	BinList = c(Rchoice,Cchoice)

	MBM_B[data] = MAXIMUM_BIP_MOD(A,Rchoice,Cchoice)
	
	NMI[data] = NormalisedMutualInformation(QuanList,BinList)	
}

NMIout = as.matrix(NMI)
rownames(NMIout) = Networks

write.csv(NMIout,paste(FRONT,"NMIoutput.csv",sep=""))


MaxModularityOut = as.matrix(cbind(MBM_B,MBM_Q))
rownames(MaxModularityOut) = Networks

write.csv(MaxModularityOut,paste(FRONT,"MAXMODoutput.csv",sep=""))

changeQ = QW-QB

NQB = QB/MBM_B #normalised binary modularity
NQW = QW/MBM_Q #normalised weighted modularity

changeQnorm = NQW - NQB


colw = "grey"

plot(0,0,type='n',xlab="NMI",ylab=expression(paste(Delta, "Q"^norm)),bty="n",ylim=c(min(changeQnorm),0.3),xlim=c(0,1))
lines(0:1,rep(0,2),lty=3,lwd=4,col="grey80")
points(NMI,changeQnorm)
axis(1,col=colw,labels=FALSE,lwd=2)
axis(2,col=colw,labels=FALSE,lwd=2)

subset = c(1,4,5,7,14,15)#6,8,9,13,14,19) # labels to put above point
subset2 = c(19,22)

text(NMI[subset], changeQnorm[subset], labels = Networks[subset], cex= 0.7, pos=3, offset = 0.4)  #above
text(NMI[subset2], changeQnorm[subset2], labels = Networks[subset2], cex= 0.7, pos=4, offset = 0.4)  #above
text(NMI[-c(subset,subset2)], changeQnorm[-c(subset,subset2)], labels = Networks[-c(subset,subset2)], cex= 0.7, pos=1, offset = 0.4)  #below



dev.copy2eps(file="changeQnorm_NMI.eps")
dev.off()



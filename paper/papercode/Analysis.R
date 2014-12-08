#Analysis script

library(bipartite) # import bipartite library for datasets

#Function to find realized modularity
RealizedMod <- function(MATRIX,redlabels,bluelabels) {
	
	E = sum(MATRIX)
	W=0
	for(aa in 1:nrow(MATRIX)) {
		for (bb in 1:ncol(MATRIX)) {
			if(redlabels[aa] == bluelabels[bb]) {
				W = W + MATRIX[aa,bb]
			}
		}
	}

	return( 2*(W/E) - 1 )
}



#DATA
DATALIST = list(Safariland,barrett1987,bezerra2009,elberling1999,inouye1988,junker2013,kato1990,kevan1970,memmott1999,mosquin1967,motten1982,olesen2002aigrettes,olesen2002flores,ollerton2003,schemske1978,small1976,vazarr,vazcer,vazllao,vazmasc,vazmasnc,vazquec,vazquenc)

NETWORKS = list("Safariland","barrett1987","bezerra2009","elberling1999","inouye1988","junker2013","kato1990","kevan1970","memmott1999","mosquin1967","motten1982","olesen2002aigrettes","olesen2002flores","ollerton2003","schemske1978","small1976","vazarr","vazcer","vazllao","vazmasc","vazmasnc","vazquec","vazquenc")


front1 = "output/summary/"
front2 = "output/"
fQua = "Qua"
fq = "Q"
fBin = "Bin"
fbin = "bin"
frow = "rows"
fcol = "cols"
fEL="EXLPAwb+"
fL ="LPAwb+"
fQ = "QBM"


source("NormalisedMutualInformation.R")

#LOAD DATA for modularity and timings

BinModularityScore_LPA_wb_plus = as.matrix(read.csv(paste(front1,fBin,"_Modscore_LPA_wb_plus.csv",sep=""),row.names=1))
BinModularityScore_QuanBiMo = as.matrix(read.csv(paste(front1,fBin,"_Modscore_QuanBiMo.csv",sep=""),row.names=1))
BinModularityScore_EXLPAwbplus = as.matrix(read.csv(paste(front1,fBin,"_Modscore_EXLPA_wb_plus.csv",sep=""),row.names=1))
BinTimeData_LPA_wb_plus = as.matrix(read.csv(paste(front1,fBin,"_Times_LPA_wb_plus.csv",sep=""),row.names=1))
BinTimeData_QuanBiMo = as.matrix(read.csv(paste(front1,fBin,"_Times_QuanBiMo.csv",sep=""),row.names=1))
BinTimeData_EXLPAwbplus = as.matrix(read.csv(paste(front1,fBin,"_Times_EXLPA_wb_plus.csv",sep=""),row.names=1))


ModularityScore_LPA_wb_plus = as.matrix(read.csv(paste(front1,fQua,"_Modscore_LPA_wb_plus.csv",sep=""),row.names=1))
ModularityScore_QuanBiMo = as.matrix(read.csv(paste(front1,fQua,"_Modscore_QuanBiMo.csv",sep=""),row.names=1))
ModularityScore_EXLPAwbplus = as.matrix(read.csv(paste(front1,fQua,"_Modscore_EXLPA_wb_plus.csv",sep=""),row.names=1))
TimeData_LPA_wb_plus = as.matrix(read.csv(paste(front1,fQua,"_Times_LPA_wb_plus.csv",sep=""),row.names=1))
TimeData_QuanBiMo = as.matrix(read.csv(paste(front1,fQua,"_Times_QuanBiMo.csv",sep=""),row.names=1))
TimeData_EXLPAwbplus = as.matrix(read.csv(paste(front1,fQua,"_Times_EXLPA_wb_plus.csv",sep=""),row.names=1))


#COMPUTE BASIC STATISTICS
#BL - binary LPAwb+
#QL - quantitative LPAwb+
#BQ - binary QuanBiMo
#QQ - quantitative QuanBiMo
#BE - binary ExhaustiveLPAwb+
#QE - quantitative ExhaustiveLPAwb+

AVG_time_BL = rowMeans(BinTimeData_LPA_wb_plus)
AVG_time_QL = rowMeans(TimeData_LPA_wb_plus)
AVG_time_BQ = rowMeans(BinTimeData_QuanBiMo)
AVG_time_QQ = rowMeans(TimeData_QuanBiMo)
AVG_time_BE = rowMeans(BinTimeData_EXLPAwbplus)
AVG_time_QE = rowMeans(TimeData_EXLPAwbplus)


MAX_QW_BL= apply(BinModularityScore_LPA_wb_plus,1,max,na.rm=TRUE)
MAX_QW_QL= apply(ModularityScore_LPA_wb_plus,1,max,na.rm=TRUE)
MAX_QW_BQ = apply(BinModularityScore_QuanBiMo,1,max,na.rm=TRUE)
MAX_QW_QQ = apply(ModularityScore_QuanBiMo,1,max,na.rm=TRUE)
MAX_QW_BE = apply(BinModularityScore_EXLPAwbplus,1,max,na.rm=TRUE)
MAX_QW_QE = apply(ModularityScore_EXLPAwbplus,1,max,na.rm=TRUE)



MEDIAN_QW_BL = apply(BinModularityScore_LPA_wb_plus,1,median,na.rm=TRUE)
MEDIAN_QW_QL = apply(ModularityScore_LPA_wb_plus,1,median,na.rm=TRUE)
MEDIAN_QW_BQ = apply(BinModularityScore_QuanBiMo,1,median,na.rm=TRUE)
MEDIAN_QW_QQ = apply(ModularityScore_QuanBiMo,1,median,na.rm=TRUE)
MEDIAN_QW_BE = apply(BinModularityScore_EXLPAwbplus,1,median,na.rm=TRUE)
MEDIAN_QW_QE = apply(ModularityScore_EXLPAwbplus,1,median,na.rm=TRUE)

MIN_QW_BL = apply(BinModularityScore_LPA_wb_plus,1,min,na.rm=TRUE)
MIN_QW_QL = apply(ModularityScore_LPA_wb_plus,1,min,na.rm=TRUE)
MIN_QW_BQ = apply(BinModularityScore_QuanBiMo,1,min,na.rm=TRUE)
MIN_QW_QQ = apply(ModularityScore_QuanBiMo,1,min,na.rm=TRUE)
MIN_QW_BE = apply(BinModularityScore_EXLPAwbplus,1,min,na.rm=TRUE)
MIN_QW_QE = apply(ModularityScore_EXLPAwbplus,1,min,na.rm=TRUE)

MEAN_QW_BL = apply(BinModularityScore_LPA_wb_plus,1,mean,na.rm=TRUE)
MEAN_QW_QL = apply(ModularityScore_LPA_wb_plus,1,mean,na.rm=TRUE)
MEAN_QW_BQ = apply(BinModularityScore_QuanBiMo,1,mean,na.rm=TRUE)
MEAN_QW_QQ = apply(ModularityScore_QuanBiMo,1,mean,na.rm=TRUE)
MEAN_QW_BE = apply(BinModularityScore_EXLPAwbplus,1,mean,na.rm=TRUE)
MEAN_QW_QE = apply(ModularityScore_EXLPAwbplus,1,mean,na.rm=TRUE)



NUMFAIL_BL = apply(BinModularityScore_LPA_wb_plus,1,function(x) sum(is.na(x)))
NUMFAIL_QL = apply(ModularityScore_LPA_wb_plus,1,function(x) sum(is.na(x)))
NUMFAIL_BQ = apply(BinModularityScore_QuanBiMo,1,function(x) sum(is.na(x)))
NUMFAIL_QQ = apply(ModularityScore_QuanBiMo,1,function(x) sum(is.na(x)))
NUMFAIL_BE = apply(BinModularityScore_EXLPAwbplus,1,function(x) sum(is.na(x)))
NUMFAIL_QE = apply(ModularityScore_EXLPAwbplus,1,function(x) sum(is.na(x)))

NUMMAX_BL = c()
NUMMAX_QL = c()
NUMMAX_BQ = c()
NUMMAX_QQ = c()
NUMMAX_BE = c()
NUMMAX_QE = c()
for( aa in 1:dim(ModularityScore_LPA_wb_plus)[1] ) {
	NUMMAX_BL[aa] = length(which(BinModularityScore_LPA_wb_plus[aa,] == MAX_QW_BL[aa]))
	NUMMAX_QL[aa] = length(which(ModularityScore_LPA_wb_plus[aa,] == MAX_QW_QL[aa]))
	NUMMAX_BQ[aa] = length(which(BinModularityScore_QuanBiMo[aa,] == MAX_QW_BQ[aa]))
	NUMMAX_QQ[aa] = length(which(ModularityScore_QuanBiMo[aa,] == MAX_QW_QQ[aa]))
	NUMMAX_BE[aa] = length(which(BinModularityScore_EXLPAwbplus[aa,] == MAX_QW_BE[aa]))
	NUMMAX_QE[aa] = length(which(ModularityScore_EXLPAwbplus[aa,] == MAX_QW_QE[aa]))
}





#ADVANCED


MODULES_BL = c()
UniqueConfigs_BL = c()
MODULES_QL = c()
UniqueConfigs_QL = c()
MODULES_BQ = c()
UniqueConfigs_BQ = c()
MODULES_QQ = c()
UniqueConfigs_QQ = c()
MODULES_BE = c()
UniqueConfigs_BE = c()
MODULES_QE = c()
UniqueConfigs_QE = c()

QR_BL =c()
QR_QL =c()
QR_BQ =c()
QR_QQ =c()
QR_BE =c()
QR_QE =c()

NMI1 =c()
NMI2 =c()
NMI3 =c()
NMI4 =c()
NMI5 =c()
NMI6 =c()

data = 1

##FIND #MODULES and #Configurations

for( aa in DATALIST ) {
	print(NETWORKS[data])
#read in network	
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

#binary LP
	B_LP_R = as.matrix(read.csv(paste(front2,fbin,frow,fL,NETWORKS[data],".csv",sep=""),row.names=1))
	B_LP_C = as.matrix(read.csv(paste(front2,fbin,fcol,fL,NETWORKS[data],".csv",sep=""),row.names=1))

	HAM=which(BinModularityScore_LPA_wb_plus[data,] == MAX_QW_BL[data])
	BLPindices = HAM
	HMM = c()
	MODULES = c()
	MODULES[1] = length(unique(B_LP_R[HAM[1],]))
	Realized=c()
	Realized[1] = RealizedMod(A,B_LP_R[HAM[1],],B_LP_C[HAM[1],])
	QR_BL[data] = Realized[1]
	HAH=1
	if(length(HAM)>1) {
		for( bb in 2:length(HAM)) {
			HMM[bb-1] = NormalisedMutualInformation(c(B_LP_R[HAM[1],],B_LP_C[HAM[1],]) , c(B_LP_R[HAM[bb],],B_LP_C[HAM[bb],])) #compare first best with rest NMI of 1 indicates same config.
			MODULES[bb] = length(unique(B_LP_R[HAM[bb],]))
			Realized[bb] = RealizedMod(A,B_LP_R[HAM[bb],],B_LP_C[HAM[bb],])
		}
		
		if(length(unique(Realized))>1) {
			print("multiple Qr values: binary LPAwb+")
			print(NETWORKS[[data]])
		}
	HAH=unique(HMM)
	}
	UniqueConfigs_BL[data] = length(HAH)
	if( length(unique(MODULES)) == 1) {
		MODULES_BL[data] = MODULES[1]
	} else {	print("multiple module numbers found: binary LPAwb+")
		print(NETWORKS[[data]])
		print(MODULES) }

	#May want to check these for extra configs...hard to code!!
#quantitative LP
	Q_LP_R = as.matrix(read.csv(paste(front2,fq,frow,fL,NETWORKS[data],".csv",sep=""),row.names=1))
	Q_LP_C = as.matrix(read.csv(paste(front2,fq,fcol,fL,NETWORKS[data],".csv",sep=""),row.names=1))

	HAM=which(ModularityScore_LPA_wb_plus[data,] == MAX_QW_QL[data])
	QLPindices = HAM
	HMM = c()
	MODULES = c()
	MODULES[1] = length(unique(Q_LP_R[HAM[1],]))
	Realized=c()
	Realized[1] = RealizedMod(Q,Q_LP_R[HAM[1],],Q_LP_C[HAM[1],])
	QR_QL[data] = Realized[1]
	HAH=1
	if(length(HAM)>1) {
		for( bb in 2:length(HAM)) {
			HMM[bb-1] = NormalisedMutualInformation(c(Q_LP_R[HAM[1],],Q_LP_C[HAM[1],]) , c(Q_LP_R[HAM[bb],],Q_LP_C[HAM[bb],])) #compare first best with rest NMI of 1 indicates same config.
			MODULES[bb] = length(unique(Q_LP_R[HAM[bb],]))
			Realized[bb] = RealizedMod(Q,Q_LP_R[HAM[bb],],Q_LP_C[HAM[bb],])
		}
		
		if(length(unique(Realized))>1) {
			print("multiple Qr values: quantitative LPAwb+")
			print(NETWORKS[[data]])
		}
	HAH=unique(HMM)
	}
	UniqueConfigs_QL[data] = length(HAH)
	if( length(unique(MODULES)) == 1) {
		MODULES_QL[data] = MODULES[1]
	} else {	print("multiple module numbers found: quantitative LPAwb+")
		print(NETWORKS[[data]])
		print(MODULES) }

	#May want to check these for extra configs...hard to code!!
#binary QBM
	B_QB_R = as.matrix(read.csv(paste(front2,fbin,frow,fQ,NETWORKS[data],".csv",sep=""),row.names=1))
	B_QB_C = as.matrix(read.csv(paste(front2,fbin,fcol,fQ,NETWORKS[data],".csv",sep=""),row.names=1))

	HAM=which(BinModularityScore_QuanBiMo[data,] == MAX_QW_BQ[data])
	BQBindices = HAM
	HMM = c()
	MODULES = c()
	MODULES[1] = length(unique(B_QB_R[HAM[1],]))
	Realized=c()
	Realized[1] = RealizedMod(A,B_QB_R[HAM[1],],B_QB_C[HAM[1],])
	QR_BQ[data] = Realized[1]
	HAH=1
	if(length(HAM)>1) {
		for( bb in 2:length(HAM)) {
			HMM[bb-1] = NormalisedMutualInformation(c(B_QB_R[HAM[1],],B_QB_C[HAM[1],]) , c(B_QB_R[HAM[bb],],B_QB_C[HAM[bb],])) #compare first best with rest NMI of 1 indicates same config.
			MODULES[bb] = length(unique(B_QB_R[HAM[bb],]))
			Realized[bb] = RealizedMod(A,B_QB_R[HAM[bb],],B_QB_C[HAM[bb],])
		}
		if(length(unique(Realized))>1) {
			print("multiple Qr values: binary QuanBiMo")
			print(NETWORKS[[data]]) }
	HAH=unique(HMM)
	}
	UniqueConfigs_BQ[data] = length(HAH)
	if( length(unique(MODULES)) == 1) {
		MODULES_BQ[data] = MODULES[1]
	} else {	print("multiple module numbers found: binary QuanBiMo")
		print(NETWORKS[[data]])
		print(MODULES) }

	#May want to check these for extra configs...hard to code!!
#quantitative QBM
	Q_QB_R = as.matrix(read.csv(paste(front2,fq,frow,fQ,NETWORKS[data],".csv",sep=""),row.names=1))
	Q_QB_C = as.matrix(read.csv(paste(front2,fq,fcol,fQ,NETWORKS[data],".csv",sep=""),row.names=1))

	HAM=which(ModularityScore_QuanBiMo[data,] == MAX_QW_QQ[data])
	QQBindices = HAM
	HMM = c()
	MODULES = c()
	MODULES[1] = length(unique(Q_QB_R[HAM[1],]))
	Realized=c()
	Realized[1] = RealizedMod(Q,Q_QB_R[HAM[1],],Q_QB_C[HAM[1],])
	QR_QQ[data] = Realized[1]
	HAH=1
	if(length(HAM)>1) {
		for( bb in 2:length(HAM)) {
			HMM[bb-1] = NormalisedMutualInformation(c(Q_QB_R[HAM[1],],Q_QB_C[HAM[1],]) , c(Q_QB_R[HAM[bb],],Q_QB_C[HAM[bb],])) #compare first best with rest NMI of 1 indicates same config.
			MODULES[bb] = length(unique(Q_QB_R[HAM[bb],]))
			Realized[bb] = RealizedMod(Q,Q_QB_R[HAM[bb],],Q_QB_C[HAM[bb],])
		}
		if(length(unique(Realized))>1) {
			print("multiple Qr values: Quantitative QuanBiMo")
			print(NETWORKS[[data]]) }
	HAH=unique(HMM)
	}
	UniqueConfigs_QQ[data] = length(HAH)
	if( length(unique(MODULES)) == 1) {
		MODULES_QQ[data] = MODULES[1]
	} else {	print("multiple module numbers found: Quantitative QuanBiMo")
		print(NETWORKS[[data]])
		print(MODULES) }

	#May want to check these for extra configs...hard to code!!
#binary ExhaustiveLPAwb+
	B_EX_R = as.matrix(read.csv(paste(front2,fbin,frow,fEL,NETWORKS[data],".csv",sep=""),row.names=1))
	B_EX_C = as.matrix(read.csv(paste(front2,fbin,fcol,fEL,NETWORKS[data],".csv",sep=""),row.names=1))

	HAM=which(BinModularityScore_EXLPAwbplus[data,] == MAX_QW_BE[data])
	BEXindices = HAM
	HMM = c()
	MODULES = c()
	MODULES[1] = length(unique(B_EX_R[HAM[1],]))
	Realized=c()
	Realized[1] = RealizedMod(A,B_EX_R[HAM[1],],B_EX_C[HAM[1],])
	QR_BE[data] = Realized[1]
	HAH=1
	if(length(HAM)>1) {
		for( bb in 2:length(HAM)) {
			HMM[bb-1] = NormalisedMutualInformation(c(B_EX_R[HAM[1],],B_EX_C[HAM[1],]) , c(B_EX_R[HAM[bb],],B_EX_C[HAM[bb],])) #compare first best with rest NMI of 1 indicates same config.
			MODULES[bb] = length(unique(B_EX_R[HAM[bb],]))
			Realized[bb] = RealizedMod(A,B_EX_R[HAM[bb],],B_EX_C[HAM[bb],])
		}
		if(length(unique(Realized))>1) {
			print("multiple Qr values: binary EXLPAwbplus")
			print(NETWORKS[[data]]) }
	HAH=unique(HMM)
	}
	UniqueConfigs_BE[data] = length(HAH)
	if( length(unique(MODULES)) == 1) {
		MODULES_BE[data] = MODULES[1]
	} else {	print("multiple module numbers found: binary EXLPAwbplus")
		print(NETWORKS[[data]])
		print(MODULES) }

	#May want to check these for extra configs...hard to code!!
#quantitative ExhaustiveLPAwb+
	Q_EX_R = as.matrix(read.csv(paste(front2,fq,frow,fEL,NETWORKS[data],".csv",sep=""),row.names=1))
	Q_EX_C = as.matrix(read.csv(paste(front2,fq,fcol,fEL,NETWORKS[data],".csv",sep=""),row.names=1))

	HAM=which(ModularityScore_EXLPAwbplus[data,] == MAX_QW_QE[data])
	QQBindices = HAM
	HMM = c()
	MODULES = c()
	MODULES[1] = length(unique(Q_EX_R[HAM[1],]))
	Realized=c()
	Realized[1] = RealizedMod(Q,Q_EX_R[HAM[1],],Q_EX_C[HAM[1],])
	QR_QE[data] = Realized[1]
	HAH=1
	if(length(HAM)>1) {
		for( bb in 2:length(HAM)) {
			HMM[bb-1] = NormalisedMutualInformation(c(Q_EX_R[HAM[1],],Q_EX_C[HAM[1],]) , c(Q_EX_R[HAM[bb],],Q_EX_C[HAM[bb],])) #compare first best with rest NMI of 1 indicates same config.
			MODULES[bb] = length(unique(Q_EX_R[HAM[bb],]))
			Realized[bb] = RealizedMod(Q,Q_EX_R[HAM[bb],],Q_EX_C[HAM[bb],])
		}
	if(length(unique(Realized))>1) {
		print("multiple Qr values: Quantitative EXLPAwbplus")
		print(NETWORKS[[data]]) }
	HAH=unique(HMM)
	}
	UniqueConfigs_QE[data] = length(HAH)
	if( length(unique(MODULES)) == 1) {
		MODULES_QE[data] = MODULES[1]
	} else {	print("multiple module numbers found: Quantitative EXLPAwbplus")
		print(NETWORKS[[data]])
		print(MODULES) }

	#May want to check these for extra configs...hard to code!!	


	data = data + 1
}



FILE1 = cbind(AVG_time_BL,MIN_QW_BL,MAX_QW_BL,MEAN_QW_BL,MEDIAN_QW_BL,NUMFAIL_BL,NUMMAX_BL,MODULES_BL,UniqueConfigs_BL,QR_BL)

FILE2 = cbind(AVG_time_QL,MIN_QW_QL,MAX_QW_QL,MEAN_QW_QL,MEDIAN_QW_QL,NUMFAIL_QL,NUMMAX_QL,MODULES_QL,UniqueConfigs_QL,QR_QL)

FILE3 = cbind(AVG_time_BQ,MIN_QW_BQ,MAX_QW_BQ,MEAN_QW_BQ,MEDIAN_QW_BQ,NUMFAIL_BQ,NUMMAX_BQ,MODULES_BQ,UniqueConfigs_BQ,QR_BQ)

FILE4 = cbind(AVG_time_QQ,MIN_QW_QQ,MAX_QW_QQ,MEAN_QW_QQ,MEDIAN_QW_QQ,NUMFAIL_QQ,NUMMAX_QQ,MODULES_QQ,UniqueConfigs_QQ,QR_QQ)

FILE5 = cbind(AVG_time_BE,MIN_QW_BE,MAX_QW_BE,MEAN_QW_BE,MEDIAN_QW_BE,NUMFAIL_BE,NUMMAX_BE,MODULES_BE,UniqueConfigs_BE,QR_BE)

FILE6 = cbind(AVG_time_QE,MIN_QW_QE,MAX_QW_QE,MEAN_QW_QE,MEDIAN_QW_QE,NUMFAIL_QE,NUMMAX_QE,MODULES_QE,UniqueConfigs_QE,QR_QE)




write.csv(FILE1,paste(front1,"summary",fBin,fL,".csv",sep=""))
write.csv(FILE2,paste(front1,"summary",fQua,fL,".csv",sep=""))
write.csv(FILE3,paste(front1,"summary",fBin,fQ,".csv",sep=""))
write.csv(FILE4,paste(front1,"summary",fQua,fQ,".csv",sep=""))
write.csv(FILE5,paste(front1,"summary",fBin,fEL,".csv",sep=""))
write.csv(FILE6,paste(front1,"summary",fQua,fEL,".csv",sep=""))



print("completed")



















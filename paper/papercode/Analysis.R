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
f3 = "configurations/"
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

source("../../code/R/convert2moduleWeb.R")

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
#BE - binary DIRTLPAwb+
#QE - quantitative DIRTLPAwb+

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



	minNMI_B_QB_LP=c()
	minNMI_B_QB_EX=c()
	minNMI_B_LP_EX=c()
	maxNMI_B_QB_LP=c()
	maxNMI_B_QB_EX=c()
	maxNMI_B_LP_EX=c()
	minNMI_Q_QB_LP=c()
	minNMI_Q_QB_EX=c()
	minNMI_Q_LP_EX=c()
	maxNMI_Q_QB_LP=c()
	maxNMI_Q_QB_EX=c()
	maxNMI_Q_LP_EX=c()


B_LP=c()
B_QB=c()
B_EX=c()
Q_LP=c()
Q_QB=c()
Q_EX=c()

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
#binary DIRTLPAwb+
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
			print("multiple Qr values: binary DIRTLPAwbplus")
			print(NETWORKS[[data]]) }
	HAH=unique(HMM)
	}
	UniqueConfigs_BE[data] = length(HAH)
	if( length(unique(MODULES)) == 1) {
		MODULES_BE[data] = MODULES[1]
	} else {	print("multiple module numbers found: binary DIRTLPAwbplus")
		print(NETWORKS[[data]])
		print(MODULES) }

	#May want to check these for extra configs...hard to code!!
#quantitative DIRTLPAwb+
	Q_EX_R = as.matrix(read.csv(paste(front2,fq,frow,fEL,NETWORKS[data],".csv",sep=""),row.names=1))
	Q_EX_C = as.matrix(read.csv(paste(front2,fq,fcol,fEL,NETWORKS[data],".csv",sep=""),row.names=1))

	HAM=which(ModularityScore_EXLPAwbplus[data,] == MAX_QW_QE[data])
	QEXindices = HAM
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


	





	
	NMI_B_QB_LP=c()
	NMI_B_QB_EX=c()
	NMI_B_LP_EX=c()
	ix=1
	ix2=1
	for(aa in 1:length(BQBindices)) {
		#Binary NMI QB vs LP
		for(bb in 1:length(BLPindices)) {
			NMI_B_QB_LP[ix]=NormalisedMutualInformation(c(B_QB_R[BQBindices[aa],],B_QB_C[BQBindices[aa],]) , c(B_LP_R[BLPindices[bb],],B_LP_C[BLPindices[bb],]))
			ix=ix+1
		}
		#Binary NMI QB vs EX
		for(bb in 1:length(BEXindices)) {
			NMI_B_QB_EX[ix2]=NormalisedMutualInformation(c(B_QB_R[BQBindices[aa],],B_QB_C[BQBindices[aa],]) , c(B_EX_R[BEXindices[bb],],B_EX_C[BEXindices[bb],]))
			ix2=ix2+1
		}
	}
	#Binary NMI LP vs EX
	ix=1
	for(aa in 1:length(BLPindices)) {
		for(bb in 1:length(BEXindices)) {
			NMI_B_LP_EX[ix]=NormalisedMutualInformation(c(B_LP_R[BLPindices[aa],],B_LP_C[BLPindices[aa],]) , c(B_EX_R[BEXindices[bb],],B_EX_C[BEXindices[bb],]))
			ix=ix+1
		}
	}
	minNMI_B_QB_LP[data] = min(NMI_B_QB_LP)
	minNMI_B_QB_EX[data] = min(NMI_B_QB_EX)
	minNMI_B_LP_EX[data] = min(NMI_B_LP_EX)
	maxNMI_B_QB_LP[data] = max(NMI_B_QB_LP)
	maxNMI_B_QB_EX[data] = max(NMI_B_QB_EX)
	maxNMI_B_LP_EX[data] = max(NMI_B_LP_EX)

#FIND CONFIGURATIONS THAT RESULT IN HIGHEST (most similar) and SMALLEST (most different) NMI

	if(minNMI_B_QB_LP[data]==maxNMI_B_QB_LP[data]) {
		#save to file single config data
		IX = which(NMI_B_QB_LP == min(NMI_B_QB_LP))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1 = floor((IX-1)/length(BLPindices))+1 #BQB config, found using division (BLP is the inside loop)
		IX_2 = ((IX-1) %% length(BLPindices)) +1 #BLP config; the remainder.

		B_LP$Row_labels = B_LP_R[BLPindices[IX_2],]
		B_LP$Col_labels = B_LP_C[BLPindices[IX_2],]
		B_LP$modularity = MAX_QW_BL[data]
		B_QB$Row_labels = B_QB_R[BQBindices[IX_1],]
		B_QB$Col_labels = B_QB_C[BQBindices[IX_1],]
		B_QB$modularity = MAX_QW_BQ[data]

		moduleWebobj = convert2moduleWeb(A,B_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBLP_BQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(A,B_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBLP_BLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	} else {
		#want min config and max config for data
		IXmin = which(NMI_B_QB_LP == min(NMI_B_QB_LP))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1n = floor((IXmin-1)/length(BLPindices))+1 #BQB config, found using division (BLP is the inside loop)
		IX_2n = ((IXmin-1) %% length(BLPindices)) +1 #BLP config; the remainder.
		IXmax = which(NMI_B_QB_LP == max(NMI_B_QB_LP))[[1]]
		IX_1x = floor((IXmax-1)/length(BLPindices))+1 #BQB config, found using division (BLP is the inside loop)
		IX_2x = ((IXmax-1) %% length(BLPindices)) +1 #BLP config; the remainder.

		B_LP$Row_labels = B_LP_R[BLPindices[IX_2x],]
		B_LP$Col_labels = B_LP_C[BLPindices[IX_2x],]
		B_LP$modularity = MAX_QW_BL[data]
		B_QB$Row_labels = B_QB_R[BQBindices[IX_1x],]
		B_QB$Col_labels = B_QB_C[BQBindices[IX_1x],]
		B_QB$modularity = MAX_QW_BQ[data]

		moduleWebobj = convert2moduleWeb(A,B_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBLP_maxBQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(A,B_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBLP_maxBLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()

		B_LP$Row_labels = B_LP_R[BLPindices[IX_2n],]
		B_LP$Col_labels = B_LP_C[BLPindices[IX_2n],]
		B_LP$modularity = MAX_QW_BL[data]
		B_QB$Row_labels = B_QB_R[BQBindices[IX_1n],]
		B_QB$Col_labels = B_QB_C[BQBindices[IX_1n],]
		B_QB$modularity = MAX_QW_BQ[data]

		moduleWebobj = convert2moduleWeb(A,B_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBLP_minBQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(A,B_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBLP_minBLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	}

	if(minNMI_B_QB_EX[data]==maxNMI_B_QB_EX[data]) {
		#save to file single config data
		IX = which(NMI_B_QB_EX == min(NMI_B_QB_EX))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1 = floor((IX-1)/length(BEXindices))+1 #BQB config, found using division (BEX is the inside loop)
		IX_2 = ((IX-1) %% length(BEXindices)) +1 #BEX config; the remainder.

		B_EX$Row_labels = B_EX_R[BEXindices[IX_2],]
		B_EX$Col_labels = B_EX_C[BEXindices[IX_2],]
		B_EX$modularity = MAX_QW_BE[data]
		B_QB$Row_labels = B_QB_R[BQBindices[IX_1],]
		B_QB$Col_labels = B_QB_C[BQBindices[IX_1],]
		B_QB$modularity = MAX_QW_BQ[data]

		moduleWebobj = convert2moduleWeb(A,B_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBEX_BQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(A,B_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBEX_BEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	} else {
		#want min config and max config for data
		IXmin = which(NMI_B_QB_EX == minNMI_B_QB_EX[data])[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1n = floor((IXmin-1)/length(BEXindices))+1 #BQB config, found using division (BEX is the inside loop)
		IX_2n = ((IXmin-1) %% length(BEXindices)) +1 #BEX config; the remainder.
		IXmax = which(NMI_B_QB_EX == maxNMI_B_QB_EX[data])[[1]]
		IX_1x = floor((IXmax-1)/length(BEXindices))+1 #BQB config, found using division (BEX is the inside loop)
		IX_2x = ((IXmax-1) %% length(BEXindices)) +1 #BEX config; the remainder.
	
		B_EX$Row_labels = B_EX_R[BEXindices[IX_2x],]
		B_EX$Col_labels = B_EX_C[BEXindices[IX_2x],]
		B_EX$modularity = MAX_QW_BE[data]
		B_QB$Row_labels = B_QB_R[BQBindices[IX_1x],]
		B_QB$Col_labels = B_QB_C[BQBindices[IX_1x],]
		B_QB$modularity = MAX_QW_BQ[data]

		moduleWebobj = convert2moduleWeb(A,B_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBEX_maxBQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(A,B_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBEX_maxBEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()

		B_EX$Row_labels = B_EX_R[BEXindices[IX_2n],]
		B_EX$Col_labels = B_EX_C[BEXindices[IX_2n],]
		B_EX$modularity = MAX_QW_BE[data]
		B_QB$Row_labels = B_QB_R[BQBindices[IX_1n],]
		B_QB$Col_labels = B_QB_C[BQBindices[IX_1n],]
		B_QB$modularity = MAX_QW_BQ[data]

		moduleWebobj = convert2moduleWeb(A,B_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBEX_minBQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(A,B_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BQBEX_minBEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	}
	
	if(minNMI_B_LP_EX[data]==maxNMI_B_LP_EX[data]) {
		#save to file single config data
		IX = which(NMI_B_LP_EX == min(NMI_B_LP_EX))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1 = floor((IX-1)/length(BEXindices))+1 #BLP config, found using division (BEX is the inside loop)
		IX_2 = ((IX-1) %% length(BEXindices)) +1 #BEX config; the remainder.

		B_EX$Row_labels = B_EX_R[BEXindices[IX_2],]
		B_EX$Col_labels = B_EX_C[BEXindices[IX_2],]
		B_EX$modularity = MAX_QW_BE[data]
		B_LP$Row_labels = B_LP_R[BLPindices[IX_1],]
		B_LP$Col_labels = B_LP_C[BLPindices[IX_1],]
		B_LP$modularity = MAX_QW_BL[data]

		moduleWebobj = convert2moduleWeb(A,B_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BLPEX_BLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(A,B_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BLPEX_BEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	} else {
		#want min config and max config for data
		IXmin = which(NMI_B_LP_EX == min(NMI_B_LP_EX))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1n = floor((IXmin-1)/length(BEXindices))+1 #BLP config, found using division (BEX is the inside loop)
		IX_2n = ((IXmin-1) %% length(BEXindices)) +1 #BEX config; the remainder.
		IXmax = which(NMI_B_LP_EX == max(NMI_B_LP_EX))[[1]]
		IX_1x = floor((IXmax-1)/length(BEXindices))+1 #BLP config, found using division (BEX is the inside loop)
		IX_2x = ((IXmax-1) %% length(BEXindices)) +1 #BEX config; the remainder.

		B_EX$Row_labels = B_EX_R[BEXindices[IX_2x],]
		B_EX$Col_labels = B_EX_C[BEXindices[IX_2x],]
		B_EX$modularity = MAX_QW_BE[data]
		B_LP$Row_labels = B_LP_R[BLPindices[IX_1x],]
		B_LP$Col_labels = B_LP_C[BLPindices[IX_1x],]
		B_LP$modularity = MAX_QW_BL[data]

		moduleWebobj = convert2moduleWeb(A,B_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BLPEX_maxBLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(A,B_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BLPEX_maxBEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()

		B_EX$Row_labels = B_EX_R[BEXindices[IX_2n],]
		B_EX$Col_labels = B_EX_C[BEXindices[IX_2n],]
		B_EX$modularity = MAX_QW_BE[data]
		B_LP$Row_labels = B_LP_R[BLPindices[IX_1n],]
		B_LP$Col_labels = B_LP_C[BLPindices[IX_1n],]
		B_LP$modularity = MAX_QW_BL[data]

		moduleWebobj = convert2moduleWeb(A,B_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BLPEX_minBLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(A,B_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_BLPEX_minBEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	}


	NMI_Q_QB_LP=c()
	NMI_Q_QB_EX=c()
	NMI_Q_LP_EX=c()
	ix=1
	ix2=1
	for(aa in 1:length(QQBindices)) {
		#Quant NMI QB vs LP
		for(bb in 1:length(QLPindices)) {
			NMI_Q_QB_LP[ix]=NormalisedMutualInformation(c(Q_QB_R[QQBindices[aa],],Q_QB_C[QQBindices[aa],]) , c(Q_LP_R[QLPindices[bb],],Q_LP_C[QLPindices[bb],]))
			ix=ix+1
		}
		#Quant NMI QB vs EX
		for(bb in 1:length(QEXindices)) {
			NMI_Q_QB_EX[ix2]=NormalisedMutualInformation(c(Q_QB_R[QQBindices[aa],],Q_QB_C[QQBindices[aa],]) , c(Q_EX_R[QEXindices[bb],],Q_EX_C[QEXindices[bb],]))
			ix2=ix2+1
		}
	}
	#Quant NMI LP vs EX
	ix=1
	for(aa in 1:length(QLPindices)) {
		for(bb in 1:length(QEXindices)) {
			NMI_Q_LP_EX[ix]=NormalisedMutualInformation(c(Q_LP_R[QLPindices[aa],],Q_LP_C[QLPindices[aa],]) , c(Q_EX_R[QEXindices[bb],],Q_EX_C[QEXindices[bb],]))
			ix=ix+1
		}
	}
	minNMI_Q_QB_LP[data] = min(NMI_Q_QB_LP)
	minNMI_Q_QB_EX[data] = min(NMI_Q_QB_EX)
	minNMI_Q_LP_EX[data] = min(NMI_Q_LP_EX)
	maxNMI_Q_QB_LP[data] = max(NMI_Q_QB_LP)
	maxNMI_Q_QB_EX[data] = max(NMI_Q_QB_EX)
	maxNMI_Q_LP_EX[data] = max(NMI_Q_LP_EX)

#FIND CONFIGURATIONS THAT RESULT IN HIGHEST (most similar) and SMALLEST (most different) NMI

	if(minNMI_Q_QB_LP[data]==maxNMI_Q_QB_LP[data]) {
		#save to file single config data
		IX = which(NMI_Q_QB_LP == min(NMI_Q_QB_LP))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1 = floor((IX-1)/length(QLPindices))+1 #QQB config, found using division (QLP is the inside loop)
		IX_2 = ((IX-1) %% length(QLPindices)) +1 #QLP config; the remainder.

		Q_LP$Row_labels = Q_LP_R[QLPindices[IX_2],]
		Q_LP$Col_labels = Q_LP_C[QLPindices[IX_2],]
		Q_LP$modularity = MAX_QW_QL[data]
		Q_QB$Row_labels = Q_QB_R[QQBindices[IX_1],]
		Q_QB$Col_labels = Q_QB_C[QQBindices[IX_1],]
		Q_QB$modularity = MAX_QW_QQ[data]

		moduleWebobj = convert2moduleWeb(Q,Q_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBLP_QQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(Q,Q_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBLP_QLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	} else {
		#want min config and max config for data
		IXmin = which(NMI_Q_QB_LP == min(NMI_Q_QB_LP))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1n = floor((IXmin-1)/length(QLPindices))+1 #QQB config, found using division (QLP is the inside loop)
		IX_2n = ((IXmin-1) %% length(QLPindices))+1 #QLP config; the remainder.
		IXmax = which(NMI_Q_QB_LP == max(NMI_Q_QB_LP))[[1]]
		IX_1x = floor((IXmax-1)/length(QLPindices))+1 #QQB config, found using division (QLP is the inside loop)
		IX_2x = ((IXmax-1) %% length(QLPindices))+1 #QLP config; the remainder.

		Q_LP$Row_labels = Q_LP_R[QLPindices[IX_2x],]
		Q_LP$Col_labels = Q_LP_C[QLPindices[IX_2x],]
		Q_LP$modularity = MAX_QW_QL[data]
		Q_QB$Row_labels = Q_QB_R[QQBindices[IX_1x],]
		Q_QB$Col_labels = Q_QB_C[QQBindices[IX_1x],]
		Q_QB$modularity = MAX_QW_QQ[data]

		moduleWebobj = convert2moduleWeb(Q,Q_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBLP_maxQQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(Q,Q_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBLP_maxQLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()

		Q_LP$Row_labels = Q_LP_R[QLPindices[IX_2n],]
		Q_LP$Col_labels = Q_LP_C[QLPindices[IX_2n],]
		Q_LP$modularity = MAX_QW_QL[data]
		Q_QB$Row_labels = Q_QB_R[QQBindices[IX_1n],]
		Q_QB$Col_labels = Q_QB_C[QQBindices[IX_1n],]
		Q_QB$modularity = MAX_QW_QQ[data]

		moduleWebobj = convert2moduleWeb(Q,B_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBLP_minQQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(Q,B_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBLP_minQLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	}

	if(minNMI_Q_QB_EX[data]==maxNMI_Q_QB_EX[data]) {
		#save to file single config data
		IX = which(NMI_Q_QB_EX == min(NMI_Q_QB_EX))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1 = floor((IX-1)/length(QEXindices))+1 #QQB config, found using division (QEX is the inside loop)
		IX_2 = ((IX-1) %% length(QEXindices))+1 #QEX config; the remainder.

		Q_EX$Row_labels = Q_EX_R[QEXindices[IX_2],]
		Q_EX$Col_labels = Q_EX_C[QEXindices[IX_2],]
		Q_EX$modularity = MAX_QW_QE[data]
		Q_QB$Row_labels = Q_QB_R[QQBindices[IX_1],]
		Q_QB$Col_labels = Q_QB_C[QQBindices[IX_1],]
		Q_QB$modularity = MAX_QW_QQ[data]

		moduleWebobj = convert2moduleWeb(Q,Q_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBEX_QQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(Q,Q_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBEX_QEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	} else {
		#want min config and max config for data
		IXmin = which(NMI_Q_QB_EX == min(NMI_Q_QB_EX))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1n = floor((IXmin-1)/length(QEXindices))+1 #QQB config, found using division (QEX is the inside loop)
		IX_2n = ((IXmin-1) %% length(QEXindices))+1 #QEX config; the remainder.
		IXmax = which(NMI_Q_QB_EX == max(NMI_Q_QB_EX))[[1]]
		IX_1x = floor((IXmax-1)/length(QEXindices))+1 #QQB config, found using division (QEX is the inside loop)
		IX_2x = ((IXmax-1) %% length(QEXindices))+1 #QEX config; the remainder.

		Q_EX$Row_labels = Q_EX_R[QEXindices[IX_2x],]
		Q_EX$Col_labels = Q_EX_C[QEXindices[IX_2x],]
		Q_QE$modularity = MAX_QW_QE[data]
		Q_QB$Row_labels = Q_QB_R[QQBindices[IX_1x],]
		Q_QB$Col_labels = Q_QB_C[QQBindices[IX_1x],]
		Q_QB$modularity = MAX_QW_QQ[data]

		moduleWebobj = convert2moduleWeb(Q,Q_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBEX_maxQQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(Q,Q_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBEX_maxQEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()

		Q_EX$Row_labels = Q_EX_R[QEXindices[IX_2n],]
		Q_EX$Col_labels = Q_EX_C[QEXindices[IX_2n],]
		Q_EX$modularity = MAX_QW_QE[data]
		Q_QB$Row_labels = Q_QB_R[QQBindices[IX_1n],]
		Q_QB$Col_labels = Q_QB_C[QQBindices[IX_1n],]
		Q_QB$modularity = MAX_QW_QQ[data]

		moduleWebobj = convert2moduleWeb(Q,Q_QB)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBEX_minQQB_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(Q,Q_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QQBEX_minQEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	}
	
	if(minNMI_Q_LP_EX[data]==maxNMI_Q_LP_EX[data]) {
		#save to file single config data
		IX = which(NMI_Q_LP_EX == min(NMI_Q_LP_EX))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1 = floor((IX-1)/length(QEXindices))+1 #QLP config, found using division (QEX is the inside loop)
		IX_2 = ((IX-1) %% length(QEXindices))+1 #QEX config; the remainder.

		Q_EX$Row_labels = Q_EX_R[QEXindices[IX_2],]
		Q_EX$Col_labels = Q_EX_C[QEXindices[IX_2],]
		Q_EX$modularity = MAX_QW_QE[data]
		Q_LP$Row_labels = Q_LP_R[QLPindices[IX_1],]
		Q_LP$Col_labels = Q_LP_C[QLPindices[IX_1],]
		Q_LP$modularity = MAX_QW_QL[data]

		moduleWebobj = convert2moduleWeb(Q,Q_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QLPEX_QLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(Q,Q_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QLPEX_QEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	} else {
		#want min config and max config for data
		IXmin = which(NMI_Q_LP_EX == min(NMI_Q_LP_EX))[[1]]
		#convert from for loop index to indices of each algorithm
		IX_1n = floor((IXmin-1)/length(QEXindices))+1 #QLP config, found using division (QEX is the inside loop)
		IX_2n = ((IXmin-1) %% length(QEXindices))+1 #QEX config; the remainder.
		IXmax = which(NMI_Q_LP_EX == max(NMI_Q_LP_EX))[[1]]
		IX_1x = floor((IXmax-1)/length(QEXindices))+1 #QLP config, found using division (QEX is the inside loop)
		IX_2x = ((IXmax-1) %% length(QEXindices))+1 #QEX config; the remainder.

		Q_EX$Row_labels = Q_EX_R[QEXindices[IX_2x],]
		Q_EX$Col_labels = Q_EX_C[QEXindices[IX_2x],]
		Q_EX$modularity = MAX_QW_QE[data]
		Q_LP$Row_labels = Q_LP_R[QLPindices[IX_1x],]
		Q_LP$Col_labels = Q_LP_C[QLPindices[IX_1x],]
		Q_LP$modularity = MAX_QW_QL[data]

		moduleWebobj = convert2moduleWeb(Q,Q_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QLPEX_maxQLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(Q,Q_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QLPEX_maxQEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()

		Q_EX$Row_labels = Q_EX_R[QEXindices[IX_2n],]
		Q_EX$Col_labels = Q_EX_C[QEXindices[IX_2n],]
		Q_EX$modularity = MAX_QW_QE[data]
		Q_LP$Row_labels = Q_LP_R[QLPindices[IX_1n],]
		Q_LP$Col_labels = Q_LP_C[QLPindices[IX_1n],]
		Q_LP$modularity = MAX_QW_QL[data]

		moduleWebobj = convert2moduleWeb(Q,Q_LP)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QLPEX_minQLP_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
		
		moduleWebobj = convert2moduleWeb(Q,Q_EX)
		FILENAME = paste(front2,f3,NETWORKS[data],"_QLPEX_minQEX_modulelisting.csv",sep="")
		sink(FILENAME)
		printoutModuleInformation(moduleWebobj)
		sink()
	}



	data = data + 1
}



FILE1 = cbind(AVG_time_BL,MIN_QW_BL,MAX_QW_BL,MEAN_QW_BL,MEDIAN_QW_BL,NUMFAIL_BL,NUMMAX_BL,MODULES_BL,UniqueConfigs_BL,QR_BL)

FILE2 = cbind(AVG_time_QL,MIN_QW_QL,MAX_QW_QL,MEAN_QW_QL,MEDIAN_QW_QL,NUMFAIL_QL,NUMMAX_QL,MODULES_QL,UniqueConfigs_QL,QR_QL)

FILE3 = cbind(AVG_time_BQ,MIN_QW_BQ,MAX_QW_BQ,MEAN_QW_BQ,MEDIAN_QW_BQ,NUMFAIL_BQ,NUMMAX_BQ,MODULES_BQ,UniqueConfigs_BQ,QR_BQ)

FILE4 = cbind(AVG_time_QQ,MIN_QW_QQ,MAX_QW_QQ,MEAN_QW_QQ,MEDIAN_QW_QQ,NUMFAIL_QQ,NUMMAX_QQ,MODULES_QQ,UniqueConfigs_QQ,QR_QQ)

FILE5 = cbind(AVG_time_BE,MIN_QW_BE,MAX_QW_BE,MEAN_QW_BE,MEDIAN_QW_BE,NUMFAIL_BE,NUMMAX_BE,MODULES_BE,UniqueConfigs_BE,QR_BE)

FILE6 = cbind(AVG_time_QE,MIN_QW_QE,MAX_QW_QE,MEAN_QW_QE,MEDIAN_QW_QE,NUMFAIL_QE,NUMMAX_QE,MODULES_QE,UniqueConfigs_QE,QR_QE)

FILE7 = cbind(minNMI_B_QB_LP,maxNMI_B_QB_LP,minNMI_B_QB_EX,maxNMI_B_QB_EX,minNMI_B_LP_EX,maxNMI_B_LP_EX)
FILE8 = cbind(minNMI_Q_QB_LP,maxNMI_Q_QB_LP,minNMI_Q_QB_EX,maxNMI_Q_QB_EX,minNMI_Q_LP_EX,maxNMI_Q_LP_EX)



write.csv(FILE1,paste(front1,"summary",fBin,fL,".csv",sep=""))
write.csv(FILE2,paste(front1,"summary",fQua,fL,".csv",sep=""))
write.csv(FILE3,paste(front1,"summary",fBin,fQ,".csv",sep=""))
write.csv(FILE4,paste(front1,"summary",fQua,fQ,".csv",sep=""))
write.csv(FILE5,paste(front1,"summary",fBin,fEL,".csv",sep=""))
write.csv(FILE6,paste(front1,"summary",fQua,fEL,".csv",sep=""))

write.csv(FILE7,paste(front1,"summary","binaryNMIdiffs.csv",sep=""))
write.csv(FILE8,paste(front1,"summary","quantNMIdiffs.csv",sep=""))



print("completed")



















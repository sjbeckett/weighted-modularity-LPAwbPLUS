

## Function to find configuration with best modularity score
FINDCONFIG <- function(NETWORK_TO_DO, Quant) {

FRONT = "output/summary/"

filestart = "Qua"
front1 = "Q"

if(Quant==0) {
	filestart = "Bin"
	front1 = "bin"
	}


A = as.matrix(read.csv(paste(FRONT,filestart,"_Modscore_QuanBiMo.csv",sep=""), row.names=1))
B = as.matrix(read.csv(paste(FRONT,filestart,"_Modscore_LPA_wb_plus.csv",sep=""), row.names=1))
C = as.matrix(read.csv(paste(FRONT,filestart,"_Modscore_EXLPA_wb_plus.csv",sep=""), row.names=1))	

NETNAME = rownames(A)[NETWORK_TO_DO]

print(paste("Network is ",NETNAME))

Decision=1

AA = max(A[NETWORK_TO_DO,],na.rm=TRUE)
BB = max(B[NETWORK_TO_DO,],na.rn=TRUE)
CC = max(C[NETWORK_TO_DO,],na.rm=TRUE)

maxes = c(AA,BB,CC)
Maximum = max(maxes)

print(paste("with Modularity of:",Maximum))

Decision = which(maxes == Maximum)[1]  # find which alg. found the best modularity score - if multiple solutions then use first found.


if(Decision == 1) {
	RUN = which(A[NETWORK_TO_DO,] == Maximum)[1]
	front2 = "QBM"
} else if(Decision ==2) {
	RUN = which(B[NETWORK_TO_DO,] == Maximum)[1]
	front2 = "LPAwb+"
} else {
	RUN = which(C[NETWORK_TO_DO,] == Maximum)[1]
	front2 = "EXLPAwb+"
}


FRONT="output/"

ROWS = as.matrix(read.csv(paste(FRONT,front1,"rows",front2,NETNAME,".csv",sep=""),row.names=1))[RUN,]
COLS = as.matrix(read.csv(paste(FRONT,front1,"cols",front2,NETNAME,".csv",sep=""),row.names=1))[RUN,]

return(list(ROWS,COLS,NETNAME))
}




## Function to plot out modules onto an incidence  matrix
PLOTMOD <- function(MATRIX,red,blue) {
	
	MODS = unique(red)
	NUMMOD = length(MODS)
	
	Rind = c()
	Cind = c()
	ModlenR = c()
	ModlenC = c()
	
	for (aa in 1:NUMMOD) {
		Rows = which(red == MODS[aa])
		Cols = which(blue == MODS[aa])
		Rind = c(Rind,Rows)
		Cind = c(Cind,Cols)
		ModlenR = c(ModlenR,length(Rows))
		ModlenC = c(ModlenC,length(Cols))
	}


	
	visweb(MATRIX[Rind,Cind],type = 'none',text = 'interaction' ,textsize =2, textcol = "honeydew4", labsize=1.8)
	R=nrow(MATRIX)
	C=ncol(MATRIX)	
	
	xl = 0
	yd = R-ModlenR[1]
	xr = ModlenC[1]
	yu = R

	for(aa in 1:NUMMOD) {
		rect(xl,yd,xr,yu,lwd=6,bor='red')	
		if(aa != NUMMOD) {
			xl=xl+ModlenC[aa]
			yd = yd - ModlenR[aa+1]
			xr = xr + ModlenC[aa+1]
			yu = yu - ModlenR[aa]
		}	
	}

}




## Function to find best configuration of a specific network and plot the modular structure for binary and quantitative forms
PLOTTINGNETWORKMODULES <- function(NETWORK_TO_DO) {


DATALIST = list(Safariland,barrett1987,bezerra2009,elberling1999,inouye1988,junker2013,kato1990,kevan1970,memmott1999,mosquin1967,motten1982,olesen2002aigrettes,olesen2002flores,ollerton2003,schemske1978,small1976,vazarr,vazcer,vazllao,vazmasc,vazmasnc,vazquec,vazquenc)

network = as.matrix(DATALIST[[NETWORK_TO_DO]])
network = network[rowSums(network)>0,colSums(network)>0] # remove empty rows and columns from the incidence matrix (as this was done when running algorithms - empty rows/cols have no module)


BinaryConfig = FINDCONFIG(NETWORK_TO_DO,0)
QuantitativeConfig = FINDCONFIG(NETWORK_TO_DO,1)



PLOTMOD(network,BinaryConfig[[1]],BinaryConfig[[2]])
dev.copy2eps(file=paste(BinaryConfig[[3]],"_mod_bin.eps",sep=""))
dev.off()

PLOTMOD(network,QuantitativeConfig[[1]],QuantitativeConfig[[2]])
dev.copy2eps(file=paste(QuantitativeConfig[[3]],"_mod_qua.eps",sep=""))
dev.off()



}


# R - run each network 100 times through each algorithm. Want to know:
#	a) modularity score
#	b) time to find
#	c) group identifications?


#Initialise

library(bipartite)
source("../../code/R/LPA_wb_plus.R")


#NB: some of the larger datasets take a very long time to run e.g. kato1990 - consider removing these if performing a reanalysis.
DATALIST = list(Safariland,barrett1987,bezerra2009,elberling1999,inouye1988,junker2013,kato1990,kevan1970,memmott1999,mosquin1967,motten1982,olesen2002aigrettes,olesen2002flores,ollerton2003,schemske1978,small1976,vazarr,vazcer,vazllao,vazmasc,vazmasnc,vazquec,vazquenc)

Networks = c("Safariland","barrett1987","bezerra2009","elberling1999","inouye1988","junker2013","kato1990","kevan1970","memmott1999","mosquin1967","motten1982","olesen2002aigrettes","olesen2002flores","ollerton2003","schemske1978","small1976","vazarr","vazcer","vazllao","vazmasc","vazmasnc","vazquec","vazquenc")


NumNetworks <- length(DATALIST)
Repititions <- 100

#Initialise data output objects
InitOut <- matrix(-99,NumNetworks,Repititions)
rownames(InitOut) = Networks


ModularityScore_LPA_wb_plus <- InitOut
ModularityScore_EXLPA_wb_plus <- InitOut
ModularityScore_QuanBiMo <- InitOut

TimeData_LPA_wb_plus <- InitOut
TimeData_EXLPA_wb_plus <- InitOut
TimeData_QuanBiMo <- InitOut

BinModularityScore_LPA_wb_plus <- InitOut
BinModularityScore_EXLPA_wb_plus <- InitOut
BinModularityScore_QuanBiMo <- InitOut

BinTimeData_LPA_wb_plus <- InitOut
BinTimeData_EXLPA_wb_plus <- InitOut
BinTimeData_QuanBiMo <- InitOut



dataset <- 1
Time = c(3.2,3.1,1.7) # initialise with some numbers
Mod = Time

for(filename in DATALIST) {
	Q <- as.matrix(filename) #network
	
	#convert to binary network
	A = 1*(Q>0)

	#Remove empty Rows/Columns from analysis
	if(length(which(rowSums(A)==0))>0) {
		Q = Q[-which(rowSums(A)==0),]
		A = A[-which(rowSums(A)==0),] }
	if(length(which(colSums(A)==0))>0) {
		Q = Q[,-which(colSums(A)==0)]
		A = A[,-which(colSums(A)==0)] }


	# Now Q is the quantitative version of the network and A is the binary version.

	dims <- dim(A)

	## Run for Quantitative networks

	#Objects to store row/column module assignments
	rLP <- matrix(0,Repititions,dims[1])
	rEX <- rLP
	rQM <- rLP
	cLP <- matrix(0,Repititions,dims[2])
	cEX <- cLP
	cQM <- cLP

	B <- Networks[dataset]
	rLPname <- paste('output/QrowsLPAwb+',B,'.csv',sep="")
	rEXname <- paste('output/QrowsEXLPAwb+',B,'.csv',sep="")
	rQMname <- paste('output/QrowsQBM',B,'.csv',sep="")
	cLPname <- paste('output/QcolsLPAwb+',B,'.csv',sep="")
	cEXname <- paste('output/QcolsEXLPAwb+',B,'.csv',sep="")
	cQMname <- paste('output/QcolsQBM',B,'.csv',sep="")

	for( aa in 1:Repititions) {
		## LPAwb+
		Time <- system.time({
					 Mod <- tryCatch( { LPA_wb_plus(Q) }, error = function(cond) {return(NA) } )
				})
		if (sum(is.na(Mod))>0) { #If modularity partitioning failed
			ModularityScore_LPA_wb_plus[dataset,aa] <- NA
		}
		else { #If modularity partitioning successful
			ModularityScore_LPA_wb_plus[dataset,aa] <- Mod[[3]]
			rLP[aa,] <- Mod[[1]]
			cLP[aa,] <- Mod[[2]]	
		}

		ModularityScore_LPA_wb_plus[dataset,aa] <- Mod[[3]]
		TimeData_LPA_wb_plus[dataset,aa] <- Time[[3]]

		## DIRTLPAwb+
		Time <- system.time({
					 Mod <- tryCatch( { DIRT_LPA_wb_plus(Q) }, error = function(cond) {return(NA) } )
				})
		if (sum(is.na(Mod))>0) { #If modularity partitioning failed
			ModularityScore_EXLPA_wb_plus[dataset,aa] <- NA
		}
		else { #If modularity partitioning successful
			ModularityScore_EXLPA_wb_plus[dataset,aa] <- Mod[[3]]
			rEX[aa,] <- Mod[[1]]
			cEX[aa,] <- Mod[[2]]	
		}

		ModularityScore_EXLPA_wb_plus[dataset,aa] <- Mod[[3]]
		TimeData_EXLPA_wb_plus[dataset,aa] <- Time[[3]]

		## QuanBiMo appears to fail when a modularity score < than the tolerance value is found. In this case to avoid error - use a tryCatch loop.
		Time <- system.time({
					 Mod <- tryCatch( { computeModules(Q) }, error = function(cond) { return(NA) })
				})


		if (isS4(Mod)) { #If modularity partitioning successful
			ModularityScore_QuanBiMo[dataset,aa] <- Mod@likelihood
			
			#findnodelabels
			JAMMY = Mod@modules
			BLAM = JAMMY[2:nrow(JAMMY),3:ncol(JAMMY)]
			nodes = ncol(BLAM)
			#convert module matrix to red and blue labels
			rownodes=c()
			colnodes=c()
			for(bb in 1:nodes) {
				TRAM = BLAM[,bb]
				if (bb <= nrow(Q))
					rownodes[bb] = which(TRAM>0)
				else
					colnodes[bb-nrow(Q)] =which(TRAM>0)
				}

			rQM[aa,] <- rownodes
			cQM[aa,] <- colnodes
		}
		else { #If modularity partitioning failed
			ModularityScore_QuanBiMo[dataset,aa] <- NA
		}
		
		TimeData_QuanBiMo[dataset,aa] <- Time[[3]]


		
	}

	write.csv(rLP,rLPname)
	write.csv(cLP,cLPname)
	write.csv(rEX,rEXname)
	write.csv(cEX,cEXname)
	write.csv(rQM,rQMname)
	write.csv(cQM,cQMname)
	

	## Repeat for binary network


	#Objects to store row/column module assignments
	rLP <- matrix(0,Repititions,dims[1])
	rEX <- rLP
	rQM <- rLP
	cLP <- matrix(0,Repititions,dims[2])
	cEX <- cLP
	cQM <- cLP

	B <- Networks[dataset]
	rLPname <- paste('output/binrowsLPAwb+',B,'.csv',sep="")
	rEXname <- paste('output/binrowsEXLPAwb+',B,'.csv',sep="")
	rQMname <- paste('output/binrowsQBM',B,'.csv',sep="")
	cLPname <- paste('output/bincolsLPAwb+',B,'.csv',sep="")
	cEXname <- paste('output/bincolsEXLPAwb+',B,'.csv',sep="")
	cQMname <- paste('output/bincolsQBM',B,'.csv',sep="")

	for( aa in 1:Repititions) {
		## LPAwb+
		Time <- system.time({
					 Mod <- tryCatch( { LPA_wb_plus(A) }, error = function(cond) {return(NA) } )
				})
		if (sum(is.na(Mod))>0) { #If modularity partitioning failed
			BinModularityScore_LPA_wb_plus[dataset,aa] <- NA
		}
		else { #If modularity partitioning successful
			BinModularityScore_LPA_wb_plus[dataset,aa] <- Mod[[3]]
			rLP[aa,] <- Mod[[1]]
			cLP[aa,] <- Mod[[2]]	
		}

		BinModularityScore_LPA_wb_plus[dataset,aa] <- Mod[[3]]
		BinTimeData_LPA_wb_plus[dataset,aa] <- Time[[3]]

		## DIRTLPAwb+
		Time <- system.time({
					 Mod <- tryCatch( { DIRT_LPA_wb_plus(A) }, error = function(cond) {return(NA) } )
				})
		if (sum(is.na(Mod))>0) { #If modularity partitioning failed
			BinModularityScore_EXLPA_wb_plus[dataset,aa] <- NA
		}
		else { #If modularity partitioning successful
			BinModularityScore_EXLPA_wb_plus[dataset,aa] <- Mod[[3]]
			rEX[aa,] <- Mod[[1]]
			cEX[aa,] <- Mod[[2]]	
		}

		BinModularityScore_EXLPA_wb_plus[dataset,aa] <- Mod[[3]]
		BinTimeData_EXLPA_wb_plus[dataset,aa] <- Time[[3]]

		## QuanBiMo appears to fail when a modularity score < than the tolerance value is found. In this case to avoid error - use a tryCatch loop.
		Time <- system.time({
					 Mod <- tryCatch( { computeModules(A) }, error = function(cond) { return(NA) })
				})


		if (isS4(Mod)) { #If modularity partitioning successful
			BinModularityScore_QuanBiMo[dataset,aa] <- Mod@likelihood
			
			#findnodelabels
			JAMMY = Mod@modules
			BLAM = JAMMY[2:nrow(JAMMY),3:ncol(JAMMY)]
			nodes = ncol(BLAM)
			#convert module matrix to red and blue labels
			rownodes=c()
			colnodes=c()
			for(bb in 1:nodes) {
				TRAM = BLAM[,bb]
				if (bb <= nrow(A))
					rownodes[bb] = which(TRAM>0)
				else
					colnodes[bb-nrow(A)] =which(TRAM>0)
				}

			rQM[aa,] <- rownodes
			cQM[aa,] <- colnodes
		}
		else { #If modularity partitioning failed
			BinModularityScore_QuanBiMo[dataset,aa] <- NA
		}
		
		BinTimeData_QuanBiMo[dataset,aa] <- Time[[3]]

		
	}

	write.csv(rLP,rLPname)
	write.csv(cLP,cLPname)
	write.csv(rEX,rEXname)
	write.csv(cEX,cEXname)
	write.csv(rQM,rQMname)
	write.csv(cQM,cQMname)
	
dataset <- dataset + 1
}



#save.image("RcomparisonQuanBiMowithLPAwbplus.RData")
write.csv(ModularityScore_LPA_wb_plus,"output/summary/Qua_Modscore_LPA_wb_plus.csv")
write.csv(TimeData_LPA_wb_plus,"output/summary/Qua_Times_LPA_wb_plus.csv")
write.csv(ModularityScore_EXLPA_wb_plus,"output/summary/Qua_Modscore_EXLPA_wb_plus.csv")
write.csv(TimeData_EXLPA_wb_plus,"output/summary/Qua_Times_EXLPA_wb_plus.csv")
write.csv(ModularityScore_QuanBiMo,"output/summary/Qua_Modscore_QuanBiMo.csv")
write.csv(TimeData_QuanBiMo,"output/summary/Qua_Times_QuanBiMo.csv")

write.csv(BinModularityScore_LPA_wb_plus,"output/summary/Bin_Modscore_LPA_wb_plus.csv")
write.csv(BinTimeData_LPA_wb_plus,"output/summary/Bin_Times_LPA_wb_plus.csv")
write.csv(BinModularityScore_EXLPA_wb_plus,"output/summary/Bin_Modscore_EXLPA_wb_plus.csv")
write.csv(BinTimeData_EXLPA_wb_plus,"output/summary/Bin_Times_EXLPA_wb_plus.csv")
write.csv(BinModularityScore_QuanBiMo,"output/summary/Bin_Modscore_QuanBiMo.csv")
write.csv(BinTimeData_QuanBiMo,"output/summary/Bin_Times_QuanBiMo.csv")







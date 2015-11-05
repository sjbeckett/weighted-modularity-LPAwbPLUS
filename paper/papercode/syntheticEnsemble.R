	
	source("syntheticEnsemblePreparation.R") #load required functions

	#Generate ensemble of matrices to test LPAwb+ and DIRTLPAwb+ against.

	Nrow = 30
	Ncol = 50
	LowMod = 2
	HighMod = 10
	LowFill = 0.5
	HighFill = 2.5
	
	Replicates=10
	
	noiseLevels = c(0,0.01,0.25,0.5)
	noiseReplicates = 5

	#storage
	LowModLowFill = c()
	LowModHighFill = c()
	HighModLowFill = c()
	HighModHighFill = c()

	LMLF=c()
	LMHF=c()
	HMLF=c()
	HMHF=c()

	DETECT_LPA_LMLF = c()
	DETECT_EXLPA_LMLF = c()
	DETECT_LPA_LMHF = c()
	DETECT_EXLPA_LMHF = c()
	DETECT_LPA_HMLF = c()
	DETECT_EXLPA_HMLF = c()
	DETECT_LPA_HMHF = c()
	DETECT_EXLPA_HMHF = c()

	DET_TO_PERF_MODULES_LPA_LMLF = c()
	DET_TO_PERF_MODULES_EXLPA_LMLF = c()
	DET_TO_PERF_MODULES_LPA_LMHF = c()
	DET_TO_PERF_MODULES_EXLPA_LMHF = c()
	DET_TO_PERF_MODULES_LPA_HMLF = c()
	DET_TO_PERF_MODULES_EXLPA_HMLF = c()
	DET_TO_PERF_MODULES_LPA_HMHF = c()
	DET_TO_PERF_MODULES_EXLPA_HMHF = c()

	DET_TO_PERF_MODULARITY_LPA_LMLF = c()
	DET_TO_PERF_MODULARITY_EXLPA_LMLF = c()
	DET_TO_PERF_MODULARITY_LPA_LMHF = c()
	DET_TO_PERF_MODULARITY_EXLPA_LMHF = c()
	DET_TO_PERF_MODULARITY_LPA_HMLF = c()
	DET_TO_PERF_MODULARITY_EXLPA_HMLF = c()
	DET_TO_PERF_MODULARITY_LPA_HMHF = c()
	DET_TO_PERF_MODULARITY_EXLPA_HMHF = c()
		
	DET_TO_PERF_NMI_LPA_LMLF = c()
	DET_TO_PERF_NMI_EXLPA_LMLF = c()
	DET_TO_PERF_NMI_LPA_LMHF = c()
	DET_TO_PERF_NMI_EXLPA_LMHF = c()
	DET_TO_PERF_NMI_LPA_HMLF = c()
	DET_TO_PERF_NMI_EXLPA_HMLF = c()
	DET_TO_PERF_NMI_LPA_HMHF = c()
	DET_TO_PERF_NMI_EXLPA_HMHF = c()



	#Starting matrices
	for(aa in 1:Replicates) {
		LowModLowFill[[aa]] = MakeModule(Nrow,Ncol,LowMod,LowFill)
		LowModHighFill[[aa]] = MakeModule(Nrow,Ncol,LowMod,HighFill)
		HighModLowFill[[aa]] = MakeModule(Nrow,Ncol,HighMod,LowFill)
		HighModHighFill[[aa]] = MakeModule(Nrow,Ncol,HighMod,HighFill)
		print(aa) #make sure stuff is happening
	}

	



	#Show replicates with noise
	index = 1

	for(aa in 1:Replicates){
		for(bb in 1:length(noiseLevels)) {
			for (cc in 1:noiseReplicates) {
				#Make networks
				LMLF[[index]] = RewireConnections(LowModLowFill[[aa]]$matrix,noiseLevels[bb])
				LMHF[[index]] = RewireConnections(LowModHighFill[[aa]]$matrix,noiseLevels[bb])
				HMLF[[index]] = RewireConnections(HighModLowFill[[aa]]$matrix,noiseLevels[bb])
				HMHF[[index]] = RewireConnections(HighModHighFill[[aa]]$matrix,noiseLevels[bb])
					
				##Measure and evaluate networks
			
				#LPA and EXLPA
				DETECT_LPA_LMLF[[index]] = LPA_wb_plus(LMLF[[index]])
				DETECT_EXLPA_LMLF[[index]] = DIRT_LPA_wb_plus(LMLF[[index]])

				DETECT_LPA_LMHF[[index]] = LPA_wb_plus(LMHF[[index]])
				DETECT_EXLPA_LMHF[[index]] = DIRT_LPA_wb_plus(LMHF[[index]])
		
				DETECT_LPA_HMLF[[index]] = LPA_wb_plus(HMLF[[index]])
				DETECT_EXLPA_HMLF[[index]] = DIRT_LPA_wb_plus(HMLF[[index]])
	
				DETECT_LPA_HMHF[[index]] = LPA_wb_plus(HMHF[[index]])
				DETECT_EXLPA_HMHF[[index]] = DIRT_LPA_wb_plus(HMHF[[index]])	
				
				#Compare Number of Modules
				DET_TO_PERF_MODULES_LPA_LMLF[[index]] = COMPAREMODS(LowModLowFill[[aa]],DETECT_LPA_LMLF[[index]])
				DET_TO_PERF_MODULES_EXLPA_LMLF[[index]] = COMPAREMODS(LowModLowFill[[aa]],DETECT_EXLPA_LMLF[[index]])

				DET_TO_PERF_MODULES_LPA_LMHF[[index]] = COMPAREMODS(LowModHighFill[[aa]],DETECT_LPA_LMHF[[index]])
				DET_TO_PERF_MODULES_EXLPA_LMHF[[index]] = COMPAREMODS(LowModHighFill[[aa]],DETECT_EXLPA_LMHF[[index]])

				DET_TO_PERF_MODULES_LPA_HMLF[[index]] = COMPAREMODS(HighModLowFill[[aa]],DETECT_LPA_HMLF[[index]])
				DET_TO_PERF_MODULES_EXLPA_HMLF[[index]] = COMPAREMODS(HighModLowFill[[aa]],DETECT_EXLPA_HMLF[[index]])

				DET_TO_PERF_MODULES_LPA_HMHF[[index]] = COMPAREMODS(HighModHighFill[[aa]],DETECT_LPA_HMHF[[index]])
				DET_TO_PERF_MODULES_EXLPA_HMHF[[index]] = COMPAREMODS(HighModHighFill[[aa]],DETECT_EXLPA_HMHF[[index]])
				
				#Compare Modularity
				DET_TO_PERF_MODULARITY_LPA_LMLF[[index]] = DETECT_LPA_LMLF[[index]]$modularity/PERFECT_MOD(LowModLowFill[[aa]])
				DET_TO_PERF_MODULARITY_EXLPA_LMLF[[index]] = DETECT_EXLPA_LMLF[[index]]$modularity/PERFECT_MOD(LowModLowFill[[aa]])

				DET_TO_PERF_MODULARITY_LPA_LMHF[[index]] = DETECT_LPA_LMHF[[index]]$modularity/PERFECT_MOD(LowModHighFill[[aa]])
				DET_TO_PERF_MODULARITY_EXLPA_LMHF[[index]] = DETECT_EXLPA_LMHF[[index]]$modularity/PERFECT_MOD(LowModHighFill[[aa]])

				DET_TO_PERF_MODULARITY_LPA_HMLF[[index]] = DETECT_LPA_HMLF[[index]]$modularity/PERFECT_MOD(HighModLowFill[[aa]])
				DET_TO_PERF_MODULARITY_EXLPA_HMLF[[index]] = DETECT_EXLPA_HMLF[[index]]$modularity/PERFECT_MOD(HighModLowFill[[aa]])

				DET_TO_PERF_MODULARITY_LPA_HMHF[[index]] = DETECT_LPA_HMHF[[index]]$modularity/PERFECT_MOD(HighModHighFill[[aa]])
				DET_TO_PERF_MODULARITY_EXLPA_HMHF[[index]] = DETECT_EXLPA_HMHF[[index]]$modularity/PERFECT_MOD(HighModHighFill[[aa]])

				#Compare NMI
				DET_TO_PERF_NMI_LPA_LMLF[[index]] = COMPARENMI(LowModLowFill[[aa]],DETECT_LPA_LMLF[[index]])
				DET_TO_PERF_NMI_EXLPA_LMLF[[index]] = COMPARENMI(LowModLowFill[[aa]],DETECT_EXLPA_LMLF[[index]])

				DET_TO_PERF_NMI_LPA_LMHF[[index]] = COMPARENMI(LowModHighFill[[aa]],DETECT_LPA_LMHF[[index]])
				DET_TO_PERF_NMI_EXLPA_LMHF[[index]] = COMPARENMI(LowModHighFill[[aa]],DETECT_EXLPA_LMHF[[index]])

				DET_TO_PERF_NMI_LPA_HMLF[[index]] = COMPARENMI(HighModLowFill[[aa]],DETECT_LPA_HMLF[[index]])
				DET_TO_PERF_NMI_EXLPA_HMLF[[index]] = COMPARENMI(HighModLowFill[[aa]],DETECT_EXLPA_HMLF[[index]])

				DET_TO_PERF_NMI_LPA_HMHF[[index]] = COMPARENMI(HighModHighFill[[aa]],DETECT_LPA_HMHF[[index]])
				DET_TO_PERF_NMI_EXLPA_HMHF[[index]] = COMPARENMI(HighModHighFill[[aa]],DETECT_EXLPA_HMHF[[index]])
				
				print(index)
			index=index+1
			}
		}
	}


	#save everything
	save.image(file = "synthetic_networks.Rdata")


	














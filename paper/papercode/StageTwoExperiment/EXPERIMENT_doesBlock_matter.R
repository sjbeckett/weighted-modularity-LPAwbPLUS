#Experiment to deduce that Stage II is important

#code without  vs. code with ( with action to notify when code block is called!)

source("LPA_wb_plusblock.R")

RunTest <- function(nrow,ncol) {

ntest = 20 #tests on matrix
ntot = 1000 # matrices to test



	modW=rep(0,ntot)
	modB=rep(0,ntot)
	modWwo=rep(0,ntot)
	modBwo=rep(0,ntot)
	countW=0
	countB=0

for( aa in 1:ntot) {

	W = matrix(runif(nrow*ncol),nrow,ncol) # random quantitative network - each element chosen uniform randomly in [0,1]
	A = 1*(W>0.45)  #  binary version of this network according to some arbitrary threshold

	switchW = 0
	switchB = 0


	for ( bb in 1:ntest) { # test algorithm with block 5 times on each network
		
		WW = LPA_wb_plusblockon(W)
		BB = LPA_wb_plusblockon(A)

		if( WW[[4]] == 1){ # if block used for weighted
			switchW = 1
			countW=countW+1
			#save maximum modularity
			if(WW[[3]]>modW[aa]){
				modW[aa]=WW[[3]]
			}
		}
		if( BB[[4]] == 1){ # if block used for binary
			switchBB = 1
			countB=countB+1
			if(BB[[3]]>modB[aa]){
				modB[aa]=BB[[3]]
			}
		}

		

		


	}
	
	if (switchW == 1){ # if block called - check result of alg without block in Weighted
		for( bb in 1:ntest) {
				WW = LPA_wb_plusblockoff(W)
				#save maximum modularity
				if(WW[[3]]>modWwo[aa]){
					modWwo[aa]=WW[[3]]
				}
		}
	} else { # if switch not called - this is the alg without block
		modWwo[aa]=modW[aa]
	}

	if (switchB == 1){# if block called - check result of alg without block in Binary
		for( bb in 1:ntest) {
			BB = LPA_wb_plusblockoff(A)
			#save max mod
			if(BB[[3]]>modBwo[aa]){
				modBwo[aa]=BB[[3]]
			}

		}
	} else { #if no switch activated - this is alg without block.
		modBwo[aa] = modB[aa]
	}
print(aa)
}


return(list(modB=modB, modBwo=modBwo, modW=modW, modWwo=modWwo))

}




nrow = 5
ncol = 10
A = RunTest(nrow,ncol)
B = RunTest(20,20)

save.image(file = "EXPERIMENT.Rdata")

#exp 1:    5 rows,10 cols
print(paste(sum(A$modB>0), " # num binary networks stage two called in "))   #156
print(paste(sum(A$modW>0), " #num weighted networks stage two called in "))  #93

print(paste(sum(A$modB>A$modBwo), " #number of improved modularity scores due to stage two (binary)"))    #0
print(paste(sum(A$modW>A$modWwo), " # number of improved modularity scores due to stage two (weighted)"))  #93

print(paste(sum(A$modB<A$modBwo), "  # number of modularity reductions due to stage two (binary)"))    #0
print(paste(sum(A$modW<A$modWwo), "  # number of modularity reductions due to stage two (weighted)"))  #0
#exp 2 : 20 rows, 20 cols
print(paste(sum(B$modB>0), "  # num binary networks stage two called in "))    #210
print(paste(sum(B$modW>0), "  #num weighted networks stage two called in "))   #15

print(paste(sum(B$modB>B$modBwo), "  #number of improved modularity scores due to stage two (binary)"))    #0
print(paste(sum(B$modW>B$modWwo), "  # number of improved modularity scores due to stage two (weighted)"))    #15

print(paste(sum(B$modB<B$modBwo), "  # number of modularity reductions due to stage two (binary)"))   #0
print(paste(sum(B$modW<B$modWwo), "  # number of modularity reductions due to stage two (weighted)"))   #0




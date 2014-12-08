#compares information between list1 and list2 and returns normalised mutual information

NormalisedMutualInformation <- function(listlabels1,listlabels2) {

# joint list of row and col labels for different partitions of the same network

#listlabels1 = c(rowlabels1,col_labels1)
#listlabels2 = c(rowlabels2,col_labels2)

Unirow1 = unique(listlabels1)
Unirow2 = unique(listlabels2)

NumNodes = length(listlabels1) # number of nodes

NumModules1 = length(Unirow1)
NumModules2 = length(Unirow2)

ConfusionMatrix = matrix(0,NumModules1,NumModules2)

for(aa in 1:NumModules1) {
	NODES_I = which(listlabels1 == Unirow1[aa])
	for(bb in 1:NumModules2) {
		NODES_J = which(listlabels2 == Unirow2[bb])
		NODES_IJ = intersect(NODES_I,NODES_J)
		ConfusionMatrix[aa,bb] = length(NODES_IJ)
	}
}

N_A = rowSums(ConfusionMatrix)
N_B = colSums(ConfusionMatrix)
N_Tot = sum(ConfusionMatrix)


NormA = sum(N_A * log(N_A/N_Tot))
NormB = sum(N_B * log(N_B/N_Tot))

MutInformation = 0
for (aa in 1:NumModules1) {
	for(bb in 1:NumModules2) {
		if(ConfusionMatrix[aa,bb]>0)
			MutInformation = MutInformation	+ ConfusionMatrix[aa,bb]* log(ConfusionMatrix[aa,bb]*N_Tot/(N_A[aa]*N_B[bb]))
	}
}

NMI = -2* MutInformation / (NormA+NormB)

return(NMI)
}




# GetModularInformation.R
# Uses the matrix, found modularity and module labels to provide more modularity related information.
# Gives normalized and realized modularity scores and row/column node names associated with each module
# Author :  Stephen Beckett ( https://github.com/sjbeckett/weighted-modularity-LPAwbPLUS )
# MIT License

GetModularInformation <- function(MATRIX,MODULE_INFO) {
	#takes two inputs, the initial matrix (MATRIX) and the output from running LPA_wb_plus or DIRT_LPA_wb_plus (MODULE_INFO)
	#returns the list of nodes associated with each module as well as normalised and realised modularity


	norm_mod = MODULE_INFO$modularity / MAXIMUM_BIP_MOD(MATRIX, MODULE_INFO$Row_labels, MODULE_INFO$Col_labels)
	realized_mod = RealizedMod(MATRIX,  MODULE_INFO$Row_labels, MODULE_INFO$Col_labels) 




	SIZE <- dim(MATRIX)
	#get row and column names
	rNodeNames <- rownames(MATRIX)
	cNodeNames <- colnames(MATRIX)
	
	# if no names assigned - assign generic names to rows
	if(  is.null(rNodeNames)==TRUE ) { 
		rNodeNames <- c()
		for(aa in 1:SIZE[1]){
			rNodeNames[aa] <- paste("R",aa,sep="")
		}

	rownames(MATRIX) <- rNodeNames 
	}

	# if no names assigned - assign generic names to columns
	if(  is.null(cNodeNames)==TRUE ) { 
		cNodeNames <- c()
		for(aa in 1:SIZE[2]){
			cNodeNames[aa] <- paste("C",aa,sep="")
		}

	colnames(MATRIX) <- cNodeNames 
	}

	
	#List of module values
	Modules = unique(MODULE_INFO$Row_labels)
	NumModules = length(Modules)

	#create lists (one for rows, one for columns) showing which nodes are within which module
	RowNodes <- list()
	ColNodes <- list()
	
	for(eachModule in 1:NumModules) {
		# find the names of the rows where the module number is equal to the current module number being looked at
		ROW_NODES_THIS_MODULE <- rNodeNames[which(MODULE_INFO$Row_labels==Modules[eachModule])] 
		
		# find the names of the cols where the module number is equal to the current module number being looked at
		COL_NODES_THIS_MODULE <- cNodeNames[which(MODULE_INFO$Col_labels==Modules[eachModule])]

		RowNodes <- c(RowNodes,list(ROW_NODES_THIS_MODULE))
		ColNodes <- c(ColNodes,list(COL_NODES_THIS_MODULE))
	}


return(list(number_of_modules = NumModules, modularity = MODULE_INFO$modularity, normalised_modularity = norm_mod, realized_modularity = realized_mod, RowNodesInModules = RowNodes, ColNodesInModules = ColNodes, ModuleLabels=Modules, Row_labels = MODULE_INFO$Row_labels, Col_labels = MODULE_INFO$Col_labels))
}




###############################################

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



###############################################

MAXIMUM_BIP_MOD <- function(MATRIX,redlabels,bluelabels) {

	A= rowSums(MATRIX)
	B= colSums(MATRIX)
	m = sum(MATRIX)

	holdsum = 0
	for( aa in 1:nrow(MATRIX)) {
		for( bb in 1:ncol(MATRIX)) {
			if(redlabels[aa] == bluelabels[bb]) {
				holdsum = holdsum + A[[aa]]*B[[bb]]/m
			}
		}
	}

	Qmax =	(m - holdsum)/m

	return(Qmax)
}

###############################################

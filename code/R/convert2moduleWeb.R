# convert2moduleWeb.R
# Use output from LPAwb+ to create a moduleWeb class object that can interface with the bipartite library in R
# Author :  Stephen Beckett ( https://github.com/sjbeckett/weighted-modularity-LPAwbPLUS )
# MIT License

#Specifically this can be used in order to make use of the plotModuleWeb function in bipartite to visualise detected modular structure
convert2moduleWeb <- function(NET,MODINFO){
	#NET - network biadjacency matrix
	#MODINFO - output from running LPA_wb_plus or Exhaustive_LPA_wb_plus , contains
	
	#if names are NULL - assign names to network
	if (is.null(rownames(NET)))
		rownames(NET) = 1:dim(NET)[1]

	if (is.null(colnames(NET)))
		colnames(NET) = 1:dim(NET)[2]

	#module sorting
	ROW_IX <- order(MODINFO$Row_labels)
    	COL_IX <- order(MODINFO$Col_labels)
	ROWS <- MODINFO$Row_labels[ROW_IX]
    	COLS <- MODINFO$Col_labels[COL_IX]
	MODS <- unique(ROWS)
	LMod <- length(MODS)

	#Creating the modules matrix
	#as LPAwb+ does not search at multiple depths can fix the first two columns
	Col1 <- c(0,rep(1,LMod)) #depth
	Col2 <- c(1,rep(0,LMod-1),1) # module markers
	Vals <- 1:(length(ROWS)+length(COLS))
	#create matrix store, with LMod rows - where within module indices are shown and others are zero.
	store <- Vals
	for(bb in MODS) {
		Ix1 <- MODINFO$Row_labels!=bb
		Ix2 <- MODINFO$Col_labels!=bb
		Ix <- c(Ix1,Ix2)
		New <- Vals
		New[Ix] <- 0
		store <- rbind(store,New)
	}
	
	modules <- cbind(Col1,Col2,store)
	rownames(modules) <- NULL
	colnames(modules) <- NULL

	#assigning values to moduleWeb object	

	ModwebObj <- new("moduleWeb")
	ModwebObj@originalWeb <- NET
	ModwebObj@moduleWeb <- NET[ROW_IX,COL_IX]
	ModwebObj@orderA <- ROW_IX
	ModwebObj@orderB <- COL_IX
	ModwebObj@modules <- modules
	ModwebObj@likelihood <- MODINFO$modularity

return(ModwebObj)
}


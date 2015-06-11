# MODULARPLOT.R
# Plots modular bipartite network configuration using row and column module labels. Row and column numbers indicate the ordering of rows and columns within the original network.
# Author :  Stephen Beckett ( https://github.com/sjbeckett/weighted-modularity-LPAwbPLUS )
# MIT License

MODULARPLOT <- function(MATRIX,OUT) {
    #MATRIX:  original bipartite network
    #OUT: output from modularity algorithm containing row and column module labels

    #1. Sort original network, so that rows are positioned according to modules
    ROW_IX <- order(OUT$Row_labels)
    COL_IX <- order(OUT$Col_labels)
    ROWS <- OUT$Row_labels[ROW_IX]
    COLS <- OUT$Col_labels[COL_IX]
    MODS <- unique(ROWS)

    #2. Plot this network   
    MATRIXPLOT(MATRIX[ROW_IX,COL_IX], COL_IX, ROW_IX)

    #Add the modular boundaries using output information - drawing from top left to bottom right
    xcurr<-0.5
    ycurr<-length(ROWS)+0.5

    for(aa in 1:length(MODS)) {
        Rsize <- sum(ROWS==MODS[aa]) #Number of rows in this module
        Csize <- sum(COLS==MODS[aa]) #Number of columns in this module
        rect(xcurr,ycurr,xcurr+Csize,ycurr-Rsize,border='red')
        xcurr <- xcurr+Csize
        ycurr <- ycurr-Rsize
    }


}


MATRIXPLOT <- function(INMATRIX,col_ix,row_ix) {


	#Graphic Control
	BACKGROUNDCOLOR <- rgb(173,235,255,max=255)
	offset <- 0.3 #max of 0.5,425

	columns=1:dim(INMATRIX)[2]
	rows=1:dim(INMATRIX)[1]

	#flip matrix rows in reverse so can draw it onto grid.
	INMATRIX <- t(INMATRIX[ nrow(INMATRIX):1, ] )

	RANGE <- range(INMATRIX)
	minR <- RANGE[1]
	scalefactor <- 1/(RANGE[2]-RANGE[1])
	
	#Draw background and label axes
	image(columns,rows,0*INMATRIX, col=BACKGROUNDCOLOR,xaxt='n',yaxt='n')
	axis(1,at=columns,labels=col_ix)
	axis(2,rev(rows),at=rows,labels=rev(row_ix))
	box()

	#plot filled positions as rectangles
	for (b in 1:length(rows)) {
		for (a in 1:length(columns)) { 
			if (INMATRIX[a,b]!=0){
				COL <- 	1-(INMATRIX[a,b]-minR)*scalefactor
				rect(a-offset,b-offset,a+offset,b+offset,col=rgb(COL,COL,COL),border='NA')
			}
		}
	}


}



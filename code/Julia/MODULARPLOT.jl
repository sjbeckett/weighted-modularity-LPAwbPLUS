# MODULARPLOT.jl
# Plots modular bipartite network configuration using row and column module labels. Row and column numbers indicate the ordering of rows and columns within the original network.
# Plotting performed using the PyPlot library.
# Author :  Stephen Beckett ( https://github.com/sjbeckett/weighted-modularity-LPAwbPLUS )
# MIT License

function MODULARPLOT(MATRIX,OUT)
    #MATRIX:  original bipartite network
    #OUT: output from modularity algorithm containing row and column module labels
    print("Attempting to import PyPlot library...")
    require("PyPlot")
    print("Plotting...")
    #1. Sort original network, so that rows are positioned according to modules
    ROW_IX = sortperm(vec(OUT[1]));
    COL_IX = sortperm(vec(OUT[2]));
    ROWS = OUT[1][ROW_IX];
    COLS = OUT[2][COL_IX];
    MODS = unique(ROWS);

    #2. Plot this network   
    MATRIXPLOT(MATRIX[ROW_IX,COL_IX], COL_IX, ROW_IX)

    #Add the modular boundaries using output information - drawing from top left to bottom right
    xcurr=0.5
    ycurr=length(ROWS)+0.5

    for aa = 1:length(MODS)
        Rsize = sum(ROWS.==MODS[aa]) #Number of rows in this module
        Csize = sum(COLS.==MODS[aa]) #Number of columns in this module
        PyPlot.fill_betweenx([ycurr-0.5,ycurr-Rsize-0.5],xcurr-0.5,xcurr+Csize-0.5,color="none",edgecolor="red")
        xcurr = xcurr+Csize
        ycurr = ycurr-Rsize
    end


end

###############################################

function MATRIXPLOT(INMATRIX,col_ix,row_ix)
	BACKGROUNDCOLOR = (173/255,235/255,255/255);
	offset = 0.3;

	xs = length(col_ix)
	ys = length(row_ix)

	#initial background
	PyPlot.pcolor(INMATRIX)
	PyPlot.xticks((1:xs)-0.5,col_ix)
	PyPlot.yticks((1:ys)-0.5,reverse(row_ix))
	PyPlot.fill_betweenx([0:ys],0,xs,color=BACKGROUNDCOLOR)

	#fill in
	minR = minimum(INMATRIX)
	RANGE = maximum(INMATRIX)-minR
	scalefactor = 1/RANGE

	for b = 1:xs
		for a = 1:ys
			if INMATRIX[1+ys-a,b]!=0
				COL = 1-(INMATRIX[1+ys-a,b]-minR)*scalefactor
				PyPlot.fill_betweenx([a-offset-0.5,a+offset-0.5],b-offset-0.5,b+offset-0.5,color=(COL,COL,COL))
			end
		end
	end
end

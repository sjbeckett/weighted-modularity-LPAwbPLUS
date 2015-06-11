function MODULARPLOT(MAT,ROW_LABS,COL_LABS)


    %1. Sort original network, so that rows are positioned according to modules
    [~,ROW_IX] = sort(ROW_LABS);
    [~,COL_IX] = sort(COL_LABS);
    ROWS = ROW_LABS(ROW_IX);
    COLS = COL_LABS(COL_IX);
    MODS = unique(ROWS);
    %. Plot this network   
    MATRIXPLOT(MAT(ROW_IX,COL_IX), COL_IX, ROW_IX);

    %Add the modular boundaries using output information - drawing from top left to bottom right
    xcurr = 0.5;
    ycurr = 0.5;

    for aa = 1:length(MODS)
        Rsize = sum(ROWS==MODS(aa)); %Number of rows in this module
        Csize = sum(COLS==MODS(aa)); %Number of columns in this module
        rect = rectangle('Position',[xcurr,ycurr,Csize,Rsize]);
        set(rect,'EdgeColor','red')
        xcurr = xcurr+Csize;
        ycurr = ycurr+Rsize;
    end



end



function MATRIXPLOT(INMATRIX,col_ix,row_ix)

	BACKGROUNDCOLOR = [173,235,255]/255;
	offset = 0.3; %max of 0.5,425
	RL = offset*2;

	columns = 1:size(INMATRIX,2);
	rows = 1:size(INMATRIX,1);

	minR = min(INMATRIX(:));
	RANGE = max(INMATRIX(:))-minR;
	scalefactor = 1/RANGE;

	%Draw background and label axes
	imagesc(columns,rows,0*INMATRIX)
	rec=rectangle('Position',[0.5,0.5,length(columns),length(rows)]);
	set(rec,'FaceColor',BACKGROUNDCOLOR);
	set(gca,'XTick',columns,'XTickLabel',col_ix,'YTick',rows,'YTickLabel',row_ix);

	%plot filled positions as rectangles
	for a = 1:length(rows)
	        for b = 1:length(columns) 
	            if INMATRIX(a,b)~=0
				COL = 	1-(INMATRIX(a,b)-minR)*scalefactor;
				rectang= rectangle('Position',[b-offset,a-offset,RL,RL]);
		                set(rectang,'FaceColor',([COL,COL,COL]),'EdgeColor','none'); 
	            end
	        end
	end

end

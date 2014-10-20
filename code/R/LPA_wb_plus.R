###############################################

BarbersMatrix <- function(MATRIX) {
	return(MATRIX - (t(t(rowSums(MATRIX)))%*%t(colSums(MATRIX)))/sum(MATRIX))
}

###############################################

WEIGHTEDMODULARITY <- function(BMatrix,Matsum,redlabels,bluelabels) {

	holdsum = 0;

	for (rr in 1:length(redlabels)) {
	    for (cc in 1:length(bluelabels)) {
	        kroneckerdelta = redlabels[rr] == bluelabels[cc];
	        holdsum = holdsum + BMatrix[rr,cc] * kroneckerdelta;
	    }
	}
	return(holdsum/Matsum)
}

###############################################

DIVISION <- function(redlabels,bluelabels) {
	divisionsFound= c();

	#for (aa in 1:min(max(redlabels),max(bluelabels))) {
	#	if( (length(which(redlabels == aa))>0) && (length(which(bluelabels==aa))>0) ) { ###
	#		divisionsFound = c(divisionsFound,aa)
	#	}
	#}
	divisionsFound <- intersect(redlabels,bluelabels)

	return(divisionsFound)
}

###############################################

StageTwo_LPAwbdash <- function(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels, Qb_now) {

	divisionsFound = DIVISION(redlabels,bluelabels);
	IterateFlag = 1
	while(IterateFlag == 1) {
		CombinedDivisionsThisTime = 0
		if(length(divisionsFound) > 1) {
			for(div1check in 1:(length(divisionsFound)-1)) {
				for(div2check in (div1check+1):length(divisionsFound)) {
					CHECK_RED = redlabels
					CHECK_RED[which(redlabels==divisionsFound[div1check])] = divisionsFound[div2check];
					CHECK_BLUE = bluelabels
					CHECK_BLUE[which(bluelabels==divisionsFound[div1check])] = divisionsFound[div2check];
					QQ = WEIGHTEDMODULARITY(BMatrix,Matsum,redlabels,bluelabels);			
					if(QQ > Qb_now) { #If this arrangement looks good
                    				FoundBetter = 0;
                    				for(aa in 1:length(divisionsFound)) {
                        				CHECK_RED2 = redlabels
                        				CHECK_RED2[which(redlabels==divisionsFound[aa])] = divisionsFound[div1check];
                        				CHECK_BLUE2 = bluelabels
                        				CHECK_BLUE2[which(bluelabels==divisionsFound[aa])] = divisionsFound[div1check];
                        				if(WEIGHTEDMODULARITY(BMatrix,Matsum,CHECK_RED2,CHECK_BLUE2) > QQ) {
                            					FoundBetter = 1;
                        				}
                        				CHECK_RED2 = redlabels
                        				CHECK_RED2[which(redlabels==divisionsFound[aa])] = divisionsFound[div2check];
                        				CHECK_BLUE2 = bluelabels
                        				CHECK_BLUE2[which(bluelabels==divisionsFound[aa])] = divisionsFound[div2check];
                        				if(WEIGHTEDMODULARITY(BMatrix,Matsum,CHECK_RED2,CHECK_BLUE2) > QQ) {
                            					FoundBetter = 1;
                        				}
                    				}
                    				if(FoundBetter == 0) { #If no better configuration found - JOIN.
                        				redlabels = CHECK_RED;
                        				bluelabels = CHECK_BLUE;
                        				CombinedDivisionsThisTime = CombinedDivisionsThisTime + 1;
                    				}
                			}
				}
			}
			if(CombinedDivisionsThisTime == 0) {#If no divisions were joined move on
		            IterateFlag = 0;
		        }
		}
		else {
		IterateFlag = 0;		
		}
		
		outlist=StageOne_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels);  ##
		redlabels = outlist[[1]]
		bluelabels = outlist[[2]]
		Qb_now = outlist[[3]]
		divisionsFound = DIVISION(redlabels,bluelabels);
	}

return(list(redlabels, bluelabels, Qb_now))
}


###############################################

StageOne_LPAwbdash <- function(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels) {
	#Create storage containers for total marginals attached to each red(row)
	#label and blue(column) label

	BLUELABELLENGTH=length(bluelabels);
	TotalRedDegrees = NA*1:max(redlabels)
	
	TotalBlueDegrees = NA*1:max(length(bluelabels),length(redlabels))

	#Fill up these containers according to current labels
	#Red
	for(aa in 1:length(redlabels)) {
	    if(is.na(TotalRedDegrees[redlabels[aa]])) {
        	TotalRedDegrees[redlabels[aa]] = row_marginals[aa];
	    }
	    else {
        	TotalRedDegrees[redlabels[aa]] = TotalRedDegrees[redlabels[aa]] + row_marginals[aa];
	    }
	}

	#Blue
	if(sum(is.na(bluelabels)) != BLUELABELLENGTH) { #occurs first time through as blue nodes unlabelled
	    for(bb in 1:BLUELABELLENGTH) {
        	if(is.na(TotalBlueDegrees[bluelabels[bb]])) {
	            TotalBlueDegrees[bluelabels[bb]] = col_marginals[bb];
		}
        	else {
        	    TotalBlueDegrees[bluelabels[bb]] = TotalBlueDegrees[bluelabels[bb]] + col_marginals[bb];
        	}
	    }
	}
	else {
	    TotalBlueDegrees = 0*1:max(length(bluelabels),length(redlabels)); #####
	}
	
	
	#locally maximise modularity!
	outlist = LOCALMAXIMISATION(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees)
	redlabels = outlist[[1]]
	bluelabels = outlist[[2]]
	Qb_now = outlist[[3]]

	return(list(redlabels, bluelabels, Qb_now))

}



###############################################

LOCALMAXIMISATION <-  function(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees) {
	
	#Find score for current partition
	QbAfter = WEIGHTEDMODULARITY(BMatrix,Matsum,redlabels,bluelabels);

	if(is.na(QbAfter)) { QbAfter = -999 }

	IterateFlag = 1;
	while(IterateFlag == 1) {
		#Save old information
		QbBefore = QbAfter;
		old_redlabels = redlabels
		old_bluelabels = bluelabels
		old_TRD = TotalRedDegrees
		old_TBD = TotalBlueDegrees

		#Update Blue Nodes (using red node info.)
		bluelabelchoices = 1:max(redlabels);

		for(bb in 1:length(bluelabels)) {
			if(is.na(bluelabels[bb]) == FALSE) {
				TotalBlueDegrees[bluelabels[bb]] = TotalBlueDegrees[bluelabels[bb]] - col_marginals[bb]; #20
			}
			changebluelabeltest = NA*1:length(bluelabelchoices)

			for(ww in 1:length(bluelabelchoices)) {
				changebluelabeltest[ww] = sum( (redlabels == bluelabelchoices[ww]) * MATRIX[,bb])  -  col_marginals[bb]*TotalRedDegrees[bluelabelchoices[ww]]/Matsum ;
			}

			#assign new label based on maximisation of above condition  

			labels = which(changebluelabeltest == max(changebluelabeltest,na.rm =TRUE));#30
			newlabelindex = labels[sample(1:length(labels),1)];
			bluelabels[bb] = bluelabelchoices[newlabelindex[1]];
			if(bluelabels[bb] > length(TotalBlueDegrees)) {
	                	TotalBlueDegrees[bluelabels[bb]] = 0;
	        	}
            
        		#Update total marginals on new labelling
            		TotalBlueDegrees[bluelabels[bb]] = TotalBlueDegrees[bluelabels[bb]] + col_marginals[bb];
        	}

		#Now update red node labels based on blue node information
		redlabelchoices = 1:max(bluelabels);

		for(aa in 1:length(redlabels)) {
			TotalRedDegrees[redlabels[aa]] = TotalRedDegrees[redlabels[aa]] - row_marginals[aa];
			changeredlabeltest = NA*1:length(redlabelchoices)

			for(ww in 1:length(redlabelchoices)) {
				changeredlabeltest[ww] = sum( (bluelabels == redlabelchoices[ww]) * MATRIX[aa,])  -  row_marginals[aa]*TotalBlueDegrees[redlabelchoices[ww]]/Matsum ; #49
			}
			
			#assign new label based on maximisation of above condition
			labels = which(changeredlabeltest == max(changeredlabeltest,na.rm = TRUE));
			newlabelindex = labels[sample(1:length(labels),1)];
 			redlabels[aa] = redlabelchoices[newlabelindex[1]];

			if(redlabels[aa] > length(TotalRedDegrees)) {
				TotalRedDegrees[redlabels[aa]] = 0;
			}
			TotalRedDegrees[redlabels[aa]] = TotalRedDegrees[redlabels[aa]] + row_marginals[aa];
		}
		
		
		#Find the new modularity score based on node label updates.
        	QbAfter = WEIGHTEDMODULARITY(BMatrix,Matsum,redlabels,bluelabels);
        
        	#If this modularity is not as good as previous stop iterating and
        	#use that previous best information

        	if(QbAfter <= QbBefore) {
            		redlabels = old_redlabels
	            	bluelabels = old_bluelabels
        	    	TotalRedDegrees = old_TRD
        	    	TotalBlueDegrees = old_TBD
        	    	IterateFlag = 0;
        	}
		
	}

	Qb_now = QbAfter;

	return(list(redlabels, bluelabels, Qb_now))
}


###############################################

LPA_wb_plus <- function(MATRIX) {

	flipped = 0;
	if(dim(MATRIX)[1] > dim(MATRIX)[2]) {
		MATRIX = t(MATRIX);
		flipped = 1;
	}

	Matsum = sum(MATRIX);
	col_marginals = colSums(MATRIX);
	row_marginals = rowSums(MATRIX);
	BMatrix = BarbersMatrix(MATRIX);

	#initiliase labels
	redlabels = 1:dim(MATRIX)[1];
	bluelabels = NA*1:dim(MATRIX)[2];

	#Run Phase 1: Locally update lables to maximise Qb
	outlist = StageOne_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels)
	redlabels = outlist[[1]]
	bluelabels = outlist[[2]]
	Qb_now = outlist[[3]]

	#Run Phase 2: Connect divisions from top-down if it improves Qb, then run
	#phase 1 again. Repeat till Qb cannot be improved.
	outlist = StageTwo_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,Qb_now);
	redlabels = outlist[[1]]
	bluelabels = outlist[[2]]
	Qb_now = outlist[[3]]

	if(flipped==1) { #If matrix was flipped, swap row and column labels
    		holder = redlabels
		redlabels = bluelabels;
		bluelabels = holder;
	}

	return(list(redlabels, bluelabels, Qb_now))
}

###############################################




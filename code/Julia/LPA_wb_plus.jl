###############################################

function BarbersMatrix(MATRIX)
	return MATRIX - (sum(MATRIX,2)*sum(MATRIX,1))/sum(MATRIX);
end

###############################################

function WEIGHTEDMODULARITY(BMatrix,Matsum,redlabels,bluelabels)

	holdsum = 0;

	for rr = 1:length(redlabels)
	    for cc = 1:length(bluelabels)
	        kroneckerdelta = redlabels[rr] == bluelabels[cc];
	        holdsum += BMatrix[rr,cc] * kroneckerdelta;
	    end
	end
	return holdsum/Matsum;
end

###############################################

function DIVISION(redlabels,bluelabels)
	divisionsFound= Float64[];

	for aa = 1:min(maximum(redlabels),maximum(bluelabels))
		if( (isempty(findin(redlabels,aa))==false) && (isempty(findin(bluelabels,aa))==false) )
			push!(divisionsFound,aa)
		end
	end

	return divisionsFound;
end

###############################################

function StageTwo_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels, Qb_now)

	divisionsFound = DIVISION(redlabels,bluelabels);
	IterateFlag = 1;
	while IterateFlag == 1
		CombinedDivisionsThisTime = 0;
		if length(divisionsFound) > 1
			for div1check = 1:length(divisionsFound)-1
				for div2check = div1check+1:length(divisionsFound)
					CHECK_RED = copy(redlabels);
					CHECK_RED[findin(redlabels,divisionsFound[div1check])] = divisionsFound[div2check];
					CHECK_BLUE = copy(bluelabels);
					CHECK_BLUE[findin(bluelabels,divisionsFound[div1check])] = divisionsFound[div2check];
					QQ = WEIGHTEDMODULARITY(BMatrix,Matsum,redlabels,bluelabels);			
					if QQ > Qb_now #If this arrangement looks good
                    				FoundBetter = 0;
                    				for aa = 1:length(divisionsFound)
                        				CHECK_RED2 = copy(redlabels);
                        				CHECK_RED2[findin(redlabels,divisionsFound[aa])] = divisionsFound[div1check];
                        				CHECK_BLUE2 = copy(bluelabels);
                        				CHECK_BLUE2[findin(bluelabels,divisionsFound[aa])] = divisionsFound[div1check];
                        				if WEIGHTEDMODULARITY(BMatrix,Matsum,CHECK_RED2,CHECK_BLUE2) > QQ
                            					FoundBetter = 1;
                        				end
                        				CHECK_RED2 = copy(redlabels);
                        				CHECK_RED2[findin(redlabels,divisionsFound[aa])] = divisionsFound[div2check];
                        				CHECK_BLUE2 = copy(bluelabels);
                        				CHECK_BLUE2[findin(bluelabels,divisionsFound[aa])] = divisionsFound[div2check];
                        				if WEIGHTEDMODULARITY(BMatrix,Matsum,CHECK_RED2,CHECK_BLUE2) > QQ
                            					FoundBetter = 1;
                        				end
                    				end
                    				if FoundBetter == 0 #If no better configuration found - JOIN.
                        				redlabels = copy(CHECK_RED);
                        				bluelabels = copy(CHECK_BLUE);
                        				CombinedDivisionsThisTime += 1;
                    				end
                			end
				end
			end
			if CombinedDivisionsThisTime == 0 #If no divisions were joined move on
		            IterateFlag = 0;
		        end
		else
		IterateFlag = 0;		
		end
		
		(redlabels,bluelabels,Qb_now) = StageOne_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels);
		divisionsFound = DIVISION(redlabels,bluelabels);
	end

return redlabels, bluelabels, Qb_now;

end


###############################################

function StageOne_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels)
	#Create storage containers for total marginals attached to each red(row)
	#label and blue(column) label

	BLUELABELLENGTH=length(bluelabels);
	TotalRedDegrees = NaN*ones(maximum(redlabels),1);
	
	TotalBlueDegrees = NaN*ones(max(length(bluelabels),length(redlabels)),1);

	#Fill up these containers according to current labels
	#Red
	for aa = 1:length(redlabels)
	    if isnan(TotalRedDegrees[redlabels[aa]])
        	TotalRedDegrees[redlabels[aa]] = row_marginals[aa];
	    else
        	TotalRedDegrees[redlabels[aa]] += row_marginals[aa];
	    end
	end

	#Blue
	if sum(isnan(bluelabels)) != BLUELABELLENGTH #occurs first time through as blue nodes unlabelled
	    for bb = 1:BLUELABELLENGTH
        	if isnan(TotalBlueDegrees[bluelabels[bb]])
	            TotalBlueDegrees[bluelabels[bb]] = col_marginals[bb];
        	else
        	    TotalBlueDegrees[bluelabels[bb]] += col_marginals[bb];
        	end
	    end
	else
	    TotalBlueDegrees = zeros(max(length(bluelabels),length(redlabels)),1); #####
	end
	
	
	#locally maximise modularity!
	(redlabels,bluelabels,Qb_now) = LOCALMAXIMISATION(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees);

	return redlabels, bluelabels, Qb_now;

end



###############################################

function LOCALMAXIMISATION(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees)

	#Find score for current partition
	QbAfter = WEIGHTEDMODULARITY(BMatrix,Matsum,redlabels,bluelabels);

	IterateFlag = 1;
	while IterateFlag == 1
		#Save old information
		QbBefore = QbAfter;
		old_redlabels = copy(redlabels);
		old_bluelabels = copy(bluelabels);
		old_TRD = copy(TotalRedDegrees);
		old_TBD = copy(TotalBlueDegrees); #13

		#Update Blue Nodes (using red node info.)
		bluelabelchoices = 1:maximum(redlabels);

		for bb = 1:length(bluelabels)
			if isnan(bluelabels[bb]) == false
				TotalBlueDegrees[bluelabels[bb]] -= col_marginals[bb]; #20
			end
			changebluelabeltest = NaN*ones(1,length(bluelabelchoices));

			for ww = 1:length(bluelabelchoices)
				changebluelabeltest[ww] = sum( (redlabels .== bluelabelchoices[ww]) .* MATRIX[:,bb])  -  col_marginals[bb]*TotalRedDegrees[bluelabelchoices[ww]]/Matsum ;
			end

			#assign new label based on maximisation of above condition  

			labels = findin(changebluelabeltest, maximum(changebluelabeltest));#30
			newlabelindex = labels[rand(1:length(labels),1)];
			bluelabels[bb] = bluelabelchoices[newlabelindex[1]];
			
			if bluelabels[bb] > length(TotalBlueDegrees)
	                	TotalBlueDegrees[bluelabels[bb]] = 0;
	        	end
            
        		#Update total marginals on new labelling
            		TotalBlueDegrees[bluelabels[bb]] += col_marginals[bb];
        	end#40

		#Now update red node labels based on blue node information
		redlabelchoices = 1:maximum(bluelabels);

		for aa = 1:length(redlabels)
			TotalRedDegrees[redlabels[aa]] -= row_marginals[aa];
			changeredlabeltest = NaN*ones(1,length(redlabelchoices));

			for ww = 1:length(redlabelchoices)
				changeredlabeltest[ww] = sum( (bluelabels .== redlabelchoices[ww]) .* MATRIX[aa,:])  -  row_marginals[aa]*TotalBlueDegrees[redlabelchoices[ww]]/Matsum ; #49
			end
			
			#assign new label based on maximisation of above condition
			labels = findin(changeredlabeltest,maximum(changeredlabeltest));
			newlabelindex = labels[rand(1:length(labels),1)];
 			redlabels[aa] = redlabelchoices[newlabelindex[1]];

			if redlabels[aa] > length(TotalRedDegrees)
				TotalRedDegrees[redlabels[aa]] = 0;
			end
			TotalRedDegrees[redlabels[aa]] +=  row_marginals[aa];
		end
		
		
		#Find the new modularity score based on node label updates.
        	QbAfter = WEIGHTEDMODULARITY(BMatrix,Matsum,redlabels,bluelabels);
        
        	#If this modularity is not as good as previous stop iterating and
        	#use that previous best information
        	if QbAfter <= QbBefore
            		redlabels = copy(old_redlabels);
	            	bluelabels = copy(old_bluelabels);
        	    	TotalRedDegrees = copy(old_TRD);
        	    	TotalBlueDegrees = copy(old_TBD);
        	    	IterateFlag = 0;
        	end
		
	end

	Qb_now = QbAfter;

	return redlabels, bluelabels, Qb_now;
end


###############################################

function LPA_wb_plus(MATRIX)

	flipped = 0;
	if size(MATRIX,1) > size(MATRIX,2)
		MATRIX = MATRIX';
		flipped = 1;
	end

	Matsum = sum(MATRIX);
	col_marginals = sum(MATRIX,1);
	row_marginals = sum(MATRIX,2);
	BMatrix = BarbersMatrix(MATRIX);

	#initiliase labels
	redlabels = [1:size(MATRIX,1)];
	bluelabels = [NaN*ones(1,size(MATRIX,2))];

	#Run Phase 1: Locally update lables to maximise Qb
	(redlabels,bluelabels,Qb_now) = StageOne_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels);

	#Run Phase 2: Connect divisions from top-down if it improves Qb, then run
	#phase 1 again. Repeat till Qb cannot be improved.
	(redlabels,bluelabels,Qb_now) = StageTwo_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,Qb_now);

	if flipped==1 #If matrix was flipped, swap row and column labels
    		holder = copy(redlabels);
		redlabels = bluelabels;
		bluelabels = holder;
	end

	return redlabels, bluelabels, Qb_now;
end

###############################################




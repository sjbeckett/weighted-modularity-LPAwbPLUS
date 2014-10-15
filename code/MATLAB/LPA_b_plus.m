%LPAb+ ALGORITHM

function [Qb_current,redlabels,bluelabels]=LPA_b_plus(MATRIX,parallel)

flipped=0;
%Make the rows the shortest dimension
if size(MATRIX,1) > size(MATRIX,2)
    MATRIX = MATRIX';
    flipped=1;
end

if parallel == 1
    matlabpool open
end
    
    


%Find properties of Matrix
info.matrix=MATRIX>0;
[info.row_num info.col_num]=size(MATRIX);
info.edge_num = sum(sum(MATRIX>0));
info.row_degrees = sum(MATRIX>0,2);
info.col_degrees = sum(MATRIX>0,1);


%Initialise red(row) and blue(column) labels
redlabels = 1:info.row_num;

bluelabels = nan(1,info.col_num);


%Run Phase 1: Locally update labels to gain local maximum of Qb
if parallel == 1
     [redlabels,bluelabels,Qb_current]=StageOne_LPAbdashPARALLEL(info,redlabels,bluelabels);
else
     [redlabels,bluelabels,Qb_current]=StageOne_LPAbdash(info,redlabels,bluelabels);
end


%Run Phase 2: Connect divisions from top-down if it improves Qb, then run
%phase 1 again. Repeat till Qb cannot be improved.
[redlabels,bluelabels,Qb_current]=StageTwo_LPAbdash(info,redlabels,bluelabels,Qb_current);



if flipped==1 %If matrix was flipped, swap row and column labels
    holder = redlabels;
    redlabels = bluelabels;
    bluelabels = holder;
end

if parallel == 1
    matlabpool close
end


end

function [BarbersModularity]=BARBERSMODULARITY(info,redlabels,bluelabels)


ProbMatrix = (info.row_degrees * info.col_degrees)./info.edge_num; %initialise probability matrix


BarbersMatrix = info.matrix - ProbMatrix;


holdsum = 0;


for rr = 1:info.row_num
    for cc = 1:info.col_num
       
        kroneckerdelta = redlabels(rr) == bluelabels(cc);
        
        holdsum = holdsum + BarbersMatrix(rr,cc) * kroneckerdelta;
        
    end
end


BarbersModularity = holdsum / info.edge_num;

end

function [redlabels,bluelabels,Qb_current]=StageOne_LPAbdash(info,redlabels,bluelabels)

%Create storage containers for total degrees attached to each red(row)
%label and blue(column) label

TotalRedDegrees = nan(max(redlabels),1);

if sum(isnan(bluelabels)) ~= length(bluelabels) %occurs first time through - as blue nodes unlabelled
    TotalBlueDegrees = nan(max(bluelabels),1);
end
    
%Fill up these containers according to current labels

%Red
for aa = 1:info.row_num
    if isnan(TotalRedDegrees(redlabels(aa)))
        TotalRedDegrees(redlabels(aa)) = info.row_degrees(aa);
    else
        TotalRedDegrees(redlabels(aa)) = TotalRedDegrees(redlabels(aa)) + info.row_degrees(aa);
    end
end

%Blue
if sum(isnan(bluelabels)) ~= length(bluelabels) %occurs first time through as blue nodes unlabelled
    for bb = 1:info.col_num
        if isnan(TotalBlueDegrees(bluelabels(bb)))
            TotalBlueDegrees(bluelabels(bb)) = info.col_degrees(bb);
        else
            TotalBlueDegrees(bluelabels(bb)) = TotalBlueDegrees(bluelabels(bb)) + info.col_degrees(bb);
        end
    end
else
    TotalBlueDegrees = [];
end



%Locally Maximise Modularity!

[redlabels,bluelabels,Qb_current]=LOCALMAXIMISATION(info,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees);

end

function [redlabels,bluelabels,Qb_current]=StageOne_LPAbdashPARALLEL(info,redlabels,bluelabels)

%Create storage containers for total degrees attached to each red(row)
%label and blue(column) label

TotalRedDegrees = nan(max(redlabels),1);

if sum(isnan(bluelabels)) ~= length(bluelabels) %occurs first time through - as blue nodes unlabelled
    TotalBlueDegrees = nan(max(bluelabels),1);
end
    
%Fill up these containers according to current labels

%Red
for aa = 1:info.row_num
    if isnan(TotalRedDegrees(redlabels(aa)))
        TotalRedDegrees(redlabels(aa)) = info.row_degrees(aa);
    else
        TotalRedDegrees(redlabels(aa)) = TotalRedDegrees(redlabels(aa)) + info.row_degrees(aa);
    end
end

%Blue
if sum(isnan(bluelabels)) ~= length(bluelabels) %occurs first time through as blue nodes unlabelled
    for bb = 1:info.col_num
        if isnan(TotalBlueDegrees(bluelabels(bb)))
            TotalBlueDegrees(bluelabels(bb)) = info.col_degrees(bb);
        else
            TotalBlueDegrees(bluelabels(bb)) = TotalBlueDegrees(bluelabels(bb)) + info.col_degrees(bb);
        end
    end
else
    TotalBlueDegrees = [];
end



%Locally Maximise Modularity!

[redlabels,bluelabels,Qb_current]=LOCALMAXIMISATIONPARALLEL(info,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees);

end

function [redlabels,bluelabels,Qb_current]=LOCALMAXIMISATION(info,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees)

%Find score for current partition
QbAfter = BARBERSMODULARITY(info,redlabels,bluelabels);



IterateFlag = 1;

while IterateFlag == 1
    
%     while QbBefore < QbAfter (this is handled by an if statement below)
        
        %Save old information
        QbBefore = QbAfter;
        old_redlabels = redlabels;
        old_bluelabels = bluelabels;
        old_TotRedDegree = TotalRedDegrees;
        old_TotBlueDegree = TotalBlueDegrees;
        
        %update blue nodes using red node information
        bluelabelchoices = 1:max(redlabels);%unique(redlabels);
                
        for bb = 1:info.col_num
            
            if isnan(bluelabels(bb))==0 %If node is labelled - reduce degree for that label so it can be reassigned
                TotalBlueDegrees(bluelabels(bb)) = TotalBlueDegrees(bluelabels(bb)) - info.col_degrees(bb);
            end
                
            changebluelabeltest=[];
            %check all possible label choices
            for ww = 1:length(bluelabelchoices)
                changebluelabeltest(ww) = sum((redlabels == bluelabelchoices(ww)).*info.matrix(:,bb)') - (info.col_degrees(bb)*TotalRedDegrees(bluelabelchoices(ww)) / info.edge_num);
            end
            
            %assign new label based on maximisation of above condition            

            labels = find(changebluelabeltest == max(changebluelabeltest)); %0's and 1's... only want 1's.
            newlabelindex = labels(randi(length(labels)));
            bluelabels(bb) = bluelabelchoices(newlabelindex);
            
            
            
            if bluelabels(bb) > length(TotalBlueDegrees)
                TotalBlueDegrees(bluelabels(bb)) = 0;
            end
            
            %Update total degrees on new labelling
            TotalBlueDegrees(bluelabels(bb)) = TotalBlueDegrees(bluelabels(bb)) + info.col_degrees(bb);

        end

        %Now update red node labels based on blue node information
        redlabelchoices = 1:max(bluelabels);%unique(bluelabels);
        
        for aa = 1:info.row_num
            
           TotalRedDegrees(redlabels(aa)) = TotalRedDegrees(redlabels(aa)) - info.row_degrees(aa);
           
           changeredlabeltest=[];
           
           for ww = 1:length(redlabelchoices)
               changeredlabeltest(ww) = sum((bluelabels == redlabelchoices(ww)).*info.matrix(aa,:)) - (info.row_degrees(aa)*TotalBlueDegrees(redlabelchoices(ww)) / info.edge_num);
           end
            
           %assign new label based on maximisation of above condition 
           
           labels = find(changeredlabeltest == max(changeredlabeltest)); %0's and 1's... only want 1's.
           newlabelindex = labels(randi(length(labels)));
           redlabels(aa) = redlabelchoices(newlabelindex);
           
           if redlabels(aa) > length(TotalRedDegrees)
                TotalRedDegrees(redlabels(aa)) = 0;
           end
           
           TotalRedDegrees(redlabels(aa)) = TotalRedDegrees(redlabels(aa)) + info.row_degrees(aa);
          
        end
    
        %Find the new modularity score based on node label updates.
        QbAfter = BARBERSMODULARITY(info,redlabels,bluelabels);
        
        %If this modularity is not as good as previous stop iterating and
        %use that previous best information
        if QbAfter <= QbBefore
            redlabels = old_redlabels;
            bluelabels = old_bluelabels;
            TotalRedDegrees = old_TotRedDegree;
            TotalBlueDegrees = old_TotBlueDegree;
            IterateFlag = 0;
        end

            
end

Qb_current = QbAfter;


end


function [redlabels,bluelabels,Qb_current]=LOCALMAXIMISATIONPARALLEL(info,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees)

%Find score for current partition
QbAfter = BARBERSMODULARITY(info,redlabels,bluelabels);



IterateFlag = 1;

while IterateFlag == 1
    
%     while QbBefore < QbAfter (this is handled by an if statement below)
        
        %Save old information
        QbBefore = QbAfter;
        old_redlabels = redlabels;
        old_bluelabels = bluelabels;
        old_TotRedDegree = TotalRedDegrees;
        old_TotBlueDegree = TotalBlueDegrees;
        
        %update blue nodes using red node information
        bluelabelchoices = 1:max(redlabels);%unique(redlabels);
                
        for bb = 1:info.col_num
            
            if isnan(bluelabels(bb))==0 %If node is labelled - reduce degree for that label so it can be reassigned
                TotalBlueDegrees(bluelabels(bb)) = TotalBlueDegrees(bluelabels(bb)) - info.col_degrees(bb);
            end
                
            changebluelabeltest=[];
            %check all possible label choices
            parfor ww = 1:length(bluelabelchoices)
                changebluelabeltest(ww) = sum((redlabels == bluelabelchoices(ww)).*info.matrix(:,bb)') - (info.col_degrees(bb)*TotalRedDegrees(bluelabelchoices(ww)) / info.edge_num);
            end
            
            %assign new label based on maximisation of above condition            

            labels = find(changebluelabeltest == max(changebluelabeltest)); %0's and 1's... only want 1's.
            newlabelindex = labels(randi(length(labels)));
            bluelabels(bb) = bluelabelchoices(newlabelindex);
            
            
            
            if bluelabels(bb) > length(TotalBlueDegrees)
                TotalBlueDegrees(bluelabels(bb)) = 0;
            end
            
            %Update total degrees on new labelling
            TotalBlueDegrees(bluelabels(bb)) = TotalBlueDegrees(bluelabels(bb)) + info.col_degrees(bb);

        end

        %Now update red node labels based on blue node information
        redlabelchoices = 1:max(bluelabels);%unique(bluelabels);
        
        for aa = 1:info.row_num
            
           TotalRedDegrees(redlabels(aa)) = TotalRedDegrees(redlabels(aa)) - info.row_degrees(aa);
           
           changeredlabeltest=[];
           

           
           parfor ww = 1:length(redlabelchoices)
               changeredlabeltest(ww) = sum((bluelabels == redlabelchoices(ww)).*info.matrix(aa,:)) - (info.row_degrees(aa)*TotalBlueDegrees(redlabelchoices(ww)) / info.edge_num);
           end
            
           %assign new label based on maximisation of above condition 
           
           labels = find(changeredlabeltest == max(changeredlabeltest)); %0's and 1's... only want 1's.
           newlabelindex = labels(randi(length(labels)));
           redlabels(aa) = redlabelchoices(newlabelindex);
           
           if redlabels(aa) > length(TotalRedDegrees)
                TotalRedDegrees(redlabels(aa)) = 0;
           end
           
           TotalRedDegrees(redlabels(aa)) = TotalRedDegrees(redlabels(aa)) + info.row_degrees(aa);
          
        end
    
        %Find the new modularity score based on node label updates.
        QbAfter = BARBERSMODULARITY(info,redlabels,bluelabels);
        
        %If this modularity is not as good as previous stop iterating and
        %use that previous best information
        if QbAfter <= QbBefore
            redlabels = old_redlabels;
            bluelabels = old_bluelabels;
            TotalRedDegrees = old_TotRedDegree;
            TotalBlueDegrees = old_TotBlueDegree;
            IterateFlag = 0;
        end

            
end

Qb_current = QbAfter;


end

function [divisionsFound]=DIVISION(redlabels,bluelabels)


%Find the set of non-empty communities
divisionsFound = [];


for ii=1:min(max(redlabels),max(bluelabels))
    
        
    if isempty(find(redlabels==ii,1))==0 && isempty(find(bluelabels==ii,1))==0
        
        divisionsFound=[divisionsFound ii];
        
    end
end
    


end

function [redlabels,bluelabels,Qb_current]=StageTwo_LPAbdash(info,redlabels,bluelabels,Qb_current)

[divisionsFound] = DIVISION(redlabels,bluelabels);

IterateFlag = 1;

while IterateFlag == 1
    CombinedDivisionsThisTime = 0;
    
    if length(divisionsFound) > 1  %Can only combine divisions if enough exist
    
        %Pairwise check the divisions found
        for div1check = 1:length(divisionsFound)-1
            for div2check = div1check+1:length(divisionsFound)
                %Make labels of division div1check equal to the labels for div2check
                CHECK_RED = redlabels;
                CHECK_RED(redlabels == divisionsFound(div1check)) = divisionsFound(div2check);
                CHECK_BLUE = bluelabels;
                CHECK_BLUE(bluelabels == divisionsFound(div1check)) = divisionsFound(div2check);
                
                
                QQ = BARBERSMODULARITY(info,CHECK_RED,CHECK_BLUE);
                
                if QQ > Qb_current %If this arrangement looks good
                    FoundBetter = 0;
                    % check against all other divisions to make sure their isn't a better arrangement
                    for aa = 1:length(divisionsFound)
                        CHECK_RED2 = redlabels;
                        CHECK_RED2(redlabels == divisionsFound(aa)) = divisionsFound(div1check);
                        CHECK_BLUE2 = bluelabels;
                        CHECK_BLUE2(bluelabels == divisionsFound(aa)) = divisionsFound(div1check);
                        
                        if BARBERSMODULARITY(info,CHECK_RED2,CHECK_BLUE2) > QQ
                            FoundBetter = 1;
                        end
                        
                        CHECK_RED2 = redlabels;
                        CHECK_RED2(redlabels == divisionsFound(aa)) = divisionsFound(div2check);
                        CHECK_BLUE2 = bluelabels;
                        CHECK_BLUE2(bluelabels == divisionsFound(aa)) = divisionsFound(div2check);
                        
                        if BARBERSMODULARITY(info,CHECK_RED2,CHECK_BLUE2) > QQ
                            FoundBetter = 1;
                        end
                        
                    end
                    
                    if FoundBetter == 0 %If no better configuration found - JOIN.
                        redlabels = CHECK_RED;
                        bluelabels = CHECK_BLUE;
                        CombinedDivisionsThisTime = CombinedDivisionsThisTime + 1;
                    end
                end
            end
        end
        
        if CombinedDivisionsThisTime == 0 %If no divisions were joined move on
            IterateFlag = 0;
        end
        
    else
        %If not enough divisions to join cannot continue
        IterateFlag = 0;
    end
    
    
    [redlabels,bluelabels,Qb_current]=StageOne_LPAbdash(info,redlabels,bluelabels);
    
    [divisionsFound]=DIVISION(redlabels,bluelabels);

    
end




end


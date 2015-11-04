% DIRT_LPA_wb_plus.R
% Label propagation algorithm for weighted bipartite networks that finds modularity.
% Contains the DIRTLPAwb+ algorithm
% Author :  Stephen Beckett ( https://github.com/sjbeckett/weighted-modularity-LPAwbPLUS )
% MIT License

function [Qb_current,redlabels,bluelabels] = DIRT_LPA_wb_plus(MATRIX,mini,reps,parallel)


    %Default settings are invoked when only some inputs are supplied
    if nargin < 4
        parallel = 0;
        if nargin < 3
            reps = 10;
            if nargin < 2
                mini = 4;
            end
        end
    end
    
    
    [A,R,C] = LPA_wb_plus(MATRIX,parallel);
    
    mods = length(unique(R));
    
    if (mods-mini) > 0
        for aa = mini:mods
            for bb = 1:reps
            
                [B,R2,C2] = LPA_wb_plus(MATRIX,parallel,aa);
                
                if B > A %If new modularity is better
                    A = B;
                    R = R2;
                    C = C2;
                end
            end
        end
    end
    
    
    Qb_current = A;
    redlabels = R;
    bluelabels = C;

end

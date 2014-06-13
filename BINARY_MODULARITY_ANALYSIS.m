%A weighted modularity algorithm for bipartite networks. Stephen Beckett. figshare. http://dx.doi.org/10.6084/m9.figshare.999114

function BINARY_MODULARITY_ANALYSIS(MATRIX, nametosave)

%Code to perform modularity analysis on a binary bipartite matrix given as MATRIX using algorithm LPAb+ (Liu & Murata,2010).
%This code finds the modularity of MATRIX and performs significance testing against two binary null models: a Bernoulli random null model and the null model of Bascompte et al., 2003.
%As LPAb+ can provide stochastic output we find the best modularity score from 1000 initialisations each time we run it.
%Output is saved to the name given as a string variable to nametosave.

%Liu X., Murata T. 2010. An Efficient Algorithm for Optimizing Bipartite Modularity in Bipartite Networks. Journal of Advanced Computational Intelligence and Intelligent Informatics (JACIII) 14(4): 408-415.

%Bascompte J., Jordano P., Melián CJ., Olesen JM. 2003. The nested assembly of plant--animal mutualistic networks. Proceedings of the National Academy of Sciences 100(16): 9383-9387.

	


	parallelise = 0; %To use (1) or not use (0) parallel processing to evaluate LPAb+ - this is recommended for large networks.

	%% 1. Find modularity of MATRIX


	Qb = -9; %Initialise Qb with a score lower than possible.
	QT = zeros(1,1000);

	for a = 1:1000 %for this number of initialisations
	    [QT(a),red,blue] = LPA_b_plus(MATRIX,parallelise);
    
	    if QT(a) > Qb
		%save parameters from setup that returns the best modularity score from MATRIX
	        Qb = QT(a);
        	redstore=red;
        	bluestore=blue;
	    end
	end

	%Now can compare modularity score Qb from MATRIX against that of null matrices under different null models.

	%% 2. Significance testing

	%Create variables to store output from ensemble for each null model
	Q1=zeros(1,1000);
	Q2=Q1;


	for c = 1:1000 %number of null models
	    c
	    NULL1 = PROBNULLMODELS(MATRIX,1); %Bernoulli - EE
	    NULL2 = PROBNULLMODELS(MATRIX,2); %DD - Bascompte et al., 2003 null model
    
	    Q1m=-9;
	    Q2m=Q1m;

    
	    for b = 1:1000 %perform a few random initialisations to try and maximise.
        
        	[Q1t,~,~] = LPA_b_plus(NULL1,parallelise);
	        [Q2t,~,~] = LPA_b_plus(NULL2,parallelise);

    
        	if Q1t > Q1m
        	    Q1m=Q1t;
        	end
    
        	if Q2t >Q2m
        	    Q2m=Q2t;
        	end
        

	    end
    
	    Q1(c) = Q1m; %record scores
	    Q2(c) = Q2m;
    
    
	end
    

	%% 3. Assemble statistics from the ensembles
    
	P1=sum(Q1>Qb)/1000;  %p-value
 	MAX1 = max(Q1);
 	MIN1 = min(Q1);
 	M1 = mean(Q1);
 	SDev1 = std(Q1);    %standard deviation
 	Z1 = (Qb - M1 )/SDev1;   % sample z-score

 
	P2=sum(Q2>Qb)/1000;  %p-value
 	MAX2 = max(Q2);
 	MIN2 = min(Q2);
 	M2 = mean(Q2);
 	SDev2 = std(Q2);  %standard deviation
 	Z2 = (Qb - M2 )/SDev2;  % sample z-score

 

	%% 4. Save ouput to file
	 save(nametosave);


end


%%Function for creating null matrices for each of the null models.

function [TEST] = PROBNULLMODELS(MATRIX,option)

	[rows,cols]=size(MATRIX);
	TEST=MATRIX.*0;


	if option ==1
	    %USE Equiprobable - Equiprobable Bernoulli null model    
	    TEST = rand(rows,cols)  < ( sum(MATRIX(:)) / (rows*cols) );

	elseif option ==2	
    
	    %USE Degreeprobable - Degreeprobable null model from Bascompte et al., 2003.
    
	    %Fill up each matrix element probabilistically depending on the matrix dimensions and
	    %degree distribution

	    coldegreesprop=sum(MATRIX>0,1)./rows;
	    rowdegreesprop=sum(MATRIX>0,2)./cols;
    
	    TEST = rand(rows,cols) < 0.5.*( coldegreesprop(ones(rows,1),:) + rowdegreesprop(:,ones(1,cols)) );

	end

end


source("../../code/R/LPA_wb_plus.R")
source("NormalisedMutualInformation.R")




MakeModule <- function(Nrow,Ncol,Nmod,Filling) {

	#Make module labels
	flag=1;
	while(flag==1) {
		Lrows = sort(sample(Nmod,Nrow,replace=T));
		Lcols = sort(sample(Nmod,Ncol,replace=T));

		if(length(unique(Lrows))==Nmod & length(unique(Lcols))==Nmod)
			flag=0;
		}


	#Create module matrix
	A = matrix(0,Nrow,Ncol);
	for(a in 1:Nrow){
		for(b in 1:Ncol){
			if(Lrows[a] == Lcols[b])
				A[a,b]=1
		}
	}
	

	#Put in Quantitative links
	DistSize=100000 # make a large negative binomial distribution
	W = A #copy A to fill with weighted links
	dist = rnbinom(DistSize,size=Filling,mu=4) #negative binomial distribution

	flag = 1
	while(flag==1) {	
		want = sample(dist,sum(A>0)) # sample from distribution to replace 1's in A
	
		W[A>0] = want
		image(W)

		if( length(which(rowSums(W)==0))==0 && length(which(colSums(W)==0))==0 ) { # check no zero rows/columns
			flag=0
		}
	}

	
	return(list(matrix = W,row_labels = Lrows,col_labels = Lcols))
}


RewireConnections <- function(Matrix,Threshold) {
	#Threshold E [0,1]

	#removal of some cells
	A = runif(sum(Matrix>0))<Threshold
	B = Matrix[Matrix>0]
	C = B[A]
	B[A] = 0
	Matrix[Matrix>0] = B

	#readdition
	D = which(Matrix==0)
	E = sample(D, length(C))
	Matrix[E] = C

	return(Matrix)
}

PERFECT_MOD <- function(AA) {
	#find modularity according to perfect module layout
	return(WEIGHTEDMODULARITY2(BarbersMatrix(AA$matrix),sum(AA$matrix),AA$row_labels,AA$col_labels))
}


COMPARENMI <- function(AA,BB) {
	#find difference in NMI
	return( NormalisedMutualInformation(c(AA$row_labels,AA$col_labels),c(BB$Row_labels,BB$Col_labels)) )
}

COMPAREMODS <- function(AA,BB) {
	#find ratio of module number FOUND/PERFECT
	return( length(unique(BB$Row_labels))/length(unique(AA$row_labels)) )
}






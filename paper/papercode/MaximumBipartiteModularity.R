MAXIMUM_BIP_MOD <- function(MAT,red,blue) {

	A= rowSums(MAT)
	B= colSums(MAT)
	m = sum(MAT)

	holdsum = 0
	for( aa in 1:nrow(MAT)) {
		for( bb in 1:ncol(MAT)) {
			if(red[aa] == blue[bb])
				holdsum = holdsum + A[[aa]]*B[[bb]]/m
		}
	}

	Qmax =	(m - holdsum)/m

	return(Qmax)
}

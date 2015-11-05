# NegativeBinomialTest

#Find out how many non-zero elements will be returned from the negative binomial distribution if size is varied for fixed mu

DistSize=100000
sizes = c(seq(0.001,2,0.01),seq(2,50,0.5))
PropNonZero=c()

for(aaa in 1:length(sizes)){

	dist = rnbinom(DistSize,size=sizes[aaa],mu=4) #negative binomial distribution
	PropNonZero[aaa] = sum(dist>0)/DistSize
}

#Large scale plot
#plot(sizes,PropNonZero,ylim=c(0,1),type='l')

#Small scale plot

plot(sizes[1:205],PropNonZero[1:205],ylim=c(0,1),type='l',xlab="Size",ylab="Proportion non-zero",main=expression(paste(mu ,"= 4")))
lines(c(2.5,2.5),c(0,1),col="red") # ~ 90% fill    HF
lines(c(0.5,0.5),c(0,1),col="red") # ~ 65% fill    LF
dev.copy2eps(file="negativeBinomialSizevariation.eps")
dev.off()


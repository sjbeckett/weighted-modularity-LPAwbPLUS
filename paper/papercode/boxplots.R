 


NOISEBOXPLOT <- function(Hm1,Hm2,Hm3,Hm4,Hm5,Hm6,Hm7,Hm8,ylabel,legon,plottype) {
colw="grey"
if(plottype==1) { #module number ratio
	YLIMITS = c(0,14)
	ATY = 0:14
}else if(plottype==2) { #modularity ratio
	YLIMITS = c(0.3,1.1)
	ATY = seq(0.3,1.1,0.1)
} else { #NMI
	YLIMITS = c(0,1)
	ATY = seq(0,1,0.2)
}


boxplot(Hm1,Hm2,Hm3,Hm4,Hm5,Hm6,Hm7,Hm8, boxwex=0.3, col=c("#99cd99","#ff99cc"),axes=FALSE,frame.plot=FALSE,at=c(0.84,1.16,1.84,2.16,2.84,3.16,3.84,4.16),xlim=c(0.7,4.3),xlab="Level of noise",ylab=ylabel,boxlwd=0.5,ylim=YLIMITS)

axis(1, at=seq(1,4,1),labels=c(0,0.01,0.25,0.5),col=colw,lwd=2)
axis(2,col=colw,lwd=2,at=ATY)

if(legon==1)
	legend(2,0.1,c("LPAwb+","DIRTLPAwb+"),fill=c("#99cd99","#ff99cc"),border=c("#99cd99","#ff99cc"),bty="n")

}



MODULESANDFILLBOXPLOT <- function(Hm1,Hm2,Hm3,Hm4,Hm5,Hm6,Hm7,Hm8,ylabel,legon,plottype) {
colw="grey"
if(plottype==1) { #module number ratio
	YLIMITS = c(0,14)
	ATY = 0:14
}else if(plottype==2) { #modularity ratio
	YLIMITS = c(0.3,1.1)
	ATY = seq(0.3,1.1,0.1)
} else { #NMI
	YLIMITS = c(0,1)
	ATY = seq(0,1,0.2)
}

boxplot(Hm1,Hm2,Hm3,Hm4,Hm5,Hm6,Hm7,Hm8, boxwex=0.3, col=c("#99cd99","#ff99cc"),axes=FALSE,frame.plot=FALSE,at=c(0.84,1.16,1.84,2.16,2.84,3.16,3.84,4.16),xlim=c(0.7,4.3),xlab="Network type",ylab=ylabel,boxlwd=0.5,ylim=YLIMITS)
axis(1, at=seq(1,4,1),labels=c("LMLF","LMHF","HMLF","HMHF"),col=colw,lwd=2)
axis(2,col=colw,lwd=2,at =ATY)


if(legon==1)
	legend(2,0.1,c("LPAwb+","DIRTLPAwb+"),fill=c("#99cd99","#ff99cc"),border=c("#99cd99","#ff99cc"),bty="n")


}






Ind_Rewire_0 = c()
Ind_Rewire_0_01 = c()
Ind_Rewire_0_25 = c()
Ind_Rewire_0_5 = c()

index = 1

	for(aa in 1:10) {
		for(bb in 1:4) {
			for(cc in 1:5) {

				if(bb == 1) {
					Ind_Rewire_0 = c(Ind_Rewire_0,index)
				} else if(bb==2) {
					Ind_Rewire_0_01 = c(Ind_Rewire_0_01,index)
				} else if(bb==3) {
					Ind_Rewire_0_25 = c(Ind_Rewire_0_25,index)
				} else {
					Ind_Rewire_0_5 = c(Ind_Rewire_0_5,index)
				}

				index=index+1
				print(index)
			}
		}
	}







setEPS()
postscript("NoiseAssessment.eps",height=5.5,width=9)

par(mfrow=c(1,3))
Hm1 = c(DET_TO_PERF_MODULES_LPA_LMLF[Ind_Rewire_0],DET_TO_PERF_MODULES_LPA_LMHF[Ind_Rewire_0],DET_TO_PERF_MODULES_LPA_HMLF[Ind_Rewire_0],DET_TO_PERF_MODULES_LPA_HMHF[Ind_Rewire_0])
Hm2 = c(DET_TO_PERF_MODULES_EXLPA_LMLF[Ind_Rewire_0],DET_TO_PERF_MODULES_EXLPA_LMHF[Ind_Rewire_0],DET_TO_PERF_MODULES_EXLPA_HMLF[Ind_Rewire_0],DET_TO_PERF_MODULES_EXLPA_HMHF[Ind_Rewire_0])
Hm3 = c(DET_TO_PERF_MODULES_LPA_LMLF[Ind_Rewire_0_01],DET_TO_PERF_MODULES_LPA_LMHF[Ind_Rewire_0_01],DET_TO_PERF_MODULES_LPA_HMLF[Ind_Rewire_0_01],DET_TO_PERF_MODULES_LPA_HMHF[Ind_Rewire_0_01])
Hm4 = c(DET_TO_PERF_MODULES_EXLPA_LMLF[Ind_Rewire_0_01],DET_TO_PERF_MODULES_EXLPA_LMHF[Ind_Rewire_0_01],DET_TO_PERF_MODULES_EXLPA_HMLF[Ind_Rewire_0_01],DET_TO_PERF_MODULES_EXLPA_HMHF[Ind_Rewire_0_01])
Hm5 = c(DET_TO_PERF_MODULES_LPA_LMLF[Ind_Rewire_0_25],DET_TO_PERF_MODULES_LPA_LMHF[Ind_Rewire_0_25],DET_TO_PERF_MODULES_LPA_HMLF[Ind_Rewire_0_25],DET_TO_PERF_MODULES_LPA_HMHF[Ind_Rewire_0_25])
Hm6 = c(DET_TO_PERF_MODULES_EXLPA_LMLF[Ind_Rewire_0_25],DET_TO_PERF_MODULES_EXLPA_LMHF[Ind_Rewire_0_25],DET_TO_PERF_MODULES_EXLPA_HMLF[Ind_Rewire_0_25],DET_TO_PERF_MODULES_EXLPA_HMHF[Ind_Rewire_0_25])
Hm7 = c(DET_TO_PERF_MODULES_LPA_LMLF[Ind_Rewire_0_5],DET_TO_PERF_MODULES_LPA_LMHF[Ind_Rewire_0_5],DET_TO_PERF_MODULES_LPA_HMLF[Ind_Rewire_0_5],DET_TO_PERF_MODULES_LPA_HMHF[Ind_Rewire_0_5])
Hm8 = c(DET_TO_PERF_MODULES_EXLPA_LMLF[Ind_Rewire_0_5],DET_TO_PERF_MODULES_EXLPA_LMHF[Ind_Rewire_0_5],DET_TO_PERF_MODULES_EXLPA_HMLF[Ind_Rewire_0_5],DET_TO_PERF_MODULES_EXLPA_HMHF[Ind_Rewire_0_5])
NOISEBOXPLOT(Hm1,Hm2,Hm3,Hm4,Hm5,Hm6,Hm7,Hm8,"Module number ratio",0,1)
abline(h=1,lty=3,lwd=4,col="grey80")
Hm1 = c(DET_TO_PERF_MODULARITY_LPA_LMLF[Ind_Rewire_0],DET_TO_PERF_MODULARITY_LPA_LMHF[Ind_Rewire_0],DET_TO_PERF_MODULARITY_LPA_HMLF[Ind_Rewire_0],DET_TO_PERF_MODULARITY_LPA_HMHF[Ind_Rewire_0])
Hm2 = c(DET_TO_PERF_MODULARITY_EXLPA_LMLF[Ind_Rewire_0],DET_TO_PERF_MODULARITY_EXLPA_LMHF[Ind_Rewire_0],DET_TO_PERF_MODULARITY_EXLPA_HMLF[Ind_Rewire_0],DET_TO_PERF_MODULARITY_EXLPA_HMHF[Ind_Rewire_0])
Hm3 = c(DET_TO_PERF_MODULARITY_LPA_LMLF[Ind_Rewire_0_01],DET_TO_PERF_MODULARITY_LPA_LMHF[Ind_Rewire_0_01],DET_TO_PERF_MODULARITY_LPA_HMLF[Ind_Rewire_0_01],DET_TO_PERF_MODULARITY_LPA_HMHF[Ind_Rewire_0_01])
Hm4 = c(DET_TO_PERF_MODULARITY_EXLPA_LMLF[Ind_Rewire_0_01],DET_TO_PERF_MODULARITY_EXLPA_LMHF[Ind_Rewire_0_01],DET_TO_PERF_MODULARITY_EXLPA_HMLF[Ind_Rewire_0_01],DET_TO_PERF_MODULARITY_EXLPA_HMHF[Ind_Rewire_0_01])
Hm5 = c(DET_TO_PERF_MODULARITY_LPA_LMLF[Ind_Rewire_0_25],DET_TO_PERF_MODULARITY_LPA_LMHF[Ind_Rewire_0_25],DET_TO_PERF_MODULARITY_LPA_HMLF[Ind_Rewire_0_25],DET_TO_PERF_MODULARITY_LPA_HMHF[Ind_Rewire_0_25])
Hm6 = c(DET_TO_PERF_MODULARITY_EXLPA_LMLF[Ind_Rewire_0_25],DET_TO_PERF_MODULARITY_EXLPA_LMHF[Ind_Rewire_0_25],DET_TO_PERF_MODULARITY_EXLPA_HMLF[Ind_Rewire_0_25],DET_TO_PERF_MODULARITY_EXLPA_HMHF[Ind_Rewire_0_25])
Hm7 = c(DET_TO_PERF_MODULARITY_LPA_LMLF[Ind_Rewire_0_5],DET_TO_PERF_MODULARITY_LPA_LMHF[Ind_Rewire_0_5],DET_TO_PERF_MODULARITY_LPA_HMLF[Ind_Rewire_0_5],DET_TO_PERF_MODULARITY_LPA_HMHF[Ind_Rewire_0_5])
Hm8 = c(DET_TO_PERF_MODULARITY_EXLPA_LMLF[Ind_Rewire_0_5],DET_TO_PERF_MODULARITY_EXLPA_LMHF[Ind_Rewire_0_5],DET_TO_PERF_MODULARITY_EXLPA_HMLF[Ind_Rewire_0_5],DET_TO_PERF_MODULARITY_EXLPA_HMHF[Ind_Rewire_0_5])
NOISEBOXPLOT(Hm1,Hm2,Hm3,Hm4,Hm5,Hm6,Hm7,Hm8,"Modularity ratio",0,2)
abline(h=1,lty=3,lwd=4,col="grey80")
Hm1 = c(DET_TO_PERF_NMI_LPA_LMLF[Ind_Rewire_0],DET_TO_PERF_NMI_LPA_LMHF[Ind_Rewire_0],DET_TO_PERF_NMI_LPA_HMLF[Ind_Rewire_0],DET_TO_PERF_NMI_LPA_HMHF[Ind_Rewire_0])
Hm2 = c(DET_TO_PERF_NMI_EXLPA_LMLF[Ind_Rewire_0],DET_TO_PERF_NMI_EXLPA_LMHF[Ind_Rewire_0],DET_TO_PERF_NMI_EXLPA_HMLF[Ind_Rewire_0],DET_TO_PERF_NMI_EXLPA_HMHF[Ind_Rewire_0])
Hm3 = c(DET_TO_PERF_NMI_LPA_LMLF[Ind_Rewire_0_01],DET_TO_PERF_NMI_LPA_LMHF[Ind_Rewire_0_01],DET_TO_PERF_NMI_LPA_HMLF[Ind_Rewire_0_01],DET_TO_PERF_NMI_LPA_HMHF[Ind_Rewire_0_01])
Hm4 = c(DET_TO_PERF_NMI_EXLPA_LMLF[Ind_Rewire_0_01],DET_TO_PERF_NMI_EXLPA_LMHF[Ind_Rewire_0_01],DET_TO_PERF_NMI_EXLPA_HMLF[Ind_Rewire_0_01],DET_TO_PERF_NMI_EXLPA_HMHF[Ind_Rewire_0_01])
Hm5 = c(DET_TO_PERF_NMI_LPA_LMLF[Ind_Rewire_0_25],DET_TO_PERF_NMI_LPA_LMHF[Ind_Rewire_0_25],DET_TO_PERF_NMI_LPA_HMLF[Ind_Rewire_0_25],DET_TO_PERF_NMI_LPA_HMHF[Ind_Rewire_0_25])
Hm6 = c(DET_TO_PERF_NMI_EXLPA_LMLF[Ind_Rewire_0_25],DET_TO_PERF_NMI_EXLPA_LMHF[Ind_Rewire_0_25],DET_TO_PERF_NMI_EXLPA_HMLF[Ind_Rewire_0_25],DET_TO_PERF_NMI_EXLPA_HMHF[Ind_Rewire_0_25])
Hm7 = c(DET_TO_PERF_NMI_LPA_LMLF[Ind_Rewire_0_5],DET_TO_PERF_NMI_LPA_LMHF[Ind_Rewire_0_5],DET_TO_PERF_NMI_LPA_HMLF[Ind_Rewire_0_5],DET_TO_PERF_NMI_LPA_HMHF[Ind_Rewire_0_5])
Hm8 = c(DET_TO_PERF_NMI_EXLPA_LMLF[Ind_Rewire_0_5],DET_TO_PERF_NMI_EXLPA_LMHF[Ind_Rewire_0_5],DET_TO_PERF_NMI_EXLPA_HMLF[Ind_Rewire_0_5],DET_TO_PERF_NMI_EXLPA_HMHF[Ind_Rewire_0_5])
NOISEBOXPLOT(Hm1,Hm2,Hm3,Hm4,Hm5,Hm6,Hm7,Hm8,"NMI",1,3)

dev.off()




setEPS()
postscript("ModuleFillingAssessment.eps",height=5.5,width=9)

par(mfrow=c(1,3))
Hm1 = DET_TO_PERF_MODULES_LPA_LMLF
Hm2 = DET_TO_PERF_MODULES_EXLPA_LMLF
Hm3 = DET_TO_PERF_MODULES_LPA_LMHF
Hm4 = DET_TO_PERF_MODULES_EXLPA_LMHF
Hm5 = DET_TO_PERF_MODULES_LPA_HMLF
Hm6 = DET_TO_PERF_MODULES_EXLPA_HMLF
Hm7 = DET_TO_PERF_MODULES_LPA_HMHF
Hm8 = DET_TO_PERF_MODULES_EXLPA_HMHF
MODULESANDFILLBOXPLOT(Hm1,Hm2,Hm3,Hm4,Hm5,Hm6,Hm7,Hm8,"Module number ratio",0,1)
abline(h=1,lty=3,lwd=4,col="grey80")
Hm1 = DET_TO_PERF_MODULARITY_LPA_LMLF
Hm2 = DET_TO_PERF_MODULARITY_EXLPA_LMLF
Hm3 = DET_TO_PERF_MODULARITY_LPA_LMHF
Hm4 = DET_TO_PERF_MODULARITY_EXLPA_LMHF
Hm5 = DET_TO_PERF_MODULARITY_LPA_HMLF
Hm6 = DET_TO_PERF_MODULARITY_EXLPA_HMLF
Hm7 = DET_TO_PERF_MODULARITY_LPA_HMHF
Hm8 = DET_TO_PERF_MODULARITY_EXLPA_HMHF
MODULESANDFILLBOXPLOT(Hm1,Hm2,Hm3,Hm4,Hm5,Hm6,Hm7,Hm8,"Modularity ratio",0,2)
abline(h=1,lty=3,lwd=4,col="grey80")
Hm1 = DET_TO_PERF_NMI_LPA_LMLF
Hm2 = DET_TO_PERF_NMI_EXLPA_LMLF
Hm3 = DET_TO_PERF_NMI_LPA_LMHF
Hm4 = DET_TO_PERF_NMI_EXLPA_LMHF
Hm5 = DET_TO_PERF_NMI_LPA_HMLF
Hm6 = DET_TO_PERF_NMI_EXLPA_HMLF
Hm7 = DET_TO_PERF_NMI_LPA_HMHF
Hm8 = DET_TO_PERF_NMI_EXLPA_HMHF
MODULESANDFILLBOXPLOT(Hm1,Hm2,Hm3,Hm4,Hm5,Hm6,Hm7,Hm8,"NMI",1,3)

dev.off()

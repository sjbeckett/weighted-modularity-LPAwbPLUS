LPAwb+ and DIRTLPAwb+ in R
=============================


```r
source("LPA_wb_plus.R")  #read in functions

MAT = matrix(sample(0:3,20*20,replace=TRUE),20,20) # create an example matrix

#example scripts

MOD1 = LPA_wb_plus(MAT) # find labels and weighted modularity using LPAwb+

MOD2 = DIRT_LPA_wb_plus(MAT) # find labels and weighted modularity using DIRTLPAwb+

MOD3 = DIRT_LPA_wb_plus(MAT>0, 2, 20) # find labels and binary modularity using DIRTLPAwb+ checking from a minimum of 2 modules and 20 replicates

#plotting

source("MODULARPLOT.R") #read in plotting function
MODULARPLOT(MAT,MOD1) # show the modular network configuration found in MOD1. Row and column numbering indicates the ordering of rows and columns in MAT. Modules are highlighted in red rectangles.



#use with R library 'bipartite'
library("bipartite")
source("convert2moduleWeb.R") # read in conversion function 

MOD1modWeb = convert2moduleWeb(MAT,MOD1) # converts the configuration found in MOD1 to a moduleWeb object
plotModuleWeb(MOD1modWeb) # plot the corresponding moduleWeb object using plotModuleWeb in library bipartite


```


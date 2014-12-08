Data and code
=============================

This folder contains the code used to create the material for the paper and the resulting figures.
The folder output contains output data files from running analyses:
	- the row and column module label assignments for each call to the algorithms for each network
	- the time it took to execute each algorithm on the networks
	- the modularity score found from each algorithm call
	- summary data based on the above
	- the normalised mutual information between the best binary and weighted partitions
	- realized modularity scores of the best partitions

The workflow for running the weighted modularity algorithms and interrogating the resulting data can be found in **workflow.R**
This is the best place to start!

The input networks were directly called from the bipartite package in R.


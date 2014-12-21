weighted-modularity-LPAwbPLUS
=============================

Two algorithms for finding modularity in bipartite networks: LPAwb+ and Exhaustive LPAwb+.

These are based on the LPAb+ algorithm of Liu & Murata, 2010 for use on bipartite/two-mode networks. The algorithm has been modified such that the weighted modularity of a network can be found (note that the result is equivalent to that found by the LPAb+ algorithm if the network is binary). Code is currently available for Julia, MATLAB/Octave and R langauges.

For details of the methods please [view the draft paper](https://github.com/sjbeckett/weighted-modularity-LPAwbPLUS/blob/master/paper/weightedModularityDraft.pdf?raw=true) (in the 'paper' directory).


Usage
---------

Language specific instructions are included with the code.
Source the relevant code in your favourite language; then run


```julia

LPA_wb_plus(MATRIX) # find labels and weighted modularity using LPAwb+


```
where MATRIX is the incidence/biadjacency matrix describing the input network. Three outputs are returned: `redlabels` - the module labels for each row in the input matrix, `bluelabels` - the module labels for each row in the input matrix and `Q` - the modularity score.




Details
---------

For details of the LPAb+ algorithm please see [Liu & Murata, 2010](https://www.fujipress.jp/finder/xslt.php?mode=present&inputfile=JACII001400040010.xml) .

Liu X., Murata T. 2010. An Efficient Algorithm for Optimizing Bipartite Modularity in Bipartite Networks. Journal of Advanced Computational Intelligence and Intelligent Informatics (JACIII) 14(4): 408-415.

Further details for LPAwb+ are being written in (Beckett, unsubmitted)


Contact
--------

This package is authored by Stephen Beckett ([@BeckettStephen](https://twitter.com/BeckettStephen)) at the University of Exeter. If you have any questions please contact me!

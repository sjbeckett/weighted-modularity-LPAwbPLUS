Change notes
---------

v1.1 (November 19th 2015) - Revision of code

* Updated Stage Two of the LPAwb+ algorithm. Previously a block of code was not being called as it was checking the modularity against itself (within R and Julia). This has been updated for all langauges.
* To test the importance of this code block an experiment to test the effect has been added to papercode/code/StageTwoExperiment
* To reflect the change in the code, all previous results have been recomputed. (The previous results computed in v1.01 are archived at http://dx.doi.org/10.5281/zenodo.19585 )
* Exhaustive LPAwb+ has been renamed as DIRTLPAwb+ (Different Initialisations with Repeated Trials of LPAwb+), the code and paper has been updated to reflect this change.
* An error was discovered in the creation of the synthetic ensembles, that resulted in matrices being computed that had columns with no positive values within them (nodes that are unconnected from the rest of the network). This was corrected and new values for the parameterisation of the negative binomial distribution are used in the new synthetic ensemble.
* Various changes have been made to the manuscript to improve its clarity

v1.01 (July 7th 2015) - update of the DOI.

v1.0 (July 6th 2015) - First release of code


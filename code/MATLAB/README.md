LPAwb+ and DIRTLPAwb+ in MATLAB
=============================


```matlab


MAT = randi(4,20,20)-1 % create an example matrix

%example scripts

[MOD1, R1, C1] = LPA_wb_plus(MAT) % find labels and weighted modularity using LPAwb+

[MOD2, R2, C2] = DIRT_LPA_wb_plus(MAT) % find labels and weighted modularity using DIRTLPAwb+

[MOD3, R3, C3] = DIRT_LPA_wb_plus(MAT>0, 2, 20) % find labels and binary modularity using DIRTLPAwb+ checking from a minimum of 2 modules and 20 replicates

%plotting

MODULARPLOT(MAT,R1,C1) % show the modular network configuration found in MOD1. Row and column numbering indicates the ordering of rows and columns in MAT. Modules are highlighted in red rectangles.


```


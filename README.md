ExaML
=====

Exascale Maximum Likelihood (ExaML) code for phylogenetic inference using MPI.

This code implements the popular RAxML search algorithm for maximum likelihood based inference 
of phylogenetic trees.

It uses a radically new MPI parallelization approach that yields improved parallel efficiency, 
in particular on partitioned multi-gene or whole-genome datasets.

When using ExaML please cite the following paper: 

Alexey M. Kozlov, Andre J. Aberer, Alexandros Stamatakis: "ExaML Version 3: A Tool for Phylogenomic Analyses on Supercomputers." Bioinformatics (2015) 31 (15): 2577-2579.

It is up to 4 times faster than RAxML-Light [1].

As RAxML-Light, ExaML also implements checkpointing, SSE3, AVX vectorization and 
memory saving techniques.

[1] A. Stamatakis,  A.J. Aberer, C. Goll, S.A. Smith, S.A. Berger, F. Izquierdo-Carrasco: 
    "RAxML-Light: A Tool for computing TeraByte Phylogenies", 
    Bioinformatics 2012; doi: 10.1093/bioinformatics/bts309.


Intel Xeon Phi
--------------

For details on running ExaML on Intel Xeon Phi, please refer to:

    Knights Landing: [README_KNL.txt]

    Knights Corner: [README_MIC.txt]



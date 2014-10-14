ExaML
=====

Exascale Maximum Likelihood (ExaML) code for phylogenetic inference using MPI.

This code implements the popular RAxML search algorithm for maximum likelihood based inference 
of phylogenetic trees.
It uses a radically new MPI parallelization approach that yields improved parallel efficiency, 
in particular on partitioned multi-gene or whole-genome datasets.

It is up to 4 times faster than RAxML-Light [1].

As RAxML-Light, ExaML also implements checkpointing, SSE3, AVX vectorization and 
memory saving techniques.

[1] A. Stamatakis,  A.J. Aberer, C. Goll, S.A. Smith, S.A. Berger, F. Izquierdo-Carrasco: 
    "RAxML-Light: A Tool for computing TeraByte Phylogenies", 
    Bioinformatics 2012; doi: 10.1093/bioinformatics/bts309.

Please refer to README_MIC.txt For details on running ExaML on Intel MIC (aka Xeon Phi).
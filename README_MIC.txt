Using ExaML on the Intel MIC/Intel Xeon Phi coprocessors

Compiling under Linux
---------------------

Please set your MPI/MIC environment (ask your sysadmin if unsure) and then run:

   make -f Makefile.AVX.gcc
   make -f Makefile.MIC.icc clean
   make -f Makefile.MIC.icc

This will create two executables for both host(=CPU) and MIC - they will be 
named examl-AVX and examl-MIC, respectively.


Running
----------------------

1. Use parse-examl to generate a binary alignment file as usual.

2. You might want to allocate MPI ranks on both host CPUs and MICs (hybrid mode)
or just on the MICs, depending on your configuration.

Sample command line for running ExaML in hybrid mode (16 CPU core + 2 MIC cards):

  mpiexec -host myhost-ib -n 16 /scratch/examl-AVX -n mictest -s /scratch/mictest.binary -t /scratch/start.tre -m GAMMA -w /scratch : \
          -host myhost-mic0 -n 30 -env OMP_NUM_THREADS 4 -env KMP_AFFINITY "granularity=fine,balanced" /scratch/examl-MIC -n mictest \
          -s /scratch/mictest.binary -t /scratch/start.tre -m GAMMA -w /scratch : \
          -host myhost-mic1 -n 30 -env OMP_NUM_THREADS 4 -env KMP_AFFINITY "granularity=fine,balanced" /scratch/examl-MIC -n mictest \
          -s /scratch/mictest.binary -t /scratch/start.tre -m GAMMA -w /scratch

Here, we use 1 MPI rank per core on the host CPUs. On each MIC, we start 30 ranks x 4 OpenMP threads, 
which gives 120 threads in total or 2 threads per MIC core. Changing the ratio of CPU:MIC ranks allows
to fine-tune load balance for the specific hardware configuration at hand.


Limitations & caveats
---------------------

1. Supported on the MIC:

   + DNA and AA alignments
   + GAMMA model of rate heterogeneity
   + multiple partitions 
   + all AA substitution matrices supported by ExaML, including LG4

2. Currently NOT supported:

   - binary and generic multi-state alignments
   - PSR model
   - memory saving for gappy alignments (-S option)

3. Memory 

  Compared to traditional CPUs, MIC cards have significantly lower memory-per-core value,
  which poses a problem for memory-intensive ML computations. Thus you should plan carefully 
  and split your run over multiple cards, if needed.

  To estimate memory requirements for your dataset, you can use the web-calculator here:

    http://sco.h-its.org/exelixis/web/software/raxml/index.html#memcalc

  A similar tool tailored for MICs is coming soon, stay tuned :)

4. Performance
  
  ExaML-MIC performs best on alignments with large number of sites and few taxa.
  The latter is due to the limited on-card memory of the MICs (s. above), so you 
  might need to use multiple cards if the number of taxa is large.

  For details, please refer to: http://www.hicomb.org/papers/HICOMB2014-04.pdf
  
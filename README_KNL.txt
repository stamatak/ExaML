Using ExaML on the Intel Xeon Phi (Knights Landing) coprocessors

Compiling under Linux
---------------------

Please set your MPI/MIC environment (ask your sysadmin if unsure) and then run:

   cd examl
   make -f Makefile.KNL.icc

This will create an executable named examl-KNL. 


Running
----------------------

1. Use parse-examl to generate a binary alignment file as usual

2. You can run examl-KNL in OpenMP-only, MPI-only or in hybrid OpenMP/MPI mode 
   (s. sample commands below). Unlike with KNC, I didn't notice any significant performance
   differences between those three configurations, so just choose whatever is easier/more convenient:

   OpenMP:

     OMP_NUM_THREADS=128 ./examl-KNL -s myTest.binary -m GAMMA -t myStart.tre -n myTest

   MPI:

     mpirun -n 128 -env OMP_NUM_THREADS 1 ./examl-KNL -s myTest.binary -m GAMMA -t myStart.tre -n myTest

   Hybrid:

     mpirun -n 8 -env OMP_NUM_THREADS 16 ./examl-KNL -s myTest.binary -m GAMMA -t myStart.tre -n myTest
   
  NOTE: although KNL has 4 logical threads/core, it usually doesn't make sense to use more then 2 threads/core
        with ExaML (since ExaML doesn't benefit for hyper-threading). Furthermore, please consider the general 
        recommendations regarding number of alignment patterns per core given in the ExaML manual.

3. IMPORTANT: KNL on-card memory can be configured in one of two modes: "Flat" or "Cache".
   You should find out which one is set on you card(s), since it has important performance implications.

   - in "Cache" mode, no 

   - in "Flat" mode, you should use numactl to explicitly bind ExaML process to the fast memory NUMA domain:

        numactl --membind=1 mpirun -n 128 -env OMP_NUM_THREADS 1 ./examl-KNL -s myTest.binary -m GAMMA -t myStart.tre -n myTest

     Obviously, it is not possible if memory requirements for your analysis exceed the on-card memory size (typically 16GB). 
     In this case, you should either switch into "Cache" mode, or let ExaML run in the (slow) main memory. The latter option
     will typically induce a huge performance penalty (up to 5x), and thus is not recommended.

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

3. Performance
  
  ExaML-MIC performs best on alignments with large number of sites and few taxa.
  The latter is due to the limited on-card memory of the MICs (s. above), so you 
  might need to use multiple cards if the number of taxa is large.

  For details, please refer to: http://www.hicomb.org/papers/HICOMB2014-04.pdf and
                                https://doi.org/10.1093/bioinformatics/btv184


Contact & Support
--------------------

Please use RAxML google group to ask questions:

https://groups.google.com/forum/?hl=en#!forum/raxml


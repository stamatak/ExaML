/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses
 *  with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef _AXML_H
#define _AXML_H


#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include "../versionHeader/version.h"

#ifdef __MIC_NATIVE
#define BYTE_ALIGNMENT 64
#define VECTOR_PADDING 8
#else
#define BYTE_ALIGNMENT 32
#define VECTOR_PADDING 1
#endif

#define GET_PADDED_WIDTH(w) w % VECTOR_PADDING == 0 ? w : w + (VECTOR_PADDING - (w % VECTOR_PADDING))

#include <mpi.h>

#ifdef _USE_OMP
#include "omp.h"
#endif

/* BEGIN: file streams */
#ifdef _GNU_SOURCE

/* notice, that the gnu source macro implies posix compliance */


/* these are posix compliant functions. They potentially work on files
   larger than 2 GB (for gcc this can be ensured using the following
   macro) */
#define exa_fseek fseeko
#define exa_ftell ftello 
#define exa_off_t off_t

/* only usefull for ftello/fseeko: ensure that we are using 64-bit
   types for representing an offset */
#define _FILE_OFFSET_BITS 64 

#else 

#define exa_fseek fseek 
#define exa_ftell ftell 
#define exa_off_t long 

#endif
/* END: file streams  */


#define MAX_TIP_EV     0.999999999 /* max tip vector value, sum of EVs needs to be smaller than 1.0, otherwise the numerics break down */
#define smoothings     32          /* maximum smoothing passes through tree */
#define iterations     10          /* maximum iterations of iterations per insert */
#define newzpercycle   1           /* iterations of makenewz per tree traversal */
#define nmlngth        256         /* number of characters in species name */
#define deltaz         0.00001     /* test of net branch length change in update */
#define defaultz       0.9         /* value of z assigned as starting point */
#define unlikely       -1.0E300    /* low likelihood for initialization */


#define AUTO_ML   0
#define AUTO_BIC  1
#define AUTO_AIC  2
#define AUTO_AICC 3

#define SUMMARIZE_LENGTH -3
#define SUMMARIZE_LH     -2
#define NO_BRANCHES      -1

#define MASK_LENGTH 32
#define GET_BITVECTOR_LENGTH(x) ((x % MASK_LENGTH) ? (x / MASK_LENGTH + 1) : (x / MASK_LENGTH))

#define zmin       1.0E-15  /* max branch prop. to -log(zmin) (= 34) */
#define zmax (1.0 - 1.0E-6) /* min branch prop. to 1.0-zmax (= 1.0E-6) */

#define twotothe256  \
  115792089237316195423570985008687907853269984665640564039457584007913129639936.0  
                                                     /*  2**256 (exactly)  */

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood




/* 18446744073709551616.0 */

/*4294967296.0*/

/* 18446744073709551616.0 */

/*  2**64 (exactly)  */
/* 4294967296 2**32 */

#define badRear         -1

#define NUM_BRANCHES     256

#define TRUE             1
#define FALSE            0



#define LIKELIHOOD_EPSILON 0.0000001

#define AA_SCALE 10.0
#define AA_SCALE_PLUS_EPSILON 10.001

/* ALPHA_MIN is critical -> numerical instability, eg for 4 discrete rate cats                    */
/* and alpha = 0.01 the lowest rate r_0 is                                                        */
/* 0.00000000000000000000000000000000000000000000000000000000000034878079110511010487             */
/* which leads to numerical problems Table for alpha settings below:                              */
/*                                                                                                */
/* 0.010000 0.00000000000000000000000000000000000000000000000000000000000034878079110511010487    */
/* 0.010000 yielded nasty numerical bugs in at least one case !                                   */
/* 0.020000 0.00000000000000000000000000000044136090435925743185910935350715027016962154188875    */
/* 0.030000 0.00000000000000000000476844846859006690412039180149775802624789852441798419292220    */
/* 0.040000 0.00000000000000049522423236954066431210260930029681736928018820007024736185030633    */
/* 0.050000 0.00000000000050625351310359203371872643495343928538368616365517027588794007897377    */
/* 0.060000 0.00000000005134625283884191118711474021861409372524676086868566926568746566772461    */
/* 0.070000 0.00000000139080650074206434685544624965062437960128249869740102440118789672851562    */
/* 0.080000 0.00000001650681201563587066858709818343436959153791576682124286890029907226562500    */
/* 0.090000 0.00000011301977332931251259273962858978301859735893231118097901344299316406250000    */
/* 0.100000 0.00000052651925834844387815526344648331402709118265192955732345581054687500000000    */


#define ALPHA_MIN    0.02
#define ALPHA_MAX    1000.0

#define RATE_MIN     0.0000001
#define RATE_MAX     1000000.0

#define INVAR_MIN    0.0001
#define INVAR_MAX    0.9999

#define TT_MIN       0.0000001
#define TT_MAX       1000000.0

#define FREQ_MIN     0.001

#define LG4X_RATE_MIN 0.0000001
#define LG4X_RATE_MAX 1000.0

/* 
   previous values between 0.001 and 0.000001

   TO AVOID NUMERICAL PROBLEMS WHEN FREQ == 0 IN PARTITIONED MODELS, ESPECIALLY WITH AA 
   previous value of FREQ_MIN was: 0.000001, but this seemed to cause problems with some 
   of the 7-state secondary structure models with some rather exotic small toy test datasets,
   on the other hand 0.001 caused problems with some of the 16-state secondary structure models

   For some reason the frequency settings seem to be repeatedly causing numerical problems
   
*/

#define ITMAX 100



#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))

#define ABS(x)    (((x)<0)   ?  (-(x)) : (x))
#define MIN(x,y)  (((x)<(y)) ?    (x)  : (y))
#define MAX(x,y)  (((x)>(y)) ?    (x)  : (y))
#define NINT(x)   ((int) ((x)>0 ? ((x)+0.5) : ((x)-0.5)))


#define LOG(x)  log(x)

#define FABS(x) fabs(x)


#define EXP(x)  exp(x)






#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))

//#define programName        "ExaML"
//#define programVersion     "2.0.3"
//#define programDate        "June 25 2014"


#define  TREE_EVALUATION            0
#define  BIG_RAPID_MODE             1
#define  QUARTET_CALCULATION        2


#define M_GTRCAT         1
#define M_GTRGAMMA       2
#define M_BINCAT         3
#define M_BINGAMMA       4
#define M_PROTCAT        5
#define M_PROTGAMMA      6
#define M_32CAT          7
#define M_32GAMMA        8
#define M_64CAT          9
#define M_64GAMMA        10


#define DAYHOFF    0
#define DCMUT      1
#define JTT        2
#define MTREV      3
#define WAG        4
#define RTREV      5
#define CPREV      6
#define VT         7
#define BLOSUM62   8
#define MTMAM      9
#define LG         10
#define MTART      11
#define MTZOA      12
#define PMB        13
#define HIVB       14
#define HIVW       15
#define JTTDCMUT   16
#define FLU        17 
#define STMTREV    18
#define AUTO       19
#define LG4M       20
#define LG4X       21
#define GTR        22  /* GTR always needs to be the last one */

#define NUM_PROT_MODELS 23

/* bipartition stuff */

#define BIPARTITIONS_ALL       0
#define GET_BIPARTITIONS_BEST  1
#define DRAW_BIPARTITIONS_BEST 2
#define BIPARTITIONS_BOOTSTOP  3
#define BIPARTITIONS_RF  4



/* bootstopping stuff */

#define BOOTSTOP_PERMUTATIONS 100
#define START_BSTOP_TEST      10

#define FC_THRESHOLD          99
#define FC_SPACING            50
#define FC_LOWER              0.99
#define FC_INIT               20

#define FREQUENCY_STOP 0
#define MR_STOP        1
#define MRE_STOP       2
#define MRE_IGN_STOP   3

#define MR_CONSENSUS 0
#define MRE_CONSENSUS 1
#define STRICT_CONSENSUS 2



/* bootstopping stuff end */


#define TIP_TIP     0
#define TIP_INNER   1
#define INNER_INNER 2

#define MIN_MODEL        -1
#define BINARY_DATA      0
#define DNA_DATA         1
#define AA_DATA          2
#define SECONDARY_DATA   3
#define SECONDARY_DATA_6 4
#define SECONDARY_DATA_7 5
#define GENERIC_32       6
#define GENERIC_64       7
#define MAX_MODEL        8

#define SEC_6_A 0
#define SEC_6_B 1
#define SEC_6_C 2
#define SEC_6_D 3
#define SEC_6_E 4

#define SEC_7_A 5
#define SEC_7_B 6
#define SEC_7_C 7
#define SEC_7_D 8
#define SEC_7_E 9
#define SEC_7_F 10

#define SEC_16   11
#define SEC_16_A 12
#define SEC_16_B 13
#define SEC_16_C 14
#define SEC_16_D 15
#define SEC_16_E 16
#define SEC_16_F 17
#define SEC_16_I 18
#define SEC_16_J 19
#define SEC_16_K 20

#define ORDERED_MULTI_STATE 0
#define MK_MULTI_STATE      1
#define GTR_MULTI_STATE     2





#define CAT         0
#define GAMMA       1
#define GAMMA_I     2



typedef  int boolean;


typedef struct {
  double lh;
  int tree;
  double weight;
} elw;

struct ent
{
  unsigned int *bitVector;
  unsigned int *treeVector;
  unsigned int amountTips;
  int *supportVector;
  unsigned int bipNumber;
  unsigned int bipNumber2;
  unsigned int supportFromTreeset[2]; 
  struct ent *next;
};

typedef struct ent entry;

typedef unsigned int hashNumberType;



/*typedef uint_fast32_t parsimonyNumber;*/

#define PCF 32

/*
  typedef uint64_t parsimonyNumber;

  #define PCF 16


typedef unsigned char parsimonyNumber;

#define PCF 2
*/

typedef struct
{
  hashNumberType tableSize;
  entry **table;
  hashNumberType entryCount;
}
  hashtable;


struct stringEnt
{
  int nodeNumber;
  char *word;
  struct stringEnt *next;
};

typedef struct stringEnt stringEntry;
 
typedef struct
{
  hashNumberType tableSize;
  stringEntry **table;
}
  stringHashtable;





typedef struct ratec
{
  double accumulatedSiteLikelihood;
  double rate;
}
  rateCategorize;


typedef struct
{
  int tipCase;
  int pNumber;
  int qNumber;
  int rNumber;
  double qz[NUM_BRANCHES];
  double rz[NUM_BRANCHES];
} traversalInfo;

typedef struct
{
  traversalInfo *ti;
  int count;
  int functionType;
  boolean traversalHasChanged;
  boolean *executeModel;
  double  *parameterValues;
} traversalData;


struct noderec;



typedef struct
{
 

  unsigned int *vector; 
  int support;   
  struct noderec *oP;
  struct noderec *oQ;
} branchInfo;








typedef struct
{
  boolean valid;
  int partitions;
  int *partitionList;
}
  linkageData;

typedef struct
{
  int entries;
  linkageData* ld;
}
  linkageList;


typedef  struct noderec
{
  double           z[NUM_BRANCHES];
#ifdef _BAYESIAN 
  double           z_tmp[NUM_BRANCHES];
#endif 
  struct noderec  *next;
  struct noderec  *back;
  hashNumberType   hash;
  int              number;
  char             x;
  char             xPars;
  char             xBips;
}
  node, *nodeptr;

typedef struct
  {
    double lh;
    int number;
  }
  info;

typedef struct bInf {
  double likelihood;
  nodeptr node;
} bestInfo;

typedef struct iL {
  bestInfo *list;
  int n;
  int valid;
} infoList;



typedef unsigned int parsimonyNumber;




typedef struct {
  int     states;
  int     maxTipStates;
  size_t     lower;
  size_t     upper;
  size_t     width;

  size_t offset; 		/* NEW: makes the data assigned to
				   this process identifiable (since we
				   now, that all data from one
				   partition must be in one contiguous
				   chunk).  */

  int     dataType;
  int     protModels;
  int     autoProtModels;
  int     protFreqs;
  boolean nonGTR;
  boolean optimizeBaseFrequencies;
  int     numberOfCategories;

  char   *partitionName;
  int    *symmetryVector;
  int    *frequencyGrouping;
    
  double *sumBuffer; 
  double gammaRates[4];
  double *EIGN;
  double *EV;
  double *EI;
  double *left;
  double *right;

   /* LG4 */

  double *rawEIGN_LG4[4];
  double *EIGN_LG4[4];
  double *EV_LG4[4];
  double *EI_LG4[4];   

  double *frequencies_LG4[4];
  double *tipVector_LG4[4];
  double *substRates_LG4[4];
  
  /* LG4X */

  double weights[4];
  double weightExponents[4];

  double weightsBuffer[4];
  double weightExponentsBuffer[4];

  /* LG4 */

  double *frequencies;
  double *freqExponents;
  double *empiricalFrequencies;
  double *tipVector; 
  double *substRates;    
  double *perSiteRates;
  int    *wgt;			
  int    *rateCategory;
  double alpha;

  double          **xVector;
  size_t           *xSpaceVector;
  unsigned char   **yVector;
  unsigned char    *yResource; 	/* contains the entire array, that is referenced in yVector */
  unsigned int     *globalScaler; 

  int               gapVectorLength;
  unsigned int     *gapVector;
  double           *gapColumn; 

  size_t parsimonyLength;
  parsimonyNumber *parsVect; 

  double *lhs;
  double *patrat;

#ifdef _USE_OMP
  /* thread-private data for OMP version */
  unsigned int **threadGlobalScaler;
  double *reductionBuffer;
  double *reductionBuffer2;
#endif

#ifdef __MIC_NATIVE
  double *mic_EV;
  double *mic_tipVector;

  /* these arrays will store the precomputed product of tipVector and left/right P-matrix */
  double *mic_umpLeft;
  double *mic_umpRight;
#endif  

} pInfo;



typedef struct 
{
  int left;
  int right;
  double likelihood;
} lhEntry;


typedef struct 
{
  int count;
  int size;
  lhEntry *entries;
} lhList;


typedef struct List_{
  void *value; 			
  struct List_ *next; 
} List;


#define REARR_SETTING 1
#define FAST_SPRS     2
#define SLOW_SPRS     3
#define MOD_OPT       4
#define QUARTETS      5

typedef struct {
  boolean useMedian;
  int saveBestTrees;
  boolean saveMemory;
  boolean searchConvergenceCriterion;
  boolean perGeneBranchLengths; //adef
  double likelihoodEpsilon; //adef
  int categories;
  int mode; //adef
  int fastTreeEvaluation;
  boolean initialSet;//adef
  int initial;//adef
  int rateHetModel;
  int autoProteinSelectionType;
  
  //quartets 
  boolean useQuartetGrouping;//adef
  unsigned long int numberRandomQuartets;//adef

} commandLine;

typedef struct {
 
  int state;

  /* search algorithm */

  unsigned int vLength;

  boolean constraintTree;
  
  int rearrangementsMax;
  int rearrangementsMin;
  int thoroughIterations;
  int fastIterations;
  int treeVectorLength;  
  int mintrav;
  int maxtrav;
  int bestTrav;
  int    Thorough;
  int    optimizeRateCategoryInvocations;
  
  double accumulatedTime;

  double startLH; 
  double lh;
  double previousLh;
  double difference;
  double epsilon;
  
  boolean impr;
  boolean cutoff;  
       
  double tr_startLH;
  double tr_endLH;
  double tr_likelihood;
  double tr_bestOfNode;
  
  double tr_lhCutoff;
  double tr_lhAVG;
  double tr_lhDEC;
  int    tr_NumberOfCategories;
  int    tr_itCount;  
  int    tr_doCutoff;

  /* modOpt */

  int catOpt;
  int treeIteration; 
  /* quartets */

  long seed;
  int flavor;   
  uint64_t quartetCounter;  
  long filePosition;
  char quartetFileName[1024];
  //FILE NAME???

  /* command line settings */

  commandLine cmd;
  
} checkPointState;


typedef struct {
  double EIGN[19] __attribute__ ((aligned (BYTE_ALIGNMENT)));             
  double EV[400] __attribute__ ((aligned (BYTE_ALIGNMENT)));                
  double EI[380] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double substRates[190];        
  double frequencies[20] ;      
  double tipVector[460] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double left[1600] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double right[1600] __attribute__ ((aligned (BYTE_ALIGNMENT)));
} siteAAModels;


typedef struct assign
{
  int partitionId;
  int procId;	     /* to which process is the partition assigned  */
  size_t offset;     /* what is the offset of this assignment */
  size_t width; 
} Assign ; 


typedef  struct  {

  int *ti;

  unsigned int randomSeed;
  boolean constraintTree;
  boolean useGappedImplementation;
  boolean saveMemory;  
  int              saveBestTrees;

  stringHashtable  *nameHash;

  pInfo            *partitionData;
  

  char             *secondaryStructureInput;

  boolean          *executeModel;

  double           *perPartitionLH;

  traversalData    td[1];

  int              maxCategories;
  int              categories;
  
  double           coreLZ[NUM_BRANCHES];
  int              numBranches;
  
  
 
  branchInfo	   *bInf;

  int              multiStateModel;


  boolean curvatOK[NUM_BRANCHES];
  /* the stuff below is shared among DNA and AA, span does
     not change depending on datatype */

  /* model stuff end */

  unsigned char             **yVector;
  int              secondaryStructureModel;
  size_t           originalCrunchedLength;
 
 
  int              *secondaryStructurePairs;


  double            *partitionContributions;
  double            *partitionWeights;

  double            lhCutoff;
  double            lhAVG;
  unsigned long     lhDEC;
  unsigned long     itCount;
  int               numberOfInvariableColumns;
  int               weightOfInvariableColumns;
  int               rateHetModel;

  double           startLH;
  double           endLH;
  double           likelihood;
  
 
  node           **nodep;
  nodeptr          nodeBaseAddress;
  node            *start;
  int              mxtips;  

  int              *constraintVector;
  int              numberOfSecondaryColumns;
  boolean          searchConvergenceCriterion;
  int              ntips;
  int              nextnode;  
  int              NumberOfModels;    

  boolean          bigCutoff;
  boolean          partitionSmoothed[NUM_BRANCHES];
  boolean          partitionConverged[NUM_BRANCHES];
  boolean          rooted;
  boolean          doCutoff;
 


  double         gapyness;

  char **nameList;
  char *tree_string;
  char *treeStrings;
  char *tree0;
  char *tree1;
  int treeStringLength;
 
  unsigned int bestParsimony;
  unsigned int *parsimonyScore;
  
  double bestOfNode;
  nodeptr removeNode;
  nodeptr insertNode;

  double zqr[NUM_BRANCHES];
  double currentZQR[NUM_BRANCHES];

  double currentLZR[NUM_BRANCHES];
  double currentLZQ[NUM_BRANCHES];
  double currentLZS[NUM_BRANCHES];
  double currentLZI[NUM_BRANCHES];
  double lzs[NUM_BRANCHES];
  double lzq[NUM_BRANCHES];
  double lzr[NUM_BRANCHES];
  double lzi[NUM_BRANCHES];

 
 


  unsigned int **bitVectors;

  unsigned int vLength;

  hashtable *h;

  char bits_in_16bits [0x1u << 16];
  
  boolean useMedian;

  int autoProteinSelectionType;

  int numberOfTrees;

  double *likelihoods;

  boolean fastTreeEvaluation;
  
  int numAssignments; 
  Assign *partAssigns;

  /** 
      IMPORTANT: 
      
      introducing a few resource pointers. All memeory needed for
      example for per-partition patrats is owned by these
      basepointers, the per-partition pointer just points at the
      contiguous block of memory.
      
      The big advantage, why I really think, this is worth it is, that
      these base pointers can be used to conveniently gather/scatter
      all data at a master. The master still has to reorder the
      gathered data, but less copying is necessary at the workers

      REQUIREMENTS: 

      * all memory necessary for a partition must be in a contiguous
      block,

      * memory for partitions is ordered by partition id (first
      * partition 1, then partition 2,... )
   */ 

  double *patrat_basePtr; 
  int *rateCategory_basePtr; 
  double *lhs_basePtr;

#ifdef _USE_OMP
  /* number of OMP threads*/
  int nThreads;

  /* maximum number of partitions assigned to a single thread */
  int maxModelsPerThread;

  /* maximum number of threads assigned to a single partition */
  int maxThreadsPerModel;

  /* partition-to-threads assignments: indexed by thread */
  Assign **threadPartAssigns;

  /* partition-to-threads assignments: indexed by partition id */
  Assign **partThreadAssigns;

#endif

} tree;


/***************************************************************/

typedef struct {
  int partitionNumber;
  int partitionLength;
} partitionType;

typedef struct
{
  double z[NUM_BRANCHES];
  nodeptr p, q;
  int cp, cq;
}
  connectRELL, *connptrRELL;

typedef  struct
{
  connectRELL     *connect; 
  int             start;
  double          likelihood;
}
  topolRELL;


typedef  struct
{
  int max;
  topolRELL **t;
}
  topolRELL_LIST;


/**************************************************************/



typedef struct conntyp {
    double           z[NUM_BRANCHES];           /* branch length */
    node            *p, *q;       /* parent and child sectors */
    void            *valptr;      /* pointer to value of subtree */
    int              descend;     /* pointer to first connect of child */
    int              sibling;     /* next connect from same parent */
    } connect, *connptr;

typedef  struct {
    double           likelihood;
  int              initialTreeNumber;
    connect         *links;       /* pointer to first connect (start) */
    node            *start;
    int              nextlink;    /* index of next available connect */
                                  /* tr->start = tpl->links->p */
    int              ntips;
    int              nextnode;
    int              scrNum;      /* position in sorted list of scores */
    int              tplNum;      /* position in sorted list of trees */

    } topol;

typedef struct {
    double           best;        /* highest score saved */
    double           worst;       /* lowest score saved */
    topol           *start;       /* starting tree for optimization */
    topol          **byScore;
    topol          **byTopol;
    int              nkeep;       /* maximum topologies to save */
    int              nvalid;      /* number of topologies saved */
    int              ninit;       /* number of topologies initialized */
    int              numtrees;    /* number of alternatives tested */
    boolean          improved;
    } bestlist;

#define randomTree    0
#define givenTree     1 
#define parsimonyTree 2

typedef  struct {
  int              bestTrav;
  int              max_rearrange;
  int              stepwidth;
  int              initial;
  boolean          initialSet;
  int              mode; 
  boolean        perGeneBranchLengths;
  boolean        permuteTreeoptimize; 
  boolean        compressPatterns;
  double         likelihoodEpsilon;
  boolean        useCheckpoint;
  boolean        useQuartetGrouping;
  unsigned long int numberRandomQuartets;
  unsigned long int quartetCkpInterval; 
 
#ifdef _BAYESIAN 
  boolean       bayesian;
#endif
} analdef;




typedef struct 
{
  int leftLength;
  int rightLength;
  int eignLength;
  int evLength;
  int eiLength;
  int substRatesLength;
  int frequenciesLength;
  int tipVectorLength;
  int symmetryVectorLength;
  int frequencyGroupingLength;

  boolean nonGTR;

  int undetermined;

  const char *inverseMeaning;

  int states;

  boolean smoothFrequencies;

  const unsigned  int *bitVector;

} partitionLengths;

/****************************** FUNCTIONS ****************************************************/

#ifdef _BAYESIAN 
extern void mcmc(tree *tr, analdef *adef);
#endif


boolean isThisMyPartition(tree *localTree, int tid, int model);

extern boolean allSmoothed(tree *tr);

extern int treeFindTipName(FILE *fp, tree *tr, boolean check);

extern void computePlacementBias(tree *tr, analdef *adef);

extern int lookupWord(char *s, stringHashtable *h);

extern void getDataTypeString(tree *tr, int model, char typeOfData[1024]);

extern unsigned int genericBitCount(unsigned int* bitVector, unsigned int bitVectorLength);
extern int countTips(nodeptr p, int numsp);
extern entry *initEntry(void);
extern void computeRogueTaxa(tree *tr, char* treeSetFileName, analdef *adef);
extern unsigned int precomputed16_bitcount(unsigned int n, char *bits_in_16bits);





extern size_t discreteRateCategories(int rateHetModel);

extern partitionLengths * getPartitionLengths(pInfo *p);
extern boolean getSmoothFreqs(int dataType);
extern const unsigned int *getBitVector(int dataType);
extern int getUndetermined(int dataType);
extern int getStates(int dataType);
extern char getInverseMeaning(int dataType, unsigned char state);
extern double gettime ( void );
extern int gettimeSrand ( void );
extern double randum ( long *seed );

extern void getxnode ( nodeptr p );
extern void hookup ( nodeptr p, nodeptr q, double *z, int numBranches);
extern void hookupDefault ( nodeptr p, nodeptr q, int numBranches);
extern boolean whitechar ( int ch );
extern void errorExit ( int e );
extern void printResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printBootstrapResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printBipartitionResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printLog ( tree *tr);
extern void printStartingTree ( tree *tr, analdef *adef, boolean finalPrint );
extern void writeInfoFile ( analdef *adef, tree *tr, double t );
extern int main ( int argc, char *argv[] );
extern void calcBipartitions ( tree *tr, analdef *adef, char *bestTreeFileName, char *bootStrapFileName );
extern void initReversibleGTR (tree *tr, int model);
extern double LnGamma ( double alpha );
extern double IncompleteGamma ( double x, double alpha, double ln_gamma_alpha );
extern double PointNormal ( double prob );
extern double PointChi2 ( double prob, double v );
extern void makeGammaCats (double alpha, double *gammaRates, int K, boolean useMedian);
extern void initModel ( tree *tr);
extern void doAllInOne ( tree *tr, analdef *adef );

extern void classifyML(tree *tr, analdef *adef);

extern void resetBranches ( tree *tr );
extern void modOpt ( tree *tr, double likelihoodEpsilon, analdef *adef, int treeIteration);



extern void computeBOOTRAPID (tree *tr, analdef *adef, long *radiusSeed);
extern void optimizeRAPID ( tree *tr, analdef *adef );
extern void thoroughOptimization ( tree *tr, analdef *adef, topolRELL_LIST *rl, int index );
extern int treeOptimizeThorough ( tree *tr, int mintrav, int maxtrav);
extern void computeQuartets(tree *tr, analdef *adef);

extern void makeRandomTree ( tree *tr);
extern void nodeRectifier ( tree *tr );
extern void makeParsimonyTreeFast(tree *tr);
extern void allocateParsimonyDataStructures(tree *tr);
extern void freeParsimonyDataStructures(tree *tr);
extern void parsimonySPR(nodeptr p, tree *tr);

extern FILE *myfopen(const char *path, const char *mode);


extern boolean initrav ( tree *tr, nodeptr p );
extern void initravPartition ( tree *tr, nodeptr p, int model );
extern boolean update ( tree *tr, nodeptr p );
extern boolean smooth ( tree *tr, nodeptr p );
extern boolean smoothTree ( tree *tr, int maxtimes );
extern boolean localSmooth ( tree *tr, nodeptr p, int maxtimes );
extern boolean localSmoothMulti(tree *tr, nodeptr p, int maxtimes, int model);
extern void initInfoList ( int n );
extern void freeInfoList ( void );
extern void insertInfoList ( nodeptr node, double likelihood );
extern boolean smoothRegion ( tree *tr, nodeptr p, int region );
extern boolean regionalSmooth ( tree *tr, nodeptr p, int maxtimes, int region );
extern nodeptr removeNodeBIG ( tree *tr, nodeptr p, int numBranches);
extern nodeptr removeNodeRestoreBIG ( tree *tr, nodeptr p );
extern boolean insertBIG ( tree *tr, nodeptr p, nodeptr q, int numBranches);
extern boolean insertRestoreBIG ( tree *tr, nodeptr p, nodeptr q );
extern boolean testInsertBIG ( tree *tr, nodeptr p, nodeptr q );
extern void addTraverseBIG ( tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav );
extern int rearrangeBIG ( tree *tr, nodeptr p, int mintrav, int maxtrav );
extern void traversalOrder ( nodeptr p, int *count, nodeptr *nodeArray );
extern double treeOptimizeRapid ( tree *tr, int mintrav, int maxtrav, analdef *adef, bestlist *bt, bestlist *bestML);
extern boolean testInsertRestoreBIG ( tree *tr, nodeptr p, nodeptr q );
extern void restoreTreeFast ( tree *tr );
extern int determineRearrangementSetting ( tree *tr, analdef *adef, bestlist *bestT, bestlist *bt, bestlist *bestML);
extern void computeBIGRAPID ( tree *tr, analdef *adef, boolean estimateModel);
extern boolean treeEvaluate ( tree *tr, double smoothFactor );
extern boolean treeEvaluatePartition ( tree *tr, double smoothFactor, int model );

extern void meshTreeSearch(tree *tr, analdef *adef, int thorough);

extern void initTL ( topolRELL_LIST *rl, tree *tr, int n );
extern void freeTL ( topolRELL_LIST *rl);
extern void restoreTL ( topolRELL_LIST *rl, tree *tr, int n );
extern void resetTL ( topolRELL_LIST *rl );
extern void saveTL ( topolRELL_LIST *rl, tree *tr, int index );

extern int  saveBestTree (bestlist *bt, tree *tr, boolean keepIdenticalTrees);
extern int  recallBestTree (bestlist *bt, int rank, tree *tr);
extern int initBestTree ( bestlist *bt, int newkeep, int numsp );
extern void resetBestTree ( bestlist *bt );
extern boolean freeBestTree ( bestlist *bt );


extern char *Tree2String ( char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, boolean printLikelihood, 
			   boolean rellTree, boolean finalPrint, int perGene, boolean branchLabelSupport, boolean printSHSupport);
extern void printTreePerGene(tree *tr, analdef *adef, char *fileName, char *permission);



extern int treeReadLen (FILE *fp, tree *tr, boolean readBranches, boolean readNodeLabels, boolean topologyOnly);
extern void treeReadTopologyString(char *treeString, tree *tr);
extern boolean treeReadLenMULT ( FILE *fp, tree *tr, int *partCount);

extern void getStartingTree ( tree *tr);

extern void computeBootStopOnly(tree *tr, char *bootStrapFileName, analdef *adef);
extern boolean bootStop(tree *tr, hashtable *h, int numberOfTrees, double *pearsonAverage, unsigned int **bitVectors, int treeVectorLength, unsigned int vectorLength);
extern void computeConsensusOnly(tree *tr, char* treeSetFileName, analdef *adef);
extern double evaluatePartialGeneric (tree *, int i, double ki, int _model);
extern void evaluateGeneric (tree *tr, nodeptr p, boolean fullTraversal);
extern void newviewGeneric (tree *tr, nodeptr p, boolean masked);
extern void newviewGenericMulti (tree *tr, nodeptr p, int model);
extern void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, boolean mask);
extern void makenewzGenericDistance(tree *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2);
extern double evaluatePartitionGeneric (tree *tr, nodeptr p, int model);
extern void newviewPartitionGeneric (tree *tr, nodeptr p, int model);
extern double evaluateGenericVector (tree *tr, nodeptr p);
extern void categorizeGeneric (tree *tr, nodeptr p);
extern double makenewzPartitionGeneric(tree *tr, nodeptr p, nodeptr q, double z0, int maxiter, int model);
extern boolean isTip(int number, int maxTips);
extern void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches, boolean partialTraversal);



extern void   newviewIterative(tree *tr, int startIndex);

extern void evaluateIterative(tree *);

extern void *malloc_aligned( size_t size);

extern void storeExecuteMaskInTraversalDescriptor(tree *tr);
extern void storeValuesInTraversalDescriptor(tree *tr, double *value);




extern void makenewzIterative(tree *);
extern void execCore(tree *, volatile double *dlnLdlz, volatile double *d2lnLdlz2);



extern void determineFullTraversal(nodeptr p, tree *tr);
/*extern void optRateCat(tree *, int i, double lower_spacing, double upper_spacing, double *lhs);*/





extern double evaluateGenericInitravPartition(tree *tr, nodeptr p, int model);
extern void evaluateGenericVectorIterative(tree *, int startIndex, int endIndex);
extern void categorizeIterative(tree *, int startIndex, int endIndex);

extern void fixModelIndices(tree *tr, int endsite, boolean fixRates);
extern void calculateModelOffsets(tree *tr);
extern void gammaToCat(tree *tr);
extern void catToGamma(tree *tr, analdef *adef);


extern nodeptr findAnyTip(nodeptr p, int numsp);

extern void parseProteinModel(analdef *adef);



extern void computeNextReplicate(tree *tr, long *seed, int *originalRateCategories, int *originalInvariant, boolean isRapid, boolean fixRates);
/*extern void computeNextReplicate(tree *tr, analdef *adef, int *originalRateCategories, int *originalInvariant);*/

extern void putWAG(double *ext_initialRates);

extern void reductionCleanup(tree *tr, int *originalRateCategories, int *originalInvariant);
extern void parseSecondaryStructure(tree *tr, analdef *adef, int sites);
extern void printPartitions(tree *tr);
extern void compareBips(tree *tr, char *bootStrapFileName, analdef *adef);
extern void computeRF(tree *tr, char *bootStrapFileName, analdef *adef);


extern  unsigned int **initBitVector(int mxtips, unsigned int *vectorLength);
extern hashtable *copyHashTable(hashtable *src, unsigned int vectorLength);
extern hashtable *initHashTable(unsigned int n);
extern void cleanupHashTable(hashtable *h, int state);
extern double convergenceCriterion(hashtable *h, int mxtips);
extern void freeBitVectors(unsigned int **v, int n);
extern void freeHashTable(hashtable *h);
extern stringHashtable *initStringHashTable(hashNumberType n);
extern void addword(char *s, stringHashtable *h, int nodeNumber);


extern void printBothOpen(const char* format, ... );
extern void initRateMatrix(tree *tr);

extern void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h, int treeNumber, int function, branchInfo *bInf,
				    int *countBranches, int treeVectorLength, boolean traverseOnly, boolean computeWRF);

extern int getIncrement(tree *tr, int model);



extern void writeBinaryModel(tree *tr);
extern void readBinaryModel(tree *tr);
extern void treeEvaluateRandom (tree *tr, double smoothFactor);
extern void treeEvaluateProgressive(tree *tr);

extern void testGapped(tree *tr);

extern boolean issubset(unsigned int* bipA, unsigned int* bipB, unsigned int vectorLen);
extern boolean compatible(entry* e1, entry* e2, unsigned int bvlen);



extern int *permutationSH(tree *tr, int nBootstrap, long _randomSeed);

extern void checkPerSiteRates(const tree * const tr ); 

extern void restart(tree *tr, analdef *adef);

extern void writeCheckpoint(tree *tr, analdef *adef);

extern boolean isGap(unsigned int *x, int pos);
extern boolean noGap(unsigned int *x, int pos);

extern void scaleLG4X_EIGN(tree *tr, int model);

extern void myBinFwrite(void *ptr, size_t size, size_t nmemb, FILE *byteFile);
extern void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile);

extern void newviewGTRGAMMAPROT_AVX_LG4(int tipCase,
					double *x1, double *x2, double *x3, double *extEV[4], double *tipVector[4],
					int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, 
					double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling);

extern void newviewGTRCAT_AVX(int tipCase,  double *EV,  int *cptr,
			      double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
			      unsigned char *tipX1, unsigned char *tipX2,
			      int n,  double *left, double *right, int *wgt, int *scalerIncrement);


extern void newviewGTRCATPROT_AVX(int tipCase, double *extEV,
				  int *cptr,
				  double *x1, double *x2, double *x3, double *tipVector,
				  unsigned char *tipX1, unsigned char *tipX2,
				  int n, double *left, double *right, int *wgt, int *scalerIncrement);


extern void newviewGTRGAMMA_AVX(int tipCase,
				double *x1_start, double *x2_start, double *x3_start,
				double *EV, double *tipVector,
				unsigned char *tipX1, unsigned char *tipX2,
				const int n, double *left, double *right, int *wgt, int *scalerIncrement
				);

extern void newviewGTRGAMMAPROT_AVX(int tipCase,
				    double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				    unsigned char *tipX1, unsigned char *tipX2, int n, 
				    double *left, double *right, int *wgt, int *scalerIncrement);

/* memory saving functions */

void newviewGTRCAT_AVX_GAPPED_SAVE(int tipCase,  double *EV,  int *cptr,
				   double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
				   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

void newviewGTRCATPROT_AVX_GAPPED_SAVE(int tipCase, double *extEV,
				       int *cptr,
				       double *x1, double *x2, double *x3, double *tipVector,
				       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				       int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
				       unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				       double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

void  newviewGTRGAMMA_AVX_GAPPED_SAVE(int tipCase,
				      double *x1_start, double *x2_start, double *x3_start,
				      double *extEV, double *tipVector,
				      int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				      const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
				      unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
				      double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
				      );

void newviewGTRGAMMAPROT_AVX_GAPPED_SAVE(int tipCase,
					 double *x1_start, double *x2_start, double *x3_start, double *extEV, double *tipVector,
					 int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, 
					 double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
					 unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
					 double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn); 



/* from communication.c */
void calculateLengthAndDisplPerProcess(tree *tr, int **length_result, int **disp_result);
void scatterDistrbutedArray(tree *tr, void *src, void *destination, MPI_Datatype type, int *countPerProc, int *displPerProc);
void gatherDistributedArray(tree *tr, void **destination, void *src, MPI_Datatype type, int* countPerProc, int *displPerProc);


#endif



#

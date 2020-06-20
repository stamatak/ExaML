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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifdef WIN32
#include <direct.h>
#endif

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <unistd.h>
#include <getopt.h>

#include <mpi.h>

#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse.h"
#if defined(SIMDE_SSE_NATIVE)
#include <xmmintrin.h>
/*
  special bug fix, enforces denormalized numbers to be flushed to zero,
  without this program is a tiny bit faster though.
  #include <emmintrin.h> 
  #define MM_DAZ_MASK    0x0040
  #define MM_DAZ_ON    0x0040
  #define MM_DAZ_OFF    0x0000
*/
#endif

#include "axml.h"


#include "globalVariables.h"

#include "byteFile.h"
#include "partitionAssignment.h"

#ifdef __MIC_NATIVE
#include "mic_native.h"
#endif

/***************** UTILITY FUNCTIONS **************************/

/*pInfo *cleanPinfoInit()
{
  pInfo *p = (pInfo*)malloc(sizeof(pInfo));

  
  return p;
  }*/


void storeExecuteMaskInTraversalDescriptor(tree *tr)
{
   int model;
      
   for(model = 0; model < tr->NumberOfModels; model++)
     tr->td[0].executeModel[model] = tr->executeModel[model];
}

void storeValuesInTraversalDescriptor(tree *tr, double *value)
{
   int model;
      
   for(model = 0; model < tr->NumberOfModels; model++)
     tr->td[0].parameterValues[model] = value[model];
}



void myBinFwrite(void *ptr, size_t size, size_t nmemb, FILE *byteFile)
{
  size_t
    bytes_read;
  
  bytes_read = fwrite(ptr, size, nmemb, byteFile);

  assert(bytes_read == nmemb);
}

void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile)
{  
  size_t
    bytes_read;
  
  bytes_read = fread(ptr, size, nmemb, byteFile);

  assert(bytes_read == nmemb);
}


static void outOfMemory(void)
{
  printf("ExaML process %d was not able to allocate enough memory.\n", processID);
  printf("Please check the approximate memory consumption of your dataset using\n");
  printf("the memory calculator at http://www.exelixis-lab.org/web/software/raxml/index.html.\n");
  printf("ExaML will exit now\n");

  
  MPI_Abort(MPI_COMM_WORLD, -1);

  exit(-1);
 }

void *malloc_aligned(size_t size) 
{
  void 
    *ptr = (void *)NULL;
 
  int 
    res;
  

#ifdef WIN32
  ptr = _aligned_malloc(size, BYTE_ALIGNMENT);;
#else
  res = posix_memalign( &ptr, BYTE_ALIGNMENT, size );

  if(res != 0)
  {
    outOfMemory();
    assert(0);
  }
#endif 
   
  return ptr;
}







static void printBoth(FILE *f, const char* format, ... )
{
  if(processID == 0)
    {
      va_list args;
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
      
      va_start(args, format);
      vprintf(format, args );
      va_end(args);
    }
}




void printBothOpen(const char* format, ... )
{
  if(processID == 0)
    {
      FILE *f = myfopen(infoFileName, "ab");
      
      va_list args;
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
      
      va_start(args, format);
      vprintf(format, args );
      va_end(args);
      
      fclose(f);
    }
}

static void printBothOpenDifferentFile(char *fileName, const char* format, ... )
{
  if(processID == 0)
    {
      FILE 
	*f = myfopen(fileName, "ab");
      
      va_list 
	args;
      
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
            
      fclose(f);
    }
}



boolean getSmoothFreqs(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}

const unsigned int *getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}


int getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].states;
}

int getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].undetermined;
}



char getInverseMeaning(int dataType, unsigned char state)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return  pLengths[dataType].inverseMeaning[state];
}

partitionLengths *getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
    states    = p->states,
    tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  pLength.leftLength = pLength.rightLength = states * states;
  pLength.eignLength = states;
  pLength.evLength   = states * states;
  pLength.eiLength   = states * states;
  pLength.substRatesLength = (states * states - states) / 2;
  pLength.frequenciesLength = states;
  pLength.tipVectorLength   = tipLength * states;
  pLength.symmetryVectorLength = (states * states - states) / 2;
  pLength.frequencyGroupingLength = states;
  pLength.nonGTR = FALSE;

  return (&pLengths[dataType]); 
}










size_t discreteRateCategories(int rateHetModel)
{
  size_t 
    result;

  switch(rateHetModel)
    {
    case CAT:
      result = 1;
      break;
    case GAMMA:
      result = 4;
      break;
    default:
      assert(0);
    }

  return result;
}



double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}

int gettimeSrand(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 24*60*60*localtm.tm_yday + 60*60*localtm.tm_hour + 60*localtm.tm_min  + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec;
#endif
}

double randum (long  *seed)
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}

static int filexists(char *filename)
{
  FILE 
    *fp;
  
  int 
    res;
  
  fp = fopen(filename,"rb");

  if(fp)
    {
      res = 1;
      fclose(fp);
    }
  else
    res = 0;

  return res;
}


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("The file %s you want to open for reading does not exist, exiting ...\n", path);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }
  else
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("The file %s ExaML wants to open for writing or appending can not be opened [mode: %s], exiting ...\n",
		   path, mode);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }


}





/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/


boolean isTip(int number, int maxTips)
{
  assert(number > 0);

  if(number <= maxTips)
    return TRUE;
  else
    return FALSE;
}









void getxnode (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->x || (s = s->next)->x)
    {
      p->x = s->x;
      s->x = 0;
    }

  assert(p->x);
}





void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}

void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;
}


/***********************reading and initializing input ******************/







boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}










static unsigned int KISS32(void)
{
  static unsigned int 
    x = 123456789, 
    y = 362436069,
    z = 21288629,
    w = 14921776,
    c = 0;

  unsigned int t;

  x += 545925293;
  y ^= (y<<13); 
  y ^= (y>>17); 
  y ^= (y<<5);
  t = z + w + c; 
  z = w; 
  c = (t>>31); 
  w = t & 2147483647;

  return (x+y+w);
}

static boolean setupTree (tree *tr)
{
  nodeptr  
    p0, 
    p, 
    q;
  
  int
    i,
    j,   
    tips,
    inter; 
  
  tr->bigCutoff = FALSE;
  
  tr->maxCategories = MAX(4, tr->categories);
  
  tr->partitionContributions = (double *)malloc(sizeof(double) * tr->NumberOfModels);
  tr->partitionWeights       = (double *)malloc(sizeof(double) * tr->NumberOfModels);

  for(i = 0; i < tr->NumberOfModels; i++)
    {
      tr->partitionContributions[i] = -1.0;
      tr->partitionWeights[i] = -1.0;
    }
  
  tr->perPartitionLH = (double *)malloc(sizeof(double) * tr->NumberOfModels);
    
  for(i = 0; i < tr->NumberOfModels; i++)    
    tr->perPartitionLH[i] = 0.0;	    
     
  tips  = tr->mxtips;
  inter = tr->mxtips - 1;

  /* printf("%d tips\t%d inner\n", tips, inter); */
   
  
  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)calloc(tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)calloc(tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)calloc(tr->treeStringLength, sizeof(char));


  /* TODO, must that be so long ? */
  /* assert(0);  */


  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * tr->mxtips);
  tr->td[0].executeModel = (boolean *)malloc(sizeof(boolean) * tr->NumberOfModels);
  tr->td[0].parameterValues = (double *)malloc(sizeof(double) * tr->NumberOfModels);  
  
  tr->constraintVector = (int *)malloc((2 * tr->mxtips) * sizeof(int));


  if (!(p0 = (nodeptr) malloc((tips + 3*inter) * sizeof(node))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }
  
  tr->nodeBaseAddress = p0;


  if (!(tr->nodep = (nodeptr *) malloc((2*tr->mxtips) * sizeof(nodeptr))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory, too\n");
      return  FALSE;
    }

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
    {
      p = p0++;

      p->hash   =  KISS32(); /* hast table stuff */
      p->x      =  0;
      p->xBips  =  0;
      p->number =  i;
      p->next   =  p;
      p->back   = (node *)NULL;     
      tr->nodep[i] = p;
    }

  for (i = tips + 1; i <= tips + inter; i++)
    {
      q = (node *) NULL;
      for (j = 1; j <= 3; j++)
	{	 
	  p = p0++;
	  if(j == 1)
	    {
	      p->xBips = 1;
	      p->x = 1;
	    }
	  else
	    {
	      p->xBips = 0;
	      p->x =  0;
	    }
	  p->number = i;
	  p->next   = q;	  
	  p->back   = (node *) NULL;
	  p->hash   = 0;       
	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  tr->likelihood  = unlikely;
  tr->start       = (node *) NULL;  

  tr->ntips       = 0;
  tr->nextnode    = 0;
 
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;

  tr->bitVectors = (unsigned int **)NULL;

  tr->vLength = 0;

  tr->h = (hashtable*)NULL;
  
  tr->nameHash = initStringHashTable(10 * tr->mxtips);

  return TRUE;
}




static void initAdef(analdef *adef)
{   
  adef->max_rearrange          = 21;
  adef->stepwidth              = 5;
  adef->initial                = 10;
  adef->bestTrav               = 10;
  adef->initialSet             = FALSE; 
  adef->mode                   = BIG_RAPID_MODE; 
  adef->likelihoodEpsilon      = 0.1;
 
  adef->permuteTreeoptimize    = FALSE; 
  adef->perGeneBranchLengths   = FALSE;  
 
  adef->useCheckpoint          = FALSE;
   
  adef->useQuartetGrouping        = FALSE;
  adef->numberRandomQuartets      = 0;

  adef->quartetCkpInterval        = 1000;

#ifdef _BAYESIAN 
  adef->bayesian               = FALSE;
#endif

}



static int modelExists(char *model, tree *tr)
{  
   if(strcmp(model, "PSR\0") == 0)
    {
      tr->rateHetModel = CAT;
      return 1;
    }

  if(strcmp(model, "GAMMA\0") == 0)
    {
      tr->rateHetModel = GAMMA;
      return 1;
    }

  
  return 0;
}






/*********************************** *********************************************************/


static void printVersionInfo(void)
{
  if(processID == 0)
    printf("\n\nThis is %s version %s released by Alexandros Stamatakis, Andre J. Aberer, and Alexey Kozlov on %s.\n\n",  programName, programVersion, programDate); 
}

static void printMinusFUsage(void)
{
  printf("\n");
 

  printf("              \"-f d\": new rapid hill-climbing \n");
  printf("                      DEFAULT: ON\n");

  printf("\n");

  printf("               \"-f e\": compute the likelihood of a bunch of trees passed via -t\n");
  printf("                this option will do a quick and dirty optimization without re-optimizng\n");
  printf("                the model parameters for each tree\n");

  printf("\n");

  printf("               \"-f E\": compute the likelihood of a bunch of trees passed via -t\n");
  printf("                this option will do a thorough optimization that re-optimizes\n");
  printf("                the model parameters for each tree\n");

  printf("\n");

  printf("              \"-f o\": old and slower rapid hill-climbing without heuristic cutoff\n");

  printf("\n");

  printf("              \"-f q\": fast quartet calculator\n");
  
  printf("\n");

  printf("              DEFAULT for \"-f\": new rapid hill climbing\n");

  printf("\n");
}


static void printREADME(void)
{
  if(processID == 0)
    {
      printVersionInfo();
      printf("\n");  
      printf("\nTo report bugs use the RAxML google group\n");
      printf("Please send me all input files, the exact invocation, details of the HW and operating system,\n");
      printf("as well as all error messages printed to screen.\n\n\n");
      
      printf("examl|examl-AVX\n");
      printf("      -s binarySequenceFileName\n");
      printf("      -n outputFileNames\n");
      printf("      -m rateHeterogeneityModel\n");
      printf("      -t userStartingTree|-R binaryCheckpointFile|-g constraintTree -p randomNumberSeed\n");
      printf("      [-a]\n");
      printf("      [-B numberOfMLtreesToSave]\n"); 
      printf("      [-c numberOfCategories]\n");
      printf("      [-D]\n");
      printf("      [-e likelihoodEpsilon] \n");
      printf("      [-f d|e|E|o|q]\n");    
      printf("      [-h] \n");
      printf("      [-i initialRearrangementSetting] \n");
      printf("      [-I quartetCheckpointInterval] \n");
      printf("      [-M]\n");
      printf("      [-r randomQuartetNumber] \n");
      printf("      [-S]\n");
      printf("      [-v]\n"); 
      printf("      [-w outputDirectory] \n"); 
      printf("      [-Y quartetGroupingFileName]\n");
      printf("      [--auto-prot=ml|bic|aic|aicc]\n");
      printf("\n");  
      printf("      -a      use the median for the discrete approximation of the GAMMA model of rate heterogeneity\n");
      printf("\n");
      printf("              DEFAULT: OFF\n");
      printf("\n");
      printf("      -B      specify the number of best ML trees to save and print to file\n");
      printf("\n");
      printf("      -c      Specify number of distinct rate catgories for ExaML when modelOfEvolution\n");
      printf("              is set to GTRPSR\n");
      printf("              Individual per-site rates are categorized into numberOfCategories rate \n");
      printf("              categories to accelerate computations. \n");
      printf("\n");
      printf("              DEFAULT: 25\n");
      printf("\n");
      printf("      -D      ML search convergence criterion. This will break off ML searches if the relative \n");
      printf("              Robinson-Foulds distance between the trees obtained from two consecutive lazy SPR cycles\n");
      printf("              is smaller or equal to 1%s. Usage recommended for very large datasets in terms of taxa.\n", "%");
      printf("              On trees with more than 500 taxa this will yield execution time improvements of approximately 50%s\n",  "%");
      printf("              While yielding only slightly worse trees.\n");
      printf("\n");
      printf("              DEFAULT: OFF\n");    
      printf("\n");
      printf("      -e      set model optimization precision in log likelihood units for final\n");
      printf("              optimization of model parameters\n");
      printf("\n");
      printf("              DEFAULT: 0.1 \n"); 
      printf("\n");
      printf("      -f      select algorithm:\n");
      
      printMinusFUsage();
 
      printf("\n");
      printf("      -g      Pass a multi-furcating constraint tree to ExaML. The tree needs to contain all taxa of the alignment!\n");
      printf("              When using this option you also need to specify a random number seed via \"-p\"\n");
      printf("\n");
      printf("      -h      Display this help message.\n");
      printf("\n");  
      printf("      -i      Initial rearrangement setting for the subsequent application of topological \n");
      printf("              changes phase\n");
      printf("\n");
      printf("      -I      Set after how many quartet evaluations a new checkpoint will be printed.\n");
      printf("\n");
      printf("              DEFAULT: 1000\n");
      printf("\n");
      printf("      -m      Model of rate heterogeneity\n");
      printf("\n"); 
      printf("              select \"-m PSR\" for the per-site rate category model (this used to be called CAT in RAxML)\n");
      printf("              select \"-m GAMMA\" for the gamma model of rate heterogeneity with 4 discrete rates\n");
      printf("\n");
      printf("      -M      Switch on estimation of individual per-partition branch lengths. Only has effect when used in combination with \"-q\"\n");
      printf("              Branch lengths for individual partitions will be printed to separate files\n");
      printf("              A weighted average of the branch lengths is computed by using the respective partition lengths\n");
      printf("\n");
      printf("              DEFAULT: OFF\n");
      printf("\n");
      printf("      -n      Specifies the name of the output file.\n"); 
      printf("\n");
      printf("      -p      Specify a random number seed, required in conjunction with the \"-g\" option for constraint trees\n");
      printf("\n");
      printf("      -R      read in a binary checkpoint file called ExaML_binaryCheckpoint.RUN_ID_number\n");
      printf("\n");
      printf("      -r      Pass the number of quartets to randomly sub-sample from the possible number of quartets for the given taxon set.\n");
      printf("              Only works in combination with -f q !\n");
      printf("\n");
      printf("      -s      Specify the name of the BINARY alignment data file generated by the parser component\n");
      printf("\n");
      printf("      -S      turn on memory saving option for gappy multi-gene alignments. For large and gappy datasets specify -S to save memory\n");
      printf("              This will produce slightly different likelihood values, may be a bit slower but can reduce memory consumption\n");
      printf("              from 70GB to 19GB on very large and gappy datasets\n");
      printf("\n");
      printf("      -t      Specify a user starting tree file name in Newick format\n");
      printf("\n");
      printf("      -v      Display version information\n");
      printf("\n");
      printf("      -w      FULL (!) path to the directory into which ExaML shall write its output files\n");
      printf("\n");
      printf("              DEFAULT: current directory\n");  
      printf("\n"); 
      printf("      -Y      Pass a quartet grouping file name defining four groups from which to draw quartets\n");
      printf("              The file input format must contain 4 groups in the following form:\n");
      printf("              (Chicken, Human, Loach), (Cow, Carp), (Mouse, Rat, Seal), (Whale, Frog);\n");
      printf("              Only works in combination with -f q !\n");
      printf("\n");
      
      printf("\n");
      printf("      --auto-prot=ml|bic|aic|aicc When using automatic protein model selection you can chose the criterion for selecting these models.\n");
      printf("              RAxML will test all available prot subst. models except for LG4M, LG4X and GTR-based models, with and without empirical base frequencies.\n");
      printf("              You can chose between ML score based selection and the BIC, AIC, and AICc criteria.\n");
      printf("\n");
      printf("              DEFAULT: ml\n");
      printf("\n\n\n\n");
    }
}




static void analyzeRunId(char id[128])
{
  int i = 0;

  while(id[i] != '\0')
    {    
      if(i >= 128)
	{
	  printf("Error: run id after \"-n\" is too long, it has %d characters please use a shorter one\n", i);
	  assert(0);
	}
      
      if(id[i] == '/')
	{
	  printf("Error character %c not allowed in run ID\n", id[i]);
	  assert(0);
	}


      i++;
    }

  if(i == 0)
    {
      printf("Error: please provide a string for the run id after \"-n\" \n");
      assert(0);
    }

}

static void get_args(int argc, char *argv[], analdef *adef, tree *tr)
{
  boolean   
    resultDirSet = FALSE;

  char
    resultDir[1024] = "",          
    //*optarg,
    model[1024] = "",       
    modelChar;

  double 
    likelihoodEpsilon;
  
  int        
    fOptionCount = 0,
    c,
    nameSet = 0,
    treeSet = 0,   
    modelSet = 0, 
    byteFileSet = 0,
    seedSet = 0;


  /*********** tr inits **************/ 
 
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->rateHetModel = GAMMA;
 
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->useGappedImplementation = FALSE;
  tr->saveMemory = FALSE;
  tr->constraintTree = FALSE;

  tr->fastTreeEvaluation = FALSE;

  /* tr->manyPartitions = FALSE; */

  tr->categories             = 25;

  tr->gapyness               = 0.0; 
  tr->saveBestTrees          = 0;

  tr->useMedian = FALSE;
  
  tr->autoProteinSelectionType = AUTO_ML;
  
  /********* tr inits end*************/
	
  //while(!bad_opt && ((c = mygetopt(argc,argv,"R:B:e:c:f:i:m:t:g:w:n:s:p:vhMSDa", &optind, &optarg))!=-1))
	
  static 
    int flag;
  
  while(1)
    {
      static struct 
	option long_options[2] =
	{	 	 
	  {"auto-prot",   required_argument, &flag, 1},	   	  	 	 
	  {0, 0, 0, 0}
	};
      
      int 
	option_index;
      
      flag = 0;        

      c = getopt_long(argc, argv, "R:B:Y:I:e:c:f:i:m:t:g:w:n:s:p:r:vhMSDa", long_options, &option_index);    
    
      if(c == -1)
	break;
      
      if(flag > 0)
	{
	  switch(option_index)
	    {
	    case 0:
	      {
		char 
		  *autoModels[4] = {"ml", "bic", "aic", "aicc"};

		int 
		  k;

		for(k = 0; k < 4; k++)		  
		  if(strcmp(optarg, autoModels[k]) == 0)
		    break;

		if(k == 4)
		  {
		    printf("\nError, unknown protein model selection type, you can specify one of the following selection criteria:\n\n");
		    for(k = 0; k < 4; k++)
		      printf("--auto-prot=%s\n", autoModels[k]);
		    printf("\n");
		    errorExit(-1);
		  }
		else
		  {
		    switch(k)
		      {
		      case 0:
			tr->autoProteinSelectionType = AUTO_ML;
			break;
		      case 1:
			tr->autoProteinSelectionType = AUTO_BIC;
			break;
		      case 2:
			tr->autoProteinSelectionType = AUTO_AIC;
			break;
		      case 3:
			tr->autoProteinSelectionType = AUTO_AICC;
			break;
		      default:
			assert(0);
		      }
		  }
	      }
	      break;
	    default:
	      assert(0);
	    }
	}
      else	
	switch(c)
	  {    
	  case 'Y':
	    adef->useQuartetGrouping = TRUE;	 
	    strcpy(quartetGroupingFileName, optarg);
	    break;
	  case 'r':
	    sscanf(optarg, "%lu", &(adef->numberRandomQuartets));	    
	    assert(adef->numberRandomQuartets > 0);
	    break; 
	  case 'p':
	    sscanf(optarg,"%u", &(tr->randomSeed));
	    seedSet = 1;
	    break;
	  case 'a':
	    tr->useMedian = TRUE;	
	    break;
	  case 'B':
	    sscanf(optarg,"%d", &(tr->saveBestTrees));
	    if(tr->saveBestTrees < 0)
	      {
		printf("Number of best trees to save must be greater than 0!\n");
		errorExit(-1);	 
	      }
	    break;       
	  case 's':	    
	    strcpy(byteFileName, optarg);	 	
	    byteFileSet = TRUE;
	    /*printf("%s \n", byteFileName);*/
	    break;      
	  case 'S':
	    tr->saveMemory = TRUE;
	    break;
	  case 'D':
	    tr->searchConvergenceCriterion = TRUE;	
	    break;
	  case 'R':
	    adef->useCheckpoint = TRUE;
	    strcpy(binaryCheckpointInputName, optarg);
	    break;          
	  case 'I':
	    sscanf(optarg, "%lu", &(adef->quartetCkpInterval));
	    break;
	  case 'M':
	    adef->perGeneBranchLengths = TRUE;
	    break;                                 
	  case 'e':
	    sscanf(optarg,"%lf", &likelihoodEpsilon);
	    adef->likelihoodEpsilon = likelihoodEpsilon;
	    break;    	    
	  case 'v':
	    printVersionInfo();
	    errorExit(0);	    
	  case 'h':
	    printREADME();
	    errorExit(0);     
	  case 'c':
	    sscanf(optarg, "%d", &tr->categories);
	    break;     
	  case 'f':
	    sscanf(optarg, "%c", &modelChar);
	    fOptionCount++;
	    if(fOptionCount > 1) 
	      {
		printf("\nError: only one of the various \"-f \" options can be used per ExaML run!\n");
		printf("They are mutually exclusive! exiting ...\n\n");
		errorExit(-1);
	      }
	    switch(modelChar)
	      {	 
	      case 'e':
		adef->mode = TREE_EVALUATION;
		tr->fastTreeEvaluation = TRUE;
		break;
	      case 'E':
		adef->mode = TREE_EVALUATION;
		tr->fastTreeEvaluation = FALSE;
		break;
	      case 'd':
		adef->mode = BIG_RAPID_MODE;
		tr->doCutoff = TRUE;
		break;	  
	      case 'o':
		adef->mode = BIG_RAPID_MODE;
		tr->doCutoff = FALSE;
		break;	    	  	  	     
	      case 'q':
		adef->mode = QUARTET_CALCULATION;
		break;	
	      default:
		{
		  if(processID == 0)
		    {
		      printf("Error select one of the following algorithms via -f :\n");
		      printMinusFUsage();
		    }
		  errorExit(-1);
		}
	      }
	    break;
	  case 'i':
	    sscanf(optarg, "%d", &adef->initial);
	    adef->initialSet = TRUE;
	    break;
	  case 'n':
	    strcpy(run_id,optarg);
	    analyzeRunId(run_id);
	    nameSet = 1;
	    break;
	  case 'w':
	    strcpy(resultDir, optarg);
	    resultDirSet = TRUE;
	    break;
	  case 't':
	    strcpy(tree_file, optarg);       
	    treeSet = 1;       
	    break;
	  case 'g':
	    strcpy(tree_file, optarg);       
	    treeSet = 1;       
	    tr->constraintTree = TRUE;
	    break;	
	  case 'm':
	    strcpy(model,optarg);
	    if(modelExists(model, tr) == 0)
	      {
		if(processID == 0)
		  {
		    printf("Rate heterogeneity Model %s does not exist\n\n", model);               
		    printf("For per site rates (called CAT in previous versions) use: PSR\n");	
		    printf("For GAMMA use: GAMMA\n");		
		  }
		errorExit(-1);
	      }
	    else
	      modelSet = 1;
	    break;     
	  default:
	    errorExit(-1);
	  }
    }
  
  if(adef->useQuartetGrouping && adef->mode != QUARTET_CALCULATION)
    {
      if(processID == 0)
	printf("\nError, you must specify \"-Y quartetGroupingFileName\" in combination with \"-f q\"\n");
      errorExit(-1);
    }

  if(adef->numberRandomQuartets > 0 && adef->mode != QUARTET_CALCULATION)
    {
       if(processID == 0)
	printf("\nError, you must specify \"-r randomQuartetNumber\" in combination with \"-f q\"\n");
      errorExit(-1);
    }

  if((adef->numberRandomQuartets > 0) && (adef->useQuartetGrouping))
    {
      if(processID == 0)
	printf("\nError, you must specify either \"-r randomQuartetNumber\" or \"-Y quartetGroupingFileName\"\n");
      errorExit(-1);
    }  

  if(tr->constraintTree)
    {
      if(!seedSet && processID == 0)
	{
	  printf("\nError, you must specify a random number seed via \"-p\" when  using a constraint\n");
	  printf("tree via \"-g\" \n");
	  errorExit(-1);
	}
    }

  if(!byteFileSet)
    {
      if(processID == 0)
	printf("\nError, you must specify a binary format data file with the \"-s\" option\n");
      errorExit(-1);
    }

  if(!modelSet)
    {
      if(processID == 0)
	printf("\nError, you must specify a model of rate heterogeneity with the \"-m\" option\n");
      errorExit(-1);
    }

  if(!nameSet)
    {
      if(processID == 0)
	printf("\nError: please specify a name for this run with -n\n");
      errorExit(-1);
    }

  if(!treeSet && !adef->useCheckpoint)
    {
      if(processID == 0)
	{
	  printf("\nError: please either specify a starting tree for this run with -t\n");
	  printf("or re-start the run from a checkpoint with -R\n");
	}
      
      errorExit(-1);
    }
  
   {

    const 
      char *separator = "/";

    if(resultDirSet)
      {
	char 
	  dir[1024] = "";
	

	if(resultDir[0] != separator[0])
	  strcat(dir, separator);
	
	strcat(dir, resultDir);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	strcpy(workdir, dir);
      }
    else
      {
	char 
	  dir[1024] = "",
	  *result = getcwd(dir, sizeof(dir));
	
	assert(result != (char*)NULL);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	
	strcpy(workdir, dir);		
      }
   }

  return;
}




void errorExit(int e)
{
  MPI_Finalize();

  exit(e);
}



static void makeFileNames(void)
{
  int 
    infoFileExists = 0;
    
  strcpy(resultFileName,       workdir);
  strcpy(logFileName,          workdir);  
  strcpy(infoFileName,         workdir);
  strcpy(treeFileName,         workdir);
  strcpy(binaryCheckpointName, workdir);
  strcpy(modelFileName, workdir);
  strcpy(quartetFileName,         workdir);

  strcat(resultFileName,       "ExaML_result.");
  strcat(logFileName,          "ExaML_log.");  
  strcat(infoFileName,         "ExaML_info.");
  strcat(binaryCheckpointName, "ExaML_binaryCheckpoint.");
  strcat(modelFileName,        "ExaML_modelFile.");
  strcat(treeFileName,         "ExaML_TreeFile.");
  strcat(quartetFileName,      "ExaML_quartets.");
  
  strcat(resultFileName,       run_id);
  strcat(logFileName,          run_id);  
  strcat(infoFileName,         run_id); 
  strcat(binaryCheckpointName, run_id);
  strcat(modelFileName,        run_id);
  strcat(treeFileName,         run_id);
  strcat(quartetFileName,         run_id);

  infoFileExists = filexists(infoFileName);

  if(infoFileExists)
    {
      if(processID == 0)
	{
	  printf("ExaML output files with the run ID <%s> already exist \n", run_id);
	  printf("in directory %s ...... exiting\n", workdir);
	}

      errorExit(-1);	
    }
}




 




/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/


static void printModelAndProgramInfo(tree *tr, analdef *adef, int argc, char *argv[])
{
  if(processID == 0)
    {
      int i, model;
      FILE *infoFile = myfopen(infoFileName, "ab");
      char modelType[128];

      
      if(tr->useMedian)
	strcpy(modelType, "GAMMA with Median");
      else
	strcpy(modelType, "GAMMA");   
     
      printBoth(infoFile, "\n\nThis is %s version %s released by Alexandros Stamatakis, Andre Aberer, and Alexey Kozlov in %s.\n\n",  programName, programVersion, programDate);
                     
      printBoth(infoFile, "\nAlignment has %zu distinct alignment patterns\n\n",  tr->originalCrunchedLength);
                 
      printBoth(infoFile, "Proportion of gaps and completely undetermined characters in this alignment: %3.2f%s\n", 100.0 * tr->gapyness, "%");
      
      switch(adef->mode)
	{	
	case  BIG_RAPID_MODE:	 
	  printBoth(infoFile, "\nExaML rapid hill-climbing mode\n\n");
	  break;
	case TREE_EVALUATION:
	  printBoth(infoFile, "\nExaML %s tree evaluation mode\n\n", (tr->fastTreeEvaluation)?"fast":"slow");
	  break;
	case QUARTET_CALCULATION:
	  printBoth(infoFile, "\nExaML quartet evaluation mode\n\n");
	  break;
	default:
	  assert(0);
	}
     	  
      if(adef->perGeneBranchLengths)
	printBoth(infoFile, "Using %d distinct models/data partitions with individual per partition branch length optimization\n\n\n", tr->NumberOfModels);
      else
	printBoth(infoFile, "Using %d distinct models/data partitions with joint branch length optimization\n\n\n", tr->NumberOfModels);	
	      
      printBoth(infoFile, "All free model parameters will be estimated by ExaML\n");
           	
      if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
	printBoth(infoFile, "%s model of rate heteorgeneity, ML estimate of alpha-parameter\n\n", modelType);
      else
	{
	  printBoth(infoFile, "ML estimate of %d per site rate categories\n\n", tr->categories);
	 
	}               
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  printBoth(infoFile, "Partition: %d\n", model);
	  printBoth(infoFile, "Alignment Patterns: %d\n", tr->partitionData[model].upper - tr->partitionData[model].lower);
	  printBoth(infoFile, "Name: %s\n", tr->partitionData[model].partitionName);
	  
	  switch(tr->partitionData[model].dataType)
	    {
	    case DNA_DATA:
	      printBoth(infoFile, "DataType: DNA\n");	     
	      printBoth(infoFile, "Substitution Matrix: GTR\n");
	      if(tr->partitionData[model].optimizeBaseFrequencies)
		printBoth(infoFile, "ML optimization of base frequencies\n");
	      break;
	    case AA_DATA:
	      assert(tr->partitionData[model].protModels >= 0 && tr->partitionData[model].protModels < NUM_PROT_MODELS);
	      printBoth(infoFile, "DataType: AA\n");	      
	      printBoth(infoFile, "Substitution Matrix: %s\n", protModels[tr->partitionData[model].protModels]);
	      if(!tr->partitionData[model].optimizeBaseFrequencies)
		printBoth(infoFile, "Using %s Base Frequencies\n", (tr->partitionData[model].protFreqs == 1)?"empirical":"fixed");	     
	      else		
		printBoth(infoFile, "ML optimization of base frequencies\n");
	      break;
	    case BINARY_DATA:
	      printBoth(infoFile, "DataType: BINARY/MORPHOLOGICAL\n");	      
	      printBoth(infoFile, "Substitution Matrix: Uncorrected\n");
	      break;
	    
	      /*
		case SECONDARY_DATA:
		printBoth(infoFile, "DataType: SECONDARY STRUCTURE\n");	     
		printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
		break;
		case SECONDARY_DATA_6:
		printBoth(infoFile, "DataType: SECONDARY STRUCTURE 6 STATE\n");	     
		printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
		break;
		case SECONDARY_DATA_7:
		printBoth(infoFile, "DataType: SECONDARY STRUCTURE 7 STATE\n");	      
		printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
		break;
		case GENERIC_32:
		printBoth(infoFile, "DataType: Multi-State with %d distinct states in use (maximum 32)\n",tr->partitionData[model].states);		  
		switch(tr->multiStateModel)
		{
		case ORDERED_MULTI_STATE:
		printBoth(infoFile, "Substitution Matrix: Ordered Likelihood\n");
		break;
		case MK_MULTI_STATE:
		printBoth(infoFile, "Substitution Matrix: MK model\n");
		break;
		case GTR_MULTI_STATE:
		printBoth(infoFile, "Substitution Matrix: GTR\n");
		break;
		default:
		assert(0);
		}
		break;
		case GENERIC_64:
		printBoth(infoFile, "DataType: Codon\n");		  
		break;	
	      */
	    default:
	      assert(0);
	    }
	  printBoth(infoFile, "\n\n\n");
	}
      
      printBoth(infoFile, "\n");

      printBoth(infoFile, "ExaML was called as follows:\n\n");
      for(i = 0; i < argc; i++)
	printBoth(infoFile,"%s ", argv[i]);
      printBoth(infoFile,"\n\n\n");

      fclose(infoFile);
    }
}

void printResult(tree *tr, analdef *adef, boolean finalPrint)
{
  if(processID == 0)
    {
      FILE *logFile;
      char temporaryFileName[1024] = "";
      
      strcpy(temporaryFileName, resultFileName);
      
      switch(adef->mode)
	{    
	case TREE_EVALUATION:
	  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, SUMMARIZE_LH, FALSE, FALSE);
	  
	  logFile = myfopen(temporaryFileName, "wb");
	  fprintf(logFile, "%s", tr->tree_string);
	  fclose(logFile);
	  
	  if(adef->perGeneBranchLengths)
	    printTreePerGene(tr, adef, temporaryFileName, "wb");
	  break;
	case BIG_RAPID_MODE:     
	  if(finalPrint)
	    {
	      switch(tr->rateHetModel)
		{
		case GAMMA:
		case GAMMA_I:
		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint,
			      SUMMARIZE_LH, FALSE, FALSE);
		  
		  logFile = myfopen(temporaryFileName, "wb");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);
		  
		  if(adef->perGeneBranchLengths)
		    printTreePerGene(tr, adef, temporaryFileName, "wb");
		  break;
		case CAT:
		  /*Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef,
		    NO_BRANCHES, FALSE, FALSE);*/
		  
		  
		  
		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE,
			      TRUE, SUMMARIZE_LH, FALSE, FALSE);
		  
		  
		  
		  
		  logFile = myfopen(temporaryFileName, "wb");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);
		  
		  break;
		default:
		  assert(0);
		}
	    }
	  else
	    {
	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint,
			  NO_BRANCHES, FALSE, FALSE);
	      logFile = myfopen(temporaryFileName, "wb");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile);
	    }    
	  break;
	default:
	  printf("FATAL ERROR call to printResult from undefined STATE %d\n", adef->mode);
	  exit(-1);
	  break;
	}
    }
}








void printLog(tree *tr)
{
  if(processID == 0)
    {
      FILE *logFile;
      double t;
      
      t = gettime() - masterTime;
      
      logFile = myfopen(logFileName, "ab");
      
      /* printf("%f %1.40f\n", t, tr->likelihood); */

      fprintf(logFile, "%f %f\n", t, tr->likelihood);
      
      fclose(logFile);
    }
	     
}









void getDataTypeString(tree *tr, int model, char typeOfData[1024])
{
  switch(tr->partitionData[model].dataType)
    {
    case AA_DATA:
      strcpy(typeOfData,"AA");
      break;
    case DNA_DATA:
      strcpy(typeOfData,"DNA");
      break;
    case BINARY_DATA:
      strcpy(typeOfData,"BINARY/MORPHOLOGICAL");
      break;
    case SECONDARY_DATA:
      strcpy(typeOfData,"SECONDARY 16 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_6:
      strcpy(typeOfData,"SECONDARY 6 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_7:
      strcpy(typeOfData,"SECONDARY 7 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case GENERIC_32:
      strcpy(typeOfData,"Multi-State");
      break;
    case GENERIC_64:
      strcpy(typeOfData,"Codon"); 
      break;
    default:
      assert(0);
    }
}
static void printRatesDNA_BIN(int n, double *r, char **names, char *fileName)
{
  int i, j, c;

  for(i = 0, c = 0; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	{
	  if(i == n - 2 && j == n - 1)
	    printBothOpenDifferentFile(fileName, "rate %s <-> %s: %f\n", names[i], names[j], 1.0);
	  else
	    printBothOpenDifferentFile(fileName, "rate %s <-> %s: %f\n", names[i], names[j], r[c]);
	  c++;
	}
    }
}

static void printRatesRest(int n, double *r, char **names, char *fileName)
{
  int i, j, c;

  for(i = 0, c = 0; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	{
	  printBothOpenDifferentFile(fileName, "rate %s <-> %s: %f\n", names[i], names[j], r[c]);
	  c++;
	}
    }
}
static double branchLength(int model, double *z, tree *tr)
{
  double x;
  
  x = z[model];
  assert(x > 0);
  if (x < zmin) 
    x = zmin;  
  
 
  assert(x <= zmax);
  
  x = -log(x);
  
  return x;
}


static double treeLengthRec(nodeptr p, tree *tr, int model)
{  
  double 
    x = branchLength(model, p->z, tr);

  if(isTip(p->number, tr->mxtips))  
    return x;    
  else
    {
      double acc = 0;
      nodeptr q;                
     
      q = p->next;      

      while(q != p)
	{
	  acc += treeLengthRec(q->back, tr, model);
	  q = q->next;
	}

      return acc + x;
    }
}

static double treeLength(tree *tr, int model)
{  
  return treeLengthRec(tr->start->back, tr, model);
}

static void printFreqs(int n, double *f, char **names, char *fileName)
{
  int k;

  for(k = 0; k < n; k++)
    printBothOpenDifferentFile(fileName, "freq pi(%s): %f\n", names[k], f[k]);
}

static void printModelParams(tree *tr, analdef *adef, int treeIteration)
{
  int
    model;

  double
    *f = (double*)NULL,
    *r = (double*)NULL;

  char 
    fileName[2048],
    buf[64];
  
  
  strcpy(fileName, modelFileName);

  if(treeIteration >= 0)
    {
      strcat(fileName, ".");
      sprintf(buf, "%d", treeIteration);
      strcat(fileName, buf);
    }

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      double tl;
      char typeOfData[1024];

      getDataTypeString(tr, model, typeOfData);      

      printBothOpenDifferentFile(fileName, "\n\n");

      printBothOpenDifferentFile(fileName, "Model Parameters of Partition %d, Name: %s, Type of Data: %s\n",
				 model, tr->partitionData[model].partitionName, typeOfData);
      
      if(tr->rateHetModel == GAMMA)
	printBothOpenDifferentFile(fileName, "alpha: %f\n", tr->partitionData[model].alpha);
     

      if(adef->perGeneBranchLengths)
	tl = treeLength(tr, model);
      else
	tl = treeLength(tr, 0);

      printBothOpenDifferentFile(fileName, "Tree-Length: %f\n", tl);

      f = tr->partitionData[model].frequencies;
      r = tr->partitionData[model].substRates;

      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  {
	    char *freqNames[20] = {"A", "R", "N ","D", "C", "Q", "E", "G",
				   "H", "I", "L", "K", "M", "F", "P", "S",
				   "T", "W", "Y", "V"};

	     if(tr->partitionData[model].protModels == LG4M || tr->partitionData[model].protModels == LG4X)
	      {
		int 
		  k;
		
		for(k = 0; k < 4; k++)
		  {
		    printBothOpenDifferentFile(fileName, "LGM %d\n", k);
		    printRatesRest(20, tr->partitionData[model].substRates_LG4[k], freqNames, fileName);
		    printBothOpenDifferentFile(fileName, "\n");
		    printFreqs(20, tr->partitionData[model].frequencies_LG4[k], freqNames, fileName);
		  }
	      }

	    printRatesRest(20, r, freqNames, fileName);
	    printBothOpenDifferentFile(fileName, "\n");
	    printFreqs(20, f, freqNames, fileName);
	  }
	  break;
	case DNA_DATA:
	  {
	    char *freqNames[4] = {"A", "C", "G", "T"};

	    printRatesDNA_BIN(4, r, freqNames, fileName);
	    printBothOpenDifferentFile(fileName, "\n");
	    printFreqs(4, f, freqNames, fileName);
	  }
	  break;
	case BINARY_DATA:
	  {
	    char *freqNames[2] = {"0", "1"};

	    printRatesDNA_BIN(2, r, freqNames, fileName);
	    printBothOpenDifferentFile(fileName, "\n");
	    printFreqs(2, f, freqNames, fileName);
	  }
	  break;
	default:
	  assert(0);
	}

      printBothOpenDifferentFile(fileName, "\n");
    }

  printBothOpenDifferentFile(fileName, "\n");
}


static void finalizeInfoFile(tree *tr, analdef *adef)
{
  if(processID == 0)
    {
      double t;

      t = gettime() - masterTime;
      accumulatedTime = accumulatedTime + t;

      switch(adef->mode)
	{	
	case  BIG_RAPID_MODE:	 
	  printBothOpen("\n\nOverall Time for 1 Inference %f\n", t);
	  printBothOpen("\nOverall accumulated Time (in case of restarts): %f\n\n", accumulatedTime);
	  printBothOpen("Likelihood   : %f\n", tr->likelihood);
	  printBothOpen("\n\n");	  	  
	  printBothOpen("Model parameters written to:           %s\n", modelFileName);
	  printBothOpen("Final tree written to:                 %s\n", resultFileName);
	  printBothOpen("Execution Log File written to:         %s\n", logFileName);
	  printBothOpen("Execution information file written to: %s\n",infoFileName);	
	  break;
	case TREE_EVALUATION:	
	  printBothOpen("\n\nOverall Time for evaluating the likelihood of %d trees: %f secs\n\n", tr->numberOfTrees, t); 
	  printBothOpen("\n\nThe model parameters of the trees have been written to files called %s.i\n", modelFileName);
	  printBothOpen("where i is the number of the tree\n\n");
	  printBothOpen("Note that, in case of a restart from a checkpoint, some tree model files will have been produced by previous runs!\n\n");
	  printBothOpen("The trees with branch lengths have been written to file: %s\n", treeFileName);
	  printBothOpen("They are in the same order as in the input file!\n\n");
	  break;
	case QUARTET_CALCULATION:
	  printBothOpen("\n\nOverall quartet computation time: %f secs\n", t);
	  printBothOpen("\nAll quartets and corresponding likelihoods written to file %s\n", quartetFileName);
	  break;
	default:
	  assert(0);
	}

	 
    }

}


/************************************************************************************/


static int iterated_bitcount(unsigned int n)
{
    int 
      count=0;    
    
    while(n)
      {
        count += n & 0x1u ;    
        n >>= 1 ;
      }
    
    return count;
}

/*static char bits_in_16bits [0x1u << 16];*/

static void compute_bits_in_16bits(char *bits_in_16bits)
{
    unsigned int i;    
    
    /* size is 65536 */

    for (i = 0; i < (0x1u<<16); i++)
        bits_in_16bits[i] = iterated_bitcount(i);       

    return ;
}

unsigned int precomputed16_bitcount (unsigned int n, char *bits_in_16bits)
{
  /* works only for 32-bit int*/
    
    return bits_in_16bits [n         & 0xffffu]
        +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}


static void clean_MPI_Exit(void)
{
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}

static void error_MPI_Exit(void)
{
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  exit(1);
}


static void initializePartitions(tree *tr)
{ 
  size_t
    i,
    len, 
    j,    
    width;

  int
    model, 
    maxCategories;

  compute_bits_in_16bits(tr->bits_in_16bits);

  maxCategories = tr->maxCategories;

  for(model = 0; model < tr->NumberOfModels; model++)
    {    
      const partitionLengths 
	*pl = getPartitionLengths(&(tr->partitionData[model])); 

      //must already be set as a consequence of alloc in function readPartitions
      //and the subsequent copy of bf->partitions into tr->partitions!
      assert(tr->partitionData[model].partitionName != (char*)NULL);

      //printf("Partition name %s\n", tr->partitionData[model].partitionName);

      width = tr->partitionData[model].width;
	
      /* 
	 globalScaler needs to be 2 * tr->mxtips such that scalers of inner AND tip nodes can be added without a case switch
	 to this end, it must also be initialized with zeros -> calloc
       */
      
      len = 2 * tr->mxtips; 
      tr->partitionData[model].globalScaler       = (unsigned int *)calloc(len, sizeof(unsigned int));

#ifdef _USE_OMP
      tr->partitionData[model].threadGlobalScaler = (unsigned int**) calloc(tr->nThreads, sizeof(unsigned int*));

      tr->partitionData[model].reductionBuffer 	  = (double*) calloc(tr->nThreads, sizeof(double));
      tr->partitionData[model].reductionBuffer2   = (double*) calloc(tr->nThreads, sizeof(double));

      int 
	t;
      
      for (t = 0; t < tr->maxThreadsPerModel; ++t)
	{
	  Assign*
	    pAss = tr->partThreadAssigns[model * tr->maxThreadsPerModel + t];

	  if (pAss)
	    {
	      int
		tid = pAss->procId;

	      tr->partitionData[model].threadGlobalScaler[tid]    = (unsigned int *)calloc(len, sizeof(unsigned int));
	    }
	}
#endif

      tr->partitionData[model].left              = (double *)malloc_aligned(pl->leftLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[model].right             = (double *)malloc_aligned(pl->rightLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[model].EIGN              = (double*)malloc(pl->eignLength * sizeof(double));
      tr->partitionData[model].EV                = (double*)malloc_aligned(pl->evLength * sizeof(double));
      tr->partitionData[model].EI                = (double*)malloc(pl->eiLength * sizeof(double));
      
      tr->partitionData[model].substRates        = (double *)malloc(pl->substRatesLength * sizeof(double));


      //must already be set as a consequence of alloc in function readPartitions
      //and the subsequent copy of bf->partitions into tr->partitions!
      assert(tr->partitionData[model].frequencies != (double*)NULL);
      //tr->partitionData[model].frequencies       = (double*)malloc(pl->frequenciesLength * sizeof(double));

     

      tr->partitionData[model].freqExponents     = (double*)malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[model].empiricalFrequencies       = (double*)malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[model].tipVector         = (double *)malloc_aligned(pl->tipVectorLength * sizeof(double));


      if(tr->partitionData[model].protModels == LG4M || tr->partitionData[model].protModels == LG4X)      
	{	  	  
	  int 
	    k;
	  
	  for(k = 0; k < 4; k++)
	    {	    
	      tr->partitionData[model].rawEIGN_LG4[k]              = (double*)malloc(pl->eignLength * sizeof(double));
	      tr->partitionData[model].EIGN_LG4[k]              = (double*)malloc(pl->eignLength * sizeof(double));
	      tr->partitionData[model].EV_LG4[k]                = (double*)malloc_aligned(pl->evLength * sizeof(double));
	      tr->partitionData[model].EI_LG4[k]                = (double*)malloc(pl->eiLength * sizeof(double));
	      tr->partitionData[model].substRates_LG4[k]        = (double *)malloc(pl->substRatesLength * sizeof(double));
	      tr->partitionData[model].frequencies_LG4[k]       = (double*)malloc(pl->frequenciesLength * sizeof(double));
	      tr->partitionData[model].tipVector_LG4[k]         = (double *)malloc_aligned(pl->tipVectorLength * sizeof(double));
	    }
	}


      tr->partitionData[model].symmetryVector    = (int *)malloc(pl->symmetryVectorLength  * sizeof(int));
      tr->partitionData[model].frequencyGrouping = (int *)malloc(pl->frequencyGroupingLength  * sizeof(int));
      
      tr->partitionData[model].perSiteRates      = (double *)malloc(sizeof(double) * tr->maxCategories);
            
      //      tr->partitionData[model].nonGTR = FALSE; 
      //      tr->partitionData[model].optimizeBaseFrequencies = FALSE; 
      

      //tr->partitionData[model].gammaRates = (double*)malloc(sizeof(double) * 4);

      tr->partitionData[model].xVector = (double **)malloc(sizeof(double*) * tr->mxtips);   
      	
      for(j = 0; j < (size_t)tr->mxtips; j++)	        	  	  	  	 
	  tr->partitionData[model].xVector[j]   = (double*)NULL;   

      tr->partitionData[model].xSpaceVector = (size_t *)calloc(tr->mxtips, sizeof(size_t));  

#ifdef __MIC_NATIVE
      tr->partitionData[model].mic_EV                = (double*)malloc_aligned(4 * pl->evLength * sizeof(double));
      tr->partitionData[model].mic_tipVector         = (double*)malloc_aligned(4 * pl->tipVectorLength * sizeof(double));
      tr->partitionData[model].mic_umpLeft           = (double*)malloc_aligned(4 * pl->tipVectorLength * sizeof(double));
      tr->partitionData[model].mic_umpRight           = (double*)malloc_aligned(4 * pl->tipVectorLength * sizeof(double));

      /* for Xeon Phi, sumBuffer must be padded to the multiple of 8 (because of site blocking in kernels) */
      const int padded_width = GET_PADDED_WIDTH(width);
      const int span = (size_t)(tr->partitionData[model].states) *
              discreteRateCategories(tr->rateHetModel);

      tr->partitionData[model].sumBuffer = (double *)malloc_aligned(padded_width *
									   span * sizeof(double));

      /* fill padding entries with 1. (will be corrected for with zero site weights in wgt) */
      {
          int k;
          for (k = width*span; k < padded_width*span; ++k)
              tr->partitionData[model].sumBuffer[k] = 1.;
      }
#else
      tr->partitionData[model].sumBuffer = (double *)malloc_aligned(width *
									   (size_t)(tr->partitionData[model].states) *
									   discreteRateCategories(tr->rateHetModel) *
									   sizeof(double));
#endif

      /* tr->partitionData[model].wgt = (int *)malloc_aligned(width * sizeof(int));	   */

      /* rateCategory must be assigned using calloc() at start up there is only one rate category 0 for all sites */

      if(width > 0 && tr->saveMemory)
	{
	  tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
	  
	  len = tr->partitionData[model].gapVectorLength * 2 * tr->mxtips; 
	  tr->partitionData[model].gapVector = (unsigned int*)calloc(len, sizeof(unsigned int));	  	    	  	  
	    
	  tr->partitionData[model].gapColumn = (double *)malloc_aligned(((size_t)tr->mxtips) *								      
									       ((size_t)(tr->partitionData[model].states)) *
									       discreteRateCategories(tr->rateHetModel) * sizeof(double));
	}
      else
	{
	   tr->partitionData[model].gapVectorLength = 0;
	    
	   tr->partitionData[model].gapVector = (unsigned int*)NULL; 	  	    	   
	    
	   tr->partitionData[model].gapColumn = (double*)NULL;	    	    	   
	}              
    }


  /* set up the averaged frac changes per partition such that no further reading accesses to aliaswgt are necessary
     and we can free the array for the GAMMA model */
 
  {      
    /* definitions: 
       sizeof(short) <= sizeof(int) <= sizeof(long)
       size_t defined by address space (here: 64 bit). 
       
       size_t + MPI is a bad idea: in the mpi2.2 standard, they do
	 not mention it once.
    */

    unsigned long 
      *modelWeights = (unsigned long*) calloc(tr->NumberOfModels, sizeof(unsigned long)); 
    
    size_t
      wgtsum = 0;  

    /* determine my weights per partition    */
    for(model = 0; model < tr->NumberOfModels; model++)      
      {
	const pInfo 
	  partition =  tr->partitionData[model] ; 
	   
	size_t 
	  i = 0; 
	   
	for(i = 0; i < partition.width; ++i)
	  modelWeights[model] += (long) partition.wgt[i]; 
      }
    MPI_Allreduce(MPI_IN_PLACE, modelWeights, tr->NumberOfModels, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD); 
       
    /* determine sum */
    for(model = 0; model < tr->NumberOfModels; ++model)
      wgtsum += modelWeights[model]; 

    for(model = 0; model < tr->NumberOfModels; model++)      	
      {
	tr->partitionWeights[model]       = (double)modelWeights[model];
	tr->partitionContributions[model] = ((double)modelWeights[model]) / ((double)wgtsum); 
      }
       
    free(modelWeights);
  }

  /* initialize gap bit vectors at tips when memory saving option is enabled */
  
  if(tr->saveMemory)
    {
      for(model = 0; model <tr->NumberOfModels; model++)
	{
	  int        
	    undetermined = getUndetermined(tr->partitionData[model].dataType);
	  	 
	  width =  tr->partitionData[model].width;
	    
	  if(width > 0)
	    {	   	    	      	    	     
	      for(j = 1; j <= (size_t)(tr->mxtips); j++)
		for(i = 0; i < width; i++)
		  if(tr->partitionData[model].yVector[j][i] == undetermined)
		    tr->partitionData[model].gapVector[tr->partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];	    
	    }     
	}
    }
}



static void initializeTree(tree *tr, analdef *adef)
{
  size_t 
    i ;

  if(adef->perGeneBranchLengths)
    tr->numBranches = tr->NumberOfModels;
  else
    tr->numBranches = 1;


  if(NUM_BRANCHES < tr->numBranches)
    {
      if(processID == 0 )
	printf("You have specified per-partition branch lengths (-M option) with %d  models. \n\
Please set #define NUM_BRANCHES in axml.h to %d (or higher) and recompile %s\n", 
	       tr->NumberOfModels,tr->NumberOfModels, programName );
      error_MPI_Exit();
    }


  /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
     let's not worry about this right now, because it is indeed ExaML-specific */

  tr->executeModel   = (boolean *)calloc( tr->NumberOfModels, sizeof(boolean));
  
  for(i = 0; i < (size_t)tr->NumberOfModels; i++)
    tr->executeModel[i] = TRUE;
  
  setupTree(tr); 
  
  if(tr->searchConvergenceCriterion && processID == 0)
    {                     
      tr->bitVectors = initBitVector(tr->mxtips, &(tr->vLength));
      tr->h = initHashTable(tr->mxtips * 4);     
    }
  
  for(i = 1; i <= (size_t)tr->mxtips; i++)
    addword(tr->nameList[i], tr->nameHash, i);
   
  initializePartitions(tr);

  initModel(tr);
}


static int getNumberOfTrees(char *fileName, boolean getOffsets, exa_off_t *treeOffsets)
{
  FILE 
    *f = myfopen(fileName, "r");

  int 
    trees = 0,
    ch;

  if(getOffsets)
    treeOffsets[trees] = 0;

  while((ch = fgetc(f)) != EOF)
    {
      if(ch == ';')
	{
	  trees++;
	  if(getOffsets)	    
	    treeOffsets[trees] = exa_ftell(f) + 1;	 	      	   
	}
    }

  assert(trees > 0);

  fclose(f);

  return trees;
}

static void optimizeTrees(tree *tr, analdef *adef)
{
  exa_off_t
    *treeOffsets;

  int 
    i;   

  tr->numberOfTrees = getNumberOfTrees(tree_file, FALSE, (exa_off_t *)NULL);
  
  if(processID == 0)
    accumulatedTime = 0.0;

  treeOffsets = (exa_off_t *)malloc(sizeof(exa_off_t) * (tr->numberOfTrees + 1));

  tr->likelihoods = (double *)malloc(sizeof(double) * tr->numberOfTrees);
  tr->treeStrings = (char   *)malloc(sizeof(char) * (size_t)tr->treeStringLength * (size_t)tr->numberOfTrees);

  getNumberOfTrees(tree_file, TRUE, treeOffsets);
  
  if(processID == 0)   
    printBothOpen("\n\nFound %d trees to evaluate\n\n", tr->numberOfTrees);
  
  i = 0;

  if(adef->useCheckpoint)
    {      
      restart(tr, adef);       		   	    
	  
      i = ckp.treeIteration;
	       
      if(tr->fastTreeEvaluation && i > 0)	
	treeEvaluate(tr, 2);	
      else
	modOpt(tr, 0.1, adef, i);
      
      tr->likelihoods[i] = tr->likelihood;
      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);
      memcpy(&(tr->treeStrings[(size_t)tr->treeStringLength * (size_t)i]), tr->tree_string, sizeof(char) * tr->treeStringLength);
      

      if(processID == 0)
	printModelParams(tr, adef, i);
      
      i++;
    }
       
  for(; i < tr->numberOfTrees; i++)
    {     
      FILE 
	*treeFile = myfopen(tree_file, "rb");
	    
      if(exa_fseek(treeFile, treeOffsets[i], SEEK_SET) != 0)
	assert(0);

      tr->likelihood = unlikely;
   
      treeReadLen(treeFile, tr, FALSE, FALSE, FALSE);
               
      fclose(treeFile);
 
      tr->start = tr->nodep[1];     
	  
      if(i > 0)
	resetBranches(tr);
      
      evaluateGeneric(tr, tr->start, TRUE);	
      	  
      if(tr->fastTreeEvaluation && i > 0)
	{
	  ckp.state = MOD_OPT;	  	 

	  ckp.treeIteration = i;
	  
	  writeCheckpoint(tr, adef);

	  treeEvaluate(tr, 2);
	}
      else
	{
	  treeEvaluate(tr, 1);      
	  modOpt(tr, 0.1, adef, i);
	}
      
      tr->likelihoods[i] = tr->likelihood;
      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);
      memcpy(&(tr->treeStrings[(size_t)tr->treeStringLength * (size_t)i]), tr->tree_string, sizeof(char) * tr->treeStringLength);

      if(processID == 0)
	printModelParams(tr, adef, i);
    }

  if(processID == 0)
    {
      FILE 
	*f = myfopen(treeFileName, "w");
      
      for(i = 0; i < tr->numberOfTrees; i++)
	{
	  printBothOpen("Likelihood tree %d: %f \n", i, tr->likelihoods[i]);    
	  fprintf(f, "%s", &(tr->treeStrings[(size_t)tr->treeStringLength * (size_t)i]));
	}
      
      fclose(f);
    }
}










static void readByteFile (tree *tr, int commRank, int commSize )
{
  /* read stuff that is cheap; do not change the order! */
  ByteFile 
    *bFile = NULL; 
  
  initializeByteFile(&bFile, byteFileName); 
  readHeader(bFile);
  readTaxa(bFile);
  readPartitions(bFile); 

  /* calculate optimal distribution of data */
  PartitionAssignment 
    *pAss = NULL; 
  
  initializePartitionAssignment(&pAss, bFile->partitions, bFile->numPartitions, commSize); 
  assign(pAss);

  if(commRank == 0 )
    {
      printf("\n"); 
      printAssignments(pAss); 
      printf("\n"); 
      printLoad(pAss); 
      printf("\n");
    }

  /* now the data of this process is in this struct */
  readMyData(bFile,pAss, commRank );

  /* carry over the information to the tree */
  initializeTreeFromByteFile(bFile, tr); 
  
  /* just fills up tr->partAssigns that contains the representation of
     the assignment that we will need */
  copyAssignmentInfoToTree(pAss, tr);

  deletePartitionAssignment(pAss);
  deleteByteFile(bFile);
}

#ifdef _USE_OMP
void allocateXVectors(tree* tr)
{
  nodeptr
    p = tr->start,
    q = p->back;

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  tr->td[0].count = 1;

  computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches, FALSE);

  traversalInfo
    *ti = tr->td[0].ti;

  int
    i,
    model;

  for(i = 1; i < tr->td[0].count; i++)
    {
      traversalInfo *tInfo = &ti[i];

      /* now loop over all partitions for nodes p, q, and r of the current traversal vector entry */

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  /* printf("new view on model %d with width %d\n", model, width);  */

	  size_t
	    width  = (size_t)tr->partitionData[model].width;

	  double
	    *x3_start = tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1];

	  size_t
	    rateHet = discreteRateCategories(tr->rateHetModel),

	    /* get the number of states in the data stored in partition model */

	    states = (size_t)tr->partitionData[model].states,

	    /* get the length of the current likelihood array stored at node p. This is
	       important mainly for the SEV-based memory saving option described in here:

	       F. Izquierdo-Carrasco, S.A. Smith, A. Stamatakis: "Algorithms, Data Structures, and Numerics for Likelihood-based Phylogenetic Inference of Huge Trees".

	       So tr->partitionData[model].xSpaceVector[i] provides the length of the allocated conditional array of partition model
	       and node i
	    */

	    availableLength = tr->partitionData[model].xSpaceVector[(tInfo->pNumber - tr->mxtips - 1)],
	    requiredLength = 0;

	  /* memory saving stuff, not important right now, but if you are interested ask Fernando */

	  if(tr->saveMemory)
	    {
	      size_t
		j,
		setBits = 0;

	      unsigned int
		*x1_gap = &(tr->partitionData[model].gapVector[tInfo->qNumber * tr->partitionData[model].gapVectorLength]),
		*x2_gap = &(tr->partitionData[model].gapVector[tInfo->rNumber * tr->partitionData[model].gapVectorLength]),
		*x3_gap = &(tr->partitionData[model].gapVector[tInfo->pNumber * tr->partitionData[model].gapVectorLength]);

	      for(j = 0; j < (size_t)tr->partitionData[model].gapVectorLength; j++)
		{
		  x3_gap[j] = x1_gap[j] & x2_gap[j];
		  setBits += (size_t)(precomputed16_bitcount(x3_gap[j], tr->bits_in_16bits));
		}

	      requiredLength = (width - setBits)  * rateHet * states * sizeof(double);
	    }
	  else
	    /* if we are not trying to save memory the space required to store an inner likelihood array
	       is the number of sites in the partition times the number of states of the data type in the partition
	       times the number of discrete GAMMA rates (1 for CAT essentially) times 8 bytes */
	    requiredLength  =  width * rateHet * states * sizeof(double);

	  /* Initially, even when not using memory saving no space is allocated for inner likelihood arrats hence
	     availableLength will be zero at the very first time we traverse the tree.
	     Hence we need to allocate something here */

	  if(requiredLength != availableLength)
	    {
	      /* if there is a vector of incorrect length assigned here i.e., x3 != NULL we must free
		 it first */
	      if(x3_start)
		free(x3_start);

	      /* allocate memory: note that here we use a byte-boundary aligned malloc, because we need the vectors
		 to be aligned at 16 BYTE (SSE3) or 32 BYTE (AVX) boundaries! */

	      x3_start = (double*)malloc_aligned(requiredLength);

	      /* update the data structures for consistent bookkeeping */
	      tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1] = x3_start;
	      tr->partitionData[model].xSpaceVector[(tInfo->pNumber - tr->mxtips - 1)] = requiredLength;
	    }
	} // for model
    } // for traversal
}

void assignPartitionsToThreads(tree *tr, int commRank)
{
  pInfo** rankPartitions = (pInfo **)calloc(tr->NumberOfModels, sizeof(pInfo*) );
  int i;
  for (i = 0; i < tr->NumberOfModels; ++i)
  {
    rankPartitions[i] = (pInfo *)calloc(1, sizeof(pInfo));
    rankPartitions[i]->lower = 0;
    rankPartitions[i]->upper = tr->partitionData[i].width;
    rankPartitions[i]->width = rankPartitions[i]->upper;
    rankPartitions[i]->states = tr->partitionData[i].states;
  }

  PartitionAssignment *pAss = NULL;
  initializePartitionAssignment(&pAss, rankPartitions, tr->NumberOfModels, tr->nThreads);

  /* */
  for(i = 0; i < pAss->numPartitions; ++i)
    {
      Partition
	*p = pAss->partitions + i;
      p->width = (int) ceil((float) p->width / (float) VECTOR_PADDING);
    }
  assign(pAss);

  /* Align partition sizes to the boundary (needed for site-blocking on the MIC) */
  int j;
  for(i = 0; i < pAss->numProc; ++i)
    {
      for(j = 0; j < pAss->numAssignPerProc[i] ; ++j)
	{
	  Assignment *a = &pAss->assignPerProc[i][j];
	  a->offset *= VECTOR_PADDING;
	  a->width *= VECTOR_PADDING;

	  /* adjust width of last chunk -> must NOT include padding */
	  size_t realWidth = rankPartitions[a->partId]->width;
	  if (a->offset + a->width > realWidth)
	    a->width = realWidth - a->offset;
	}
    }

  printf("Partition assignments to threads: \n");
  printAssignments(pAss);
  printf("\n");
  printLoad(pAss);
  printf("\n");

  copyThreadAssignmentInfoToTree(pAss, tr);

  deletePartitionAssignment(pAss);
  for (i = 0; i < tr->NumberOfModels; ++i)
    free(rankPartitions[i]);
  free(rankPartitions);
}
#endif


int main (int argc, char *argv[])
{ 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);
  printf("\nThis is ExaML FINE-GRAIN MPI Process Number: %d\n", processID);   
  MPI_Barrier(MPI_COMM_WORLD);
  
  {
    tree  
      *tr = (tree*)malloc(sizeof(tree));
  
    analdef 
      *adef = (analdef*)malloc(sizeof(analdef));   

    /* 
       tell the CPU to ignore exceptions generated by denormalized floating point values.
       If this is not done, depending on the input data, the likelihood functions can exhibit 
       substantial run-time differences for vectors of equal length.
    */
    
#if defined(SIMDE_SSE_NATIVE)
# include <xmmintrin.h>
    _mm_setcsr( _mm_getcsr() | _MM_FLUSH_ZERO_ON);
#endif   

  /* get the start time */
   
    masterTime = gettime();         
    
  /* initialize the analysis parameters in struct adef to default values */
    
    initAdef(adef);

  /* parse command line arguments: this has a side effect on tr struct and adef struct variables */
  
    get_args(argc, argv, adef, tr); 
  
  /* generate the ExaML output file names and store them in strings */
    
    makeFileNames();

#ifdef _USE_OMP
    if(tr->saveMemory)
      {
	printBothOpen("\nError: Memory saving option \"-S\" is not supported by the OpenMP version of ExaML!\n\n");
	error_MPI_Exit();
      }
#endif

    readByteFile(tr, processID, processes );

#ifdef _USE_OMP
    tr->nThreads = omp_get_max_threads();
    assignPartitionsToThreads(tr, processID);
#endif

    initializeTree(tr, adef);

    if(processID == 0)  
      {	
	printModelAndProgramInfo(tr, adef, argc, argv);  
	printBothOpen("Memory Saving Option: %s\n", (tr->saveMemory == TRUE)?"ENABLED":"DISABLED");   	             
      }
	
    /* do some error checks for the LG4 model and the binary models and the MIC and exit gracefully */

    {
      int 
	countBinary = 0,
	countLG4 = 0,
	model;
	
#ifdef __MIC_NATIVE
      if(tr->saveMemory)
	{
	  printBothOpen("Error: There is no MIC support yet for the memory saving option \"-S\"!\n\n");	  
	  error_MPI_Exit();  	      
	}
      
      if(tr->rateHetModel == CAT)
	{
	  printBothOpen("Error: There is no MIC support yet for the PSR model!\n\n");	  
	  error_MPI_Exit(); 
	}
#endif


      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->partitionData[model].protModels == LG4M ||  tr->partitionData[model].protModels == LG4X)
	    countLG4++;
	  if(tr->partitionData[model].states == 2)
	    countBinary++;
	}

      if(countLG4 > 0)
	{
	  if(tr->saveMemory == TRUE)
	    {
	      printBothOpen("Error: the LG4 substitution model does not work in combination with the \"-S\" memory saving flag!\n\n");	  
	      error_MPI_Exit();
	    }

	  if(tr->rateHetModel == CAT)
	    {
	      printBothOpen("Error: the LG4 substitution model does not work for proportion of invariavble sites estimates!\n\n");
	      error_MPI_Exit();
	    }
	}

      if(countBinary > 0)
	{
	  if(tr->saveMemory == TRUE)
	    {
	      printBothOpen("Error: Binary data partitions can not be used in combination with the \"-S\" memory saving flag!\n\n");	  
	      error_MPI_Exit();
	    }
	  
#ifdef __MIC_NATIVE
	  printBothOpen("Error: There is no MIC support yet for binary data partitions!\n\n");	  
	  error_MPI_Exit();  	      
#endif
	}

      if(countBinary > 0)
	{
	  if(tr->saveMemory == TRUE)
	    {
	      printBothOpen("Error: Binary data partitions can not be used in combination with the \"-S\" memory saving flag!\n\n");	  
	      error_MPI_Exit();
	    }
	  
#ifdef __MIC_NATIVE
	  printBothOpen("Error: There is no MIC support yet for binary data partitions!\n\n");	  
	  error_MPI_Exit();  	      
#endif
	}
    }
	             
    /* 
       this will re-start ExaML exactly where it has left off from a checkpoint file,
       while checkpointing is important and has to be implemented for the library we should not worry about this right now 
    */
  
   

    switch(adef->mode)
      {
      case TREE_EVALUATION:
	optimizeTrees(tr, adef);	
	break;
      case BIG_RAPID_MODE:
	if(adef->useCheckpoint)
	  {      
	    /* read checkpoint file */
	    restart(tr, adef);       	
	    
	    /* continue tree search where we left it off */
	    computeBIGRAPID(tr, adef, TRUE); 

	    /* now print the model parameters to file */
	    if(processID == 0)
	      printModelParams(tr, adef, -1);
	  }
	else
	  {
	    /* not important, only used to keep track of total accumulated exec time 
	       when checkpointing and restarts were used */
	    
	    if(processID == 0)
	      accumulatedTime = 0.0;
	    
	    /* get the starting tree: here we just parse the tree passed via the command line 
	   and do an initial likelihood computation traversal 
	   which we maybe should skip, TODO */
	    
	    getStartingTree(tr); 
	    
#ifdef _USE_OMP
	    allocateXVectors(tr);
#endif

	    /* 
	       here we do an initial full tree traversal on the starting tree using the Felsenstein pruning algorithm 
	       This should basically be the first call to the library that actually computes something :-)
	    */
	    
	    evaluateGeneric(tr, tr->start, TRUE);	       	  

	    /* the treeEvaluate() function repeatedly iterates over the entire tree to optimize branch lengths until convergence */
	    
	    treeEvaluate(tr, 1); 

	    /* now start the ML search algorithm */
	    
	    computeBIGRAPID(tr, adef, TRUE); 			     

	    /* now print the model parameters to file */
	    if(processID == 0)
	      printModelParams(tr, adef, -1);
	    
	  }         
	break;
      case QUARTET_CALCULATION:	 
	computeQuartets(tr, adef);          
	break;
      default:
	assert(0);
      }
      
    /* print some more nonsense into the ExaML_info file */
  
    if(processID == 0)
      finalizeInfoFile(tr, adef);
  }
  
  /* return 0 which means that our unix program terminated correctly, the return value is not 1 here */

  clean_MPI_Exit();

  return 0;
}



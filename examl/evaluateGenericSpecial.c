/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein. 
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

#ifndef WIN32 
#include <unistd.h>
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "axml.h"

/* the set of functions in here computes the log likelihood at a given branch (the virtual root of a tree) */

/* includes for using SSE3 intrinsics */

#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse3.h"

#ifdef __MIC_NATIVE
#include "mic_native.h"
#endif


/* 
   global variables of pthreads version, reductionBuffer is the global array 
   that is used for implementing deterministic reduction operations, that is,
   the total log likelihood over the partial log lieklihoods for the sites that each thread has computed 
   
   NumberOfThreads is just the number of threads.

   Note the volatile modifier here, that guarantees that the compiler will not do weird optimizations 
   rearraengements of the code accessing those variables, because it does not know that several concurrent threads 
   will access those variables simulatenously 
*/


extern const char inverseMeaningDNA[16];
extern int processID;

/* a pre-computed 32-bit integer mask */

extern const unsigned int mask32[32];

/* the function below computes the P matrix from the decomposition of the Q matrix and the respective rate categories for a single partition */
   

static void calcDiagptable(const double z, const int states, const int numberOfCategories, const double *rptr, const double *EIGN, double *diagptable)
{
  int 
    i, 
    l;
  
  double 
    lz,
    *lza = (double *)malloc(sizeof(double) * states);

  /* transform the root branch length to the log and check if it is not too small */

  if (z < zmin) 
    lz = log(zmin);
  else
    lz = log(z);
  
  /* do some pre-computations to avoid redundant computations further below */

  for(i = 0; i < states; i++)      
    lza[i] = EIGN[i] * lz; 

  /* loop over the number of per-site or discrete gamma rate categories */

  for(i = 0; i < numberOfCategories; i++)
    {	      	       
      /* 
	 diagptable is a pre-allocated array of doubles that stores the P-Matrix 
	 the first entry is always 1.0 
       */
      diagptable[i * states] = 1.0;

      /* compute the P matrix for all remaining states of the model */

      for(l = 1; l < states; l++)
	diagptable[i * states + l] = EXP(rptr[i] * lza[l]);
    }
  
  free(lza);
}


static void calcDiagptableFlex_LG4(double z, int numberOfCategories, double *rptr, double *EIGN[4], double *diagptable, const int numStates)
{
  int 
    i, 
    l;
  
  double 
    lz;
  
  assert(numStates <= 64);
  
  if (z < zmin) 
    lz = log(zmin);
  else
    lz = log(z);

  for(i = 0; i <  numberOfCategories; i++)
    {	      	       
      diagptable[i * numStates + 0] = 1.0;

      for(l = 1; l < numStates; l++)
	diagptable[i * numStates + l] = EXP(rptr[i] * EIGN[i][l] * lz);     	          
    }        
}





/* below are the function headers for unreadeble highly optimized versions of the above functions 
   for DNA and protein data that also use SSE3 intrinsics and implement some memory saving tricks.
   The actual functions can be found at the end of this source file. 
   All other likelihood function implementation files:

   newviewGenericSpacial.c
   makenewzSpecial.c
   evaluatePartialGenericSpecial.c

   are also structured like this 

   To decide which set of function implementations to use you will have to undefine or define _OPTIMIZED_FUNCTIONS 
   in the Makefile 
*/
   

static double evaluateGTRGAMMA_BINARY(int *ex1, int *ex2, int *wptr,
                                      double *x1_start, double *x2_start, 
                                      double *tipVector, 
                                      unsigned char *tipX1, const int n, double *diagptable, const boolean fastScaling);

static double evaluateGTRCAT_BINARY (int *ex1, int *ex2, int *cptr, int *wptr,
                                     double *x1_start, double *x2_start, double *tipVector,                   
                                     unsigned char *tipX1, int n, double *diagptable_start, const boolean fastScaling);

static double evaluateGTRGAMMAPROT_LG4(int *ex1, int *ex2, int *wptr,
				       double *x1, double *x2,  
				       double *tipVector[4], 
				       unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling, double *weights);

/* GAMMA for proteins with memory saving */

static double evaluateGTRGAMMAPROT_GAPPED_SAVE (int *wptr,
						double *x1, double *x2,  
						double *tipVector, 
						unsigned char *tipX1, int n, double *diagptable, 
						double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);


/* GAMMA for proteins */

static double evaluateGTRGAMMAPROT (int *wptr,
				    double *x1, double *x2,  
				    double *tipVector, 
				    unsigned char *tipX1, int n, double *diagptable);

/* CAT for proteins */

static double evaluateGTRCATPROT (int *cptr, int *wptr,
				  double *x1, double *x2, double *tipVector,
				  unsigned char *tipX1, int n, double *diagptable_start);


/* CAT for proteins with memory saving */

static double evaluateGTRCATPROT_SAVE (int *cptr, int *wptr,
				       double *x1, double *x2, double *tipVector,
				       unsigned char *tipX1, int n, double *diagptable_start, 
				       double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

/* analogous DNA fuctions */

static double evaluateGTRCAT_SAVE (int *cptr, int *wptr,
				   double *x1_start, double *x2_start, double *tipVector, 		      
				   unsigned char *tipX1, int n, double *diagptable_start,
				   double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

static double evaluateGTRGAMMA_GAPPED_SAVE(int *wptr,
					   double *x1_start, double *x2_start, 
					   double *tipVector, 
					   unsigned char *tipX1, const int n, double *diagptable,
					   double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

static double evaluateGTRGAMMA(int *wptr,
			       double *x1_start, double *x2_start, 
			       double *tipVector, 
			       unsigned char *tipX1, const int n, double *diagptable);


static double evaluateGTRCAT (int *cptr, int *wptr,
			      double *x1_start, double *x2_start, double *tipVector, 		      
			      unsigned char *tipX1, int n, double *diagptable_start);




/* This is the core function for computing the log likelihood at a branch */

void evaluateIterative(tree *tr)
{
  /* the branch lengths and node indices of the virtual root branch are always the first one that 
     are stored in the very important traversal array data structure that describes a partial or full tree traversal */

  /* get the branch length at the root */
  double 
    *pz = tr->td[0].ti[0].qz;   

  /* get the node number of the node to the left and right of the branch that defines the virtual rooting */

  int    
    pNumber = tr->td[0].ti[0].pNumber, 
    qNumber = tr->td[0].ti[0].qNumber;
 
  /* before we can compute the likelihood at the virtual root, we need to do a partial or full tree traversal to compute 
     the conditional likelihoods of the vectors as specified in the traversal descriptor. Maintaining this tarversal descriptor consistent 
     will unfortunately be the responsibility of users. This is tricky, if as planned for here, we use a rooted view (described somewhere in Felsenstein's book)
     for the conditional vectors with respect to the tree
  */
     
  /* iterate over all valid entries in the traversal descriptor */
  newviewIterative(tr, 1);

  int 
    m;

#ifdef _USE_OMP
#pragma omp parallel for
#endif
  for(m = 0; m < tr->NumberOfModels; m++)
    {
      /* check if this partition has to be processed now - otherwise no need to compute P matrix */
	if(!tr->td[0].executeModel[m] || tr->partitionData[m].width == 0)
	  continue;

	int
	  categories,
	  states = tr->partitionData[m].states;

	double
	  z,
	  *rateCategories,
	  *diagptable = tr->partitionData[m].left;

	/* if we are using a per-partition branch length estimate, the branch has an index, otherwise, for a joint branch length
	   estimate over all partitions we just use the branch length value with index 0 */
	if(tr->numBranches > 1)
	  z = pz[m];
	else
	  z = pz[0];


	  /*
	     figure out if we are using the CAT or GAMMA model of rate heterogeneity
	     and set pointers to the rate heterogeneity rate arrays and also set the
	     number of distinct rate categories appropriately.

	     Under GAMMA this is constant and hard-coded as 4, weheras under CAT
	     the number of site-wise rate categories can vary in the course of computations
	     up to a user defined maximum value of site categories (default: 25)
	   */
	if(tr->rateHetModel == CAT)
	  {
	    rateCategories = tr->partitionData[m].perSiteRates;
	    categories = tr->partitionData[m].numberOfCategories;
	  }
	else
	  {
	    rateCategories = tr->partitionData[m].gammaRates;
	    categories = 4;
	  }

	if(tr->partitionData[m].protModels == LG4M || tr->partitionData[m].protModels == LG4X)
	  calcDiagptableFlex_LG4(z, 4, tr->partitionData[m].gammaRates, tr->partitionData[m].EIGN_LG4, diagptable, 20);
	else
	  calcDiagptable(z, states, categories, rateCategories, tr->partitionData[m].EIGN, diagptable);
    }

  /* after the above call we are sure that we have properly and consistently computed the 
     conditionals to the right and left of the virtual root and we can now invoke the 
     the log likelihood computation */

  /* we need to loop over all partitions. Note that we may have a mix of DNA, protein binary data etc partitions */
#ifdef _USE_OMP
#pragma omp parallel
#endif
  {
    int
      m,
      model,
      maxModel;

#ifdef _USE_OMP
    maxModel = tr->maxModelsPerThread;
#else
    maxModel = tr->NumberOfModels;
#endif

  for(m = 0; m < maxModel; m++)
    {    
      /* just defaults -> if partion wasn't assigned to this thread, it will be ignored later on */
      size_t
	width = 0,
	offset = 0;

      double
	*diagptable     = (double*)NULL,
	*perPartitionLH = (double*)NULL;

      unsigned int
	*globalScaler = (unsigned int*)NULL;


#ifdef _USE_OMP
    	  int
    	    tid = omp_get_thread_num();

    	  /* check if this thread should process this partition */
    	  Assign* 
	    pAss = tr->threadPartAssigns[tid * tr->maxModelsPerThread + m];

    	  if(pAss)
	    {
	      model  = pAss->partitionId;
	      width  = pAss->width;
	      offset = pAss->offset;
	      
	      assert(model < tr->NumberOfModels);
	      
	      diagptable = tr->partitionData[model].left;
	      globalScaler = tr->partitionData[model].threadGlobalScaler[tid];
	      perPartitionLH = &tr->partitionData[model].reductionBuffer[tid];
	    }
    	  else
    	    break;
	  
#else
    	  model = m;

    	  /* number of sites in this partition */
	  width  = (size_t)tr->partitionData[model].width;
	  offset = 0;

	  /* set this pointer to the memory area where space has been reserved a priori for storing the
	     P matrix at the root */
	  diagptable = tr->partitionData[model].left;
	  globalScaler = tr->partitionData[model].globalScaler;
	  perPartitionLH = &tr->perPartitionLH[model];
#endif

        
      /* 
	 Important part of the tarversal descriptor: 
	 figure out if we need to recalculate the likelihood of this 
	 partition: 

	 The reasons why this is important in terms of performance are given in this paper 
	 here which you should actually read:
	 
	 A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009, accepted for publication, Vienna, Austria, September 2009
	 
	 The width > 0 check is for checking if under the cyclic data distribution of per-partition sites to threads this thread does indeed have a site 
	 of the current partition.

       */

      if(tr->td[0].executeModel[model] && width > 0)
	{	
	  int 
	    rateHet = (int)discreteRateCategories(tr->rateHetModel),
	    
	    /* get the number of states in the partition, e.g.: 4 = DNA, 20 = Protein */
	    states = tr->partitionData[model].states,

	    /* span for single alignment site (in doubles!) */
	    span = rateHet * states;

	  size_t
	    /* offset for current thread's data in global xVector (in doubles!) */
	    x_offset = offset * (size_t)span;

	  int
	    /* integer weight vector with pattern compression weights */
	    *wgt = tr->partitionData[model].wgt + offset,

	    /* integer rate category vector (for each pattern, _number_ of PSR category assigned to it, NOT actual rate!) */
	    *rateCategory = tr->partitionData[model].rateCategory + offset;
	  
	  double 
	    partitionLikelihood = 0.0, 	 
	    *weights = tr->partitionData[model].weights,
	    *x1_start   = (double*)NULL, 
	    *x2_start   = (double*)NULL,
	    *x1_gapColumn = (double*)NULL,
	    *x2_gapColumn = (double*)NULL;
	  	    	 	  
	  unsigned int
	    *x1_gap = (unsigned int*)NULL,
	    *x2_gap = (unsigned int*)NULL;	 
	  
	  unsigned char 
	    *tip = (unsigned char*)NULL;	  

	  /* figure out if we need to address tip vectors (a char array that indexes into a precomputed tip likelihood 
	     value array or if we need to address inner vectors */

	  /* either node p or node q is a tip */
	  
	  if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
	    {	        	    
	      /* q is a tip */

	      if(isTip(qNumber, tr->mxtips))
		{	
		  /* get the start address of the inner likelihood vector x2 for partition model,
		     note that inner nodes are enumerated/indexed starting at 0 to save allocating some 
		     space for additional pointers */
		  		 
		  x2_start = tr->partitionData[model].xVector[pNumber - tr->mxtips -1] + x_offset;

		  /* get the corresponding tip vector */

		  tip      = tr->partitionData[model].yVector[qNumber] + offset;

		  /* memory saving stuff, let's deal with this later or ask Fernando ;-) */

		  if(tr->saveMemory)
		    {
		      x2_gap         = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
		      x2_gapColumn   = &(tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet]);
		    }
		}           
	      else
		{	
		  /* p is a tip, same as above */
	 
		  x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1] + x_offset;
		  tip = tr->partitionData[model].yVector[pNumber] + offset;

		  if(tr->saveMemory)
		    {
		      x2_gap         = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
		      x2_gapColumn   = &(tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet]);
		    }

		}
	    }
	  else
	    {  
	 
	      /* neither p nor q are tips, hence we need to get the addresses of two inner vectors */
    
	      x1_start = tr->partitionData[model].xVector[pNumber - tr->mxtips - 1] + x_offset;
	      x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1] + x_offset;

	      /* memory saving option */

	      if(tr->saveMemory)
		{
		  x1_gap = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
		  x2_gap = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
		  x1_gapColumn   = &tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet];
		  x2_gapColumn   = &tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet];
		}
	
	    }
	  

	  /* for the optimized functions we have a dedicated, optimized function implementation 
	     for each rate heterogeneity and data type combination, we switch over the number of states 
	     and the rate heterogeneity model */
	  
	  switch(states)
	    { 	  
	    case 2:
#ifdef __MIC_NATIVE
 	      assert(0 && "Binary data model is not implemented on Intel MIC");
#else
	      assert(!tr->saveMemory);
	      if(tr->rateHetModel == CAT)
		partitionLikelihood = evaluateGTRCAT_BINARY((int *)NULL, (int *)NULL, rateCategory, wgt,
				      x1_start, x2_start, tr->partitionData[model].tipVector, 
				      tip, width, diagptable, TRUE);
	      else				  
		partitionLikelihood = evaluateGTRGAMMA_BINARY((int *)NULL, (int *)NULL, wgt,
							     x1_start, x2_start, 
							     tr->partitionData[model].tipVector,
							     tip, width, diagptable, TRUE);	      	      
#endif
	      break;
	    case 4: /* DNA */
	      {
		if(tr->rateHetModel == CAT)
		  {		  		  
		    if(tr->saveMemory)
#ifdef __MIC_NATIVE
		     assert(0 && "Neither CAT model of rate heterogeneity nor memory saving are implemented on Intel MIC");
#else
		      partitionLikelihood =  evaluateGTRCAT_SAVE(rateCategory, wgt,
								 x1_start, x2_start, tr->partitionData[model].tipVector, 
								 tip, width, diagptable, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
#endif
		    else
#ifdef __MIC_NATIVE
		     assert(0 && "CAT model of rate heterogeneity is not implemented on Intel MIC");
#else
		      partitionLikelihood =  evaluateGTRCAT(rateCategory, wgt,
							    x1_start, x2_start, tr->partitionData[model].tipVector, 
							    tip, width, diagptable);
#endif
		  }
		else
		  {		
		    if(tr->saveMemory)		   
#ifdef __MIC_NATIVE
 		      assert(0 && "Memory saving is not implemented on Intel MIC");
#else
		      partitionLikelihood =  evaluateGTRGAMMA_GAPPED_SAVE(wgt,
									  x1_start, x2_start, tr->partitionData[model].tipVector,
									  tip, width, diagptable,
									  x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
#endif
		    else
#ifdef __MIC_NATIVE
              partitionLikelihood = evaluateGAMMA_MIC(wgt,
	                                 x1_start, x2_start, tr->partitionData[model].mic_tipVector,
	                                 tip, width, diagptable);
#else
		      partitionLikelihood =  evaluateGTRGAMMA(wgt,
							      x1_start, x2_start, tr->partitionData[model].tipVector,
							      tip, width, diagptable);
#endif
		  }
	      }
	      break;	  	   		   
	    case 20: /* proteins */
	      {
		if(tr->rateHetModel == CAT)
		  {		   		  
		    if(tr->saveMemory)
#ifdef __MIC_NATIVE
		     assert(0 && "Neither CAT model of rate heterogeneity nor memory saving are implemented on Intel MIC");
#else
		      partitionLikelihood = evaluateGTRCATPROT_SAVE(rateCategory, wgt,
								    x1_start, x2_start, tr->partitionData[model].tipVector,
								    tip, width, diagptable,  x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
#endif
		    else
#ifdef __MIC_NATIVE
		     assert(0 && "CAT model of rate heterogeneity is not implemented on Intel MIC");
#else
		      partitionLikelihood = evaluateGTRCATPROT(rateCategory, wgt,
							       x1_start, x2_start, tr->partitionData[model].tipVector,
							       tip, width, diagptable);
#endif
		  }
		else
		  {		    		    		      
		    if(tr->saveMemory)
#ifdef __MIC_NATIVE
 		      assert(0 && "Memory saving is not implemented on Intel MIC");
#else
		      partitionLikelihood = evaluateGTRGAMMAPROT_GAPPED_SAVE(wgt,
									     x1_start, x2_start, tr->partitionData[model].tipVector,
									     tip, width, diagptable,
									     x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
#endif
		    else
		      {
			if(tr->partitionData[model].protModels == LG4M || tr->partitionData[model].protModels == LG4X)
#ifdef __MIC_NATIVE
			 partitionLikelihood = evaluateGAMMAPROT_LG4_MIC(wgt,
                               x1_start, x2_start, tr->partitionData[model].mic_tipVector,
                               tip, width, diagptable, weights);
#else
			  partitionLikelihood =  evaluateGTRGAMMAPROT_LG4((int *)NULL, (int *)NULL, wgt,
									  x1_start, x2_start, tr->partitionData[model].tipVector_LG4,
									  tip, width, diagptable, TRUE, weights);
#endif
			else
#ifdef __MIC_NATIVE
	            partitionLikelihood = evaluateGAMMAPROT_MIC(wgt,
	                               x1_start, x2_start, tr->partitionData[model].mic_tipVector,
	                               tip, width, diagptable);
#else
			  partitionLikelihood = evaluateGTRGAMMAPROT(wgt,
								     x1_start, x2_start, tr->partitionData[model].tipVector,
								     tip, width, diagptable);
#endif
		      }
		  }
	      }
	      break;	      		    
	    default:
	      assert(0);	    
	    }	
	  
	  /* now here is a nasty part, for each partition and each node we maintain an integer counter to count how often 
	     how many entries per node were scaled by a constant factor. Here we use this information generated during Felsenstein's 
	     pruning algorithm by the newview() functions to undo the preceding scaling multiplications at the root, for mathematical details 
	     you should actually read:

	     A. Stamatakis: "Orchestrating the Phylogenetic Likelihood Function on Emerging Parallel Architectures". 
	     In B. Schmidt, editor, Bioinformatics: High Performance Parallel Computer Architectures, 85-115, CRC Press, Taylor & Francis, 2010.

	     There's a copy of this book in my office 
	  */

	  partitionLikelihood += (globalScaler[pNumber] + globalScaler[qNumber]) * LOG(minlikelihood);

	  /* check that there was no major numerical screw-up, the log likelihood should be < 0.0 always */

	  

	  assert(partitionLikelihood < 0.0);

	  /* now we have the correct log likelihood for the current partition after undoing scaling multiplications */	  	 
	  
	  /* finally, we also store the per partition log likelihood which is important for optimizing the alpha parameter 
	     of this partition for example */

	  *perPartitionLH = partitionLikelihood;
	}
      else
	{
	  /* if the current thread does not have a single site of this partition
	     it is important to set the per partition log like to 0.0 because 
	     of the reduction operation that will take place later-on.
	     That is, the values of tr->perPartitionLH across all threads 
	     need to be in a consistent state, always !
	  */

	  if(width == 0)	    
	    *perPartitionLH = 0.0;
	  else
	    {
	      assert(tr->td[0].executeModel[model] == FALSE && *perPartitionLH < 0.0);
	    }
	}
    }  /* for model */
  }  /* OMP parallel */


#ifdef _USE_OMP
  /* perform reduction of per-partition LH scores */
  int
    model,
    t;

  for(model = 0; model < tr->NumberOfModels; model++)
    {
     if (!tr->td[0].executeModel[model])
       continue;

      tr->perPartitionLH[model] = 0.0;
      for(t = 0; t < tr->maxThreadsPerModel; t++)
	{
	  Assign*
	    pAss = tr->partThreadAssigns[model * tr->maxThreadsPerModel + t];

	  if (pAss)
	    {
	      int
		tid = pAss->procId;

	      tr->perPartitionLH[model] += tr->partitionData[model].reductionBuffer[tid];
	    }
	}
    }
#endif
}




void evaluateGeneric (tree *tr, nodeptr p, boolean fullTraversal)
{
  /* now this may be the entry point of the library to compute 
     the log like at a branch defined by p and p->back == q */

  volatile double 
    result = 0.0;
  
  nodeptr 
    q = p->back; 
  
  int 
    i,
    model;

 
  /* set the first entry of the traversal descriptor to contain the indices
     of nodes p and q */

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;          
  
  /* copy the branch lengths of the tree into the first entry of the traversal descriptor.
     if -M is not used tr->numBranches must be 1 */

  for(i = 0; i < tr->numBranches; i++)    
    tr->td[0].ti[0].qz[i] =  q->z[i];
  
  /* now compute how many conditionals must be re-computed/re-oriented by newview
     to be able to calculate the likelihood at the root defined by p and q.
  */

  /* one entry in the traversal descriptor is already used, hence set the tarversal length counter to 1 */
  tr->td[0].count = 1;

  /* do we need to recompute any of the vectors at or below p ? */
  
  if(fullTraversal)
    { 
      assert(isTip(p->number, tr->mxtips));
      computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches, FALSE);     
    }
  else
    {
      if(!p->x)
	computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches, TRUE);

      /* recompute/reorient any descriptors at or below q ? 
	 computeTraversalInfo computes and stores the newview() to be executed for the traversal descriptor */
      
      if(!q->x)
	computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches, TRUE);  
    }
   
      /* now we copy this partition execute mask into the traversal descriptor which must come from the 
	 calling program, the logic of this should not form part of the library */

  storeExecuteMaskInTraversalDescriptor(tr);  
  
  /* also store in the traversal descriptor that something has changed i.e., in the parallel case that the 
     traversal descriptor list of nodes needs to be broadcast once again */
  
  tr->td[0].traversalHasChanged = TRUE;


  evaluateIterative(tr);  
  
  {
    double 
      *recv = (double *)malloc(sizeof(double) * tr->NumberOfModels);
    
#ifdef _USE_ALLREDUCE   
    MPI_Allreduce(tr->perPartitionLH, recv, tr->NumberOfModels, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    MPI_Reduce(tr->perPartitionLH, recv, tr->NumberOfModels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(recv, tr->NumberOfModels, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    memcpy(tr->perPartitionLH, recv, tr->NumberOfModels * sizeof(double));

    for(model = 0; model < tr->NumberOfModels; model++)        
      result += tr->perPartitionLH[model];
         
    free(recv);
  }


  /* set the tree data structure likelihood value to the total likelihood */

  tr->likelihood = result;    
  
  /* 
     MPI_Barrier(MPI_COMM_WORLD);
     printf("Process %d likelihood: %f\n", processID, tr->likelihood);
     MPI_Barrier(MPI_COMM_WORLD);
  */

  /* do some bookkeeping to have traversalHasChanged in a consistent state */

  tr->td[0].traversalHasChanged = FALSE;  

  
  
 
}







/* below are the optimized function versions with geeky intrinsics */

/* binary data */

static double evaluateGTRCAT_BINARY (int *ex1, int *ex2, int *cptr, int *wptr,
                                     double *x1_start, double *x2_start, double *tipVector,                   
                                     unsigned char *tipX1, int n, double *diagptable_start, const boolean fastScaling)
{
  double  sum = 0.0, term;       
  int     i;
  double  *diagptable, *x1, *x2;                            
 
  if(tipX1)
    {          
      for (i = 0; i < n; i++) 
        {
          double 
	    t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));

          x1 = &(tipVector[2 * tipX1[i]]);
          x2 = &(x2_start[2 * i]);
          
          diagptable = &(diagptable_start[2 * cptr[i]]);                          
        

          _mm_store_pd(t, _mm_mul_pd(_mm_load_pd(x1), _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(diagptable))));
          
          if(fastScaling)
            term = log(fabs(t[0] + t[1]));
          else
            term = log(fabs(t[0] + t[1])) + (ex2[i] * log(minlikelihood));                           

          sum += wptr[i] * term;
        }       
    }               
  else
    {
      for (i = 0; i < n; i++) 
        {       
          double 
	    t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));                                 
            
          x1 = &x1_start[2 * i];
          x2 = &x2_start[2 * i];
          
          diagptable = &diagptable_start[2 * cptr[i]];            

          _mm_store_pd(t, _mm_mul_pd(_mm_load_pd(x1), _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(diagptable))));
          
          if(fastScaling)
            term = log(fabs(t[0] + t[1]));
          else
            term = log(fabs(t[0] + t[1])) + ((ex1[i] + ex2[i]) * log(minlikelihood));                        

          
          sum += wptr[i] * term;
        }          
    }
       
  return  sum;         
} 


static double evaluateGTRGAMMA_BINARY(int *ex1, int *ex2, int *wptr,
                                      double *x1_start, double *x2_start, 
                                      double *tipVector, 
                                      unsigned char *tipX1, const int n, double *diagptable, const boolean fastScaling)
{
  double   sum = 0.0, term;    
  int     i, j;  
  double  *x1, *x2;             

  if(tipX1)
    {          
      for (i = 0; i < n; i++)
        {
          double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
          __m128d termv, x1v, x2v, dv;
	  
          x1 = &(tipVector[2 * tipX1[i]]);       
          x2 = &x2_start[8 * i];                                

          termv = _mm_set1_pd(0.0);                
          
          for(j = 0; j < 4; j++)
            {
              x1v = _mm_load_pd(&x1[0]);
              x2v = _mm_load_pd(&x2[j * 2]);
              dv   = _mm_load_pd(&diagptable[j * 2]);
              
              x1v = _mm_mul_pd(x1v, x2v);
              x1v = _mm_mul_pd(x1v, dv);
              
              termv = _mm_add_pd(termv, x1v);                 
            }
          
          _mm_store_pd(t, termv);               
          
          if(fastScaling)
            term = log(0.25 * (fabs(t[0] + t[1])));
          else
            term = log(0.25 * (fabs(t[0] + t[1]))) + (ex2[i] * log(minlikelihood));       
 
          
          sum += wptr[i] * term;
        }         
    }
  else
    {         
      for (i = 0; i < n; i++) 
        {

          double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
          __m128d termv, x1v, x2v, dv;
                        
          x1 = &x1_start[8 * i];
          x2 = &x2_start[8 * i];
                  

          termv = _mm_set1_pd(0.0);                
          
          for(j = 0; j < 4; j++)
            {
              x1v = _mm_load_pd(&x1[j * 2]);
              x2v = _mm_load_pd(&x2[j * 2]);
              dv   = _mm_load_pd(&diagptable[j * 2]);
              
              x1v = _mm_mul_pd(x1v, x2v);
              x1v = _mm_mul_pd(x1v, dv);
              
              termv = _mm_add_pd(termv, x1v);                 
            }
          
          _mm_store_pd(t, termv);
          
          
          if(fastScaling)
            term = log(0.25 * (fabs(t[0] + t[1])));
          else
            term = log(0.25 * (fabs(t[0] + t[1]))) + ((ex1[i] +ex2[i]) * log(minlikelihood));     


          sum += wptr[i] * term;
        }                       
    }

  return sum;
} 


/* binary data end */


static double evaluateGTRGAMMAPROT_LG4(int *ex1, int *ex2, int *wptr,
				       double *x1, double *x2,  
				       double *tipVector[4], 
				       unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling, double *weights)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
	{
	  __m128d 
	    tv = _mm_setzero_pd();
	      	 	  	 
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double 
		*d = &diagptable[j * 20];

	      __m128d 
		t = _mm_setzero_pd(),
		w = _mm_set1_pd(weights[j]);
	      
	      
	      left = &(tipVector[j][20 * tipX1[i]]);
	      right = &(x2[80 * i + 20 * j]);
	      
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  t = _mm_add_pd(t, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}
	      
	      tv = _mm_add_pd(tv, _mm_mul_pd(t, w));	      	      	     
	    }
	  
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	  
	  sum += wptr[i] * term;
	}    	        
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{	  	 	             
	  __m128d 
	    tv = _mm_setzero_pd();	 	  	  
	      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double 
		*d = &diagptable[j * 20];

	      __m128d 
		t = _mm_setzero_pd(),
		w = _mm_set1_pd(weights[j]);
	      
	      left  = &(x1[80 * i + 20 * j]);
	      right = &(x2[80 * i + 20 * j]);
	      
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  t = _mm_add_pd(t, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 

	       tv = _mm_add_pd(tv, _mm_mul_pd(t, w));
	    }
	  
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);	  
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}



static double evaluateGTRGAMMAPROT_GAPPED_SAVE (int *wptr,
						double *x1, double *x2,  
						double *tipVector, 
						unsigned char *tipX1, int n, double *diagptable, 
						double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)					   
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  
    *left, 
    *right,
    *x1_ptr = x1,
    *x2_ptr = x2,
    *x1v,
    *x2v;              
  
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
	{
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2v = x2_gapColumn;
	  else
	    {
	      x2v = x2_ptr;
	      x2_ptr += 80;
	    }

	  __m128d tv = _mm_setzero_pd();
	  left = &(tipVector[20 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double *d = &diagptable[j * 20];
	      right = &(x2v[20 * j]);
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 		
	    }

	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
	  

	  
	  term = LOG(0.25 * FABS(term));	  
	  
	  sum += wptr[i] * term;
	}    	        
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{
	  if(x1_gap[i / 32] & mask32[i % 32])
	    x1v = x1_gapColumn;
	  else
	    {
	      x1v = x1_ptr;
	      x1_ptr += 80;
	    }

	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2v = x2_gapColumn;
	  else
	    {
	      x2v = x2_ptr;
	      x2_ptr += 80;
	    }
	  	 	             
	  __m128d tv = _mm_setzero_pd();	 	  	  
	      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double *d = &diagptable[j * 20];
	      left  = &(x1v[20 * j]);
	      right = &(x2v[20 * j]);
	      
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 		
	    }
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);	  
	  
	 
	  term = LOG(0.25 * FABS(term));
	
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}



static double evaluateGTRGAMMAPROT (int *wptr,
				    double *x1, double *x2,  
				    double *tipVector, 
				    unsigned char *tipX1, int n, double *diagptable)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
	{

	  __m128d tv = _mm_setzero_pd();
	  left = &(tipVector[20 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double *d = &diagptable[j * 20];
	      right = &(x2[80 * i + 20 * j]);
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 		
	    }
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
	  
	  
	 
	  term = LOG(0.25 * FABS(term));
		 
	  
	  sum += wptr[i] * term;
	}    	        
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{	  	 	             
	  __m128d tv = _mm_setzero_pd();	 	  	  
	      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double *d = &diagptable[j * 20];
	      left  = &(x1[80 * i + 20 * j]);
	      right = &(x2[80 * i + 20 * j]);
	      
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 		
	    }
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);	  
	  
	
	  term = LOG(0.25 * FABS(term));
	  
	  
	  sum += wptr[i] * term;
	}
    }
       
  return  sum;
}


static double evaluateGTRCATPROT (int *cptr, int *wptr,
				  double *x1, double *x2, double *tipVector,
				  unsigned char *tipX1, int n, double *diagptable_start)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           
  
  if(tipX1)
    {                 
      for (i = 0; i < n; i++) 
	{	       	
	  left = &(tipVector[20 * tipX1[i]]);
	  right = &(x2[20 * i]);
	  
	  diagptable = &diagptable_start[20 * cptr[i]];	           	 

	  __m128d tv = _mm_setzero_pd();	    
	  
	  for(l = 0; l < 20; l+=2)
	    {
	      __m128d lv = _mm_load_pd(&left[l]);
	      __m128d rv = _mm_load_pd(&right[l]);
	      __m128d mul = _mm_mul_pd(lv, rv);
	      __m128d dv = _mm_load_pd(&diagptable[l]);
	      
	      tv = _mm_add_pd(tv, _mm_mul_pd(mul, dv));		   
	    }		 		
	  
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
  
	  
	  term = LOG(FABS(term));
	  	  
	  sum += wptr[i] * term;
	}      
    }    
  else
    {
    
      for (i = 0; i < n; i++) 
	{		       	      	      
	  left  = &x1[20 * i];
	  right = &x2[20 * i];
	  
	  diagptable = &diagptable_start[20 * cptr[i]];	  	

	  __m128d tv = _mm_setzero_pd();	    
	      	    
	  for(l = 0; l < 20; l+=2)
	    {
	      __m128d lv = _mm_load_pd(&left[l]);
	      __m128d rv = _mm_load_pd(&right[l]);
	      __m128d mul = _mm_mul_pd(lv, rv);
	      __m128d dv = _mm_load_pd(&diagptable[l]);
	      
	      tv = _mm_add_pd(tv, _mm_mul_pd(mul, dv));		   
	    }		 		
	  
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
	  	  
	  term = LOG(FABS(term));	 
	  
	  sum += wptr[i] * term;      
	}
    }
             
  return  sum;         
} 


static double evaluateGTRCATPROT_SAVE (int *cptr, int *wptr,
				       double *x1, double *x2, double *tipVector,
				       unsigned char *tipX1, int n, double *diagptable_start, 
				       double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   
    sum = 0.0, 
    term,
    *diagptable,  
    *left, 
    *right,
    *left_ptr = x1,
    *right_ptr = x2;
  
  int     
    i, 
    l;                           
  
  if(tipX1)
    {                 
      for (i = 0; i < n; i++) 
	{	       	
	  left = &(tipVector[20 * tipX1[i]]);

	  if(isGap(x2_gap, i))
	    right = x2_gapColumn;
	  else
	    {
	      right = right_ptr;
	      right_ptr += 20;
	    }	  	 
	  
	  diagptable = &diagptable_start[20 * cptr[i]];	           	 

	  __m128d tv = _mm_setzero_pd();	    
	  
	  for(l = 0; l < 20; l+=2)
	    {
	      __m128d lv = _mm_load_pd(&left[l]);
	      __m128d rv = _mm_load_pd(&right[l]);
	      __m128d mul = _mm_mul_pd(lv, rv);
	      __m128d dv = _mm_load_pd(&diagptable[l]);
	      
	      tv = _mm_add_pd(tv, _mm_mul_pd(mul, dv));		   
	    }		 		
	  
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
    
	  
	  term = LOG(FABS(term));
	  	  
	  sum += wptr[i] * term;
	}      
    }    
  else
    {
    
      for (i = 0; i < n; i++) 
	{		       	      	      	  
	  if(isGap(x1_gap, i))
	    left = x1_gapColumn;
	  else
	    {
	      left = left_ptr;
	      left_ptr += 20;
	    }
	  
	  if(isGap(x2_gap, i))
	    right = x2_gapColumn;
	  else
	    {
	      right = right_ptr;
	      right_ptr += 20;
	    }
	  
	  diagptable = &diagptable_start[20 * cptr[i]];	  	

	  __m128d tv = _mm_setzero_pd();	    
	  
	  for(l = 0; l < 20; l+=2)
	    {
	      __m128d lv = _mm_load_pd(&left[l]);
	      __m128d rv = _mm_load_pd(&right[l]);
	      __m128d mul = _mm_mul_pd(lv, rv);
	      __m128d dv = _mm_load_pd(&diagptable[l]);
	      
	      tv = _mm_add_pd(tv, _mm_mul_pd(mul, dv));		   
	    }		 		
	  
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
	  	  
	  term = LOG(FABS(term));	 
	  
	  sum += wptr[i] * term;      
	}
    }
             
  return  sum;         
} 


static double evaluateGTRCAT_SAVE (int *cptr, int *wptr,
				   double *x1_start, double *x2_start, double *tipVector, 		      
				   unsigned char *tipX1, int n, double *diagptable_start,
				   double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double  sum = 0.0, term;       
  int     i;

  double  *diagptable, 
    *x1, 
    *x2,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start;
 
  if(tipX1)
    {           
      for (i = 0; i < n; i++) 
	{	
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	  __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;

	  x1 = &(tipVector[4 * tipX1[i]]);

	  if(isGap(x2_gap, i))
	    x2 = x2_gapColumn;
	  else
	    {
	      x2 = x2_ptr;
	      x2_ptr += 4;
	    }
	  
	  diagptable = &diagptable_start[4 * cptr[i]];
	  	    	  
	  x1v1 =  _mm_load_pd(&x1[0]);
	  x1v2 =  _mm_load_pd(&x1[2]);
	  x2v1 =  _mm_load_pd(&x2[0]);
	  x2v2 =  _mm_load_pd(&x2[2]);
	  dv1  =  _mm_load_pd(&diagptable[0]);
	  dv2  =  _mm_load_pd(&diagptable[2]);
	  
	  x1v1 = _mm_mul_pd(x1v1, x2v1);
	  x1v1 = _mm_mul_pd(x1v1, dv1);
	  
	  x1v2 = _mm_mul_pd(x1v2, x2v2);
	  x1v2 = _mm_mul_pd(x1v2, dv2);
	  
	  x1v1 = _mm_add_pd(x1v1, x1v2);
	  
	  _mm_store_pd(t, x1v1);
	  	  
	  term = LOG(FABS(t[0] + t[1]));
	      
	 

	  sum += wptr[i] * term;
	}	
    }               
  else
    {
      for (i = 0; i < n; i++) 
	{ 
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	  __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;
	   
	  if(isGap(x1_gap, i))
	    x1 = x1_gapColumn;
	  else
	    {
	      x1 = x1_ptr;
	      x1_ptr += 4;
	    }
	  
	  if(isGap(x2_gap, i))
	    x2 = x2_gapColumn;
	  else
	    {
	      x2 = x2_ptr;
	      x2_ptr += 4;
	    }
	  
	  diagptable = &diagptable_start[4 * cptr[i]];	
	  
	  x1v1 =  _mm_load_pd(&x1[0]);
	  x1v2 =  _mm_load_pd(&x1[2]);
	  x2v1 =  _mm_load_pd(&x2[0]);
	  x2v2 =  _mm_load_pd(&x2[2]);
	  dv1  =  _mm_load_pd(&diagptable[0]);
	  dv2  =  _mm_load_pd(&diagptable[2]);
	  
	  x1v1 = _mm_mul_pd(x1v1, x2v1);
	  x1v1 = _mm_mul_pd(x1v1, dv1);
	  
	  x1v2 = _mm_mul_pd(x1v2, x2v2);
	  x1v2 = _mm_mul_pd(x1v2, dv2);
	  
	  x1v1 = _mm_add_pd(x1v1, x1v2);
	  
	  _mm_store_pd(t, x1v1);
	  
	 
	  term = LOG(FABS(t[0] + t[1]));
	  
	  sum += wptr[i] * term;
	}    
    }
       
  return  sum;         
} 


static double evaluateGTRGAMMA_GAPPED_SAVE(int *wptr,
					   double *x1_start, double *x2_start, 
					   double *tipVector, 
					   unsigned char *tipX1, const int n, double *diagptable,
					   double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   sum = 0.0, term;    
  int     i, j;
  double  
    *x1, 
    *x2,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start;

 

  if(tipX1)
    {        
     
      
      for (i = 0; i < n; i++)
	{
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	  __m128d termv, x1v, x2v, dv;

	  x1 = &(tipVector[4 * tipX1[i]]);	 
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    {
	      x2 = x2_ptr;	 
	      x2_ptr += 16;
	    }
	  
	
	  termv = _mm_set1_pd(0.0);	    	   
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[0]);
	      x2v = _mm_load_pd(&x2[j * 4]);
	      dv   = _mm_load_pd(&diagptable[j * 4]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	      
	      x1v = _mm_load_pd(&x1[2]);
	      x2v = _mm_load_pd(&x2[j * 4 + 2]);
	      dv   = _mm_load_pd(&diagptable[j * 4 + 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	    }
	  
	  _mm_store_pd(t, termv);	  	 

	 
	  term = LOG(0.25 * FABS(t[0] + t[1]));
	   
	  
	  sum += wptr[i] * term;
	}     
    }
  else
    {        
      
      for (i = 0; i < n; i++) 
	{

	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	  __m128d termv, x1v, x2v, dv;

	  if(x1_gap[i / 32] & mask32[i % 32])
	    x1 = x1_gapColumn;
	  else
	    {
	      x1 = x1_ptr; 	  	  
	      x1_ptr += 16;
	    }
	 	      
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    {
	      x2 = x2_ptr;
	      x2_ptr += 16;
	    }
	
	  termv = _mm_set1_pd(0.0);	  	 
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[j * 4]);
	      x2v = _mm_load_pd(&x2[j * 4]);
	      dv   = _mm_load_pd(&diagptable[j * 4]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	      
	      x1v = _mm_load_pd(&x1[j * 4 + 2]);
	      x2v = _mm_load_pd(&x2[j * 4 + 2]);
	      dv   = _mm_load_pd(&diagptable[j * 4 + 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	    }
	  
	  _mm_store_pd(t, termv);

	 
	  term = LOG(0.25 * FABS(t[0] + t[1]));
	 	  
	  
	  sum += wptr[i] * term;
	}                      	
    }

  return sum;
} 


static double evaluateGTRGAMMA(int *wptr,
			       double *x1_start, double *x2_start, 
			       double *tipVector, 
			       unsigned char *tipX1, const int n, double *diagptable)
{
  double   sum = 0.0, term;    
  int     i, j;

  double  *x1, *x2;             

 

  if(tipX1)
    {          	
      for (i = 0; i < n; i++)
	{

	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	  __m128d termv, x1v, x2v, dv;

	  x1 = &(tipVector[4 * tipX1[i]]);	 
	  x2 = &x2_start[16 * i];	 
	  
	
	  termv = _mm_set1_pd(0.0);	    	   
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[0]);
	      x2v = _mm_load_pd(&x2[j * 4]);
	      dv   = _mm_load_pd(&diagptable[j * 4]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	      
	      x1v = _mm_load_pd(&x1[2]);
	      x2v = _mm_load_pd(&x2[j * 4 + 2]);
	      dv   = _mm_load_pd(&diagptable[j * 4 + 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	    }
	  
	  _mm_store_pd(t, termv);
	  
	  
	
	  term = LOG(0.25 * FABS(t[0] + t[1]));
	  
	 
	  
	  sum += wptr[i] * term;
	}     
    }
  else
    {        
      for (i = 0; i < n; i++) 
	{

	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	  __m128d termv, x1v, x2v, dv;

	  	 	  	  
	  x1 = &x1_start[16 * i];
	  x2 = &x2_start[16 * i];	  	  
	
	
	  termv = _mm_set1_pd(0.0);	  	 
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[j * 4]);
	      x2v = _mm_load_pd(&x2[j * 4]);
	      dv   = _mm_load_pd(&diagptable[j * 4]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	      
	      x1v = _mm_load_pd(&x1[j * 4 + 2]);
	      x2v = _mm_load_pd(&x2[j * 4 + 2]);
	      dv   = _mm_load_pd(&diagptable[j * 4 + 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	    }
	  
	  _mm_store_pd(t, termv);

	  
	    term = LOG(0.25 * FABS(t[0] + t[1]));
	 	  

	  
	  sum += wptr[i] * term;
	}                      	
    }

  return sum;
} 


static double evaluateGTRCAT (int *cptr, int *wptr,
			      double *x1_start, double *x2_start, double *tipVector, 		      
			      unsigned char *tipX1, int n, double *diagptable_start)
{
  double  sum = 0.0, term;       
  int     i;

  double  *diagptable, *x1, *x2;                      	    
 
  if(tipX1)
    {           
      for (i = 0; i < n; i++) 
	{	

	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	  __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;

	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];
	  
	  diagptable = &diagptable_start[4 * cptr[i]];
	  
	    	  
	  x1v1 =  _mm_load_pd(&x1[0]);
	  x1v2 =  _mm_load_pd(&x1[2]);
	  x2v1 =  _mm_load_pd(&x2[0]);
	  x2v2 =  _mm_load_pd(&x2[2]);
	  dv1  =  _mm_load_pd(&diagptable[0]);
	  dv2  =  _mm_load_pd(&diagptable[2]);
	  
	  x1v1 = _mm_mul_pd(x1v1, x2v1);
	  x1v1 = _mm_mul_pd(x1v1, dv1);
	  
	  x1v2 = _mm_mul_pd(x1v2, x2v2);
	  x1v2 = _mm_mul_pd(x1v2, dv2);
	  
	  x1v1 = _mm_add_pd(x1v1, x1v2);
	  
	  _mm_store_pd(t, x1v1);
	  
	  
	  term = LOG(FABS(t[0] + t[1]));
	  

	  sum += wptr[i] * term;
	}	
    }               
  else
    {
      for (i = 0; i < n; i++) 
	{ 

	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	   __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;

	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  
	  diagptable = &diagptable_start[4 * cptr[i]];	
	  
  
	  x1v1 =  _mm_load_pd(&x1[0]);
	  x1v2 =  _mm_load_pd(&x1[2]);
	  x2v1 =  _mm_load_pd(&x2[0]);
	  x2v2 =  _mm_load_pd(&x2[2]);
	  dv1  =  _mm_load_pd(&diagptable[0]);
	  dv2  =  _mm_load_pd(&diagptable[2]);
	  
	  x1v1 = _mm_mul_pd(x1v1, x2v1);
	  x1v1 = _mm_mul_pd(x1v1, dv1);
	  
	  x1v2 = _mm_mul_pd(x1v2, x2v2);
	  x1v2 = _mm_mul_pd(x1v2, dv2);
	  
	  x1v1 = _mm_add_pd(x1v1, x1v2);
	  
	  _mm_store_pd(t, x1v1);
	  
	 
	  term = LOG(FABS(t[0] + t[1]));
	  

	  sum += wptr[i] * term;
	}    
    }
       
  return  sum;         
} 

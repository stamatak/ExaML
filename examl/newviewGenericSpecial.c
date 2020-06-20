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
#include <stdint.h>
#include <limits.h>
#include "axml.h"


#include <stdint.h>
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/avx2.h"

/* required to compute the absoliute values of double precision numbers with SSE3 */

const union __attribute__ ((aligned (BYTE_ALIGNMENT)))
{
       uint64_t i[2];
       __m128d m;
} absMask = {{0x7fffffffffffffffULL , 0x7fffffffffffffffULL }};




/* includes MIC-optimized functions */

#ifdef __MIC_NATIVE
#include "mic_native.h"
#endif

extern int processID;

/* bit mask */

extern const unsigned int mask32[32];


/* generic function for computing the P matrices, for computing the conditional likelihood at a node p, given child nodes q and r 
   we compute P(z1) and P(z2) here */

static void makeP(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, double *right, boolean saveMem, int maxCat, const int states)
{
  int 
    i, 
    j, 
    k,
    /* square of the number of states = P-matrix size */
    statesSquare = states * states;
  
  /* assign some space for pre-computing and later re-using functions */

  double 
    *lz1 = (double*)malloc(sizeof(double) * states),
    *lz2 = (double*)malloc(sizeof(double) * states),
    *d1 = (double*)malloc(sizeof(double) * states),
    *d2 = (double*)malloc(sizeof(double) * states);

  /* multiply branch lengths with eigenvalues */

  for(i = 1; i < states; i++)
    {
      lz1[i] = EIGN[i] * z1;
      lz2[i] = EIGN[i] * z2;
    }


  /* loop over the number of rate categories, this will be 4 for the GAMMA model and 
     variable for the CAT model */

  for(i = 0; i < numberOfCategories; i++)
    {
      /* exponentiate the rate multiplied by the branch */

      for(j = 1; j < states; j++)
	{
	  d1[j] = EXP(rptr[i] * lz1[j]);
	  d2[j] = EXP(rptr[i] * lz2[j]);
	}

      /* now fill the P matrices for the two branch length values */

      for(j = 0; j < states; j++)
	{
	  /* left and right are pre-allocated arrays */

	  left[statesSquare * i  + states * j] = 1.0;
	  right[statesSquare * i + states * j] = 1.0;	  

	  for(k = 1; k < states; k++)
	    {
	      left[statesSquare * i + states * j + k]  = d1[k] * EI[states * j + k];
	      right[statesSquare * i + states * j + k] = d2[k] * EI[states * j + k];
	    }
	}
    }


  /* if memory saving is enabled and we are using CAT we need to do one additional P matrix 
     calculation for a rate of 1.0 to compute the entries of a column/tree site comprising only gaps */


  if(saveMem)
    {
      i = maxCat;
      
      for(j = 1; j < states; j++)
	{
	  d1[j] = EXP (lz1[j]);
	  d2[j] = EXP (lz2[j]);
	}

      for(j = 0; j < states; j++)
	{
	  left[statesSquare * i  + states * j] = 1.0;
	  right[statesSquare * i + states * j] = 1.0;

	  for(k = 1; k < states; k++)
	    {
	      left[statesSquare * i + states * j + k]  = d1[k] * EI[states * j + k];
	      right[statesSquare * i + states * j + k] = d2[k] * EI[states * j + k];
	    }
	}
    }
  
  /* free the temporary buffers */

  free(lz1);
  free(lz2);
  free(d1);
  free(d2);
}

static void makeP_FlexLG4(double z1, double z2, double *rptr, double *EI[4],  double *EIGN[4], int numberOfCategories, double *left, double *right, const int numStates)
{
  int 
    i,
    j,
    k;
  
  const int
    statesSquare = numStates * numStates;

  double    
    d1[64],  
    d2[64];

  assert(numStates <= 64);
       
  for(i = 0; i < numberOfCategories; i++)
    {
      for(j = 1; j < numStates; j++)
	{
	  d1[j] = EXP (rptr[i] * EIGN[i][j] * z1);
	  d2[j] = EXP (rptr[i] * EIGN[i][j] * z2);
	}

      for(j = 0; j < numStates; j++)
	{
	  left[statesSquare * i  + numStates * j] = 1.0;
	  right[statesSquare * i + numStates * j] = 1.0;

	  for(k = 1; k < numStates; k++)
	    {
	      left[statesSquare * i + numStates * j + k]  = d1[k] * EI[i][numStates * j + k];
	      right[statesSquare * i + numStates * j + k] = d2[k] * EI[i][numStates * j + k];
	    }
	}
    }  
}



/* The functions here are organized in a similar way as in evaluateGenericSpecial.c 
   I provide generic, slow but readable function implementations for computing the 
   conditional likelihood arrays at p, given child nodes q and r. Once again we need 
   two generic function implementations, one for CAT and one for GAMMA */



    
/* The function below computes partial traversals only down to the point/node in the tree where the 
   conditional likelihhod vector summarizing a subtree is already oriented in the correct direction */

void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches, boolean partialTraversal)
{
  /* if it's a tip we don't do anything */

  if(isTip(p->number, maxTips))
    return;

  {
    int 
      i;
    
    /* get the left and right descendants */

    nodeptr 
      q = p->next->back,
      r = p->next->next->back;   

    /* if the left and right children are tips there is not that much to do */

    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
      {
	/* fix the orientation of p->x */
	
	if (! p->x)
	  getxnode(p);	
	assert(p->x);
	  
	/* add the current node triplet p,q,r to the traversal descriptor */

	ti[*counter].tipCase = TIP_TIP;
	ti[*counter].pNumber = p->number;
	ti[*counter].qNumber = q->number;
	ti[*counter].rNumber = r->number;

	/* copy branches to traversal descriptor */

	for(i = 0; i < numBranches; i++)
	  {	    
	    ti[*counter].qz[i] = q->z[i];
	    ti[*counter].rz[i] = r->z[i];
	  }

	/* increment length counter */

	*counter = *counter + 1;
      }
    else
      {
	/* if either r or q are tips, flip them to make sure that the tip data is stored 
	   for q */

	if(isTip(r->number, maxTips) || isTip(q->number, maxTips))
	  {	    
	    if(isTip(r->number, maxTips))
	      {
		nodeptr 
		  tmp = r;
		r = q;
		q = tmp;
	      }
	   
	    /* if the orientation of the liklihood vector at r is not correct we need to re-compute it 
	       and descend into its subtree to figure out if there are more vrctors in there to re-compute and 
	       re-orient */

	    if(!r->x || !partialTraversal)
	      computeTraversalInfo(r, ti, counter, maxTips, numBranches, partialTraversal);
	    if(! p->x)
	      getxnode(p);	 
	    
	    /* make sure that everything is consistent now */

	    assert(p->x && r->x);

	    /* store data for p, q, r in the traversal descriptor */

	    ti[*counter].tipCase = TIP_INNER;
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;

	    for(i = 0; i < numBranches; i++)
	      {	
		ti[*counter].qz[i] = q->z[i];
		ti[*counter].rz[i] = r->z[i];
	      }

	    *counter = *counter + 1;
	  }
	else
	  {
	    /* same as above, only now q and r are inner nodes. Hence if they are not 
	       oriented correctly they will need to be recomputed and we need to descend into the 
	       respective subtrees to check if everything is consistent in there, potentially expanding 
	       the traversal descriptor */
	   
	    if(! q->x || !partialTraversal)
	      computeTraversalInfo(q, ti, counter, maxTips, numBranches, partialTraversal);
	    if(! r->x || !partialTraversal)
	      computeTraversalInfo(r, ti, counter, maxTips, numBranches, partialTraversal);
	    if(! p->x)
	      getxnode(p);
	     
	    /* check that the vector orientations are consistent now */

	    assert(p->x && r->x && q->x);

	    ti[*counter].tipCase = INNER_INNER;
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;

	    for(i = 0; i < numBranches; i++)
	      {	
		ti[*counter].qz[i] = q->z[i];
		ti[*counter].rz[i] = r->z[i];
	      }

	    *counter = *counter + 1;
	  }
      }
  }
}

/* below are the optimized unrolled, and vectorized versions of the above generi cfunctions 
   for computing the conditional likelihood at p given child nodes q and r. The actual implementation is located at the end/bottom of this 
   file.
*/

static void newviewGTRCAT_BINARY( int tipCase,  double *EV,  int *cptr,
                                  double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
                                  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                  int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling);

static void newviewGTRGAMMA_BINARY(int tipCase,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *EV, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling
				   );

boolean isGap(unsigned int *x, int pos)
{
  return (x[pos / 32] & mask32[pos % 32]);
}

boolean noGap(unsigned int *x, int pos)
{
  return (!(x[pos / 32] & mask32[pos % 32]));
}

/* now this is the function that just iterates over the length of the traversal descriptor and 
   just computes the conditional likelihhod arrays in the order given by the descriptor.
   So in a sense, this function has no clue that there is any tree-like structure 
   in the traversal descriptor, it just operates on an array of structs of given length */ 


extern const char inverseMeaningDNA[16]; 

void newviewIterative (tree *tr, int startIndex)
{
  traversalInfo 
    *ti   = tr->td[0].ti;

  int 
    i;

    /* loop over traversal descriptor length. Note that on average we only re-compute the conditionals on 3 -4
       nodes in RAxML */

  for(i = startIndex; i < tr->td[0].count; i++)
    {
      traversalInfo 
	*tInfo = &ti[i];
      
      int 
	model;
      
#ifdef _USE_OMP
#pragma omp parallel for
#endif
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  /* check if this partition has to be processed now - otherwise no need to compute P matrix */
	  if(!tr->td[0].executeModel[model] || tr->partitionData[model].width == 0)
	    continue;

	  int
	    categories,
	    states = tr->partitionData[model].states;
	  
	  double
	    qz,
	    rz,
	    *rateCategories,
	    *left = tr->partitionData[model].left,
	    *right = tr->partitionData[model].right;
	  
	  /* figure out what kind of rate heterogeneity approach we are using */
	  if(tr->rateHetModel == CAT)
	    {
	      rateCategories = tr->partitionData[model].perSiteRates;
	      categories = tr->partitionData[model].numberOfCategories;
	    }
	  else
	    {
	      rateCategories = tr->partitionData[model].gammaRates;
	      categories = 4;
	    }
	  
	  /* if we use per-partition branch length optimization
	     get the branch length of partition model and take the log otherwise
	     use the joint branch length among all partitions that is always stored
	     at index [0] */
	  if(tr->numBranches > 1)
	    {
	      qz = tInfo->qz[model];
	      rz = tInfo->rz[model];
	    }
	  else
	    {
	      qz = tInfo->qz[0];
	      rz = tInfo->rz[0];
	    }
	  
	  qz = (qz > zmin) ? log(qz) : log(zmin);
	  rz = (rz > zmin) ? log(rz) : log(zmin);

	  /* compute the left and right P matrices */
#ifdef __MIC_NATIVE
	  switch (tr->partitionData[model].states)
	    {
	    case 2: /* BINARY data */
	      assert(0 && "Binary data model is not implemented on Intel MIC");
	      break;
	    case 4: /* DNA data */
	      {
		makeP_DNA_MIC(qz, rz, rateCategories,   tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN, categories,
			      left, right, tr->saveMemory, tr->maxCategories);
		
		precomputeTips_DNA_MIC(tInfo->tipCase, tr->partitionData[model].tipVector,
				       left, right,
				       tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight,
				       categories);
	      } 
	      break;
	    case 20: /* AA data */
	      {
		if(tr->partitionData[model].protModels == LG4M || tr->partitionData[model].protModels == LG4X)
		  {
		    makeP_PROT_LG4_MIC(qz, rz, tr->partitionData[model].gammaRates,
				       tr->partitionData[model].EI_LG4, tr->partitionData[model].EIGN_LG4,
				       4, left, right);
		    
		    precomputeTips_PROT_LG4_MIC(tInfo->tipCase, tr->partitionData[model].tipVector_LG4,
						left, right,
						tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight,
						categories);
		  }
		else
		  {
		    makeP_PROT_MIC(qz, rz, rateCategories, tr->partitionData[model].EI,
				   tr->partitionData[model].EIGN, categories,
				   left, right, tr->saveMemory, tr->maxCategories);
		    
		    precomputeTips_PROT_MIC(tInfo->tipCase, tr->partitionData[model].tipVector,
					    left, right,
					    tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight,
					    categories);
		  }
	      } 
	      break;
	    default:
	      assert(0);
	    }
#else
	  if(tr->partitionData[model].protModels == LG4M || tr->partitionData[model].protModels == LG4X)
	    makeP_FlexLG4(qz, rz, tr->partitionData[model].gammaRates,
			  tr->partitionData[model].EI_LG4,
			  tr->partitionData[model].EIGN_LG4,
			  4, left, right, 20);
	  else
	    makeP(qz, rz, rateCategories,   tr->partitionData[model].EI,
		  tr->partitionData[model].EIGN, categories,
		  left, right, tr->saveMemory, tr->maxCategories, states);
#endif
      } // for model


      /* now loop over all partitions for nodes p, q, and r of the current traversal vector entry */
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
	    size_t
	      width  = 0,
	      offset = 0;
	    
	    double
	      *left     = (double*)NULL,
	      *right    = (double*)NULL;
	    
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
		assert(tid == pAss->procId);
		
		model  = pAss->partitionId;
		width  = pAss->width;
		offset = pAss->offset;
		
		left  = tr->partitionData[model].left;
		right = tr->partitionData[model].right;
		globalScaler = tr->partitionData[model].threadGlobalScaler[tid];
	      }
	    else
	      break;
#else
	    model = m;	    

	    /* number of sites in this partition */
	    width  = (size_t)tr->partitionData[model].width;
	    offset = 0;

	    /* set the pointers to the left and right P matrices to the pre-allocated memory space for storing them */
	    
	    left  = tr->partitionData[model].left;
	    right = tr->partitionData[model].right;
	    globalScaler = tr->partitionData[model].globalScaler;
#endif

	    /* this conditional statement is exactly identical to what we do in evaluateIterative */
	    if(tr->td[0].executeModel[model] && width > 0)
	      {	      
		double
		  *x1_start = (double*)NULL,
		  *x2_start = (double*)NULL,
		  *x3_start = (double*)NULL, //tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1],
		  *x1_gapColumn = (double*)NULL,
		  *x2_gapColumn = (double*)NULL,
		  *x3_gapColumn = (double*)NULL;

		int
		  scalerIncrement = 0,
		
		  /* integer wieght vector with pattern compression weights */
		  *wgt = tr->partitionData[model].wgt + offset,

		  /* integer rate category vector (for each pattern, _number_ of PSR category assigned to it, NOT actual rate!) */
		  *rateCategory = tr->partitionData[model].rateCategory + offset;

		unsigned int
		  *x1_gap = (unsigned int*)NULL,
		  *x2_gap = (unsigned int*)NULL,
		  *x3_gap = (unsigned int*)NULL;

		unsigned char
		  *tipX1 = (unsigned char *)NULL,
		  *tipX2 = (unsigned char *)NULL;	
	      
		size_t
		  gapOffset = 0,
		  rateHet = discreteRateCategories(tr->rateHetModel),
		  
		  /* get the number of states in the data stored in partition model */
		  
		  states = (size_t)tr->partitionData[model].states,	
		  
		  /* span for single alignment site (in doubles!) */
		  span = rateHet * states,
		  x_offset = offset * (size_t)span,
		  
		  
		  /* get the length of the current likelihood array stored at node p. This is 
		     important mainly for the SEV-based memory saving option described in here:
		     
		     F. Izquierdo-Carrasco, S.A. Smith, A. Stamatakis: "Algorithms, Data Structures, and Numerics for Likelihood-based Phylogenetic Inference of Huge Trees".
		     
		     So tr->partitionData[model].xSpaceVector[i] provides the length of the allocated conditional array of partition model 
		     and node i 
		  */
		
		  availableLength = tr->partitionData[model].xSpaceVector[(tInfo->pNumber - tr->mxtips - 1)],
		  requiredLength = 0;	     

		x3_start = tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1] + x_offset;

		/* memory saving stuff, not important right now, but if you are interested ask Fernando */
		if(tr->saveMemory)
		  {
		    size_t
		      j,
		      setBits = 0;		  
		    
		    gapOffset = states * (size_t)getUndetermined(tr->partitionData[model].dataType);
		    
		    x1_gap = &(tr->partitionData[model].gapVector[tInfo->qNumber * tr->partitionData[model].gapVectorLength]);
		    x2_gap = &(tr->partitionData[model].gapVector[tInfo->rNumber * tr->partitionData[model].gapVectorLength]);
		    x3_gap = &(tr->partitionData[model].gapVector[tInfo->pNumber * tr->partitionData[model].gapVectorLength]);		      		  
		    
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
#ifndef _USE_OMP
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
#endif

		/* now just set the pointers for data accesses in the newview() implementations above to the corresponding values 
		   according to the tip case */
		
		switch(tInfo->tipCase)
		  {
		  case TIP_TIP:		  
		    tipX1    = tr->partitionData[model].yVector[tInfo->qNumber] + offset;
		    tipX2    = tr->partitionData[model].yVector[tInfo->rNumber] + offset;
		    
		    if(tr->saveMemory)
		      {
			x1_gapColumn   = &(tr->partitionData[model].tipVector[gapOffset]);
			x2_gapColumn   = &(tr->partitionData[model].tipVector[gapOffset]);		    
			x3_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->pNumber - tr->mxtips - 1) * states * rateHet];		    
		      }
		    
		    break;
		  case TIP_INNER:		 
		    tipX1    =  tr->partitionData[model].yVector[tInfo->qNumber] + offset;
		    x2_start = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1] + x_offset;
		    
		    if(tr->saveMemory)
		      {	
			x1_gapColumn   = &(tr->partitionData[model].tipVector[gapOffset]);	     
			x2_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->rNumber - tr->mxtips - 1) * states * rateHet];
			x3_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->pNumber - tr->mxtips - 1) * states * rateHet];
		      }
		    
		    break;
		  case INNER_INNER:		 		 
		    x1_start       = tr->partitionData[model].xVector[tInfo->qNumber - tr->mxtips - 1] + x_offset;
		    x2_start       = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1] + x_offset;
		    
		    if(tr->saveMemory)
		      {
			x1_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->qNumber - tr->mxtips - 1) * states * rateHet];
			x2_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->rNumber - tr->mxtips - 1) * states * rateHet];
			x3_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->pNumber - tr->mxtips - 1) * states * rateHet];
		      }
		    
		    break;
		  default:
		    assert(0);
		  }
		
	      /* dedicated highly optimized functions. Analogously to the functions in evaluateGeneric() 
		 we also siwtch over the state number */

	      switch(states)
		{		
		case 2:
#ifdef __MIC_NATIVE
 	      assert(0 && "Binary data model is not implemented on Intel MIC");
#else
		  assert(!tr->saveMemory);
		  if(tr->rateHetModel == CAT)
		    newviewGTRCAT_BINARY(tInfo->tipCase,  tr->partitionData[model].EV,  rateCategory,
					 x1_start,  x2_start,  x3_start, tr->partitionData[model].tipVector,
					 (int*)NULL, tipX1, tipX2,
					 width, left, right, wgt, &scalerIncrement, TRUE);
		  else
		    newviewGTRGAMMA_BINARY(tInfo->tipCase,
					   x1_start, x2_start, x3_start,
					   tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					   (int *)NULL, tipX1, tipX2,
					   width, left, right, wgt, &scalerIncrement, TRUE);		 
#endif
		  break;
		case 4:	/* DNA */
		  if(tr->rateHetModel == CAT)
		    {		    		     
		      if(tr->saveMemory)
#ifdef __MIC_NATIVE
		     assert(0 && "Neither CAT model of rate heterogeneity nor memory saving are implemented on Intel MIC");
#else
			newviewGTRCAT_AVX_GAPPED_SAVE(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
						      x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						      (int*)NULL, tipX1, tipX2,
						      width, left, right, wgt, &scalerIncrement, TRUE, x1_gap, x2_gap, x3_gap,
						      x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
#endif
		      else
#ifdef __MIC_NATIVE
		     assert(0 && "CAT model of rate heterogeneity is not implemented on Intel MIC");
#else
			newviewGTRCAT_AVX(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
					  x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					  tipX1, tipX2,
					  width, left, right, wgt, &scalerIncrement);
#endif
		    }
		  else
		    {
		      
		       
		       if(tr->saveMemory)
#ifdef __MIC_NATIVE
		     assert(0 && "Memory saving is not implemented on Intel MIC");
#else
			 newviewGTRGAMMA_AVX_GAPPED_SAVE(tInfo->tipCase,
							 x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector, (int*)NULL,
							 tipX1, tipX2,
							 width, left, right, wgt, &scalerIncrement, TRUE,
							 x1_gap, x2_gap, x3_gap, 
							 x1_gapColumn, x2_gapColumn, x3_gapColumn);
#endif
		       else
#ifdef __MIC_NATIVE
			 newviewGTRGAMMA_MIC(tInfo->tipCase,
				  x1_start, x2_start, x3_start, tr->partitionData[model].mic_EV, tr->partitionData[model].tipVector,
				  tipX1, tipX2,
				  width, left, right, wgt, &scalerIncrement,
				  tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight);
#else
			 newviewGTRGAMMA_AVX(tInfo->tipCase,
					     x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					     tipX1, tipX2,
					     width, left, right, wgt, &scalerIncrement);
#endif
		    }
		
		  break;		    
		case 20: /* proteins */

		  if(tr->rateHetModel == CAT)
		    {		     
		      if(tr->saveMemory)
			{
#ifdef __MIC_NATIVE
		     assert(0 && "Neither CAT model of rate heterogeneity nor memory saving are implemented on Intel MIC");
#else
			  newviewGTRCATPROT_AVX_GAPPED_SAVE(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
							    x1_start, x2_start, x3_start, tr->partitionData[model].tipVector, (int*)NULL,
							    tipX1, tipX2, width, left, right, wgt, &scalerIncrement, TRUE, x1_gap, x2_gap, x3_gap,
							    x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
#endif
			}
		      else
			{			 			
#ifdef __MIC_NATIVE
		     assert(0 && "CAT model of rate heterogeneity is not implemented on Intel MIC");
#else
			  newviewGTRCATPROT_AVX(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
						x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						tipX1, tipX2, width, left, right, wgt, &scalerIncrement);
#endif
			}
		    }
		  else
		    {		    			 			  
		      if(tr->saveMemory)
			{
#ifdef __MIC_NATIVE
		     assert(0 && "Memory saving is not implemented on Intel MIC");
#else
			  newviewGTRGAMMAPROT_AVX_GAPPED_SAVE(tInfo->tipCase,
							      x1_start, x2_start, x3_start,
							      tr->partitionData[model].EV,
							      tr->partitionData[model].tipVector, (int*)NULL,
							      tipX1, tipX2,
							      width, left, right, wgt, &scalerIncrement, TRUE,
							      x1_gap, x2_gap, x3_gap,
							      x1_gapColumn, x2_gapColumn, x3_gapColumn);
#endif
			}
		      else
			{
			  if(tr->partitionData[model].protModels == LG4M || tr->partitionData[model].protModels == LG4X)
			    {
#ifdef __MIC_NATIVE
			      newviewGTRGAMMAPROT_LG4_MIC(tInfo->tipCase,
							x1_start, x2_start, x3_start, tr->partitionData[model].mic_EV, tr->partitionData[model].mic_tipVector,
							tipX1, tipX2,
							width, left, right, wgt, &scalerIncrement,
							tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight);
#else
			      newviewGTRGAMMAPROT_AVX_LG4(tInfo->tipCase,
							  x1_start, x2_start, x3_start,
							  tr->partitionData[model].EV_LG4,
							  tr->partitionData[model].tipVector_LG4,
							  (int*)NULL, tipX1, tipX2,
							  width, left, right, wgt, &scalerIncrement, TRUE);
#endif			    
			    }
			  else
			    {
#ifdef __MIC_NATIVE
			      newviewGTRGAMMAPROT_MIC(tInfo->tipCase,
							x1_start, x2_start, x3_start, tr->partitionData[model].mic_EV, tr->partitionData[model].mic_tipVector,
							tipX1, tipX2,
							width, left, right, wgt, &scalerIncrement,
							tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight);
#else
			      newviewGTRGAMMAPROT_AVX(tInfo->tipCase,
						      x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
						      tipX1, tipX2,
						      width, left, right, wgt, &scalerIncrement);
#endif
			    }
			}
		    }	
		  break;	
		default:
		  assert(0);
		}

	      /* important step, here we essentiallt recursively compute the number of scaling multiplications 
		 at node p: it's the sum of the number of scaling multiplications already conducted 
		 for computing nodes q and r plus the scaling multiplications done at node p */

	      globalScaler[tInfo->pNumber] =
		globalScaler[tInfo->qNumber] +
		globalScaler[tInfo->rNumber] +
		(unsigned int)scalerIncrement;

	      /* check that we are not getting an integer overflow ! */

	      assert(globalScaler[tInfo->pNumber] < INT_MAX);
	    }	
	} // for model
    }  // omp parallel block
  }  // for traversal
}


/* here is the generic function that could be called from the user program 
   it re-computes the vector at node p (regardless of whether it's orientation is 
   correct and then it also re-computes reciursively the likelihood arrays 
   in the subtrees of p as needed and if needed */

void newviewGeneric (tree *tr, nodeptr p, boolean masked)
{  
  /* if it's a tip there is nothing to do */

  if(isTip(p->number, tr->mxtips))
    return;
  
  /* the first entry of the traversal descriptor is always reserved for evaluate or branch length optimization calls,
     hence we start filling the array at the second entry with index one. This is not very nice and should be fixed 
     at some point */

  tr->td[0].count = 0;

  /* compute the traversal descriptor */
  computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches, TRUE);

  /* the traversal descriptor has been recomputed -> not sure if it really always changes, something to 
     optimize in the future */
  tr->td[0].traversalHasChanged = TRUE;
  
  /* We do a masked newview, i.e., do not execute newvies for each partition, when for example 
     doing a branch length optimization on the entire tree when branches are estimated on a per partition basis.

     you may imagine that for partition 5 the branch length optimization has already converged whereas 
     for partition 6 we still need to go over the tree again.

     This is explained in more detail in:

     A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009

     The external boolean array tr->partitionConverged[] contains exactly that information and is copied 
     to executeModel and subsequently to the executeMask of the traversal descriptor 

  */


  if(masked)
    {
      int model;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->partitionConverged[model])
	    tr->executeModel[model] = FALSE;
	  else
	    tr->executeModel[model] = TRUE;
	}
    }

  /* if there is something to re-compute */

  if(tr->td[0].count > 0)
    {
      /* store execute mask in traversal descriptor */

      storeExecuteMaskInTraversalDescriptor(tr);           
      newviewIterative(tr, 0);
    }

  /* clean up */

  if(masked)
    {
      int model;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	tr->executeModel[model] = TRUE;
    }

  tr->td[0].traversalHasChanged = FALSE;
}


/* optimized function implementations */


/*** BINARY DATA functions *****/

static void newviewGTRCAT_BINARY( int tipCase,  double *EV,  int *cptr,
                                  double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
                                  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                  int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1, *x2, *x3;
  int i, l, scale, addScale = 0;

  switch(tipCase)
    {
    case TIP_TIP:
      {
        for(i = 0; i < n; i++)
          {
            x1 = &(tipVector[2 * tipX1[i]]);
            x2 = &(tipVector[2 * tipX2[i]]);
            x3 = &x3_start[2 * i];         

            le =  &left[cptr[i] * 4];
            ri =  &right[cptr[i] * 4];

            _mm_store_pd(x3, _mm_setzero_pd());     
                     
            for(l = 0; l < 2; l++)
              {                                                                                                                          
                __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
                __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
                
                al = _mm_hadd_pd(al, al);
                ar = _mm_hadd_pd(ar, ar);
                
                al = _mm_mul_pd(al, ar);
                
                __m128d vv  = _mm_load_pd(x3);
                __m128d EVV = _mm_load_pd(&EV[2 * l]);
                
                vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
                
                _mm_store_pd(x3, vv);                                                     
              }            
          }
      }
      break;
    case TIP_INNER:
      {
        for (i = 0; i < n; i++)
          {
            x1 = &(tipVector[2 * tipX1[i]]);
            x2 = &x2_start[2 * i];
            x3 = &x3_start[2 * i];
            
            le =  &left[cptr[i] * 4];
            ri =  &right[cptr[i] * 4];

            _mm_store_pd(x3, _mm_setzero_pd());     
                     
            for(l = 0; l < 2; l++)
              {                                                                                                                          
                __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
                __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
                
                al = _mm_hadd_pd(al, al);
                ar = _mm_hadd_pd(ar, ar);
                
                al = _mm_mul_pd(al, ar);
                
                __m128d vv  = _mm_load_pd(x3);
                __m128d EVV = _mm_load_pd(&EV[2 * l]);
                
                vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
                
                _mm_store_pd(x3, vv);                                                     
              }  
            
            __m128d minlikelihood_sse = _mm_set1_pd(minlikelihood);
         
            scale = 1;
            
            __m128d v1 = _mm_and_pd(_mm_load_pd(x3), absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;                         
            
            if(scale)
              {
                __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
                
                __m128d ex3v = _mm_load_pd(x3);           
                _mm_store_pd(x3, _mm_mul_pd(ex3v,twoto));                                                 
                
                if(useFastScaling)
                  addScale += wgt[i];
                else
                  ex3[i]  += 1;   
              }                    
          }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
        {
          x1 = &x1_start[2 * i];
          x2 = &x2_start[2 * i];
          x3 = &x3_start[2 * i];

          le = &left[cptr[i] * 4];
          ri = &right[cptr[i] * 4];

          _mm_store_pd(x3, _mm_setzero_pd());       
          
          for(l = 0; l < 2; l++)
            {                                                                                                                            
              __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
              __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
              
              al = _mm_hadd_pd(al, al);
              ar = _mm_hadd_pd(ar, ar);
              
              al = _mm_mul_pd(al, ar);
              
              __m128d vv  = _mm_load_pd(x3);
              __m128d EVV = _mm_load_pd(&EV[2 * l]);
              
              vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
              
              _mm_store_pd(x3, vv);                                                       
            }                             

          __m128d minlikelihood_sse = _mm_set1_pd(minlikelihood);
         
          scale = 1;
                  
          __m128d v1 = _mm_and_pd(_mm_load_pd(x3), absMask.m);
          v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
          if(_mm_movemask_pd( v1 ) != 3)
            scale = 0;                   
         
          if(scale)
            {
              __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
                    
              __m128d ex3v = _mm_load_pd(x3);             
              _mm_store_pd(x3, _mm_mul_pd(ex3v,twoto));                                           
             
              if(useFastScaling)
                addScale += wgt[i];
              else
                ex3[i]  += 1;     
           }             
        }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

static void newviewGTRGAMMA_BINARY(int tipCase,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *EV, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling
				   )
{
  double
    *x1, *x2, *x3;
 
  int i, k, l, scale, addScale = 0; 

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
       {
	 x1  = &(tipVector[2 * tipX1[i]]);
	 x2  = &(tipVector[2 * tipX2[i]]);
	 
	 for(k = 0; k < 4; k++)
	   {	     	     	    
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
       }
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
       {
	 x1  = &(tipVector[2 * tipX1[i]]);
	 
	 for(k = 0; k < 4; k++)
	   {	     	     
	     x2 = &(x2_start[8 * i + 2 * k]);
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
	
	 x3 = &(x3_start[8 * i]);
	 __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	 
	 scale = 1;
	 for(l = 0; scale && (l < 8); l += 2)
	   {
	     __m128d vv = _mm_load_pd(&x3[l]);
	     __m128d v1 = _mm_and_pd(vv, absMask.m);
	     v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	     if(_mm_movemask_pd( v1 ) != 3)
	       scale = 0;
	   }	    	         
	 
	 if(scale)
	   {
	     __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	     
	     for(l = 0; l < 8; l+=2)
	       {
		 __m128d ex3v = _mm_load_pd(&x3[l]);		  
		 _mm_store_pd(&x3[l], _mm_mul_pd(ex3v,twoto));	
	       }		   		  
	     
	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }	 
       }      
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {	 
	 for(k = 0; k < 4; k++)
	   {	     
	     x1 = &(x1_start[8 * i + 2 * k]);
	     x2 = &(x2_start[8 * i + 2 * k]);
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
	
	 x3 = &(x3_start[8 * i]);
	 __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	 
	 scale = 1;
	 for(l = 0; scale && (l < 8); l += 2)
	   {
	     __m128d vv = _mm_load_pd(&x3[l]);
	     __m128d v1 = _mm_and_pd(vv, absMask.m);
	     v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	     if(_mm_movemask_pd( v1 ) != 3)
	       scale = 0;
	   }	    	         
	 
	 if(scale)
	   {
	     __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	     
	     for(l = 0; l < 8; l+=2)
	       {
		 __m128d ex3v = _mm_load_pd(&x3[l]);		  
		 _mm_store_pd(&x3[l], _mm_mul_pd(ex3v,twoto));	
	       }		   		  
	     
	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }	 
       }
      break;

    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}


/**** BINARY DATA functions end ****/

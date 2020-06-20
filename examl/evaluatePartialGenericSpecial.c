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

#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse3.h"

#if !defined(__MIC_NATIVE)
static inline void computeVectorGTRCATPROT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
					   traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					   unsigned  char **yVector, int mxtips);

static double evaluatePartialGTRCATPROT(int i, double ki, int counter,  traversalInfo *ti, double qz,
					int w, double *EIGN, double *EI, double *EV,
					double *tipVector, unsigned char **yVector, 
					int branchReference, int mxtips);

static inline void computeVectorGTRGAMMAPROT(double *lVector, int *eVector, double *gammaRates, int i, double qz, double rz,
					     traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					     unsigned  char **yVector, int mxtips);

static double evaluatePartialGTRGAMMAPROT(int i, int counter,  traversalInfo *ti, double qz,
					  int w, double *EIGN, double *EI, double *EV,
					  double *tipVector, unsigned char **yVector, 
					  double *gammaRates,
					  int branchReference, int mxtips);

static inline void computeVectorGTRCAT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
				       traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
				       unsigned char **yVector, int mxtips);

static double evaluatePartialGTRCAT(int i, double ki, int counter,  traversalInfo *ti, double qz,
				    int w, double *EIGN, double *EI, double *EV,
				    double *tipVector, unsigned  char **yVector, 
				    int branchReference, int mxtips);


static inline void computeVectorGTRCAT_BINARY(double *lVector, int *eVector, double ki, int i, double qz, double rz,
					      traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					      unsigned char **yVector, int mxtips);

static double evaluatePartialGTRCAT_BINARY(int i, double ki, int counter,  traversalInfo *ti, double qz,
					   int w, double *EIGN, double *EI, double *EV,
					   double *tipVector, unsigned  char **yVector, 
					   int branchReference, int mxtips);
#endif

double evaluatePartialGeneric (tree *tr, int i, double ki, int _model)
{
  double result;
  int 
    branchReference,
    states = tr->partitionData[_model].states;
    
  int 
    index;

  index = i;

  
  if(tr->numBranches > 1)
    branchReference = _model;
  else
    branchReference = 0;


#if defined(__MIC_NATIVE)
if (tr->rateHetModel == CAT)
    result = evaluatePartialCAT_FLEX(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference],
                     tr->partitionData[_model].wgt[index],
                     tr->partitionData[_model].EIGN,
                     tr->partitionData[_model].EI,
                     tr->partitionData[_model].EV,
                     tr->partitionData[_model].tipVector,
                     tr->partitionData[_model].yVector, branchReference, tr->mxtips, states);
else
    assert(0);

#else
  switch(states)
    {
    case 2:
      assert(!tr->saveMemory);
      assert(tr->rateHetModel == CAT);

       result = evaluatePartialGTRCAT_BINARY(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					    tr->partitionData[_model].wgt[index],
					    tr->partitionData[_model].EIGN, 
					    tr->partitionData[_model].EI, 
					    tr->partitionData[_model].EV,
					    tr->partitionData[_model].tipVector,
					    tr->partitionData[_model].yVector, branchReference, tr->mxtips);

      

      break;
    case 4:   /* DNA */
      assert(tr->rateHetModel == CAT);  
      
      result = evaluatePartialGTRCAT(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
				     tr->partitionData[_model].wgt[index],
				     tr->partitionData[_model].EIGN, 
				     tr->partitionData[_model].EI, 
				     tr->partitionData[_model].EV,
				     tr->partitionData[_model].tipVector,
				     tr->partitionData[_model].yVector, branchReference, tr->mxtips);
      break;
    case 20: /* proteins */
      if(tr->rateHetModel == CAT)
	result = evaluatePartialGTRCATPROT(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					   tr->partitionData[_model].wgt[index],
					   tr->partitionData[_model].EIGN, 
					   tr->partitionData[_model].EI, 
					   tr->partitionData[_model].EV,
					   tr->partitionData[_model].tipVector, 
					   tr->partitionData[_model].yVector, branchReference, tr->mxtips);
      else
	result =  evaluatePartialGTRGAMMAPROT(index, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					      tr->partitionData[_model].wgt[index],
					      tr->partitionData[_model].EIGN, 
					      tr->partitionData[_model].EI, 
					      tr->partitionData[_model].EV,
					      tr->partitionData[_model].tipVector, 
					      tr->partitionData[_model].yVector, 
					      tr->partitionData[_model].gammaRates,
					      branchReference, tr->mxtips);
      break;   
    default:
      assert(0);
    }
  #endif
 

  return result;
}


static inline void computeVectorGTRCAT_BINARY(double *lVector, int *eVector, double ki, int i, double qz, double rz,
					      traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					      unsigned char **yVector, int mxtips)
{       
  double  d1, d2,  ump_x1, ump_x2, x1px2[2], lz1, lz2; 
  double *x1, *x2, *x3;
  int 
    j, k,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &lVector[2 * (pNumber  - mxtips)];  

  switch(ti->tipCase)
    {
    case TIP_TIP:     
      x1 = &(tipVector[2 * yVector[qNumber][i]]);
      x2 = &(tipVector[2 * yVector[rNumber][i]]);   
      break;
    case TIP_INNER:     
      x1 = &(tipVector[2 * yVector[qNumber][i]]);
      x2 = &lVector[2 * (rNumber - mxtips)];                    
      break;
    case INNER_INNER:            
      x1 = &lVector[2 * (qNumber - mxtips)];
      x2 = &lVector[2 * (rNumber - mxtips)];               
      break;
    default:
      assert(0);
    }
     
  lz1 = qz * ki;  
  lz2 = rz * ki;
  
 
  d1 = x1[1] * EXP(EIGN[1] * lz1);
  d2 = x2[1] * EXP(EIGN[1] * lz2);	        
 
  for(j = 0; j < 2; j++)
    {     
      ump_x1 = x1[0];
      ump_x2 = x2[0];
      
      ump_x1 += d1 * EI[j * 2 + 1];
      ump_x2 += d2 * EI[j * 2 + 1];
	
      x1px2[j] = ump_x1 * ump_x2;
    }
  
  for(j = 0; j < 2; j++)
    x3[j] = 0.0;

  for(j = 0; j < 2; j++)          
    for(k = 0; k < 2; k++)	
      x3[k] +=  x1px2[j] *  EV[2 * j + k];	   
      
  
  if (x3[0] < minlikelihood && x3[0] > minusminlikelihood &&
      x3[1] < minlikelihood && x3[1] > minusminlikelihood
      )
    {	     
      x3[0]   *= twotothe256;
      x3[1]   *= twotothe256;     
      *eVector = *eVector + 1;
    }	              

  return;
}

static double evaluatePartialGTRCAT_BINARY(int i, double ki, int counter,  traversalInfo *ti, double qz,
					   int w, double *EIGN, double *EI, double *EV,
					   double *tipVector, unsigned  char **yVector, 
					   int branchReference, int mxtips)
{
  double lz, term;       
  double  d;
  double   *x1, *x2; 
  int scale = 0, k;
  double *lVector = (double *)malloc(sizeof(double) * 2 * mxtips);  
  traversalInfo *trav = &ti[0];
 
  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[2 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)  
    {
      double 
	qz = ti[k].qz[branchReference],
	rz = ti[k].rz[branchReference];
      
      qz = (qz > zmin) ? log(qz) : log(zmin);
      rz = (rz > zmin) ? log(rz) : log(zmin);

      computeVectorGTRCAT_BINARY(lVector, &scale, ki, i, qz, rz, &ti[k], 
				 EIGN, EI, EV, 
				 tipVector, yVector, mxtips);       
    }
   
  x2 = &lVector[2 * (trav->qNumber - mxtips)];
     
  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
       
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;  
  
  d = EXP(EIGN[1] * lz);
  
  term =  x1[0] * x2[0];
  term += x1[1] * x2[1] * d; 

  term = LOG(FABS(term)) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
  
  return  term;
}



static inline void computeVectorGTRCATPROT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
					   traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					   unsigned  char **yVector, int mxtips)
{       
  double   *x1, *x2, *x3;  
  int
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &(lVector[20 * (pNumber  - mxtips)]);     

  switch(ti->tipCase)
    {
    case TIP_TIP:    
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(tipVector[20 * yVector[rNumber][i]]);     
      break;
    case TIP_INNER:     
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(  lVector[20 * (rNumber - mxtips)]);                    
      break;
    case INNER_INNER:            
      x1 = &(lVector[20 * (qNumber - mxtips)]);
      x2 = &(lVector[20 * (rNumber - mxtips)]);                 
      break;    
    default:
      assert(0);
    }
     
  {
    double  
      e1[20] __attribute__ ((aligned (BYTE_ALIGNMENT))),
      e2[20] __attribute__ ((aligned (BYTE_ALIGNMENT))),
      d1[20] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
      d2[20] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
      lz1, 
      lz2;  
    
    int 
      l, 
      k, 
      scale;
     
    lz1 = qz * ki;            
    lz2 = rz * ki;        

    e1[0] = 1.0;
    e2[0] = 1.0;
    
    for(l = 1; l < 20; l++)
      {
	e1[l] = EXP(EIGN[l] * lz1);
	e2[l] = EXP(EIGN[l] * lz2);
      }

    for(l = 0; l < 20; l+=2)
      {
	__m128d d1v = _mm_mul_pd(_mm_load_pd(&x1[l]), _mm_load_pd(&e1[l]));
	__m128d d2v = _mm_mul_pd(_mm_load_pd(&x2[l]), _mm_load_pd(&e2[l]));
	
	_mm_store_pd(&d1[l], d1v);
	_mm_store_pd(&d2[l], d2v);	
      }

    __m128d zero = _mm_setzero_pd();

    for(l = 0; l < 20; l+=2)
      _mm_store_pd(&x3[l], zero);
                
    for(l = 0; l < 20; l++)
      { 	      
	double *ev = &EV[l * 20];
	__m128d ump_x1v = _mm_setzero_pd();
	__m128d ump_x2v = _mm_setzero_pd();
	__m128d x1px2v;

	for(k = 0; k < 20; k+=2)
	  {       
	    __m128d eiv = _mm_load_pd(&EI[20 * l + k]);
	    __m128d d1v = _mm_load_pd(&d1[k]);
	    __m128d d2v = _mm_load_pd(&d2[k]);
	    ump_x1v = _mm_add_pd(ump_x1v, _mm_mul_pd(d1v, eiv));
	    ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(d2v, eiv));	  
	  }

	ump_x1v = _mm_hadd_pd(ump_x1v, ump_x1v);
	ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);

	x1px2v = _mm_mul_pd(ump_x1v, ump_x2v);

	for(k = 0; k < 20; k+=2)
	  {
	    __m128d ex3v = _mm_load_pd(&x3[k]);
	    __m128d EVV  = _mm_load_pd(&ev[k]);
	    ex3v = _mm_add_pd(ex3v, _mm_mul_pd(x1px2v, EVV));
	    
	    _mm_store_pd(&x3[k], ex3v);	   	   
	  }
      }                      
    
    scale = 1;
    for(l = 0; scale && (l < 20); l++)
      scale = ((x3[l] < minlikelihood) && (x3[l] > minusminlikelihood));	       	      	      	       	       
    
    if(scale)
      {	      
	__m128d twoto = _mm_set_pd(twotothe256, twotothe256);

	for(l = 0; l < 20; l+=2)
	  {
	    __m128d ex3v = _mm_mul_pd(_mm_load_pd(&x3[l]),twoto);
	    _mm_store_pd(&x3[l], ex3v);	
	  }
 	


	*eVector = *eVector + 1;
      }
    
    return;      
  }
}

static double evaluatePartialGTRCATPROT(int i, double ki, int counter,  traversalInfo *ti, double qz,
					int w, double *EIGN, double *EI, double *EV,
					double *tipVector, unsigned char **yVector, 
					int branchReference, int mxtips)
{
  double lz, term;       
  double  d[20];
  double   *x1, *x2; 
  int scale = 0, k, l;
  double 
    *lVector = (double *)malloc_aligned(sizeof(double) * 20 * mxtips),
    myEI[400]  __attribute__ ((aligned (BYTE_ALIGNMENT)));

  traversalInfo *trav = &ti[0];

  

  for(k = 0; k < 20; k++)
    {     
      for(l = 0; l < 20; l++)
	myEI[k * 20 + l] = EI[k * 20 + l];
    }

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[20 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    {
      double 
	qz = ti[k].qz[branchReference],
	rz = ti[k].rz[branchReference];
      
      qz = (qz > zmin) ? log(qz) : log(zmin);
      rz = (rz > zmin) ? log(rz) : log(zmin);

      computeVectorGTRCATPROT(lVector, &scale, ki, i, qz, rz, 
			      &ti[k], EIGN, myEI, EV, 
			      tipVector, yVector, mxtips);       
    }
   
  x2 = &lVector[20 * (trav->qNumber - mxtips)];

       

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;
  
  d[0] = 1.0;
  for(l = 1; l < 20; l++)
    d[l] = EXP (EIGN[l] * lz);

  term = 0.0;
  
  for(l = 0; l < 20; l++)
    term += x1[l] * x2[l] * d[l];   

  term = LOG(FABS(term)) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
  
 
  return  term;
}
static inline void computeVectorGTRGAMMAPROT(double *lVector, int *eVector, double *gammaRates, int i, double qz, double rz,
					     traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					     unsigned  char **yVector, int mxtips)
{       
  double   
    *x1, 
    *x2, 
    *x3;  
  
  int
    s,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber,
    index1[4],
    index2[4];
  
 
  x3  = &(lVector[80 * (pNumber  - mxtips)]);     

  switch(ti->tipCase)
    {
    case TIP_TIP:    
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(tipVector[20 * yVector[rNumber][i]]);     
      for(s = 0; s < 4; s++)
	{
	  index1[s] = 0;
	  index2[s] = 0;
	}
      break;
    case TIP_INNER:     
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(  lVector[80 * (rNumber - mxtips)]);   
      for(s = 0; s < 4; s++)       
	index1[s] = 0;
      for(s = 0; s < 4; s++)     
	index2[s] = s;                     
      break;
    case INNER_INNER:            
      x1 = &(lVector[80 * (qNumber - mxtips)]);
      x2 = &(lVector[80 * (rNumber - mxtips)]); 
      for(s = 0; s < 4; s++)
	{
	  index1[s] = s;
	  index2[s] = s;
	}                
      break;    
    default:
      assert(0);
    }
     
  {
    double  
      e1[20] __attribute__ ((aligned (BYTE_ALIGNMENT))),
      e2[20] __attribute__ ((aligned (BYTE_ALIGNMENT))),
      d1[20] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
      d2[20] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
      lz1, lz2;  
    
    int 
      l, 
      k, 
      scale, 
      j;
     
    for(j = 0; j < 4; j++)
      {
	lz1 = qz * gammaRates[j];            
	lz2 = rz * gammaRates[j];        

	e1[0] = 1.0;
	e2[0] = 1.0;
    
	for(l = 1; l < 20; l++)
	  {
	    e1[l] = EXP(EIGN[l] * lz1);
	    e2[l] = EXP(EIGN[l] * lz2);
	  }

	for(l = 0; l < 20; l+=2)
	  {
	    __m128d d1v = _mm_mul_pd(_mm_load_pd(&x1[20 * index1[j] + l]), _mm_load_pd(&e1[l]));
	    __m128d d2v = _mm_mul_pd(_mm_load_pd(&x2[20 * index2[j] + l]), _mm_load_pd(&e2[l]));
	    
	    _mm_store_pd(&d1[l], d1v);
	    _mm_store_pd(&d2[l], d2v);	
	  }

	__m128d zero = _mm_setzero_pd();

	for(l = 0; l < 20; l+=2)
	  _mm_store_pd(&x3[j * 20 + l], zero);
                
	for(l = 0; l < 20; l++)
	  { 	      
	    double *ev = &EV[l * 20];
	    __m128d ump_x1v = _mm_setzero_pd();
	    __m128d ump_x2v = _mm_setzero_pd();
	    __m128d x1px2v;
	    
	    for(k = 0; k < 20; k+=2)
	      {       
		__m128d eiv = _mm_load_pd(&EI[20 * l + k]);
		__m128d d1v = _mm_load_pd(&d1[k]);
		__m128d d2v = _mm_load_pd(&d2[k]);
		ump_x1v = _mm_add_pd(ump_x1v, _mm_mul_pd(d1v, eiv));
		ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(d2v, eiv));	  
	      }

	    ump_x1v = _mm_hadd_pd(ump_x1v, ump_x1v);
	    ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);

	    x1px2v = _mm_mul_pd(ump_x1v, ump_x2v);

	    for(k = 0; k < 20; k+=2)
	      {
		__m128d ex3v = _mm_load_pd(&x3[j * 20 + k]);
		__m128d EVV  = _mm_load_pd(&ev[k]);
		ex3v = _mm_add_pd(ex3v, _mm_mul_pd(x1px2v, EVV));
		
		_mm_store_pd(&x3[j * 20 + k], ex3v);	   	   
	      }
	  }        
      }
    
    scale = 1;
    for(l = 0; scale && (l < 80); l++)
      scale = ((x3[l] < minlikelihood) && (x3[l] > minusminlikelihood));	       	      	      	       	       
    
    if(scale)
      {	      
	__m128d twoto = _mm_set_pd(twotothe256, twotothe256);

	for(l = 0; l < 80; l+=2)
	  {
	    __m128d ex3v = _mm_mul_pd(_mm_load_pd(&x3[l]),twoto);
	    _mm_store_pd(&x3[l], ex3v);	
	  }

	*eVector = *eVector + 1;
      }
    
    return;      
  }
}


static double evaluatePartialGTRGAMMAPROT(int i, int counter,  traversalInfo *ti, double qz,
					  int w, double *EIGN, double *EI, double *EV,
					  double *tipVector, unsigned char **yVector, 
					  double *gammaRates,
					  int branchReference, int mxtips)
{
  double lz, term;       
  double  d[80];
  double   *x1, *x2; 
  int scale = 0, k, l, j;
  double 
    *lVector = (double *)malloc_aligned(sizeof(double) * 80 * mxtips),
    myEI[400]  __attribute__ ((aligned (BYTE_ALIGNMENT)));

  traversalInfo 
    *trav = &ti[0];

  for(k = 0; k < 20; k++)
    {     
      for(l = 0; l < 20; l++)
	myEI[k * 20 + l] = EI[k * 20 + l];
    }

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[20 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    {
      double 
	qz = ti[k].qz[branchReference],
	rz = ti[k].rz[branchReference];
      
      qz = (qz > zmin) ? log(qz) : log(zmin);
      rz = (rz > zmin) ? log(rz) : log(zmin);

      computeVectorGTRGAMMAPROT(lVector, &scale, gammaRates, i, qz, rz, 
				&ti[k], EIGN, myEI, EV, 
				tipVector, yVector, mxtips);
    }
   
  x2 = &lVector[80 * (trav->qNumber - mxtips)];       

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  
  for(j = 0; j < 4; j++)
    {
      d[20 * j] = 1.0;
      for(l = 1; l < 20; l++)
	d[20 * j + l] = EXP(EIGN[l] * lz * gammaRates[j]);
    }

 
  for(j = 0, term = 0.0; j < 4; j++)
    {
      for(l = 0; l < 20; l++)
	term += x1[l] * x2[20 * j + l] * d[j * 20 + l];	      
    }
  
  term = LOG(0.25 * FABS(term)) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
  
 
  return  term;
}





static inline void computeVectorGTRCAT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
				       traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
				       unsigned char **yVector, int mxtips)
{       
  double  d1[3], d2[3],  ump_x1, ump_x2, x1px2[4], lz1, lz2; 
  double *x1, *x2, *x3;
  int j, k,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &lVector[4 * (pNumber  - mxtips)];  
 

  switch(ti->tipCase)
    {
    case TIP_TIP:     
      x1 = &(tipVector[4 * yVector[qNumber][i]]);
      x2 = &(tipVector[4 * yVector[rNumber][i]]);    
      break;
    case TIP_INNER:     
      x1 = &(tipVector[4 * yVector[qNumber][i]]);
      x2 = &lVector[4 * (rNumber - mxtips)];           
      break;
    case INNER_INNER:            
      x1 = &lVector[4 * (qNumber - mxtips)];
      x2 = &lVector[4 * (rNumber - mxtips)];     
      break;
    default:
      assert(0);
    }
     
  lz1 = qz * ki;  
  lz2 = rz * ki;
  
  for(j = 0; j < 3; j++)
    {
      d1[j] = 
	x1[j + 1] * 
	EXP(EIGN[j + 1] * lz1);
      d2[j] = x2[j + 1] * EXP(EIGN[j + 1] * lz2);	    
    }
 
 
  for(j = 0; j < 4; j++)
    {     
      ump_x1 = x1[0];
      ump_x2 = x2[0];
      for(k = 0; k < 3; k++)
	{
	  ump_x1 += d1[k] * EI[j * 4 + k + 1];
	  ump_x2 += d2[k] * EI[j * 4 + k + 1];
	}
      x1px2[j] = ump_x1 * ump_x2;
    }
  
  for(j = 0; j < 4; j++)
    x3[j] = 0.0;

  for(j = 0; j < 4; j++)          
    for(k = 0; k < 4; k++)	
      x3[k] +=  x1px2[j] *  EV[4 * j + k];	   
      
  
  if (x3[0] < minlikelihood && x3[0] > minusminlikelihood &&
      x3[1] < minlikelihood && x3[1] > minusminlikelihood &&
      x3[2] < minlikelihood && x3[2] > minusminlikelihood &&
      x3[3] < minlikelihood && x3[3] > minusminlikelihood)
    {	     
      x3[0]   *= twotothe256;
      x3[1]   *= twotothe256;
      x3[2]   *= twotothe256;     
      x3[3]   *= twotothe256;     
      *eVector = *eVector + 1;
    }	              

  return;
}








static double evaluatePartialGTRCAT(int i, double ki, int counter,  traversalInfo *ti, double qz,
				    int w, double *EIGN, double *EI, double *EV,
				    double *tipVector, unsigned  char **yVector, 
				    int branchReference, int mxtips)
{
  double lz, term;       
  double  d[3];
  double   *x1, *x2; 
  int scale = 0, k;
  double *lVector = (double *)malloc_aligned(sizeof(double) * 4 * mxtips);    

  traversalInfo *trav = &ti[0];
 
  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[4 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)    
    {
      double 
	qz = ti[k].qz[branchReference],
	rz = ti[k].rz[branchReference];
      
      qz = (qz > zmin) ? log(qz) : log(zmin);
      rz = (rz > zmin) ? log(rz) : log(zmin);

      computeVectorGTRCAT(lVector, &scale, ki, i, qz, rz, &ti[k], 
			  EIGN, EI, EV, 
			  tipVector, yVector, mxtips);       
    }
   
  x2 = &lVector[4 * (trav->qNumber - mxtips)]; 

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
       
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;  
  
  d[0] = EXP (EIGN[1] * lz);
  d[1] = EXP (EIGN[2] * lz);
  d[2] = EXP (EIGN[3] * lz);       	   
  
  term =  x1[0] * x2[0];
  term += x1[1] * x2[1] * d[0];
  term += x1[2] * x2[2] * d[1];
  term += x1[3] * x2[3] * d[2];     

  term = LOG(FABS(term)) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);  

  return  term;
}

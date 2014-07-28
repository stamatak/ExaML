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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands 
 *  of taxa and mixed models". 
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


static const double MNBRAK_GOLD =    1.618034;
static const double MNBRAK_TINY =      1.e-20;
static const double MNBRAK_GLIMIT =     100.0;
static const double BRENT_ZEPS  =      1.e-5;
static const double BRENT_CGOLD =   0.3819660;

extern int optimizeRatesInvocations;
extern int optimizeRateCategoryInvocations;
extern int optimizeAlphaInvocations;
extern int optimizeInvarInvocations;
extern double masterTime;
extern char ratesFileName[1024];
extern char workdir[1024];
extern char run_id[128];
extern char lengthFileName[1024];
extern char lengthFileNameModel[1024];
extern char *protModels[NUM_PROT_MODELS];

extern checkPointState ckp;

extern int processes;
extern int processID;

static void optParamGeneric(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int rateNumber, double lim_inf, double lim_sup, int whichParameterType);

// FLAG for easier debugging of model parameter optimization routines 

//#define _DEBUG_MOD_OPT


/*********************FUNCTIONS FOOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/


static void setRateModel(tree *tr, int model, double rate, int position)
{
  int
    states   = tr->partitionData[model].states,
    numRates = (states * states - states) / 2;

  if(tr->partitionData[model].dataType == DNA_DATA)
    assert(position >= 0 && position < (numRates - 1));
  else
    assert(position >= 0 && position < numRates);

  assert(tr->partitionData[model].dataType != BINARY_DATA); 

  assert(rate >= RATE_MIN && rate <= RATE_MAX);

  if(tr->partitionData[model].nonGTR)
    {    
      int 
	i, 
	index = tr->partitionData[model].symmetryVector[position],
	lastRate = tr->partitionData[model].symmetryVector[numRates - 1];

      
           
      for(i = 0; i < numRates; i++)
	{	
	  if(tr->partitionData[model].symmetryVector[i] == index)
	    {
	      if(index == lastRate)
		tr->partitionData[model].substRates[i] = 1.0;
	      else
		tr->partitionData[model].substRates[i] = rate;      
	    }
	  
	  //printf("%f ", tr->partitionData[model].substRates[i]);
	}
      //printf("\n");
    }
  else
    tr->partitionData[model].substRates[position] = rate;
}


//LIBRARY: the only thing that we will need to do here is to 
//replace linkList by a string and also add some error correction 
//code


static linkageList* initLinkageList(int *linkList, tree *tr)
{
  int 
    k,
    partitions,
    numberOfModels = 0,
    i,
    pos;
  
  linkageList
    *ll = (linkageList*)malloc(sizeof(linkageList));
      
  for(i = 0; i < tr->NumberOfModels; i++)    
    {
      assert(linkList[i] >= 0 && linkList[i] < tr->NumberOfModels);

      if(linkList[i] > numberOfModels)
	numberOfModels = linkList[i];
    }

  numberOfModels++;
  
  ll->entries = numberOfModels;
  ll->ld      = (linkageData*)malloc(sizeof(linkageData) * numberOfModels);

  for(i = 0; i < numberOfModels; i++)
    {
      ll->ld[i].valid = TRUE;
      
      partitions = 0;

      for(k = 0; k < tr->NumberOfModels; k++)	
	if(linkList[k] == i)
	  partitions++;	    

      ll->ld[i].partitions = partitions;
      ll->ld[i].partitionList = (int*)malloc(sizeof(int) * partitions);
      
      for(k = 0, pos = 0; k < tr->NumberOfModels; k++)	
	if(linkList[k] == i)
	  ll->ld[i].partitionList[pos++] = k;
    }

  return ll;
}

static linkageList* initLinkageListString(char *linkageString, tree *tr)
{
  int 
    *list = (int*)malloc(sizeof(int) * tr->NumberOfModels),
    j;

  linkageList 
    *l;

  char
    *str1,
    *saveptr,
    *ch = (char *)calloc(strlen(linkageString), sizeof(char)),
    *token;
  strncpy(ch, linkageString, strlen(linkageString));

  for(j = 0, str1 = ch; ;j++, str1 = (char *)NULL) 
    {
      token = strtok_r(str1, ",", &saveptr);
      if(token == (char *)NULL)
	break;
      assert(j < tr->NumberOfModels);
      list[j] = atoi(token);
      //printf("%d: %s\n", j, token);
    }
  
  free(ch);

  l = initLinkageList(list, tr);
  
  free(list);

  return l;
}

static void init_Q_MatrixSymmetries(char *linkageString, tree *tr, int model)
{
  int 
    states = tr->partitionData[model].states,
    numberOfRates = ((states * states - states) / 2), 
    *list = (int *)malloc(sizeof(int) * numberOfRates),
    j,
    max = -1;

  char
    *str1,
    *saveptr,
    *ch = (char*)calloc(strlen(linkageString), sizeof(char)), 
    *token;
  
  strncpy(ch, linkageString, strlen(linkageString)); 

  for(j = 0, str1 = ch; ;j++, str1 = (char *)NULL) 
    {
      token = strtok_r(str1, ",", &saveptr);
      if(token == (char *)NULL)
	break;
      assert(j < numberOfRates);
      list[j] = atoi(token);     
    }
  
  free(ch);

  for(j = 0; j < numberOfRates; j++)
    {
      assert(list[j] <= j);
      assert(list[j] <= max + 1);
      
      if(list[j] > max)
	max = list[j];
    }  

  assert(numberOfRates == 6);
  
  for(j = 0; j < numberOfRates; j++)  
    tr->partitionData[model].symmetryVector[j] = list[j];    

  //less than the maximum possible number of rate parameters

  if(max < numberOfRates - 1)    
    tr->partitionData[model].nonGTR = TRUE;

  free(list);
}



static linkageList* initLinkageListGTR(tree *tr)
{
  int
    i,
    *links = (int*)malloc(sizeof(int) * tr->NumberOfModels),
    firstAA = tr->NumberOfModels + 2,
    countGTR = 0,
    countOtherModel = 0;
  linkageList* ll;

  for(i = 0; i < tr->NumberOfModels; i++)
    {     
      if(tr->partitionData[i].dataType == AA_DATA)
	{
	  if(tr->partitionData[i].protModels == GTR)
	    {
	      if(i < firstAA)
		firstAA = i;
	      countGTR++;
	    }
	  else
	    countOtherModel++;
	}
    }
  
  assert((countGTR > 0 && countOtherModel == 0) || (countGTR == 0 && countOtherModel > 0) ||  (countGTR == 0 && countOtherModel == 0));

  if(countGTR == 0)
    {
      for(i = 0; i < tr->NumberOfModels; i++)
	links[i] = i;
    }
  else
    {
      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  switch(tr->partitionData[i].dataType)
	    {	   
	    case DNA_DATA:
	    case BINARY_DATA:
	    case GENERIC_32:
	    case GENERIC_64:
	    case SECONDARY_DATA:
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7: 
	      links[i] = i;
	      break;
	    case AA_DATA:	  
	      links[i] = firstAA;
	      break;
	    default:
	      assert(0);
	    }
	}
    }
  

  ll = initLinkageList(links, tr);

  free(links);
  
  return ll;
}



static void freeLinkageList( linkageList* ll)
{
  int i;    

  for(i = 0; i < ll->entries; i++)    
    free(ll->ld[i].partitionList);         

  free(ll->ld);
  free(ll);   
}

#define ALPHA_F 0
#define RATE_F  1
#define FREQ_F  2

static void changeModelParameters(int index, int rateNumber, double value, int whichParameterType, tree *tr)
{
  switch(whichParameterType)
    {
    case RATE_F:
      setRateModel(tr, index, value, rateNumber);  
      initReversibleGTR(tr, index);		 
      break;
    case ALPHA_F:
      tr->partitionData[index].alpha = value;
      makeGammaCats(tr->partitionData[index].alpha, tr->partitionData[index].gammaRates, 4, tr->useMedian);
      break;
    case FREQ_F:
      {
	int
	  states = tr->partitionData[index].states,
	  j;

	double 
	  w = 0.0;

	tr->partitionData[index].freqExponents[rateNumber] = value;

	for(j = 0; j < states; j++)
	  w += exp(tr->partitionData[index].freqExponents[j]);

	for(j = 0; j < states; j++)	    	    
	  tr->partitionData[index].frequencies[j] = exp(tr->partitionData[index].freqExponents[j]) / w;
	
	initReversibleGTR(tr, index);
      }
      break;
    default:
      assert(0);
    }
}

static void evaluateChange(tree *tr, int rateNumber, double *value, double *result, boolean* converged, int whichFunction, int numberOfModels, linkageList *ll, double modelEpsilon)
{ 
  int 
    i, 
    k, 
    pos;

   
  for(i = 0, pos = 0; i < ll->entries; i++)
    {
      if(ll->ld[i].valid)
	{
	  if(converged[pos])
	    {
	      //if parameter optimizations for this specific model have converged 
	      //set executeModel to FALSE 

	      for(k = 0; k < ll->ld[i].partitions; k++)
		tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;
	    }
	  else
	    {
	      for(k = 0; k < ll->ld[i].partitions; k++)
		{
		  int 
		    index = ll->ld[i].partitionList[k];

		  changeModelParameters(index, rateNumber, value[pos], whichFunction, tr);		    		  
		}
	    }
	  pos++;
	}
      else
	{
	  // if this partition is not being optimized anyway (e.g., we may be optimizing GTR rates for all DNA partitions,
	  // but there are also a couple of Protein partitions with fixed models like WAG, JTT, etc.) set executeModel to FALSE
	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;	     
	}      
    }

  assert(pos == numberOfModels);


  //LIBRARY: need to switch over parallel regions here either call 
  //the one for the rates or for alpha!

  evaluateGeneric(tr, tr->start, TRUE);  
               
  for(i = 0, pos = 0; i < ll->entries; i++)	
    {
      if(ll->ld[i].valid)
	{
	  result[pos] = 0.0;
	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    {
	      int 
		index = ll->ld[i].partitionList[k];

	      assert(tr->perPartitionLH[index] <= 0.0);
	      
	      result[pos] -= tr->perPartitionLH[index];
	      
	    }
	  pos++;
	}

      //set execute model for ALL partitions to true again 
      //for consistency 

      for(k = 0; k < ll->ld[i].partitions; k++)
	{
	  int 
	    index = ll->ld[i].partitionList[k];	  
	  tr->executeModel[index] = TRUE;
	}	  
    }
  
  assert(pos == numberOfModels);   
}



static void brentGeneric(double *ax, double *bx, double *cx, double *fb, double tol, double *xmin, double *result, int numberOfModels, 
			 int whichFunction, int rateNumber, tree *tr, linkageList *ll, double *lim_inf, double *lim_sup)
{
  int iter, i;
  double 
    *a     = (double *)malloc(sizeof(double) * numberOfModels),
    *b     = (double *)malloc(sizeof(double) * numberOfModels),
    *d     = (double *)malloc(sizeof(double) * numberOfModels),
    *etemp = (double *)malloc(sizeof(double) * numberOfModels),
    *fu    = (double *)malloc(sizeof(double) * numberOfModels),
    *fv    = (double *)malloc(sizeof(double) * numberOfModels),
    *fw    = (double *)malloc(sizeof(double) * numberOfModels),
    *fx    = (double *)malloc(sizeof(double) * numberOfModels),
    *p     = (double *)malloc(sizeof(double) * numberOfModels),
    *q     = (double *)malloc(sizeof(double) * numberOfModels),
    *r     = (double *)malloc(sizeof(double) * numberOfModels),
    *tol1  = (double *)malloc(sizeof(double) * numberOfModels),
    *tol2  = (double *)malloc(sizeof(double) * numberOfModels),
    *u     = (double *)malloc(sizeof(double) * numberOfModels),
    *v     = (double *)malloc(sizeof(double) * numberOfModels),
    *w     = (double *)malloc(sizeof(double) * numberOfModels),
    *x     = (double *)malloc(sizeof(double) * numberOfModels),
    *xm    = (double *)malloc(sizeof(double) * numberOfModels),
    *e     = (double *)malloc(sizeof(double) * numberOfModels);
  boolean *converged = (boolean *)malloc(sizeof(boolean) * numberOfModels);
  boolean allConverged;
  
  for(i = 0; i < numberOfModels; i++)    
    converged[i] = FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      e[i] = 0.0;
      d[i] = 0.0;
    }

  for(i = 0; i < numberOfModels; i++)
    {
      a[i]=((ax[i] < cx[i]) ? ax[i] : cx[i]);
      b[i]=((ax[i] > cx[i]) ? ax[i] : cx[i]);
      x[i] = w[i] = v[i] = bx[i];
      fw[i] = fv[i] = fx[i] = fb[i];
    }

  for(i = 0; i < numberOfModels; i++)
    {      
      assert(a[i] >= lim_inf[i] && a[i] <= lim_sup[i]);
      assert(b[i] >= lim_inf[i] && b[i] <= lim_sup[i]);
      assert(x[i] >= lim_inf[i] && x[i] <= lim_sup[i]);
      assert(v[i] >= lim_inf[i] && v[i] <= lim_sup[i]);
      assert(w[i] >= lim_inf[i] && w[i] <= lim_sup[i]);
    }
  
  

  for(iter = 1; iter <= ITMAX; iter++)
    {
      allConverged = TRUE;

      for(i = 0; i < numberOfModels && allConverged; i++)
	allConverged = allConverged && converged[i];

      if(allConverged)
	{
	  free(converged);
	  free(a);
	  free(b);
	  free(d);
	  free(etemp);
	  free(fu);
	  free(fv);
	  free(fw);
	  free(fx);
	  free(p);
	  free(q);
	  free(r);
	  free(tol1);
	  free(tol2);
	  free(u);
	  free(v);
	  free(w);
	  free(x);
	  free(xm);
	  free(e);
	  return;
	}     

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {	     	      
	      assert(a[i] >= lim_inf[i] && a[i] <= lim_sup[i]);
	      assert(b[i] >= lim_inf[i] && b[i] <= lim_sup[i]);
	      assert(x[i] >= lim_inf[i] && x[i] <= lim_sup[i]);
	      assert(v[i] >= lim_inf[i] && v[i] <= lim_sup[i]);
	      assert(w[i] >= lim_inf[i] && w[i] <= lim_sup[i]);
  
	      xm[i] = 0.5 * (a[i] + b[i]);
	      tol2[i] = 2.0 * (tol1[i] = tol * fabs(x[i]) + BRENT_ZEPS);
	  
	      if(fabs(x[i] - xm[i]) <= (tol2[i] - 0.5 * (b[i] - a[i])))
		{		 
		  result[i] =  -fx[i];
		  xmin[i]   = x[i];
		  converged[i] = TRUE;		  
		}
	      else
		{
		  if(fabs(e[i]) > tol1[i])
		    {		     
		      r[i] = (x[i] - w[i]) * (fx[i] - fv[i]);
		      q[i] = (x[i] - v[i]) * (fx[i] - fw[i]);
		      p[i] = (x[i] - v[i]) * q[i] - (x[i] - w[i]) * r[i];
		      q[i] = 2.0 * (q[i] - r[i]);
		      if(q[i] > 0.0)
			p[i] = -p[i];
		      q[i] = fabs(q[i]);
		      etemp[i] = e[i];
		      e[i] = d[i];
		      if((fabs(p[i]) >= fabs(0.5 * q[i] * etemp[i])) || (p[i] <= q[i] * (a[i]-x[i])) || (p[i] >= q[i] * (b[i] - x[i])))
			d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i] : b[i] - x[i]));
		      else
			{
			  d[i] = p[i] / q[i];
			  u[i] = x[i] + d[i];
			  if( u[i] - a[i] < tol2[i] || b[i] - u[i] < tol2[i])
			    d[i] = SIGN(tol1[i], xm[i] - x[i]);
			}
		    }
		  else
		    {		     
		      d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i]: b[i] - x[i]));
		    }
		  u[i] = ((fabs(d[i]) >= tol1[i]) ? (x[i] + d[i]): (x[i] +SIGN(tol1[i], d[i])));
		}

	      if(!converged[i])
		assert(u[i] >= lim_inf[i] && u[i] <= lim_sup[i]);
	    }
	}
                 
      evaluateChange(tr, rateNumber, u, fu, converged, whichFunction, numberOfModels, ll, tol);

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {
	      if(fu[i] <= fx[i])
		{
		  if(u[i] >= x[i])
		    a[i] = x[i];
		  else
		    b[i] = x[i];
		  
		  SHFT(v[i],w[i],x[i],u[i]);
		  SHFT(fv[i],fw[i],fx[i],fu[i]);
		}
	      else
		{
		  if(u[i] < x[i])
		    a[i] = u[i];
		  else
		    b[i] = u[i];
		  
		  if(fu[i] <= fw[i] || w[i] == x[i])
		    {
		      v[i] = w[i];
		      w[i] = u[i];
		      fv[i] = fw[i];
		      fw[i] = fu[i];
		    }
		  else
		    {
		      if(fu[i] <= fv[i] || v[i] == x[i] || v[i] == w[i])
			{
			  v[i] = u[i];
			  fv[i] = fu[i];
			}
		    }	    
		}
	      
	      assert(a[i] >= lim_inf[i] && a[i] <= lim_sup[i]);
	      assert(b[i] >= lim_inf[i] && b[i] <= lim_sup[i]);
	      assert(x[i] >= lim_inf[i] && x[i] <= lim_sup[i]);
	      assert(v[i] >= lim_inf[i] && v[i] <= lim_sup[i]);
	      assert(w[i] >= lim_inf[i] && w[i] <= lim_sup[i]);
	      assert(u[i] >= lim_inf[i] && u[i] <= lim_sup[i]);
	    }
	}
    }

  free(converged);
  free(a);
  free(b);
  free(d);
  free(etemp);
  free(fu);
  free(fv);
  free(fw);
  free(fx);
  free(p);
  free(q);
  free(r);
  free(tol1);
  free(tol2);
  free(u);
  free(v);
  free(w);
  free(x);
  free(xm);
  free(e);

  printf("\n. Too many iterations in BRENT !");
  assert(0);
}



static int brakGeneric(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, 
		       double *fc, double *lim_inf, double *lim_sup, 
		       int numberOfModels, int rateNumber, int whichFunction, tree *tr, linkageList *ll, double modelEpsilon)
{
  double 
    *ulim = (double *)malloc(sizeof(double) * numberOfModels),
    *u    = (double *)malloc(sizeof(double) * numberOfModels),
    *r    = (double *)malloc(sizeof(double) * numberOfModels),
    *q    = (double *)malloc(sizeof(double) * numberOfModels),
    *fu   = (double *)malloc(sizeof(double) * numberOfModels),
    *dum  = (double *)malloc(sizeof(double) * numberOfModels), 
    *temp = (double *)malloc(sizeof(double) * numberOfModels);
  
  int 
    i,
    *state    = (int *)malloc(sizeof(int) * numberOfModels),
    *endState = (int *)malloc(sizeof(int) * numberOfModels);

  boolean *converged = (boolean *)malloc(sizeof(boolean) * numberOfModels);
  boolean allConverged;

  for(i = 0; i < numberOfModels; i++)
    converged[i] = FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      state[i] = 0;
      endState[i] = 0;

      u[i] = 0.0;

      param[i] = ax[i];

      if(param[i] > lim_sup[i]) 	
	param[i] = ax[i] = lim_sup[i];
      
      if(param[i] < lim_inf[i]) 
	param[i] = ax[i] = lim_inf[i];

      assert(param[i] >= lim_inf[i] && param[i] <= lim_sup[i]);
    }
   
  
  evaluateChange(tr, rateNumber, param, fa, converged, whichFunction, numberOfModels, ll, modelEpsilon);


  for(i = 0; i < numberOfModels; i++)
    {
      param[i] = bx[i];
      if(param[i] > lim_sup[i]) 
	param[i] = bx[i] = lim_sup[i];
      if(param[i] < lim_inf[i]) 
	param[i] = bx[i] = lim_inf[i];

      assert(param[i] >= lim_inf[i] && param[i] <= lim_sup[i]);
    }
  
  evaluateChange(tr, rateNumber, param, fb, converged, whichFunction, numberOfModels, ll, modelEpsilon);

  for(i = 0; i < numberOfModels; i++)  
    {
      if (fb[i] > fa[i]) 
	{	  
	  SHFT(dum[i],ax[i],bx[i],dum[i]);
	  SHFT(dum[i],fa[i],fb[i],dum[i]);
	}
      
      cx[i] = bx[i] + MNBRAK_GOLD * (bx[i] - ax[i]);
      
      param[i] = cx[i];
      
      if(param[i] > lim_sup[i]) 
	param[i] = cx[i] = lim_sup[i];
      if(param[i] < lim_inf[i]) 
	param[i] = cx[i] = lim_inf[i];

      assert(param[i] >= lim_inf[i] && param[i] <= lim_sup[i]);
    }
  
 
  evaluateChange(tr, rateNumber, param, fc, converged, whichFunction, numberOfModels,  ll, modelEpsilon);

  while(1) 
    {       
      allConverged = TRUE;

      for(i = 0; i < numberOfModels && allConverged; i++)
	allConverged = allConverged && converged[i];

      if(allConverged)
	{
	  for(i = 0; i < numberOfModels; i++)
	    {	       
	      if(ax[i] > lim_sup[i]) 
		ax[i] = lim_sup[i];
	      if(ax[i] < lim_inf[i]) 
		ax[i] = lim_inf[i];

	      if(bx[i] > lim_sup[i]) 
		bx[i] = lim_sup[i];
	      if(bx[i] < lim_inf[i]) 
		bx[i] = lim_inf[i];
	       
	      if(cx[i] > lim_sup[i]) 
		cx[i] = lim_sup[i];
	      if(cx[i] < lim_inf[i]) 
		cx[i] = lim_inf[i];
	    }

	  free(converged);
	  free(ulim);
	  free(u);
	  free(r);
	  free(q);
	  free(fu);
	  free(dum); 
	  free(temp);
	  free(state);   
	  free(endState);
	  return 0;
	   
	}

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {
	      switch(state[i])
		{
		case 0:
		  endState[i] = 0;
		  if(!(fb[i] > fc[i]))		         
		    converged[i] = TRUE;		       		     
		  else
		    {
		   
		      if(ax[i] > lim_sup[i]) 
			ax[i] = lim_sup[i];
		      if(ax[i] < lim_inf[i]) 
			ax[i] = lim_inf[i];
		      if(bx[i] > lim_sup[i]) 
			bx[i] = lim_sup[i];
		      if(bx[i] < lim_inf[i]) 
			bx[i] = lim_inf[i];
		      if(cx[i] > lim_sup[i]) 
			cx[i] = lim_sup[i];
		      if(cx[i] < lim_inf[i]) 
			cx[i] = lim_inf[i];
		       
		      r[i]=(bx[i]-ax[i])*(fb[i]-fc[i]);
		      q[i]=(bx[i]-cx[i])*(fb[i]-fa[i]);
		      u[i]=(bx[i])-((bx[i]-cx[i])*q[i]-(bx[i]-ax[i])*r[i])/
			(2.0*SIGN(MAX(fabs(q[i]-r[i]),MNBRAK_TINY),q[i]-r[i]));
		       
		      ulim[i]=(bx[i])+MNBRAK_GLIMIT*(cx[i]-bx[i]);
		       
		      if(u[i] > lim_sup[i]) 
			u[i] = lim_sup[i];
		      if(u[i] < lim_inf[i]) 
			u[i] = lim_inf[i];
		      if(ulim[i] > lim_sup[i]) 
			ulim[i] = lim_sup[i];
		      if(ulim[i] < lim_inf[i]) 
			ulim[i] = lim_inf[i];
		       
		      if ((bx[i]-u[i])*(u[i]-cx[i]) > 0.0)
			{
			  param[i] = u[i];
			  if(param[i] > lim_sup[i]) 			     
			    param[i] = u[i] = lim_sup[i];
			  if(param[i] < lim_inf[i])
			    param[i] = u[i] = lim_inf[i];
			  endState[i] = 1;
			}
		      else 
			{
			  if ((cx[i]-u[i])*(u[i]-ulim[i]) > 0.0) 
			    {
			      param[i] = u[i];
			      if(param[i] > lim_sup[i]) 
				param[i] = u[i] = lim_sup[i];
			      if(param[i] < lim_inf[i]) 
				param[i] = u[i] = lim_inf[i];
			      endState[i] = 2;
			    }		  	       
			  else
			    {
			      if ((u[i]-ulim[i])*(ulim[i]-cx[i]) >= 0.0) 
				{
				  u[i] = ulim[i];
				  param[i] = u[i];	
				  if(param[i] > lim_sup[i]) 
				    param[i] = u[i] = ulim[i] = lim_sup[i];
				  if(param[i] < lim_inf[i]) 
				    param[i] = u[i] = ulim[i] = lim_inf[i];
				  endState[i] = 0;
				}		  		
			      else 
				{		  
				  u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
				  param[i] = u[i];
				  endState[i] = 0;
				  if(param[i] > lim_sup[i]) 
				    param[i] = u[i] = lim_sup[i];
				  if(param[i] < lim_inf[i]) 
				    param[i] = u[i] = lim_inf[i];
				}
			    }	  
			}
		    }
		  break;
		case 1:
		  endState[i] = 0;
		  break;
		case 2:
		  endState[i] = 3;
		  break;
		default:
		  assert(0);
		}
	      assert(param[i] >= lim_inf[i] && param[i] <= lim_sup[i]);
	    }
	}
             
      evaluateChange(tr, rateNumber, param, temp, converged, whichFunction, numberOfModels, ll, modelEpsilon);

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {	       
	      switch(endState[i])
		{
		case 0:
		  fu[i] = temp[i];
		  SHFT(ax[i],bx[i],cx[i],u[i]);
		  SHFT(fa[i],fb[i],fc[i],fu[i]);
		  state[i] = 0;
		  break;
		case 1:
		  fu[i] = temp[i];
		  if (fu[i] < fc[i]) 
		    {
		      ax[i]=(bx[i]);
		      bx[i]=u[i];
		      fa[i]=(fb[i]);
		      fb[i]=fu[i]; 
		      converged[i] = TRUE;		      
		    } 
		  else 
		    {
		      if (fu[i] > fb[i]) 
			{
			  assert(u[i] >= lim_inf[i] && u[i] <= lim_sup[i]);
			  cx[i]=u[i];
			  fc[i]=fu[i];
			  converged[i] = TRUE;			  
			}
		      else
			{		   
			  u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
			  param[i] = u[i];
			  if(param[i] > lim_sup[i]) {param[i] = u[i] = lim_sup[i];}
			  if(param[i] < lim_inf[i]) {param[i] = u[i] = lim_inf[i];}	  
			  state[i] = 1;		 
			}		  
		    }
		  break;
		case 2: 
		  fu[i] = temp[i];
		  if (fu[i] < fc[i]) 
		    {		     
		      SHFT(bx[i],cx[i],u[i], cx[i]+MNBRAK_GOLD*(cx[i]-bx[i]));
		      state[i] = 2;
		    }	   
		  else
		    {
		      state[i] = 0;
		      SHFT(ax[i],bx[i],cx[i],u[i]);
		      SHFT(fa[i],fb[i],fc[i],fu[i]);
		    }
		  break;	   
		case 3:		  
		  SHFT(fb[i],fc[i],fu[i], temp[i]);
		  SHFT(ax[i],bx[i],cx[i],u[i]);
		  SHFT(fa[i],fb[i],fc[i],fu[i]);
		  state[i] = 0;
		  break;
		default:
		  assert(0);
		}
	    }
	}
    }
   

  assert(0);
  free(converged);
  free(ulim);
  free(u);
  free(r);
  free(q);
  free(fu);
  free(dum); 
  free(temp);
  free(state);   
  free(endState);

  

  return(0);
}


/**********************************************************************************************************/
/* ALPHA PARAM ********************************************************************************************/



//this function is required for implementing the LG4X model later-on 

static void optAlphasGeneric(tree *tr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    non_LG4X_Partitions = 0,
    LG4X_Partitions  = 0;

  /* assumes homogeneous super-partitions, that either contain DNA or AA partitions !*/
  /* does not check whether AA are all linked */

  /* first do non-LG4X partitions */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case DNA_DATA:			  	
	case BINARY_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = TRUE;
	  non_LG4X_Partitions++;
	  break;
	case AA_DATA:	  
	  //to be implemented later-on 
	  /*if(tr->partitionData[ll->ld[i].partitionList[0]].protModels == LG4X)
	    {
	    LG4X_Partitions++;	      
	    ll->ld[i].valid = FALSE;
	    }
	    else*/
	  {
	    ll->ld[i].valid = TRUE;
	    non_LG4X_Partitions++;
	  }
	  break;
	default:
	  assert(0);
	}      
    }   

 

  if(non_LG4X_Partitions > 0)    
    optParamGeneric(tr, modelEpsilon, ll, non_LG4X_Partitions, -1, ALPHA_MIN, ALPHA_MAX, ALPHA_F);
  
  //right now this assertion shouldn't fail, undo when implementing LG4X  
  assert(LG4X_Partitions == 0);
 

  /* then LG4x partitions */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case DNA_DATA:			  	
	case BINARY_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;	  
	  break;
	case AA_DATA:	  
	  //deal with this later-on
	  /*if(tr->partitionData[ll->ld[i].partitionList[0]].protModels == LG4X)	      
	    ll->ld[i].valid = TRUE;	   
	    else*/
	  ll->ld[i].valid = FALSE;	   	    
	  break;
	default:
	  assert(0);
	}      
    }   
  
  //if(LG4X_Partitions > 0)
  //  optLG4X(tr, modelEpsilon, ll, LG4X_Partitions);

  for(i = 0; i < ll->entries; i++)
    ll->ld[i].valid = TRUE;
}


static double minFreq(int index, int whichFreq, tree *tr, double absoluteMin)
{
  double 
    min = 0.0,
    *w = tr->partitionData[index].freqExponents,
    c = 0.0;

  int
    states = tr->partitionData[index].states,
    i;

  for(i = 0; i < states; i++)
    if(i != whichFreq)
      c += exp(w[i]);

  min = log(FREQ_MIN) + log(c) - log (1.0 - FREQ_MIN);

  if(0)
    {
      double
	check = exp(min) / (exp(min) + c);
      
      printf("check %f\n", check);    

      printf("min: %f \n", min);
    }
  
  return MAX(min, absoluteMin);
}

static double maxFreq(int index, int whichFreq, tree *tr, double absoluteMax)
{
  double 
    max = 0.0,
    *w = tr->partitionData[index].freqExponents,
    c = 0.0;

  int
    states = tr->partitionData[index].states,
    i;

  for(i = 0; i < states; i++)
    if(i != whichFreq)
      c += exp(w[i]);

  max = log(1.0 - ((double)(states - 1) * FREQ_MIN)) + log(c) - log ((double)(states - 1) * FREQ_MIN);

  if(0)
    {
      double
	check = exp(max) / (exp(max) + c);
      
      printf("check max %f\n", check);    
      
      printf("max: %f \n", max);
    }
  
  return MIN(max, absoluteMax);
}


static void optParamGeneric(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int rateNumber, double _lim_inf, double _lim_sup, int whichParameterType)
{
  int
    l,
    k, 
    j, 
    pos;
    
  double 
    *startValues = (double *)malloc(sizeof(double) * numberOfModels),
    *startLH    = (double *)malloc(sizeof(double) * numberOfModels),
    *endLH      = (double *)malloc(sizeof(double) * numberOfModels),
    *_a         = (double *)malloc(sizeof(double) * numberOfModels),
    *_b         = (double *)malloc(sizeof(double) * numberOfModels),
    *_c         = (double *)malloc(sizeof(double) * numberOfModels),
    *_fa        = (double *)malloc(sizeof(double) * numberOfModels),
    *_fb        = (double *)malloc(sizeof(double) * numberOfModels),
    *_fc        = (double *)malloc(sizeof(double) * numberOfModels),
    *_param     = (double *)malloc(sizeof(double) * numberOfModels),    
    *_x         = (double *)malloc(sizeof(double) * numberOfModels),
    *lim_inf    = (double *)malloc(sizeof(double) * numberOfModels),
    *lim_sup    = (double *)malloc(sizeof(double) * numberOfModels);

  evaluateGeneric(tr, tr->start, TRUE);
  
#ifdef  _DEBUG_MOD_OPT
  double
    initialLH = tr->likelihood;
#endif

  /* 
     at this point here every worker has the traversal data it needs for the 
     search 
  */

  for(l = 0, pos = 0; l < ll->entries; l++)
    {
      if(ll->ld[l].valid)
	{
	  endLH[pos] = unlikely;
	  startLH[pos] = 0.0;

	  for(j = 0; j < ll->ld[l].partitions; j++)
	    {
	      int 
		index = ll->ld[l].partitionList[j];
	      
	      startLH[pos] += tr->perPartitionLH[index];
	      
	      switch(whichParameterType)
		{
		case ALPHA_F:
		  lim_inf[pos] = _lim_inf;
		  lim_sup[pos] = _lim_sup;
		  startValues[pos] = tr->partitionData[index].alpha;
		  break;
		case RATE_F:
		  lim_inf[pos] = _lim_inf;
		  lim_sup[pos] = _lim_sup;
		  startValues[pos] = tr->partitionData[index].substRates[rateNumber];      
		  break;
		case FREQ_F:
		  lim_inf[pos] = minFreq(index, rateNumber, tr, _lim_inf);
		  lim_sup[pos] = maxFreq(index, rateNumber, tr, _lim_sup);
		  startValues[pos] = tr->partitionData[index].freqExponents[rateNumber];
		  break;
		default:
		  assert(0);
		}
		  
	    }
	  pos++;
	}
    }  

  assert(pos == numberOfModels);
   
  for(k = 0, pos = 0; k < ll->entries; k++)
    {
      if(ll->ld[k].valid)
	{	 	 	  
	  _a[pos] = startValues[pos] + 0.1;
	  _b[pos] = startValues[pos] - 0.1;
	      
	  if(_a[pos] < lim_inf[pos]) 
	    _a[pos] = lim_inf[pos];
	  
	  if(_a[pos] > lim_sup[pos]) 
	    _a[pos] = lim_sup[pos];
	      
	  if(_b[pos] < lim_inf[pos]) 
	    _b[pos] = lim_inf[pos];
	  
	  if(_b[pos] > lim_sup[pos]) 
	    _b[pos] = lim_sup[pos];    

	  pos++;
	}
    }                    	     

  assert(pos == numberOfModels);

  brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, rateNumber, whichParameterType, tr, ll, modelEpsilon);
      
  for(k = 0; k < numberOfModels; k++)
    {
      assert(_a[k] >= lim_inf[k] && _a[k] <= lim_sup[k]);
      assert(_b[k] >= lim_inf[k] && _b[k] <= lim_sup[k]);	  
      assert(_c[k] >= lim_inf[k] && _c[k] <= lim_sup[k]);	    
    }      

  brentGeneric(_a, _b, _c, _fb, modelEpsilon, _x, endLH, numberOfModels, whichParameterType, rateNumber, tr,  ll, lim_inf, lim_sup);
		      
  for(k = 0, pos = 0; k < ll->entries; k++)
    {
      if(ll->ld[k].valid)
	{ 
	  if(startLH[pos] > endLH[pos])
	    {
	      //if the initial likelihood was better than the likelihodo after optimization, we set the values back 
	      //to their original values 

	      for(j = 0; j < ll->ld[k].partitions; j++)
		{
		  int 
		    index = ll->ld[k].partitionList[j];
		  
		  changeModelParameters(index, rateNumber, startValues[pos], whichParameterType, tr);		 
		}
	    }
	  else
	    {
	      //otherwise we set the value to the optimized value 
	      //this used to be a bug in standard RAxML, before I fixed it 
	      //I was not using _x[pos] as value that needs to be set 

	      for(j = 0; j < ll->ld[k].partitions; j++)
		{
		  int 
		    index = ll->ld[k].partitionList[j];

		  changeModelParameters(index, rateNumber, _x[pos], whichParameterType, tr);		  		  
		}
	    }
	  pos++;
	}
    }


  //LIBRARY call the barrier here in the LIBRARY to update model params at all threads/processes !
    
  assert(pos == numberOfModels);

  free(startLH);
  free(endLH);
  free(_a);
  free(_b);
  free(_c);
  free(_fa);
  free(_fb);
  free(_fc);
  free(_param);
  free(_x);  
  free(startValues);
  free(lim_inf);
  free(lim_sup);

#ifdef _DEBUG_MOD_OPT
  evaluateGeneric(tr, tr->start, TRUE);

  if(tr->likelihood < initialLH)
    printf("%f %f\n", tr->likelihood, initialLH);
  assert(tr->likelihood >= initialLH);
#endif

}



//******************** rate optimization functions ***************************************************/

static void optFreqs(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{ 
  int 
    rateNumber;

  double
    freqMin = -1000000.0,
    freqMax = 200.0;
  
  for(rateNumber = 0; rateNumber < states; rateNumber++)
    optParamGeneric(tr, modelEpsilon, ll, numberOfModels, rateNumber, freqMin, freqMax, FREQ_F);   
}

static void optBaseFreqs(tree *tr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    states,
    dnaPartitions = 0,
    aaPartitions  = 0;

  /* first do DNA */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case DNA_DATA:	
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;	 
	  if(tr->partitionData[ll->ld[i].partitionList[0]].optimizeBaseFrequencies)
	    {
	      ll->ld[i].valid = TRUE;
	      dnaPartitions++;  	    
	    }
	  else
	    ll->ld[i].valid = FALSE;
	  break;       
	case AA_DATA:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }   

  if(dnaPartitions > 0)
    optFreqs(tr, modelEpsilon, ll, dnaPartitions, states);
  
  /* then AA */

  
  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case AA_DATA:	  
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states; 	      
	  if(tr->partitionData[ll->ld[i].partitionList[0]].optimizeBaseFrequencies)
	    {
	      ll->ld[i].valid = TRUE;
	      aaPartitions++;		
	    }
	  else
	    ll->ld[i].valid = FALSE; 
	  break;
	case DNA_DATA:	    
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}	 
    }
  
  if(aaPartitions > 0)      
    optFreqs(tr, modelEpsilon, ll, aaPartitions, states);

  for(i = 0; i < ll->entries; i++)
    ll->ld[i].valid = TRUE;
}


//new version for optimizing rates, an external loop that iterates over the rates 

static void optRates(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{
  int
    rateNumber,
    numberOfRates = ((states * states - states) / 2) - 1;

  for(rateNumber = 0; rateNumber < numberOfRates; rateNumber++)
    optParamGeneric(tr, modelEpsilon, ll, numberOfModels, rateNumber, RATE_MIN, RATE_MAX, RATE_F);   
}


static boolean AAisGTR(tree *tr)
{
  int i, count = 0;

  for(i = 0; i < tr->NumberOfModels; i++)   
    {
      if(tr->partitionData[i].dataType == AA_DATA)
	{
	  count++;
	  if(tr->partitionData[i].protModels != GTR)
	    return FALSE;
	}
    }

  if(count == 0)
    return FALSE;

  return TRUE;
}

static void optRatesGeneric(tree *tr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    dnaPartitions = 0,
    aaPartitions  = 0,
    states = -1;

  /* assumes homogeneous super-partitions, that either contain DNA or AA partitions !*/
  /* does not check whether AA are all linked */

  /* first do DNA */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case DNA_DATA:	
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;	 
	  ll->ld[i].valid = TRUE;
	  dnaPartitions++;  	    
	  break;
	case BINARY_DATA:
	case AA_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }   

  if(dnaPartitions > 0)
    optRates(tr, modelEpsilon, ll, dnaPartitions, states);
  
  /* then AA for GTR */

  if(AAisGTR(tr))
    {
      for(i = 0; i < ll->entries; i++)
	{
	  switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	    {
	    case AA_DATA:
	      states = tr->partitionData[ll->ld[i].partitionList[0]].states; 	      
	      ll->ld[i].valid = TRUE;
	      aaPartitions++;		
	      break;
	    case DNA_DATA:	    
	    case BINARY_DATA:
	    case SECONDARY_DATA:	
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7:
	      ll->ld[i].valid = FALSE;
	      break;
	    default:
	      assert(0);
	    }	 
	}

      assert(aaPartitions == 1);     
      
      optRates(tr, modelEpsilon, ll, aaPartitions, states);
    }  

  for(i = 0; i < ll->entries; i++)
    ll->ld[i].valid = TRUE;
}





/*********************FUNCTIONS FOR APPROXIMATE MODEL OPTIMIZATION ***************************************/






static int catCompare(const void *p1, const void *p2)
{
  rateCategorize *rc1 = (rateCategorize *)p1;
  rateCategorize *rc2 = (rateCategorize *)p2;

  double i = rc1->accumulatedSiteLikelihood;
  double j = rc2->accumulatedSiteLikelihood;
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}


static void categorizePartition(tree *tr, rateCategorize *rc, int model, int lower, int upper, double *patrat, 
				int *rateCategory /* temporary; used to be tr->rateCategory */ 
				)
{
  

  int
    zeroCounter,
    i, 
    k;
  
  double 
    diff, 
    min;

  for (i = lower, zeroCounter = 0; i < upper; i++, zeroCounter++) 
    {
      double
	temp = patrat[i];

      int
	found = 0;
	
      for(k = 0; k < tr->partitionData[model].numberOfCategories; k++)
	{
	  if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
	    {
	      found = 1;
	      rateCategory[i] = k;
	      break;
	    }
	}
	
      if(!found)
	{
	  min = fabs(temp - rc[0].rate);
	  rateCategory[i] = 0;

	  for(k = 1; k < tr->partitionData[model].numberOfCategories; k++)
	    {
	      diff = fabs(temp - rc[k].rate);

	      if(diff < min)
		{
		  min = diff;
		  rateCategory[i] = k;
		}
	    }
	}
    }

  for(k = 0; k < tr->partitionData[model].numberOfCategories; k++)
    tr->partitionData[model].perSiteRates[k] = rc[k].rate; 
}




static void optRateCatPthreads(tree *tr, double lower_spacing, double upper_spacing)
{
  int 
    model;

  size_t
    i;

  for(model = 0; model < tr->NumberOfModels; model++)
    {   
      pInfo 
	*partition = &(tr->partitionData[model]); 
      
      for( i = 0; i < partition->width; ++i)
	{
	  double 
	    initialRate, 
	    initialLikelihood, 
	    leftLH, 
	    rightLH, 
	    leftRate, 
	    rightRate, 
	    v;
	      
	  const double 
	    epsilon = 0.00001;
	      
	  int 
	    k;	      

	  initialRate = partition->patrat[i]; 
	      
	  initialLikelihood = evaluatePartialGeneric(tr, i, initialRate, model); /* i is real i ??? */	      	      	      	      

	  leftLH = rightLH = initialLikelihood;
	  leftRate = rightRate = initialRate;
	      
	  k = 1;
	      
	  while((initialRate - k * lower_spacing > 0.0001) && 
		((v = evaluatePartialGeneric(tr, i, initialRate - k * lower_spacing, model)) 
		 > leftLH) && 
		(fabs(leftLH - v) > epsilon))  
	    {	  
#ifndef WIN32
	      if(isnan(v))
		assert(0);
#endif
		  
	      leftLH = v;
	      leftRate = initialRate - k * lower_spacing;
	      k++;	  
	    }      
	      
	  k = 1;
	      
	  while(((v = evaluatePartialGeneric(tr, i, initialRate + k * upper_spacing, model)) > rightLH) &&
		(fabs(rightLH - v) > epsilon))    	
	    {
#ifndef WIN32
	      if(isnan(v))
		assert(0);
#endif     
	      rightLH = v;
	      rightRate = initialRate + k * upper_spacing;	 
	      k++;
	    }           
	      
	  if(rightLH > initialLikelihood || leftLH > initialLikelihood)
	    {
	      if(rightLH > leftLH)	    
		{	     
		  partition->patrat[i] = rightRate;
		  partition->lhs[i]  = rightLH; 
		}
	      else
		{	      
		  partition->patrat[i] = leftRate; 
		  partition->lhs[i] = leftLH;
		}
	    }
	  else	    
	    partition->lhs[i] = initialLikelihood;	    
	}
    }
}





/** 
    determines the weighted rates for each partition. Intended for use
    with normalization of the CAT model rates.
    
    Since information about rates and weights is distributed (each
    process only has the respective info for the data assigned to it),
    we have to communicate with peer processes. Notice, that
    weightPerPart could actually be stored in a variable, since the
    result does not change...

    output:
    weightPerPart_result  -- the sum of weights per partition
    weightedRates_result  -- sum of rates per partition weighted by site weight

*/ 
static void getWeightsAndWeightedRates(const tree * const tr, int **weightPerPart_result, double **weightedRates_result )
{
  int 
    i,
    *weightPerPart = (int *)NULL;
  
  double 
    *weightedRates = (double *)NULL;
   
  *weightPerPart_result = (int*)calloc((size_t)tr->NumberOfModels, sizeof(int)); 
  *weightedRates_result = (double*) calloc((size_t)tr->NumberOfModels, sizeof(double)); 
  
  
  weightedRates = *weightedRates_result; 
  weightPerPart = *weightPerPart_result; 
  
  for(i = 0; i < tr->NumberOfModels; ++i)
    {
      size_t 
	j; 

      pInfo 
	*partition = &(tr->partitionData[i]); 
      
      for(j = 0; j < partition->width; ++j)
	{ 
	  int 
	    c = partition->rateCategory[j];
	  
	  weightPerPart[i] += partition->wgt[j];	  
	  assert(0 <= c && c < tr->maxCategories); 
	  weightedRates[i] += ((double)partition->wgt[j]) * partition->perSiteRates[c]; 
	}
    }
  
  MPI_Allreduce(MPI_IN_PLACE, weightPerPart, tr->NumberOfModels, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, weightedRates, tr->NumberOfModels, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

  for( i = 0; i < tr->NumberOfModels; ++i)
    {
      assert(weightPerPart[i] > 0 ); 
      assert(weightedRates[i] > 0.0 );
    }
}


/* 
   this used to be updatePerSiteRates without scaling. Previously,
   updatePerSiteRates without scaling only conducted a check about
   whether sites are scaled correctly.
 */

//Andre but isn't this checking that the rates have been scaled correctly?
//shouldn't the assertions fail in this case, i.e., without scaling ?
void checkPerSiteRates(const tree *const tr )
{
  int 
    i,
    *weightPerPart =  (int *)NULL; 
  
  double 
    *weightedRates = (double *)NULL; 
  
  /*
    determine the sum of weights (weightPerPart) and the sum of all
    rates of a partition weighted by site weights
   */ 
  getWeightsAndWeightedRates(tr, &weightPerPart, &weightedRates); 

  if(tr->numBranches > 1 )
    {
      /* check if the mean of rates of each partition is 1 */
      for(i = 0; i < tr->NumberOfModels; ++i)
	{
	  double accRat = weightedRates[i] / (double)weightPerPart[i]; 
	  assert(fabs(accRat - 1.0) < 1e-5); 
	}
    }
  else 
    {
      /* check, if the overall mean of rates is 1 */

      double 
	accRat = 0.0,
	accWgt = 0.0; 
      
      for(i = 0; i < tr->NumberOfModels; ++i)
	{
	  accRat += weightedRates[i]; 
	  accWgt += weightPerPart[i]; 
	}
      accRat  /= (double)accWgt; 
      
      assert(fabs(accRat - 1.0) < 1e-5); 
    }

  free(weightedRates); 
  free(weightPerPart); 
}


/** 
    updatePerSiteRates is called after the master has categorized
    rates into several categories and every process has obtained the
    categorization for only the data assigned to it. Now, we still
    have to scale the rates, s.t. they are 1 on average.

    Thus, some communication is still needed to determine the total
    weight and the weighted rates (because this information is
    destributed).

    Notice that this function previously had two modes (scaleRates =
    {TRUE,FALSE}). Previously, scaleRates = FALSE, only performed a
    check on whether rates are scaled correctly such that the average
    rate is 1. For clarity, this functionality is now in a separate
    function called checkPerSiteRates.
*/ 
static void updatePerSiteRates(tree *tr)
{
  int 
    i, 
    *weightPerPart =  (int *)NULL; 
  
  double 
    *weightedRates = (double *)NULL; 

  getWeightsAndWeightedRates(tr, &weightPerPart, &weightedRates); 

  if(tr->numBranches > 1  )	
    {
      /* scale each partition, s.t. average rate within the partition is 1 */
      for(i = 0; i < tr->NumberOfModels; ++i)
	{
	  int j; 
	  double scaler = weightedRates[i] / (double)weightPerPart[i]; 
	  scaler = 1.0 / scaler; 
	  for(j = 0; j < tr->partitionData[i].numberOfCategories; ++j)
	    tr->partitionData[i].perSiteRates[j] *= scaler; 
	}
    }
  else
    {
      /* scale, s.t. average rate is 1 */ 

      double 
	scaler = 0.0,
	accWgt = 0.0; 
      
      for(i = 0; i < tr->NumberOfModels; ++i)
	{
	  scaler += weightedRates[i]; 
	  accWgt += weightPerPart[i]; 
	}
      scaler /= (double)accWgt; 
      scaler = 1.0 / scaler; 

      for(i = 0; i < tr->NumberOfModels; ++i)
	{
	  pInfo 
	    *partition = &(tr->partitionData[i]); 
	  
	  int 
	    j; 
	  
	  for(j = 0; j < partition->numberOfCategories; ++j)
	    partition->perSiteRates[j] *= scaler; 
	}
    }

  free(weightedRates); 
  free(weightPerPart); 
  
  /* 
     finally check, whether the rates are scaled correctly, s.t. their
     mean is 1
   */ 
  checkPerSiteRates(tr);
}



/*
  gathers optimized rates and the associated persite-lnls from all
  processes at the master.

  Notice that for instance tr->patrat_basePtr already contain all rate
  data of a single process.

  Output: 
  optRates_result  -- (only at master) a pointer to an array of optimized rates (corresponds to what used to be  tr->patratStored)
  lnls_result -- (only at master) a pointer to an array with persite-lnls that correspond to the newly proposed rate  (used to be tr->lhs) 
 */ 
static void gatherOptimizedRates(tree *tr, double **optRates_result, double **lnls_result)
{
  /* determine counts and displacement for data for each processor  */
  int 
    *numPerProc = (int *)NULL, 
    *displPerProc = (int *)NULL; 
  
  calculateLengthAndDisplPerProcess(tr, &numPerProc, &displPerProc); 
  
  gatherDistributedArray( tr, (void**) optRates_result, tr->patrat_basePtr, MPI_DOUBLE, numPerProc, displPerProc); 
  gatherDistributedArray(tr , (void**) lnls_result, tr->lhs_basePtr, MPI_DOUBLE, numPerProc, displPerProc); 

  free(numPerProc);
  free(displPerProc); 
}


/* 
    The master creates rate categories and assigns the rate categories
   to processes. Only executed by the master to assure consistent
   categorization.

   This code used to be the first part of of optimizeRateCategories()
   and has only slightly been modified.

   Input: 
   patrat  -- a global array of  optimized rates (used to be tr->patratStored)  
   lnls -- a global array of per-site lnls (used to be tr->lhs)

   Output: 
   rateCategory_result  -- a pointer to a global array of rate categories (used to be tr->rateCategory) 
   
   side effect:
   tr->partitionData[i].perSiteRates gets computed in categorizePartition.

 */ 
static void categorizeTheRates(tree *tr, double *patrat, double *lnls, int maxCategories, int **rateCategory_result)
{
  int 
    model, i; 

  *rateCategory_result = (int*) calloc((size_t)tr->originalCrunchedLength, sizeof(int)); 
  
  for(model = 0; model < tr->NumberOfModels; model++)
    {    
      double 
	temp = 0.0; 
      
      int 
	where = 1,
	found = 0,
	width = tr->partitionData[model].upper -  tr->partitionData[model].lower,
	upper = tr->partitionData[model].upper,
	lower = tr->partitionData[model].lower;
	      
      rateCategorize 
	*rc = (rateCategorize *)malloc(sizeof(rateCategorize) * width);		 
	      
      for (i = 0; i < width; i++)
	{
	  rc[i].accumulatedSiteLikelihood = 0.0;
	  rc[i].rate = 0.0;
	}  
	      
      rc[0].accumulatedSiteLikelihood = lnls[lower];
      rc[0].rate = patrat[lower];

      for (i = lower + 1; i < upper; i++) 
	{
	  int k; 

	  temp = patrat[i];
	  found = 0;
		  
	  for(k = 0; k < where; k++)
	    {
	      if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
		{
		  found = 1;						
		  rc[k].accumulatedSiteLikelihood += lnls[i];	
		  break;
		}
	    }
		  
	  if(!found)
	    {	    
	      rc[where].rate = temp;	    
	      rc[where].accumulatedSiteLikelihood += lnls[i];	    
	      where++;
	    }
	}
	      
      qsort(rc, where, sizeof(rateCategorize), catCompare);
	      
      if(where < maxCategories)
	{
	  tr->partitionData[model].numberOfCategories = where;
	  categorizePartition(tr, rc, model, lower, upper, patrat, *rateCategory_result);
	}
      else
	{
	  tr->partitionData[model].numberOfCategories = maxCategories;	
	  categorizePartition(tr, rc, model, lower, upper, patrat, *rateCategory_result);
	}
	      
      free(rc);
    } 
}


/* #define PRINT_RAT_CAT */

/** 
    informs all peer processes about 
    * rateCategory
    * numberOfCategories
    * perSiteRates
    of their data
*/
static void scatterProcessedRates(tree *tr, int *rateCategory)
{
  int 
    i,
    *countPerProc = (int *)NULL, 
    *displPerProc = (int *)NULL,
    *numCatPerPart = (int*) calloc((size_t)tr->NumberOfModels, sizeof(int)); 
  
  if(processID == 0)
    {
      for(i = 0; i < tr->NumberOfModels; ++i)
	numCatPerPart[i] = tr->partitionData[i].numberOfCategories; 
    }
  MPI_Bcast(numCatPerPart, tr->NumberOfModels,  MPI_INT, 0,MPI_COMM_WORLD); 
  for(i = 0; i < tr->NumberOfModels; ++i)
    tr->partitionData[i].numberOfCategories = numCatPerPart[i]; 
  free(numCatPerPart); 
    
  /* for simplicity, broad cast all peSiteRates */
  for(i = 0; i < tr->NumberOfModels; ++i)
    MPI_Bcast(tr->partitionData[i].perSiteRates, tr->maxCategories, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  /* prepare for scattering */
  calculateLengthAndDisplPerProcess(tr,  &countPerProc, &displPerProc); 


#ifdef PRINT_RAT_CAT  
  if(processID == 0)
    {
      printf("rates BEFORE: "); 
      for(i = 0; i < tr->originalCrunchedLength; ++i)
	printf("%d,", rateCategory[i]); 
      printf("\n"); 
    }
#endif

  scatterDistrbutedArray(tr, rateCategory, tr->rateCategory_basePtr, MPI_INT, countPerProc, displPerProc); 

#ifdef PRINT_RAT_CAT
  int len = getMyCharacterLength(tr); 
  printf("basepointer AFTER: "); 
  for(i = 0; i < len ; ++i)
    printf("%d,", tr->rateCategory_basePtr[i]);
  printf("\n"); 
#endif
  
  free(countPerProc); 
  free(displPerProc); 
}



/* backup for one partition */
typedef struct 
{
  double *patrat; 
  int *rateCategory; 
  double *perSiteRates; 
  int numberOfCategories; 
} RateBackup; 


/** 
    This function creates a backup of all data relevant for CAT-rate
    assignment.
    
    Previously, the backup info has been stored in patratStored. Or
    maybe it is the otherway around and patrat was the backup, while
    patratStored contained the actual optimized rates.

    Output: 
    resultPtr -- contains the backup  
 */ 
static void backupRates(tree *tr, RateBackup** resultPtr)
{
  int 
    i,
    numCat = tr->maxCategories;
  
  RateBackup
    *backup;

  *resultPtr = (RateBackup* ) calloc((size_t)tr->NumberOfModels, sizeof(RateBackup));
  
  backup = *resultPtr; 
 
  for(i = 0; i < tr->NumberOfModels; ++i)
    {
      pInfo
	*partition = &(tr->partitionData[i]); 
      RateBackup 
	*bk = backup + i; 

      bk->patrat = (double*)calloc((size_t)partition->width, sizeof(double)); 
      bk->perSiteRates = (double*) calloc((size_t)numCat, sizeof(double)); 
      bk->rateCategory = (int*) calloc((size_t)partition->width, sizeof(int)) ;
      bk->numberOfCategories = partition->numberOfCategories; 

      memcpy(bk->patrat, partition->patrat, sizeof(double) * (size_t)partition->width); 
      memcpy(bk->perSiteRates, partition->perSiteRates, sizeof(double) * (size_t)numCat); 
      memcpy(bk->rateCategory, partition->rateCategory, sizeof(int) * (size_t)partition->width); 
    } 
}


static void restoreBackupRates(tree *tr , RateBackup *rb)
{
  int 
    numCat = tr->maxCategories,
    i; 
  
  for(i = 0; i < tr->NumberOfModels; ++i)
    {
      pInfo 
	*partition = &(tr->partitionData[i]); 
      
      RateBackup 
	*bk = rb + i; 

      partition->numberOfCategories = bk->numberOfCategories; 

      memcpy(partition->patrat, bk->patrat, sizeof(double) * (size_t)partition->width); 
      memcpy(partition->perSiteRates, bk->perSiteRates, sizeof(double) * (size_t)numCat); 
      memcpy(partition->rateCategory, bk->rateCategory, sizeof(int) * (size_t)partition->width); 
    }
}



static void deleteBackupRates(tree *tr, RateBackup** rbPtr)
{
  int i ;

  for(i = 0; i< tr->NumberOfModels; ++i)
    {
      RateBackup 
	*rb =  &((*rbPtr)[i]);
      
      free(rb->patrat);
      free(rb->perSiteRates);
      free(rb->rateCategory);
    }

  free(*rbPtr);
  rbPtr = (RateBackup **)NULL;
}


static void optimizeRateCategories(tree *tr, int _maxCategories)
{
  assert(_maxCategories > 0);  

  if(_maxCategories == 1)
    return; 


  double  
    lower_spacing, 
    upper_spacing,
    initialLH = tr->likelihood,
    *optRates = (double*)NULL,
    *lnls = (double*)NULL; 

  int
    *rateCategory = (int *)NULL,
    maxCategories = _maxCategories ; 

  RateBackup 
    *rateBackup = (RateBackup *)NULL; 

  assert(isTip(tr->start->number, tr->mxtips));         

  evaluateGeneric(tr, tr->start, TRUE);     

  if(optimizeRateCategoryInvocations == 1)
    {
      lower_spacing = 0.5 / ((double)optimizeRateCategoryInvocations);
      upper_spacing = 1.0 / ((double)optimizeRateCategoryInvocations);
    }
  else
    {
      lower_spacing = 0.05 / ((double)optimizeRateCategoryInvocations);
      upper_spacing = 0.1 / ((double)optimizeRateCategoryInvocations);
    }

  if(lower_spacing < 0.001)
    lower_spacing = 0.001;

  if(upper_spacing < 0.001)
    upper_spacing = 0.001;

  optimizeRateCategoryInvocations++;

  //store old rate category assignment 
  backupRates(tr, &rateBackup); 

  /* process specific: each process optimizes rates for data
     assigned to it */
  optRateCatPthreads(tr, lower_spacing, upper_spacing);
  
  /* gather rates and lnls at the master */
  gatherOptimizedRates(tr, &optRates, &lnls); 
  
  /* master has all necessary info now and can categorize the rates */
  if(processID == 0)
    {
      categorizeTheRates(tr, optRates, lnls, maxCategories, &rateCategory ); 
  
      /* only allocated at master  */
      free(optRates); 
      free(lnls); 
    }

  scatterProcessedRates(tr, rateCategory );
  if(processID == 0)
    free(rateCategory); 

  /* every process has now new rates and a new category
     assignment. However, we still have to scale the rates, such their
     weighted mean rate is 1.  */
  updatePerSiteRates(tr); 

  evaluateGeneric(tr, tr->start, TRUE);

  if(tr->likelihood < initialLH)
    {	 		  
      restoreBackupRates(tr, rateBackup); 

      //Andre I don't understand the comment below ... 
      //can per-site rate scaling still be dis-abled in this version of the code?

      /* 
	 => Andre: I am afraid neither do I. Comparing it to the
	 original code, I think everything should be fine: we restore
	 the previous state and check, whether rates are scaled
	 correctly.
      */
      
      /* cannot do that any more here  */
      checkPerSiteRates(tr); 
      
      evaluateGeneric(tr, tr->start, TRUE);	 

      assert(initialLH == tr->likelihood);
    }

  deleteBackupRates(tr,&rateBackup); 
}




/*****************************************************************************************************/

void resetBranches(tree *tr)
{
  nodeptr  p, q;
  int  nodes, i;

  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {   
      for(i = 0; i < tr->numBranches; i++)
	p->z[i] = defaultz;
	
      q = p->next;
      while(q != p)
	{	
	  for(i = 0; i < tr->numBranches; i++)
	    q->z[i] = defaultz;	    
	  q = q->next;
	}
      p++;
    }
}


static void printAAmatrix(tree *tr, double epsilon)
{
  if(AAisGTR(tr))
    {
      int model;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->partitionData[model].dataType == AA_DATA) 
	    {
	      char gtrFileName[1024];
	      char epsilonStr[1024];
	      FILE *gtrFile;
	      double *rates = tr->partitionData[model].substRates;
	      double *f     = tr->partitionData[model].frequencies;
	      double q[20][20];
	      int    r = 0;
	      int i, j;

	      assert(tr->partitionData[model].protModels == GTR);

	      sprintf(epsilonStr, "%f", epsilon);

	      strcpy(gtrFileName, workdir);
	      strcat(gtrFileName, "RAxML_proteinGTRmodel.");
	      strcat(gtrFileName, run_id);
	      strcat(gtrFileName, "_");
	      strcat(gtrFileName, epsilonStr);

	      gtrFile = myfopen(gtrFileName, "wb");

	      for(i = 0; i < 20; i++)
		for(j = 0; j < 20; j++)
		  q[i][j] = 0.0;

	      for(i = 0; i < 19; i++)
		for(j = i + 1; j < 20; j++)
		  q[i][j] = rates[r++];

	      for(i = 0; i < 20; i++)
		for(j = 0; j <= i; j++)
		  {
		    if(i == j)
		      q[i][j] = 0.0;
		    else
		      q[i][j] = q[j][i];
		  }
	   
	      for(i = 0; i < 20; i++)
		{
		  for(j = 0; j < 20; j++)		
		    fprintf(gtrFile, "%1.80f ", q[i][j]);
		
		  fprintf(gtrFile, "\n");
		}
	      for(i = 0; i < 20; i++)
		fprintf(gtrFile, "%1.80f ", f[i]);
	      fprintf(gtrFile, "\n");

	      fclose(gtrFile);

	      printBothOpen("\nPrinted intermediate AA substitution matrix to file %s\n\n", gtrFileName);
	      
	      break;
	    }

	}	  
    }
}


static void autoProtein(tree *tr)
{
  int 
    countAutos = 0,   
    model;  
  
  for(model = 0; model < tr->NumberOfModels; model++)	      
    if(tr->partitionData[model].protModels == AUTO)
      countAutos++;

  if(countAutos > 0)
    {
      int 
	i,
	numProteinModels = AUTO,
	*bestIndex = (int*)malloc(sizeof(int) * tr->NumberOfModels),
	*oldIndex  = (int*)malloc(sizeof(int) * tr->NumberOfModels);

      double
	startLH,
	*bestScores = (double*)malloc(sizeof(double) * tr->NumberOfModels);    

      topolRELL_LIST 
	*rl = (topolRELL_LIST *)malloc(sizeof(topolRELL_LIST));

      initTL(rl, tr, 1);
      saveTL(rl, tr, 0);

      evaluateGeneric(tr, tr->start, TRUE); 

      startLH = tr->likelihood;

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  oldIndex[model] = tr->partitionData[model].autoProtModels;
	  bestIndex[model] = -1;
	  bestScores[model] = unlikely;
	}
      
      for(i = 0; i < numProteinModels; i++)
	{
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {	   
	      if(tr->partitionData[model].protModels == AUTO)
		{
		  tr->partitionData[model].autoProtModels = i;
		  initReversibleGTR(tr, model);  
		}
	    }

	  //this barrier needs to be called in the library 
	  //#ifdef _USE_PTHREADS	
	  //masterBarrier(THREAD_COPY_RATES, tr);	   
	  //#endif
      
	  resetBranches(tr);
	  evaluateGeneric(tr, tr->start, TRUE);  
	  treeEvaluate(tr, 0.5);     
	  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->partitionData[model].protModels == AUTO)
		{		  
		  if(tr->perPartitionLH[model] > bestScores[model])
		    {
		      bestScores[model] = tr->perPartitionLH[model];
		      bestIndex[model] = i;		      
		    }
		}	      
	    }       
	}           
      
      if(processID == 0)
	printBothOpen("Automatic protein model assignment algorithm:\n\n");

      for(model = 0; model < tr->NumberOfModels; model++)
	{	   
	  if(tr->partitionData[model].protModels == AUTO)
	    {
	      tr->partitionData[model].autoProtModels = bestIndex[model];
	      initReversibleGTR(tr, model);  
	      if(processID == 0) 
		printBothOpen("\tPartition: %d best-scoring AA model: %s likelihood %f\n", model, protModels[tr->partitionData[model].autoProtModels], bestScores[model]);
	    }	 
	}
      
      if(processID == 0)
	printBothOpen("\n\n");

      //this barrier needs to be called in the library 
      //#ifdef _USE_PTHREADS	
      //masterBarrier(THREAD_COPY_RATES, tr);	   
      //#endif
          
      resetBranches(tr);
      evaluateGeneric(tr, tr->start, TRUE); 
      treeEvaluate(tr, 2.0);    
      
      if(tr->likelihood < startLH)
	{	
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->partitionData[model].protModels == AUTO)
		{
		  tr->partitionData[model].autoProtModels = oldIndex[model];
		  initReversibleGTR(tr, model);
		}
	    }
	  
	  //this barrier needs to be called in the library 	  
	  //#ifdef _USE_PTHREADS	
	  //masterBarrier(THREAD_COPY_RATES, tr);	   
	  //#endif 

	  restoreTL(rl, tr, 0);	
	  evaluateGeneric(tr, tr->start, TRUE);              
	}
      
      assert(tr->likelihood >= startLH);
      
      freeTL(rl);   
      free(rl); 
      
      free(oldIndex);
      free(bestIndex);
      free(bestScores);
    }
}

static void checkMatrixSymnmetriesAndLinkage(tree *tr, linkageList *ll)
{
  int 
    i;
  
  for(i = 0; i < ll->entries; i++)
    {
      int
	partitions = ll->ld[i].partitions;

      if(partitions > 1)
	{
	  int
	    k, 
	    reference = ll->ld[i].partitionList[0];

	  for(k = 1; k < partitions; k++)
	    {
	      int 
		index = ll->ld[i].partitionList[k];

	      int
		states = tr->partitionData[index].states,
		rates = ((states * states - states) / 2);
	      
	      if(tr->partitionData[reference].nonGTR != tr->partitionData[index].nonGTR)
		assert(0);
	      
	      if(tr->partitionData[reference].nonGTR)
		{
		  int 
		    j;
		  
		  for(j = 0; j < rates; j++)
		    {
		      if(tr->partitionData[reference].symmetryVector[j] != tr->partitionData[index].symmetryVector[j])
			assert(0);
		    }
		}
	    }	    
	}
    }
}


static void checkTolerance(double l1, double l2)
{
  if(l1 < l2)
    {   
      double 
	tolerance = MAX(l1, l2) * 0.000000000001;

      if(fabs(l1 - l2) > MIN(0.1, tolerance))
	{
	  printf("Likelihood problem in model optimization l1: %1.40f l2: %1.40f tolerance: %1.40f\n", l1, l2, tolerance);
	  assert(0);	
	}
    }
}

void modOpt(tree *tr, double likelihoodEpsilon, analdef *adef, int treeIteration)
{ 
  int 
    i, 
    catOpt = 0,
    *unlinked = (int *)malloc(sizeof(int) * tr->NumberOfModels);
  
  double 
    inputLikelihood,
    currentLikelihood,
    modelEpsilon = 0.0001;
  
  linkageList 
    *alphaList,
    *rateList,
    *freqList;      

  for(i = 0; i < tr->NumberOfModels; i++)
    unlinked[i] = i;

  //test code for library 
  if(0)
    {
      //assuming that we have three partitions for testing here 

      alphaList = initLinkageListString("0,1,2", tr);
      rateList  = initLinkageListString("0,1,1", tr);
    
      init_Q_MatrixSymmetries("0,1,2,3,4,5", tr, 0);
      init_Q_MatrixSymmetries("0,1,2,3,4,4", tr, 1);
      init_Q_MatrixSymmetries("0,1,1,2,3,4", tr, 2);
      
      //function that checks that partitions that have linked Q matrices as in our example above
      //will not have different configurations of the Q matrix as set by the init_Q_MatrixSymmetries() function
      //e.g., on would have HKY and one would have GTR, while the user claimes that they are linked
      //in our example, the Q matrices of partitions 1 and 2 are linked 
      //but we set different matrix symmetries via 
      // init_Q_MatrixSymmetries("0,1,2,3,4,4", tr, 1);
      // and
      // init_Q_MatrixSymmetries("0,1,1,2,3,4", tr, 2);
      //
      //the function just let's assertions fail for the time being .....

      checkMatrixSymnmetriesAndLinkage(tr, rateList);
    }
  else
    {
      alphaList = initLinkageList(unlinked, tr);
      freqList  = initLinkageList(unlinked, tr);
      rateList  = initLinkageListGTR(tr);
    }
   
  tr->start = tr->nodep[1];
                 
  if(adef->useCheckpoint && adef->mode == TREE_EVALUATION)
    {
      assert(ckp.state == MOD_OPT);
          	
      catOpt = ckp.catOpt;             
    }

  inputLikelihood = tr->likelihood;

  evaluateGeneric(tr, tr->start, TRUE); 

  assert(inputLikelihood == tr->likelihood);

  do
    {    
      if(adef->mode == TREE_EVALUATION)
	{
	  ckp.state = MOD_OPT;
	  
	  ckp.catOpt = catOpt;

	  ckp.treeIteration = treeIteration;
	  
	  writeCheckpoint(tr );
	}   
      
      currentLikelihood = tr->likelihood;     
          
#ifdef _DEBUG_MOD_OPT
      printf("start: %f\n", currentLikelihood);
#endif

      optRatesGeneric(tr, modelEpsilon, rateList);
                        
      evaluateGeneric(tr, tr->start, TRUE);    

#ifdef _DEBUG_MOD_OPT
      printf("after rates %f\n", tr->likelihood);
#endif                                                  

      autoProtein(tr);

      treeEvaluate(tr, 0.0625);    
      
#ifdef _DEBUG_MOD_OPT
      evaluateGeneric(tr, tr->start, TRUE); 
      printf("after br-len 1 %f\n", tr->likelihood);
#endif     

      evaluateGeneric(tr, tr->start, TRUE);

      optBaseFreqs(tr, modelEpsilon, freqList);
      
      evaluateGeneric(tr, tr->start, TRUE);
      
      treeEvaluate(tr, 0.0625);

#ifdef _DEBUG_MOD_OPT
      evaluateGeneric(tr, tr->start, TRUE); 
      printf("after optBaseFreqs 1 %f\n", tr->likelihood);
#endif 

      switch(tr->rateHetModel)
	{
	case GAMMA:      	  
	  optAlphasGeneric(tr, modelEpsilon, alphaList); 
	  
	  evaluateGeneric(tr, tr->start, TRUE); 
	 	 
#ifdef _DEBUG_MOD_OPT	 
	  printf("after alphas %f\n", tr->likelihood);
#endif	  
	  treeEvaluate(tr, 0.1);	  	 

#ifdef _DEBUG_MOD_OPT
	  evaluateGeneric(tr, tr->start, TRUE); 
	  printf("after br-len 2 %f\n", tr->likelihood);
#endif
	 
	  break;
	case CAT:
	  if(catOpt < 3)
	    {	      	     	     	     
	      evaluateGeneric(tr, tr->start, TRUE);
	      optimizeRateCategories(tr, tr->categories);	      	     	      	      	     

#ifdef _DEBUG_MOD_OPT
	      evaluateGeneric(tr, tr->start, TRUE); 
	      printf("after cat-opt %f\n", tr->likelihood);
#endif

	      catOpt++;
	    }
	  break;	  
	default:
	  assert(0);
	}                                

      checkTolerance(tr->likelihood, currentLikelihood);

      /*
	if(tr->likelihood < currentLikelihood)
	printf("%f %f\n", tr->likelihood, currentLikelihood);
	assert(tr->likelihood >= currentLikelihood);
      */
      
      printAAmatrix(tr, fabs(currentLikelihood - tr->likelihood));            
    }
  while(fabs(currentLikelihood - tr->likelihood) > likelihoodEpsilon);  
  
  free(unlinked);
  freeLinkageList(freqList);
  freeLinkageList(alphaList);
  freeLinkageList(rateList);
}


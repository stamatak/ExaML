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

#include "axml.h"

extern int processes; 

extern int Thorough;
extern int optimizeRateCategoryInvocations;
extern infoList iList;
extern char seq_file[1024];
extern char resultFileName[1024];
extern char tree_file[1024];
extern char workdir[1024];
extern char run_id[128];
extern double masterTime;
extern double accumulatedTime;

extern checkPointState ckp;
extern partitionLengths pLengths[MAX_MODEL];
extern char binaryCheckpointName[1024];
extern char binaryCheckpointInputName[1024];

extern int processID;



static int checker(tree *tr, nodeptr p)
{
  int group = tr->constraintVector[p->number];

  if(isTip(p->number, tr->mxtips))
    {
      group = tr->constraintVector[p->number];
      return group;
    }
  else
    {
      if(group != -9) 
	return group;

      group = checker(tr, p->next->back);
      if(group != -9) 
	return group;

      group = checker(tr, p->next->next->back);
      if(group != -9) 
	return group;

      return -9;
    }
}

boolean initrav (tree *tr, nodeptr p)
{ 
  nodeptr  q;
  
  if (!isTip(p->number, tr->mxtips)) 
    {      
      q = p->next;
      
      do 
	{	   
	  if (! initrav(tr, q->back))  return FALSE;		   
	  q = q->next;	
	} 
      while (q != p);  
      
      newviewGeneric(tr, p, FALSE);
    }
  
  return TRUE;
} 










/* #define _DEBUG_UPDATE */ 

boolean update(tree *tr, nodeptr p)
{       
  nodeptr  q; 
  boolean smoothedPartitions[NUM_BRANCHES];
  int i;
  double   z[NUM_BRANCHES], z0[NUM_BRANCHES];
  double _deltaz;

#ifdef _DEBUG_UPDATE
  double 
    startLH;

  evaluateGeneric(tr, p, FALSE);

  startLH = tr->likelihood;
#endif

  q = p->back;   

  for(i = 0; i < tr->numBranches; i++)
    z0[i] = q->z[i];    

  if(tr->numBranches > 1)
    makenewzGeneric(tr, p, q, z0, newzpercycle, z, TRUE);  
  else
    makenewzGeneric(tr, p, q, z0, newzpercycle, z, FALSE);
  
  for(i = 0; i < tr->numBranches; i++)    
    smoothedPartitions[i]  = tr->partitionSmoothed[i];
      
  for(i = 0; i < tr->numBranches; i++)
    {         
      if(!tr->partitionConverged[i])
	{	  
	  _deltaz = deltaz;
	    
	  if(ABS(z[i] - z0[i]) > _deltaz)  
	    {	      
	      smoothedPartitions[i] = FALSE;       
	    }	 

	  
	  
	  p->z[i] = q->z[i] = z[i];	 
	}
    }

#ifdef _DEBUG_UPDATE
  evaluateGeneric(tr, p, FALSE);

  if(tr->likelihood <= startLH)
    {
      if(fabs(tr->likelihood - startLH) > 0.01)
	{
	  printf("%f %f\n", startLH, tr->likelihood);
	  assert(0);      
	}
    }
#endif

  for(i = 0; i < tr->numBranches; i++)    
    tr->partitionSmoothed[i]  = smoothedPartitions[i];
  
  return TRUE;
}




boolean smooth (tree *tr, nodeptr p)
{
  nodeptr  q;
  
  if (! update(tr, p))               return FALSE; /*  Adjust branch */
  if (! isTip(p->number, tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  if (! smooth(tr, q->back))   return FALSE;
	  q = q->next;
	}
      
      if(tr->numBranches > 1)
	newviewGeneric(tr, p, TRUE);     
      else
	newviewGeneric(tr, p, FALSE);
    }
  
  return TRUE;
} 

boolean allSmoothed(tree *tr)
{
  int i;
  boolean result = TRUE;
  
  for(i = 0; i < tr->numBranches; i++)
    {
      if(tr->partitionSmoothed[i] == FALSE)
	result = FALSE;
      else
	tr->partitionConverged[i] = TRUE;
    }

  return result;
}



boolean smoothTree (tree *tr, int maxtimes)
{
  nodeptr  p, q;   
  int i, count = 0;
   
  p = tr->start;
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;

  while (--maxtimes >= 0) 
    {    
      for(i = 0; i < tr->numBranches; i++)	
	tr->partitionSmoothed[i] = TRUE;		

      if (! smooth(tr, p->back))       return FALSE;
      if (!isTip(p->number, tr->mxtips)) 
	{
	  q = p->next;
	  while (q != p) 
	    {
	      if (! smooth(tr, q->back))   return FALSE;
	      q = q->next;
	    }
	}
         
      count++;

      if (allSmoothed(tr)) 
	break;      
    }

  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;



  return TRUE;
} 



boolean localSmooth (tree *tr, nodeptr p, int maxtimes)
{ 
  nodeptr  q;
  int i;
  
  if (isTip(p->number, tr->mxtips)) return FALSE;
  
   for(i = 0; i < tr->numBranches; i++)	
     tr->partitionConverged[i] = FALSE;	

  while (--maxtimes >= 0) 
    {     
      for(i = 0; i < tr->numBranches; i++)	
	tr->partitionSmoothed[i] = TRUE;
	 	
      q = p;
      do 
	{
	  if (! update(tr, q)) return FALSE;
	  q = q->next;
        } 
      while (q != p);
      
      if (allSmoothed(tr)) 
	break;
    }

  for(i = 0; i < tr->numBranches; i++)
    {
      tr->partitionSmoothed[i] = FALSE; 
      tr->partitionConverged[i] = FALSE;
    }

  return TRUE;
}





static void resetInfoList(void)
{
  int i;

  iList.valid = 0;

  for(i = 0; i < iList.n; i++)    
    {
      iList.list[i].node = (nodeptr)NULL;
      iList.list[i].likelihood = unlikely;
    }    
}

void initInfoList(int n)
{
  int i;

  iList.n = n;
  iList.valid = 0;
  iList.list = (bestInfo *)malloc(sizeof(bestInfo) * n);

  for(i = 0; i < n; i++)
    {
      iList.list[i].node = (nodeptr)NULL;
      iList.list[i].likelihood = unlikely;
    }
}

void freeInfoList(void)
{ 
  free(iList.list);   
}


void insertInfoList(nodeptr node, double likelihood)
{
  int i;
  int min = 0;
  double min_l =  iList.list[0].likelihood;

  for(i = 1; i < iList.n; i++)
    {
      if(iList.list[i].likelihood < min_l)
	{
	  min = i;
	  min_l = iList.list[i].likelihood;
	}
    }

  if(likelihood > min_l)
    {
      iList.list[min].likelihood = likelihood;
      iList.list[min].node = node;
      iList.valid += 1;
    }

  if(iList.valid > iList.n)
    iList.valid = iList.n;
}


boolean smoothRegion (tree *tr, nodeptr p, int region)
{ 
  nodeptr  q;
  
  if (! update(tr, p))               return FALSE; /*  Adjust branch */

  if(region > 0)
    {
      if (!isTip(p->number, tr->mxtips)) 
	{                                 
	  q = p->next;
	  while (q != p) 
	    {
	      if (! smoothRegion(tr, q->back, --region))   return FALSE;
	      q = q->next;
	    }	
	  
	  newviewGeneric(tr, p, FALSE);
	}
    }
  
  return TRUE;
}

boolean regionalSmooth (tree *tr, nodeptr p, int maxtimes, int region)
  {
    nodeptr  q;
    int i;

    if (isTip(p->number, tr->mxtips)) return FALSE;            /* Should be an error */

    for(i = 0; i < tr->numBranches; i++)
      tr->partitionConverged[i] = FALSE;

    while (--maxtimes >= 0) 
      {	
	for(i = 0; i < tr->numBranches; i++)	  
	  tr->partitionSmoothed[i] = TRUE;
	  
	q = p;
	do 
	  {
	    if (! smoothRegion(tr, q, region)) return FALSE;
	    q = q->next;
	  } 
	while (q != p);
	
	if (allSmoothed(tr)) 
	  break;
      }

    for(i = 0; i < tr->numBranches; i++)
      tr->partitionSmoothed[i] = FALSE;
    for(i = 0; i < tr->numBranches; i++)
      tr->partitionConverged[i] = FALSE;
   
    return TRUE;
  } /* localSmooth */





nodeptr  removeNodeBIG (tree *tr, nodeptr p, int numBranches)
{  
  double   zqr[NUM_BRANCHES], result[NUM_BRANCHES];
  nodeptr  q, r;
  int i;
        
  q = p->next->back;
  r = p->next->next->back;
  
  for(i = 0; i < numBranches; i++)
    zqr[i] = q->z[i] * r->z[i];        
   
  makenewzGeneric(tr, q, r, zqr, iterations, result, FALSE);   

  for(i = 0; i < numBranches; i++)        
    tr->zqr[i] = result[i];

  hookup(q, r, result, numBranches); 
      
  p->next->next->back = p->next->back = (node *) NULL;

  return  q; 
}

nodeptr  removeNodeRestoreBIG (tree *tr, nodeptr p)
{
  nodeptr  q, r;
        
  q = p->next->back;
  r = p->next->next->back;  

  newviewGeneric(tr, q, FALSE);
  newviewGeneric(tr, r, FALSE);
  
  hookup(q, r, tr->currentZQR, tr->numBranches);

  p->next->next->back = p->next->back = (node *) NULL;
     
  return  q;
}


boolean insertBIG (tree *tr, nodeptr p, nodeptr q, int numBranches)
{
  nodeptr  r, s;
  int i;
  
  r = q->back;
  s = p->back;
      
  for(i = 0; i < numBranches; i++)
    tr->lzi[i] = q->z[i];
  
  if(Thorough)
    { 
      double  zqr[NUM_BRANCHES], zqs[NUM_BRANCHES], zrs[NUM_BRANCHES], lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax;      
      double defaultArray[NUM_BRANCHES];	
      double e1[NUM_BRANCHES], e2[NUM_BRANCHES], e3[NUM_BRANCHES];
      double *qz;
      
      qz = q->z;
      
      for(i = 0; i < numBranches; i++)
	defaultArray[i] = defaultz;
      
      makenewzGeneric(tr, q, r, qz, iterations, zqr, FALSE);           
      makenewzGeneric(tr, q, s, defaultArray, iterations, zqs, FALSE);                  
      makenewzGeneric(tr, r, s, defaultArray, iterations, zrs, FALSE);
      
      
      for(i = 0; i < numBranches; i++)
	{
	  lzqr = (zqr[i] > zmin) ? log(zqr[i]) : log(zmin); 
	  lzqs = (zqs[i] > zmin) ? log(zqs[i]) : log(zmin);
	  lzrs = (zrs[i] > zmin) ? log(zrs[i]) : log(zmin);
	  lzsum = 0.5 * (lzqr + lzqs + lzrs);
	  
	  lzq = lzsum - lzrs;
	  lzr = lzsum - lzqs;
	  lzs = lzsum - lzqr;
	  lzmax = log(zmax);
	  
	  if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
	  else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
	  else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          
	  
	  e1[i] = exp(lzq);
	  e2[i] = exp(lzr);
	  e3[i] = exp(lzs);
	}
      hookup(p->next,       q, e1, numBranches);
      hookup(p->next->next, r, e2, numBranches);
      hookup(p,             s, e3, numBranches);      		  
    }
  else
    {       
      double  z[NUM_BRANCHES]; 
      
      for(i = 0; i < numBranches; i++)
	{
	  z[i] = sqrt(q->z[i]);      
	  
	  if(z[i] < zmin) 
	    z[i] = zmin;
	  if(z[i] > zmax)
	    z[i] = zmax;
	}
      
      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);	                         
    }
  
  newviewGeneric(tr, p, FALSE);
  
  if(Thorough)
    {     
      localSmooth(tr, p, smoothings);  
      
      for(i = 0; i < numBranches; i++)
	{
	  tr->lzq[i] = p->next->z[i];
	  tr->lzr[i] = p->next->next->z[i];
	  tr->lzs[i] = p->z[i];            
	}
    }           
  
  return  TRUE;
}

boolean insertRestoreBIG (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;
  
  r = q->back;
  s = p->back;

  if(Thorough)
    {                        
      hookup(p->next,       q, tr->currentLZQ, tr->numBranches);
      hookup(p->next->next, r, tr->currentLZR, tr->numBranches);
      hookup(p,             s, tr->currentLZS, tr->numBranches);      		  
    }
  else
    {       
      double  z[NUM_BRANCHES];
      int i;
      
      for(i = 0; i < tr->numBranches; i++)
	{
	  double zz;
	  zz = sqrt(q->z[i]);     
	  if(zz < zmin) 
	    zz = zmin;
	  if(zz > zmax)
	    zz = zmax;
  	  z[i] = zz;
	}

      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);
    }   
    
  newviewGeneric(tr, p, FALSE);
       
  return  TRUE;
}


static void restoreTopologyOnly(tree *tr, bestlist *bt, bestlist *bestML)
{ 
  nodeptr p = tr->removeNode;
  nodeptr q = tr->insertNode;
  double qz[NUM_BRANCHES], pz[NUM_BRANCHES], p1z[NUM_BRANCHES], p2z[NUM_BRANCHES];
  nodeptr p1, p2, r, s;
  double currentLH = tr->likelihood;
  int i;
      
  p1 = p->next->back;
  p2 = p->next->next->back;
  
  for(i = 0; i < tr->numBranches; i++)
    {
      p1z[i] = p1->z[i];
      p2z[i] = p2->z[i];
    }
  
  hookup(p1, p2, tr->currentZQR, tr->numBranches);
  
  p->next->next->back = p->next->back = (node *) NULL;             
  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];           
    }
  
  r = q->back;
  s = p->back;
  
  if(Thorough)
    {                        
      hookup(p->next,       q, tr->currentLZQ, tr->numBranches);
      hookup(p->next->next, r, tr->currentLZR, tr->numBranches);
      hookup(p,             s, tr->currentLZS, tr->numBranches);      		  
    }
  else
    { 	
      double  z[NUM_BRANCHES];	
      for(i = 0; i < tr->numBranches; i++)
	{
	  z[i] = sqrt(q->z[i]);      
	  if(z[i] < zmin)
	    z[i] = zmin;
	  if(z[i] > zmax)
	    z[i] = zmax;
	}
      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);
    }     
  
  tr->likelihood = tr->bestOfNode;
    
  saveBestTree(bt, tr, TRUE);
  if(tr->saveBestTrees)
    saveBestTree(bestML, tr, FALSE);
  
  tr->likelihood = currentLH;
  
  hookup(q, r, qz, tr->numBranches);
  
  p->next->next->back = p->next->back = (nodeptr) NULL;
  
  if(Thorough)    
    hookup(p, s, pz, tr->numBranches);          
      
  hookup(p->next,       p1, p1z, tr->numBranches); 
  hookup(p->next->next, p2, p2z, tr->numBranches);      
}



boolean testInsertBIG (tree *tr, nodeptr p, nodeptr q)
{
  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r;
  boolean doIt = TRUE;
  double startLH = tr->endLH;
  int i;
  
  r = q->back; 
  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];
    }
  
  if(tr->constraintTree)
    {
      int rNumber, qNumber, pNumber;
      
      doIt = FALSE;
      
      rNumber = tr->constraintVector[r->number];
      qNumber = tr->constraintVector[q->number];
      pNumber = tr->constraintVector[p->number];
      
      if(pNumber == -9)
	pNumber = checker(tr, p->back);
      if(pNumber == -9)
	doIt = TRUE;
      else
	{
	  if(qNumber == -9)
	    qNumber = checker(tr, q);
	  
	  if(rNumber == -9)
	    rNumber = checker(tr, r);
	  
	  if(pNumber == rNumber || pNumber == qNumber)
	    doIt = TRUE;    	  
	}
    }
  
  if(doIt)
    {     
      if (! insertBIG(tr, p, q, tr->numBranches))       return FALSE;         
      
      evaluateGeneric(tr, p->next->next, FALSE);   
       
      if(tr->likelihood > tr->bestOfNode)
	{
	  tr->bestOfNode = tr->likelihood;
	  tr->insertNode = q;
	  tr->removeNode = p;   
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      tr->currentZQR[i] = tr->zqr[i];           
	      tr->currentLZR[i] = tr->lzr[i];
	      tr->currentLZQ[i] = tr->lzq[i];
	      tr->currentLZS[i] = tr->lzs[i];      
	    }
	}
      
      if(tr->likelihood > tr->endLH)
	{			  
	  tr->insertNode = q;
	  tr->removeNode = p;   
	  for(i = 0; i < tr->numBranches; i++)
	    tr->currentZQR[i] = tr->zqr[i];      
	  tr->endLH = tr->likelihood;                      
	}        
      
      hookup(q, r, qz, tr->numBranches);
      
      p->next->next->back = p->next->back = (nodeptr) NULL;
      
      if(Thorough)
	{
	  nodeptr s = p->back;
	  hookup(p, s, pz, tr->numBranches);      
	} 
      
      if((tr->doCutoff) && (tr->likelihood < startLH))
	{
	  tr->lhAVG += (startLH - tr->likelihood);
	  tr->lhDEC++;
	  if((startLH - tr->likelihood) >= tr->lhCutoff)
	    return FALSE;	    
	  else
	    return TRUE;
	}
      else
	return TRUE;
    }
  else
    return TRUE;  
}






 
void addTraverseBIG(tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav)
{  
  if (--mintrav <= 0) 
    {              
      if (! testInsertBIG(tr, p, q))  return;

    }
  
  if ((!isTip(q->number, tr->mxtips)) && (--maxtrav > 0)) 
    {    
      addTraverseBIG(tr, p, q->next->back, mintrav, maxtrav);
      addTraverseBIG(tr, p, q->next->next->back, mintrav, maxtrav);    
    }
} 





int rearrangeBIG(tree *tr, nodeptr p, int mintrav, int maxtrav)   
{  
  double   p1z[NUM_BRANCHES], p2z[NUM_BRANCHES], q1z[NUM_BRANCHES], q2z[NUM_BRANCHES];
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2, i;  
  boolean doP = TRUE, doQ = TRUE;
  
  if (maxtrav < 1 || mintrav > maxtrav)  return 0;
  q = p->back;
  
 
  
  if (!isTip(p->number, tr->mxtips) && doP) 
    {     
      p1 = p->next->back;
      p2 = p->next->next->back;
      
     
      if(!isTip(p1->number, tr->mxtips) || !isTip(p2->number, tr->mxtips))
	{
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      p1z[i] = p1->z[i];
	      p2z[i] = p2->z[i];	   	   
	    }
	  
	  if (! removeNodeBIG(tr, p,  tr->numBranches)) return badRear;
	  
	  if (!isTip(p1->number, tr->mxtips)) 
	    {
	      addTraverseBIG(tr, p, p1->next->back,
			     mintrav, maxtrav);         
	      addTraverseBIG(tr, p, p1->next->next->back,
			     mintrav, maxtrav);          
	    }
	  
	  if (!isTip(p2->number, tr->mxtips)) 
	    {
	      addTraverseBIG(tr, p, p2->next->back,
			     mintrav, maxtrav);
	      addTraverseBIG(tr, p, p2->next->next->back,
			     mintrav, maxtrav);          
	    }
	  	  
	  hookup(p->next,       p1, p1z, tr->numBranches); 
	  hookup(p->next->next, p2, p2z, tr->numBranches);	   	    	    
	  newviewGeneric(tr, p, FALSE);	   	    
	}
    }  
  
  if (!isTip(q->number, tr->mxtips) && maxtrav > 0 && doQ) 
    {
      q1 = q->next->back;
      q2 = q->next->next->back;
      
      /*if (((!q1->tip) && (!q1->next->back->tip || !q1->next->next->back->tip)) ||
	((!q2->tip) && (!q2->next->back->tip || !q2->next->next->back->tip))) */
      if (
	  (
	   ! isTip(q1->number, tr->mxtips) && 
	   (! isTip(q1->next->back->number, tr->mxtips) || ! isTip(q1->next->next->back->number, tr->mxtips))
	   )
	  ||
	  (
	   ! isTip(q2->number, tr->mxtips) && 
	   (! isTip(q2->next->back->number, tr->mxtips) || ! isTip(q2->next->next->back->number, tr->mxtips))
	   )
	  )
	{
	  
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      q1z[i] = q1->z[i];
	      q2z[i] = q2->z[i];
	    }
	  
	  if (! removeNodeBIG(tr, q, tr->numBranches)) return badRear;
	  
	  mintrav2 = mintrav > 2 ? mintrav : 2;
	  
	  if (/*! q1->tip*/ !isTip(q1->number, tr->mxtips)) 
	    {
	      addTraverseBIG(tr, q, q1->next->back,
			     mintrav2 , maxtrav);
	      addTraverseBIG(tr, q, q1->next->next->back,
			     mintrav2 , maxtrav);         
	    }
	  
	  if (/*! q2->tip*/ ! isTip(q2->number, tr->mxtips)) 
	    {
	      addTraverseBIG(tr, q, q2->next->back,
			     mintrav2 , maxtrav);
	      addTraverseBIG(tr, q, q2->next->next->back,
			     mintrav2 , maxtrav);          
	    }	   
	  
	  hookup(q->next,       q1, q1z, tr->numBranches); 
	  hookup(q->next->next, q2, q2z, tr->numBranches);
	  
	  newviewGeneric(tr, q, FALSE); 	   
	}
    } 
  
  return  1;
} 





double treeOptimizeRapid(tree *tr, int mintrav, int maxtrav, analdef *adef, bestlist *bt, bestlist *bestML)
{
  int 
    i, 
    index,
    *perm = (int*)NULL;   

  nodeRectifier(tr);

  if (maxtrav > tr->mxtips - 3)  
    maxtrav = tr->mxtips - 3;  
    
  resetInfoList();
  
  resetBestTree(bt);
 
  tr->startLH = tr->endLH = tr->likelihood;
 
  if(tr->doCutoff)
    {
      if(tr->bigCutoff)
	{	  
	  if(tr->itCount == 0)    
	    tr->lhCutoff = 0.5 * (tr->likelihood / -1000.0);    
	  else    		 
	    tr->lhCutoff = 0.5 * ((tr->lhAVG) / ((double)(tr->lhDEC))); 	  
	}
      else
	{
	  if(tr->itCount == 0)    
	    tr->lhCutoff = tr->likelihood / -1000.0;    
	  else    		 
	    tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));   
	}    

      tr->itCount = tr->itCount + 1;
      tr->lhAVG = 0;
      tr->lhDEC = 0;
    }
  
  /*
    printf("DoCutoff: %d\n", tr->doCutoff);
    printf("%d %f %f %f\n", tr->itCount, tr->lhAVG, tr->lhDEC, tr->lhCutoff);

    printf("%d %d\n", mintrav, maxtrav);
  */

  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {           
      tr->bestOfNode = unlikely;          

      if(adef->permuteTreeoptimize)
	index = perm[i];
      else
	index = i;     

      if(rearrangeBIG(tr, tr->nodep[index], mintrav, maxtrav))
	{    
	  if(Thorough)
	    {
	      if(tr->endLH > tr->startLH)                 	
		{			   	     
		  restoreTreeFast(tr);	 	 
		  tr->startLH = tr->endLH = tr->likelihood;	 
		  saveBestTree(bt, tr, TRUE); 
		  if(tr->saveBestTrees)
		    saveBestTree(bestML, tr, FALSE);
		}
	      else
		{ 		  
		  if(tr->bestOfNode != unlikely)		    	     
		    restoreTopologyOnly(tr, bt, bestML);		    
		}	   
	    }
	  else
	    {
	      insertInfoList(tr->nodep[index], tr->bestOfNode);	    
	      if(tr->endLH > tr->startLH)                 	
		{		      
		  restoreTreeFast(tr);	  	      
		  tr->startLH = tr->endLH = tr->likelihood;	  	 	  	  	  	  	  	  
		}	    	  
	    }
	}     
    }     

  if(!Thorough)
    {           
      Thorough = 1;  
      
      for(i = 0; i < iList.valid; i++)
	{ 	  
	  tr->bestOfNode = unlikely;
	  
	  if(rearrangeBIG(tr, iList.list[i].node, mintrav, maxtrav))
	    {	  
	      if(tr->endLH > tr->startLH)                 	
		{	 	     
		  restoreTreeFast(tr);	 	 
		  tr->startLH = tr->endLH = tr->likelihood;	 
		  saveBestTree(bt, tr, TRUE);
		  if(tr->saveBestTrees)
		    saveBestTree(bestML, tr, FALSE);
		}
	      else
		{ 
	      
		  if(tr->bestOfNode != unlikely)
		    {	     
		      restoreTopologyOnly(tr, bt, bestML);
		    }	
		}      
	    }
	}       
          
      Thorough = 0;
    }

  if(adef->permuteTreeoptimize)
    free(perm);

  return tr->startLH;     
}




boolean testInsertRestoreBIG (tree *tr, nodeptr p, nodeptr q)
{    
  if(Thorough)
    {
      if (! insertBIG(tr, p, q, tr->numBranches))       return FALSE;    
      
      evaluateGeneric(tr, p->next->next, FALSE);               
    }
  else
    {
      if (! insertRestoreBIG(tr, p, q))       return FALSE;
      
      {
	nodeptr x, y;
	x = p->next->next;
	y = p->back;
			
	if(! isTip(x->number, tr->mxtips) && isTip(y->number, tr->mxtips))
	  {
	    while ((! x->x)) 
	      {
		if (! (x->x))
		  newviewGeneric(tr, x, FALSE);		     
	      }
	  }
	
	if(isTip(x->number, tr->mxtips) && !isTip(y->number, tr->mxtips))
	  {
	    while ((! y->x)) 
	      {		  
		if (! (y->x))
		  newviewGeneric(tr, y, FALSE);
	      }
	  }
	
	if(!isTip(x->number, tr->mxtips) && !isTip(y->number, tr->mxtips))
	  {
	    while ((! x->x) || (! y->x)) 
	      {
		if (! (x->x))
		  newviewGeneric(tr, x, FALSE);
		if (! (y->x))
		  newviewGeneric(tr, y, FALSE);
	      }
	  }				      	
	
      }
	
      tr->likelihood = tr->endLH;
    }
     
  return TRUE;
} 

void restoreTreeFast(tree *tr)
{
  removeNodeRestoreBIG(tr, tr->removeNode);    
  testInsertRestoreBIG(tr, tr->removeNode, tr->insertNode);
}


static void writeTree(tree *tr, FILE *f)
{
  int 
    x = tr->mxtips + 3 * (tr->mxtips - 1);

  nodeptr
    base = tr->nodeBaseAddress;

  myBinFwrite(&(tr->start->number), sizeof(int), 1, f);
  myBinFwrite(&base, sizeof(nodeptr), 1, f);
  myBinFwrite(tr->nodeBaseAddress, sizeof(node), x, f);

}

int ckpCount = 0;


/**
    gathers patrat and rateCategory
 */
static void gatherDistributedCatInfos(tree *tr, int **rateCategory_result, double **patrat_result)
{
  /*
    countPerProc and displPerProc must be int, since the MPI functino
    signatures demand so
   */ 
  
  int 
    *countPerProc = (int*)NULL, 
    *displPerProc = (int*)NULL;

  calculateLengthAndDisplPerProcess(tr,  &countPerProc, &displPerProc);
  
  if(processID == 0)
    {
      *rateCategory_result = (int*)calloc((size_t)tr->originalCrunchedLength , sizeof(int));
      *patrat_result       = (double*)calloc((size_t)tr->originalCrunchedLength, sizeof(double));
    }
  
  gatherDistributedArray(tr, (void**) patrat_result,  tr->patrat_basePtr, MPI_DOUBLE , countPerProc, displPerProc); 
  gatherDistributedArray(tr, (void**) rateCategory_result, tr->rateCategory_basePtr, MPI_INT, countPerProc, displPerProc ); 

  free(countPerProc);
  free(displPerProc);
}


/** 
    added parameters patrat and rateCategory. The checkpoint writer
    has to gather this distributed information first. 
 */ 
static void writeCheckpointInner(tree *tr, int *rateCategory, double *patrat, analdef *adef)
{
  int   
    model; 
  
  char 
    extendedName[2048],
    buf[64];

  FILE 
    *f;

  /* only master should write the checkpoint */
  assert(processID == 0); 

  strcpy(extendedName,  binaryCheckpointName);
  strcat(extendedName, "_");
  sprintf(buf, "%d", ckpCount);
  strcat(extendedName, buf);  

  ckpCount++;

  f = myfopen(extendedName, "w"); 
  

  ckp.cmd.useMedian = tr->useMedian;
  ckp.cmd.saveBestTrees = tr->saveBestTrees;
  ckp.cmd.saveMemory = tr->saveMemory;
  ckp.cmd.searchConvergenceCriterion = tr->searchConvergenceCriterion;
  ckp.cmd.perGeneBranchLengths = adef->perGeneBranchLengths; //adef
  ckp.cmd.likelihoodEpsilon = adef->likelihoodEpsilon; //adef
  ckp.cmd.categories =  tr->categories;
  ckp.cmd.mode = adef->mode; //adef
  ckp.cmd.fastTreeEvaluation =  tr->fastTreeEvaluation;
  ckp.cmd.initialSet = adef->initialSet;//adef
  ckp.cmd.initial = adef->initial;//adef
  ckp.cmd.rateHetModel = tr->rateHetModel;
  ckp.cmd.autoProteinSelectionType = tr->autoProteinSelectionType;
  

  /* cdta */   
  
  ckp.accumulatedTime = accumulatedTime + (gettime() - masterTime);
  ckp.constraintTree = tr->constraintTree;

  /* printf("Acc time: %f\n", ckp.accumulatedTime); */

  myBinFwrite(&ckp, sizeof(checkPointState), 1, f);
  
  if(tr->constraintTree)
    myBinFwrite(tr->constraintVector, sizeof(int), 2 * tr->mxtips, f);  

  myBinFwrite(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myBinFwrite(tr->tree1, sizeof(char), tr->treeStringLength, f);


  if(tr->rateHetModel == CAT)
    {
      myBinFwrite(rateCategory, sizeof(int), tr->originalCrunchedLength, f);
      myBinFwrite(patrat, sizeof(double), tr->originalCrunchedLength, f);
    }  

  //end

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	dataType = tr->partitionData[model].dataType;
            
      myBinFwrite(&(tr->partitionData[model].numberOfCategories), sizeof(int), 1, f);
      myBinFwrite(tr->partitionData[model].perSiteRates, sizeof(double), tr->maxCategories, f);
      myBinFwrite(tr->partitionData[model].EIGN, sizeof(double), pLengths[dataType].eignLength, f);
      myBinFwrite(tr->partitionData[model].EV, sizeof(double),  pLengths[dataType].evLength, f);
      myBinFwrite(tr->partitionData[model].EI, sizeof(double),  pLengths[dataType].eiLength, f);  

      myBinFwrite(tr->partitionData[model].freqExponents, sizeof(double),  pLengths[dataType].frequenciesLength, f);
      myBinFwrite(tr->partitionData[model].frequencies,   sizeof(double),  pLengths[dataType].frequenciesLength, f);
      myBinFwrite(tr->partitionData[model].tipVector,     sizeof(double),  pLengths[dataType].tipVectorLength, f);       
      myBinFwrite(tr->partitionData[model].substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);

      //LG4X related variables 

      myBinFwrite(tr->partitionData[model].weights , sizeof(double), 4, f);
      myBinFwrite(tr->partitionData[model].weightExponents , sizeof(double), 4, f);
      //myBinFwrite(tr->partitionData[model].weightsBuffer , sizeof(double), 4, f);
      //myBinFwrite(tr->partitionData[model].weightExponentsBuffer , sizeof(double), 4, f);

      //LG4X end 

      if(tr->partitionData[model].protModels == LG4M || tr->partitionData[model].protModels == LG4X)
	{
	  int 
	    k;
	  
	  for(k = 0; k < 4; k++)
	    {
	      myBinFwrite(tr->partitionData[model].rawEIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	      myBinFwrite(tr->partitionData[model].EIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	      myBinFwrite(tr->partitionData[model].EV_LG4[k], sizeof(double),  pLengths[dataType].evLength, f);
	      myBinFwrite(tr->partitionData[model].EI_LG4[k], sizeof(double),  pLengths[dataType].eiLength, f);    
	      myBinFwrite(tr->partitionData[model].frequencies_LG4[k], sizeof(double),  pLengths[dataType].frequenciesLength, f);
	      myBinFwrite(tr->partitionData[model].tipVector_LG4[k], sizeof(double),  pLengths[dataType].tipVectorLength, f);  
	      myBinFwrite(tr->partitionData[model].substRates_LG4[k], sizeof(double),  pLengths[dataType].substRatesLength, f);    
	    }
	}
    
      myBinFwrite(&(tr->partitionData[model].alpha), sizeof(double), 1, f);
      myBinFwrite(&(tr->partitionData[model].gammaRates), sizeof(double), 4, f);
      
      myBinFwrite(&(tr->partitionData[model].protModels), sizeof(int), 1, f);
      myBinFwrite(&(tr->partitionData[model].autoProtModels), sizeof(int), 1, f);
    }
    
  if(ckp.state == MOD_OPT)
    {
      myBinFwrite(tr->likelihoods, sizeof(double), tr->numberOfTrees, f);
      myBinFwrite(tr->treeStrings, sizeof(char), (size_t)tr->treeStringLength * (size_t)tr->numberOfTrees, f);
    }

  writeTree(tr, f);

  fclose(f); 

  /* printBothOpen("\nCheckpoint written to: %s likelihood: %f\n", extendedName, tr->likelihood); */
}


void writeCheckpoint(tree *tr, analdef *adef)
{
  int 
    *rateCategory = (int *)NULL; 
  
  double 
    *patrat = (double *)NULL; 

  if(tr->rateHetModel == CAT)
    gatherDistributedCatInfos(tr, &rateCategory, &patrat); 

  if(processID == 0)
    {
      writeCheckpointInner(tr, rateCategory, patrat, adef); 

      if(tr->rateHetModel == CAT)
	{
	  free(rateCategory); 
	  free(patrat); 
	}
    }
}




static void readTree(tree *tr, FILE *f)
{
  int 
    nodeNumber,   
    x = tr->mxtips + 3 * (tr->mxtips - 1);

 

  
  
  nodeptr
    startAddress;

  myBinFread(&nodeNumber, sizeof(int), 1, f);

  tr->start = tr->nodep[nodeNumber];

  /*printf("Start: %d %d\n", tr->start->number, nodeNumber);*/

  myBinFread(&startAddress, sizeof(nodeptr), 1, f);

  /*printf("%u %u\n", (size_t)startAddress, (size_t)tr->nodeBaseAddress);*/



  myBinFread(tr->nodeBaseAddress, sizeof(node), x, f);
    
  {
    int i;    

    size_t         
      offset;

    boolean 
      addIt;

    if(startAddress > tr->nodeBaseAddress)
      {
	addIt = FALSE;
	offset = (size_t)startAddress - (size_t)tr->nodeBaseAddress;
      }
    else
      {
	addIt = TRUE;
	offset = (size_t)tr->nodeBaseAddress - (size_t)startAddress;
      }       

    for(i = 0; i < x; i++)
      {      	
	if(addIt)
	  {	    
	    tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next + offset);	
	    tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back + offset);
	  }
	else
	 {
	  
	   tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next - offset);	
	   tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back - offset);	   
	 } 
      }

  }
  
  evaluateGeneric(tr, tr->start, TRUE);  

  printBothOpen("ExaML Restart with likelihood: %1.50f\n", tr->likelihood);
}

static void genericError(void)
{
  printBothOpen("\nError: command lines used in initial run and re-start from checkpoint do not match!\n");
}

static void checkCommandLineArguments(tree *tr, analdef *adef)
{
  boolean
    match = TRUE;

  if(ckp.cmd.useMedian != tr->useMedian)
    {
      genericError();
      printBothOpen("\nDisagreement in median for gamma option: -a\n");
      match = FALSE;
    }
  
  if(ckp.cmd.saveBestTrees != tr->saveBestTrees)
    {
      genericError();
      printBothOpen("\nDisagreement in tree saving option: -B\n");
      match = FALSE;
    }
  
  if(ckp.cmd.saveMemory != tr->saveMemory)
    {
      genericError();
      printBothOpen("\nDisagreement in memory saving option: -S\n");
      match = FALSE;
    }
  
  if(ckp.cmd.searchConvergenceCriterion != tr->searchConvergenceCriterion)
    {
      genericError();
      printBothOpen("\nDisagreement in search convergence criterion: -D\n");
      match = FALSE;
    }
  
  if(ckp.cmd.perGeneBranchLengths != adef->perGeneBranchLengths)
    {
      genericError();
      printBothOpen("\nDisagreement in using per-partition branch lengths: -M\n");
      match = FALSE;
    }
  
  if(ckp.cmd.likelihoodEpsilon != adef->likelihoodEpsilon)
     {
      genericError();
      printBothOpen("\nDisagreement in likelihood epsilon value: -e\n");
      match = FALSE;
    }
  
  if(ckp.cmd.categories !=  tr->categories)
    {
      genericError();
      printBothOpen("\nDisagreement in number of PSR rate categories: -c\n");
      match = FALSE;
    }

  if(ckp.cmd.mode != adef->mode)
    {
      genericError();
      printBothOpen("\nDisagreement in tree search or evaluation mode\n");
      match = FALSE;
    }
  
  if(ckp.cmd.fastTreeEvaluation !=  tr->fastTreeEvaluation)
    {
      genericError();
      printBothOpen("\nDisagreement in fast tree evaluation: -e|-E\n");
      match = FALSE;
    }
  
  

  if(ckp.cmd.initialSet != adef->initialSet)
     {
      genericError();
      printBothOpen("\nDisagreement in rearrangement radius limitation setting: -i\n");
      match = FALSE;
    }
  
  if(ckp.cmd.initial != adef->initial)
     {
      genericError();
      printBothOpen("\nDisagreement in rearrangement radius value: -i\n");
      match = FALSE;
    }
  
  if(ckp.cmd.rateHetModel != tr->rateHetModel)
     {
      genericError();
      printBothOpen("\nDisagreement in rate heterogeneity model: -m\n");
      match = FALSE;
    }

   if(ckp.cmd.autoProteinSelectionType != tr->autoProteinSelectionType)
     {
      genericError();
      printBothOpen("\nDisagreement in protein model selection criterion: --auto-prot\n");
      match = FALSE;
    }

  if(!match)
    {
      printBothOpen("\nExaML will exit now ...\n\n");
      errorExit(-1);
    }
}

static void readCheckpoint(tree *tr, analdef *adef)
{
  int   
    model; 

  FILE 
    *f = myfopen(binaryCheckpointInputName, "rb");

  /* cdta */   

  myBinFread(&ckp, sizeof(checkPointState), 1, f);

  checkCommandLineArguments(tr, adef);

  tr->constraintTree = ckp.constraintTree;

  if(tr->constraintTree)
    myBinFread(tr->constraintVector, sizeof(int), 2 * tr->mxtips, f);  

  tr->ntips = tr->mxtips;

  

  tr->startLH    = ckp.tr_startLH;
  tr->endLH      = ckp.tr_endLH;
  tr->likelihood = ckp.tr_likelihood;
  tr->bestOfNode = ckp.tr_bestOfNode;
  
  tr->lhCutoff   = ckp.tr_lhCutoff;
  tr->lhAVG      = ckp.tr_lhAVG;
  tr->lhDEC      = ckp.tr_lhDEC;
  tr->itCount    = ckp.tr_itCount;
  Thorough       = ckp.Thorough;
  
  accumulatedTime = ckp.accumulatedTime;

  /* printf("Accumulated time so far: %f\n", accumulatedTime); */

  optimizeRateCategoryInvocations = ckp.optimizeRateCategoryInvocations;


  myBinFread(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myBinFread(tr->tree1, sizeof(char), tr->treeStringLength, f);

  if(tr->searchConvergenceCriterion && processID == 0)
    {
      int bCounter = 0;
      
      if((ckp.state == FAST_SPRS && ckp.fastIterations > 0) ||
	 (ckp.state == SLOW_SPRS && ckp.thoroughIterations > 0))
	{ 

#ifdef _DEBUG_CHECKPOINTING    
	  printf("parsing Tree 0\n");
#endif

	  treeReadTopologyString(tr->tree0, tr);   
	  
	  bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, 0, BIPARTITIONS_RF, (branchInfo *)NULL,
				  &bCounter, 1, FALSE, FALSE);
	  
	  assert(bCounter == tr->mxtips - 3);
	}
      
      bCounter = 0;
      
      if((ckp.state == FAST_SPRS && ckp.fastIterations > 1) ||
	 (ckp.state == SLOW_SPRS && ckp.thoroughIterations > 1))
	{

#ifdef _DEBUG_CHECKPOINTING
	  printf("parsing Tree 1\n");
#endif

	  treeReadTopologyString(tr->tree1, tr); 
	  
	  bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, 1, BIPARTITIONS_RF, (branchInfo *)NULL,
				  &bCounter, 1, FALSE, FALSE);
	  
	  assert(bCounter == tr->mxtips - 3);
	}
    }

  
  if(tr->rateHetModel == CAT )
    {
      /* every process reads its data */

      /* Andre will this also work if we re-start with a different
	 number of processors? have you tested? 
	 
	 => Andre: yes that works: before writing the checkpoint, we
	 gather all lhs/patrat with gatherDistributedCatInfos. This
	 function calls gatherDisributedArray in
	 communication.c. gatherDistributedArray takes care of
	 reordering the data it obtained from the various processes,
	 such the correct global array (i.e., indexing consistent with
	 character position) is obtained.  Thus, the indexing below
	 (for reading in the patrat/lhs again) works correctly.  */
      

      /* Andre I think tr->originalCrunchedLength is of type size_t???
	 -> casting required ?  
	 
	 Andre: in the very worst case, pPos overflows. There is not
	 much one can do here. See explanation about fseek/fseeko at
	 other location. But I have added an assert, in case something
	 goes wrong */
      exa_off_t
	rPos = exa_ftell(f),      
	pPos  = rPos + sizeof(int) * tr->originalCrunchedLength; 

      /* fails, in case reading failed (ftello returns -1) or an overflow happened */
      assert( ! ( rPos < 0 || pPos < 0 ) && rPos <= pPos ) ;
      
      /* first patrat then rateCategory */
      
      Assign *aIter =  tr->partAssigns,
	*aEnd = &(tr->partAssigns [ tr->numAssignments ]) ; 
      
      /* Andre coould you maybe add a drawing (scanned drawn by hand if you like) documenting this layout ? => Andre: TODO  */

      while(aIter != aEnd)
	{
	  if(aIter->procId == processID)
	    {
	      pInfo
		*partition = &(tr->partitionData[aIter->partitionId]); 
	      exa_off_t
		theOffset = pPos + (partition->lower + aIter->offset)  * sizeof(double); 
	      assert(pPos <= theOffset); 

	      exa_fseek(f,  theOffset, SEEK_SET); 
	      
	      myBinFread(partition->patrat, sizeof(double), aIter->width, f);  

	      theOffset = rPos + (partition->lower + aIter->offset) * sizeof(int); 
	      assert(rPos <= theOffset); 
	      exa_fseek(f, theOffset, SEEK_SET); 
	      myBinFread(partition->rateCategory, sizeof(int), aIter->width, f); 
	    } 
	  ++aIter; 
	}

      /* Set file pointer to the end of both of the arrays    */
      exa_fseek(f, pPos + tr->originalCrunchedLength * sizeof(double) , SEEK_SET); 
    }


  
  
 

  //end

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	dataType = tr->partitionData[model].dataType;
            
      myBinFread(&(tr->partitionData[model].numberOfCategories), sizeof(int), 1, f);
      myBinFread(tr->partitionData[model].perSiteRates, sizeof(double), tr->maxCategories, f);
      myBinFread(tr->partitionData[model].EIGN, sizeof(double), pLengths[dataType].eignLength, f);
      myBinFread(tr->partitionData[model].EV, sizeof(double),  pLengths[dataType].evLength, f);
      myBinFread(tr->partitionData[model].EI, sizeof(double),  pLengths[dataType].eiLength, f);  

      myBinFread(tr->partitionData[model].freqExponents, sizeof(double),  pLengths[dataType].frequenciesLength, f);
      myBinFread(tr->partitionData[model].frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
      myBinFread(tr->partitionData[model].tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);  
      myBinFread(tr->partitionData[model].substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);  

      //LG4X related variables 

      myBinFread(tr->partitionData[model].weights , sizeof(double), 4, f);
      myBinFread(tr->partitionData[model].weightExponents , sizeof(double), 4, f);
      //myBinFread(tr->partitionData[model].weightsBuffer , sizeof(double), 4, f);
      //myBinFread(tr->partitionData[model].weightExponentsBuffer , sizeof(double), 4, f);

      //LG4X end 

      if(tr->partitionData[model].protModels == LG4X || tr->partitionData[model].protModels == LG4M)
	{
	  int 
	    k;
	  
	  for(k = 0; k < 4; k++)
	    {
	       myBinFread(tr->partitionData[model].rawEIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	      myBinFread(tr->partitionData[model].EIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	      myBinFread(tr->partitionData[model].EV_LG4[k], sizeof(double),  pLengths[dataType].evLength, f);
	      myBinFread(tr->partitionData[model].EI_LG4[k], sizeof(double),  pLengths[dataType].eiLength, f);    
	      myBinFread(tr->partitionData[model].frequencies_LG4[k], sizeof(double),  pLengths[dataType].frequenciesLength, f);
	      myBinFread(tr->partitionData[model].tipVector_LG4[k], sizeof(double),  pLengths[dataType].tipVectorLength, f);  
	      myBinFread(tr->partitionData[model].substRates_LG4[k], sizeof(double),  pLengths[dataType].substRatesLength, f);    
	    }
	}


      myBinFread(&(tr->partitionData[model].alpha), sizeof(double), 1, f);      
      myBinFread(&(tr->partitionData[model].gammaRates), sizeof(double), 4, f);
      //conditional added by Andre modified by me
      //only overwrite values of discrete gamma cats by calling makeGammaCats if not using 
      //LG4X!
      if(tr->rateHetModel != CAT && !(tr->partitionData[model].protModels == LG4X))
	makeGammaCats(tr->partitionData[model].alpha, tr->partitionData[model].gammaRates, 4, tr->useMedian); 

      myBinFread(&(tr->partitionData[model].protModels), sizeof(int), 1, f);
      myBinFread(&(tr->partitionData[model].autoProtModels), sizeof(int), 1, f);
    }
    
  if(ckp.state == MOD_OPT)
    {
      myBinFread(tr->likelihoods, sizeof(double), tr->numberOfTrees, f);
      myBinFread(tr->treeStrings, sizeof(char), (size_t)tr->treeStringLength * (size_t)tr->numberOfTrees, f);
    }

  if(tr->rateHetModel == CAT)
    checkPerSiteRates(tr); 

  readTree(tr, f);
  fclose(f); 
}


void restart(tree *tr, analdef *adef)
{  
  readCheckpoint(tr, adef);

  switch(ckp.state)
    {
    case REARR_SETTING:      
      assert(adef->mode == BIG_RAPID_MODE);
      break;
    case FAST_SPRS:
      assert(adef->mode == BIG_RAPID_MODE);
      break;
    case SLOW_SPRS:
      assert(adef->mode == BIG_RAPID_MODE);
      break;
    case MOD_OPT:
      assert(adef->mode == TREE_EVALUATION);
      break;
    default:
      assert(0);
    }
}

int determineRearrangementSetting(tree *tr,  analdef *adef, bestlist *bestT, bestlist *bt, bestlist *bestML)
{
  const 
    int MaxFast = 26;
  
  int 
    i,   
    maxtrav = 5, 
    bestTrav = 5;

  double 
    startLH = tr->likelihood; 
  
  boolean 
    impr   = TRUE,
    cutoff = tr->doCutoff;
   
  if(adef->useCheckpoint)
    {
      assert(ckp.state == REARR_SETTING);
         
      maxtrav = ckp.maxtrav;
      bestTrav = ckp.bestTrav;
      startLH  = ckp.startLH;
      impr     = ckp.impr;
      
      cutoff = ckp.cutoff;

      adef->useCheckpoint = FALSE;
    }
  
  tr->doCutoff = FALSE;      

  resetBestTree(bt);    
 
#ifdef _DEBUG_CHECKPOINTING
  printBothOpen("MAXTRAV: %d\n", maxtrav);
#endif

  assert(Thorough == 0);

  while(impr && maxtrav < MaxFast)
    {	
      recallBestTree(bestT, 1, tr);     
      nodeRectifier(tr);            
      
      /* Andre I believe that the code below, except for
	 writeCheckpoint cann still only be executed by process 0 =>
	 Andre: all other processes need to enter writeCheckpoint,
	 because of the gather that happens there. But the assignments
	 to the checkpoint state are not necessary for all processes;
	 does it matter? */
      {
	ckp.optimizeRateCategoryInvocations = optimizeRateCategoryInvocations;
	  
	ckp.cutoff = cutoff;
	ckp.state = REARR_SETTING;     
	ckp.maxtrav = maxtrav;
	ckp.bestTrav = bestTrav;
	ckp.startLH  = startLH;
	ckp.impr = impr;
	  
	ckp.tr_startLH  = tr->startLH;
	ckp.tr_endLH    = tr->endLH;
	ckp.tr_likelihood = tr->likelihood;
	ckp.tr_bestOfNode = tr->bestOfNode;
	  
	ckp.tr_lhCutoff = tr->lhCutoff;
	ckp.tr_lhAVG    = tr->lhAVG;
	ckp.tr_lhDEC    = tr->lhDEC;      
	ckp.tr_itCount  = tr->itCount;
	  
	  
	writeCheckpoint(tr, adef);    
      }

      if (maxtrav > tr->mxtips - 3)  
	maxtrav = tr->mxtips - 3;    
 
      tr->startLH = tr->endLH = tr->likelihood;
      
      /* printBothOpen("TRAV: %d lh %f MNZC %d\n", maxtrav, tr->likelihood, mnzc); */

      {
	int changes = 0;
	
	for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	  {                	         
	    tr->bestOfNode = unlikely;
	    
	    if(rearrangeBIG(tr, tr->nodep[i], 1, maxtrav))
	      {	     
		if(tr->endLH > tr->startLH)                 	
		  {		 	 	      
		    restoreTreeFast(tr);	        	  	 	  	      
		    tr->startLH = tr->endLH = tr->likelihood;			  
		    changes++;
		  }	         	       	
	      }
	  }
	
      
	/*
	  evaluateGeneric(tr, tr->start, TRUE);	
	  
	  printBothOpen("Changes: %d TRAV: %d lh %f MNZC %d\n", changes, maxtrav, tr->likelihood, mnzc);
	*/      
      }
      
      treeEvaluate(tr, 0.25);

      /* printBothOpen("TRAV: %d lh %f MNZC %d\n", maxtrav, tr->likelihood, mnzc); */

      saveBestTree(bt, tr, TRUE); 
      if(tr->saveBestTrees)
	saveBestTree(bestML, tr, FALSE);           
                                         
#ifdef _DEBUG_CHECKPOINTING
      printBothOpen("TRAV: %d lh %f MNZC %d\n", maxtrav, tr->likelihood, mnzc);
#endif

      if(tr->likelihood > startLH)
	{	 
	  startLH = tr->likelihood; 	  	  	  
	  printLog(tr);	  
	  bestTrav = maxtrav;	 
	  impr = TRUE;
	}
      else	
	impr = FALSE;	
      
      
      
      if(tr->doCutoff)
	{
	  tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));       
  
	  tr->itCount =  tr->itCount + 1;
	  tr->lhAVG = 0;
	  tr->lhDEC = 0;
	}
      
      maxtrav += 5;
      
             
    }

  recallBestTree(bt, 1, tr);
  
  tr->doCutoff = cutoff; 
  
#ifdef _DEBUG_CHECKPOINTING
  printBothOpen("BestTrav %d\n", bestTrav);
#endif

  return bestTrav;     
}





void computeBIGRAPID (tree *tr, analdef *adef, boolean estimateModel) 
{   
  int
    i,
    impr, 
    bestTrav = 0,
    treeVectorLength = 0,
    rearrangementsMax = 0, 
    rearrangementsMin = 0,    
    thoroughIterations = 0,
    fastIterations = 0;
   
  double 
    lh = unlikely, 
    previousLh = unlikely, 
    difference, 
    epsilon;              
  
  bestlist 
    *bestML,
    *bestT, 
    *bt;        
 
  /* now here is the RAxML hill climbing search algorithm */
  
  tr->lhAVG = 0.0;
  tr->lhDEC = 0.0;

  /* initialization for the hash table to compute RF distances */

  if(tr->searchConvergenceCriterion && processID == 0)   
    treeVectorLength = 1;
     
  /* initialize two lists of size 1 and size 20 that will keep track of the best 
     and 20 best tree topologies respectively */

  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);
      
  bt = (bestlist *) malloc(sizeof(bestlist));      
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips);    



  if(tr->saveBestTrees > 0)
    { 
      bestML = (bestlist *) malloc(sizeof(bestlist));      
      bestML->ninit = 0;
      initBestTree(bestML, tr->saveBestTrees, tr->mxtips);  
    }
  else
    bestML = (bestlist *)NULL;
  
  
  /* initialize an additional data structure used by the search algo, all of this is pretty 
     RAxML-specific and should probably not be in the library */

  initInfoList(50);
 
  /* some pretty atbitrary thresholds */

  difference = 10.0;
  epsilon = 0.01;    
    
  /* Thorough = 0 means that we will do fast SPR inbsertions without optimizing the 
     three branches adjacent to the subtree insertion position via Newton-Raphson 
  */

  Thorough = 0;     
  
  /* if we are not using a checkpoint and estimateModel is set to TRUE we call the function 
     that optimizes model parameters, such as the CAT model assignment, the alpha paremeter
     or the rates in the GTR matrix. Otherwise we just optimize the branch lengths. Note that 
     the second parameter of treeEvaluate() controls how many times we will iterate over all branches 
     of the tree until we give up, provided that, the br-len opt. has not converged before.
  */

  if(!adef->useCheckpoint)
    {
      if(estimateModel)
	modOpt(tr, 10.0, adef, 0);
      else
	treeEvaluate(tr, 2);  
    }

  /* print some stuff to the RAxML_log file */

  printLog(tr); 

  /* save the current tree (which is the input tree parsed via -t in the bestT list */

  saveBestTree(bestT, tr, TRUE);
  
  /* if the rearrangmenet radius has been set by the user ie. adef->initailSet == TRUE 
     then just set the apppropriate parameter.
     Otherwise, call the function  determineRearrangementSetting() that seeks 
     for the best radius by executing SPR moves on the initial tree with different radii
     and returns the smallest radius that yields the best log likelihood score after 
     applying one cycle of SPR moves to the tree 
  */

  if(!adef->initialSet)   
    {
      if((!adef->useCheckpoint) || (adef->useCheckpoint && ckp.state == REARR_SETTING))
	{
	  bestTrav = adef->bestTrav = determineRearrangementSetting(tr, adef, bestT, bt, bestML);     	  
	  printBothOpen("\nBest rearrangement radius: %d\n", bestTrav);
	}
    }
  else
    {
      bestTrav = adef->bestTrav = adef->initial;       
      printBothOpen("\nUser-defined rearrangement radius: %d\n", bestTrav);
    }

  
  /* some checkpointing noise */
  if(!(adef->useCheckpoint && (ckp.state == FAST_SPRS || ckp.state == SLOW_SPRS)))
    {      

      /* optimize model params more thoroughly or just optimize branch lengths */
      if(estimateModel)
	modOpt(tr, 5.0, adef, 0);
      else
	treeEvaluate(tr, 1);   
    }
  
  /* save the current tree again, while the topology has not changed, the branch lengths have changed in the meantime, hence
     we need to store them again */

  saveBestTree(bestT, tr, TRUE); 

  /* set the loop variable to TRUE */

  impr = 1;

  /* this is for the additional RAxML heuristics described imn this paper here:

     A. Stamatakis,  F. Blagojevic, C.D. Antonopoulos, D.S. Nikolopoulos: "Exploring new Search Algorithms and Hardware for Phylogenetics: RAxML meets the IBM Cell". 
     In Journal of VLSI Signal Processing Systems, 48(3):271-286, 2007.

     This is turned on by default 
  */
     
  
  if(tr->doCutoff)
    tr->itCount = 0;

  /* figure out where to continue computations if we restarted from a checkpoint */

  if(adef->useCheckpoint && ckp.state == FAST_SPRS)
    goto START_FAST_SPRS;

  if(adef->useCheckpoint && ckp.state == SLOW_SPRS)
    goto START_SLOW_SPRS;
  
  while(impr)
    {              
    START_FAST_SPRS:
      /* if re-starting from checkpoint set the required variable values to the 
	 values that they had when the checkpoint was written */

      if(adef->useCheckpoint && ckp.state == FAST_SPRS)
	{
	  optimizeRateCategoryInvocations = ckp.optimizeRateCategoryInvocations;   	

  
	  impr = ckp.impr;
	  Thorough = ckp.Thorough;
	  bestTrav = ckp.bestTrav;
	  treeVectorLength = ckp.treeVectorLength;
	  rearrangementsMax = ckp.rearrangementsMax;
	  rearrangementsMin = ckp.rearrangementsMin;
	  thoroughIterations = ckp.thoroughIterations;
	  fastIterations = ckp.fastIterations;
   
  
	  lh = ckp.lh;
	  previousLh = ckp.previousLh;
	  difference = ckp.difference;
	  epsilon    = ckp.epsilon;                    
           
	
	  tr->likelihood = ckp.tr_likelihood;
              
	  tr->lhCutoff = ckp.tr_lhCutoff;
	  tr->lhAVG    = ckp.tr_lhAVG;
	  tr->lhDEC    = ckp.tr_lhDEC;   	 
	  tr->itCount = ckp.tr_itCount;	  

	  adef->useCheckpoint = FALSE;
	}
      else
	/* otherwise, restore the currently best tree */
	recallBestTree(bestT, 1, tr); 

      /* save states of algorithmic/heuristic variables for printing the next checkpoint */

      /* 
	 Andre I believe that the code below, except for
	 writeCheckpoint cann still only be executed by process 0 =>
	 Andre: see above
       */ 
      {              
	ckp.state = FAST_SPRS;  
	ckp.optimizeRateCategoryInvocations = optimizeRateCategoryInvocations;              
	  
	  
	ckp.impr = impr;
	ckp.Thorough = Thorough;
	ckp.bestTrav = bestTrav;
	ckp.treeVectorLength = treeVectorLength;
	ckp.rearrangementsMax = rearrangementsMax;
	ckp.rearrangementsMin = rearrangementsMin;
	ckp.thoroughIterations = thoroughIterations;
	ckp.fastIterations = fastIterations;
	  
	  
	ckp.lh = lh;
	ckp.previousLh = previousLh;
	ckp.difference = difference;
	ckp.epsilon    = epsilon; 
	  
	  
	ckp.bestTrav = bestTrav;       
	ckp.impr = impr;
	  
	ckp.tr_startLH  = tr->startLH;
	ckp.tr_endLH    = tr->endLH;
	ckp.tr_likelihood = tr->likelihood;
	ckp.tr_bestOfNode = tr->bestOfNode;
	  
	ckp.tr_lhCutoff = tr->lhCutoff;
	ckp.tr_lhAVG    = tr->lhAVG;
	ckp.tr_lhDEC    = tr->lhDEC;       
	ckp.tr_itCount  = tr->itCount;       
	  
	/* write a binary checkpoint */
	writeCheckpoint(tr, adef); 
      }	

      /* this is the aforementioned convergence criterion that requires computing the RF,
	 let's not worry about this right now */

      if(tr->searchConvergenceCriterion && processID == 0)
	{
	  int 
	    bCounter = 0; 

	  char 
	    *buffer = (char*)calloc(tr->treeStringLength, sizeof(char));	  	      	 	  	  	

	  if(fastIterations > 1)	    	      
	    cleanupHashTable(tr->h, (fastIterations % 2));		
	    	  	 

	  bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, fastIterations % 2, BIPARTITIONS_RF, (branchInfo *)NULL,
				  &bCounter, 1, FALSE, FALSE);	    
	  	  
	   
#ifdef _DEBUG_CHECKPOINTING
	  printf("Storing tree in slot %d\n", fastIterations % 2);
#endif

	  Tree2String(buffer, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);

	  if(fastIterations % 2 == 0)	      
	    memcpy(tr->tree0, buffer, tr->treeStringLength * sizeof(char));
	  else
	    memcpy(tr->tree1, buffer, tr->treeStringLength * sizeof(char));	    
	  
	  free(buffer);	  

	  assert(bCounter == tr->mxtips - 3);	    	   	  	 

	  if(fastIterations > 0)
	    {
	      double 
		rrf = convergenceCriterion(tr->h, tr->mxtips);
	      
	      MPI_Bcast(&rrf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	      
	      if(rrf <= 0.01) /* 1% cutoff */
		{
		  printBothOpen("ML fast search converged at fast SPR cycle %d with stopping criterion\n", fastIterations);
		  printBothOpen("Relative Robinson-Foulds (RF) distance between respective best trees after one succseful SPR cycle: %f%s\n", rrf, "%");
		  cleanupHashTable(tr->h, 0);
		  cleanupHashTable(tr->h, 1);
		  goto cleanup_fast;
		}
	      else		    
		printBothOpen("ML search convergence criterion fast cycle %d->%d Relative Robinson-Foulds %f\n", fastIterations - 1, fastIterations, rrf);
	    }
	}

      if(tr->searchConvergenceCriterion && processID != 0 && fastIterations > 0)
	{
	  double 
	    rrf;
	  
	  MPI_Bcast(&rrf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	 
	  if(rrf <= 0.01) /* 1% cutoff */		   
	    goto cleanup_fast;	      
	}

      
      /* count how many fast iterations with so-called fast SPR moves we have executed */

      fastIterations++;	

      /* optimize branch lengths */
     
     
      treeEvaluate(tr, 1.0);    
     
      /* save the tree with those branch lengths again */
      
      saveBestTree(bestT, tr, TRUE);           

      /* print the log likelihood */

      printLog(tr);    

      /* print this intermediate tree to file */
     
      printResult(tr, adef, FALSE);    

      /* update the current best likelihood */

      lh = previousLh = tr->likelihood;
            
      /* in here we actually do a cycle of SPR moves */

      treeOptimizeRapid(tr, 1, bestTrav, adef, bt, bestML);   
          
      /* set impr to 0 since in the immediately following for loop we check if the SPR moves above have generated 
	 a better tree */

      impr = 0;
	  
      /* loop over the 20 best trees generated by the fast SPR moves, and check if they improve the likelihood after all of their branch lengths
	 have been optimized */

      for(i = 1; i <= bt->nvalid; i++)
	{	    	
	  /* restore tree i from list generated by treeOptimizeRapid */
	  	   
	  recallBestTree(bt, i, tr);
	  
	  /* optimize branch lengths of this tree */

	  treeEvaluate(tr, 0.25);

	  /* calc. the likelihood improvement */

	  difference = ((tr->likelihood > previousLh)? 
			tr->likelihood - previousLh: 
			previousLh - tr->likelihood); 	    

	  /* if the likelihood has improved save the current tree as best tree and continue */
	  /* note that we always compre this tree to the likelihood of the previous best tree */
	  
	  if(tr->likelihood > lh && difference > epsilon)
	    {
	      impr = 1;	       
	      lh = tr->likelihood;	       	     
	      saveBestTree(bestT, tr, TRUE);
	      
	    }	   	   
	}
#ifdef _DEBUG_CHECKPOINTING
      printBothOpen("FAST LH: %f\n", lh);
#endif

	
    }
  
  /* needed for this RF-based convergence criterion that I actually describe in here:

     A. Stamatakis: "Phylogenetic Search Algorithms for Maximum Likelihood". In M. Elloumi, A.Y. Zomaya, editors. 
     Algorithms in Computational Biology: techniques, Approaches and Applications, John Wiley and Sons

     a copy of this book is in my office */

  if(tr->searchConvergenceCriterion && processID == 0)
    {
      cleanupHashTable(tr->h, 0);
      cleanupHashTable(tr->h, 1);
    }
  
 cleanup_fast:  
  /*
    now we have jumped out of the loop that executes 
     fast SPRs, and next we will execute a loop that executes throough SPR cycles (with SPR moves 
     that optimize via newton-Raphson all adjacent branches to the insertion point) 
     until no through SPR move can be found that improves the likelihood further. A classic 
     hill climbing algo.
  */

  Thorough = 1;
  impr = 1;
  
  /* restore the currently best tree. this si actually required, because we do not know which tree
     is actually stored in the tree data structure when the above loop exits */

  recallBestTree(bestT, 1, tr); 
  
  /* RE-TRAVERSE THE ENTIRE TREE */
  
  evaluateGeneric(tr, tr->start, TRUE);
#ifdef _DEBUG_CHECKPOINTING
  printBothOpen("After Fast SPRs Final %f\n", tr->likelihood);   
#endif
    
  /* optimize model params (including branch lengths) or just 
     optimize branch lengths and leave the other model parameters (GTR rates, alhpa) 
     alone */

  if(estimateModel)
    modOpt(tr, 1.0, adef, 0);
  else
    treeEvaluate(tr, 1.0);

  /* start loop that executes thorough SPR cycles */

  while(1)
    {	 
      /* once again if we want to restart from a checkpoint that was written during this loop we need
	 to restore the values of the variables appropriately */
    START_SLOW_SPRS:
      if(adef->useCheckpoint && ckp.state == SLOW_SPRS)
	{
	  optimizeRateCategoryInvocations = ckp.optimizeRateCategoryInvocations;   
      
	

  
	  impr = ckp.impr;
	  Thorough = ckp.Thorough;
	  bestTrav = ckp.bestTrav;
	  treeVectorLength = ckp.treeVectorLength;
	  rearrangementsMax = ckp.rearrangementsMax;
	  rearrangementsMin = ckp.rearrangementsMin;
	  thoroughIterations = ckp.thoroughIterations;
	  fastIterations = ckp.fastIterations;
   
  
	  lh = ckp.lh;
	  previousLh = ckp.previousLh;
	  difference = ckp.difference;
	  epsilon    = ckp.epsilon;                    
           
	
	  tr->likelihood = ckp.tr_likelihood;
              
	  tr->lhCutoff = ckp.tr_lhCutoff;
	  tr->lhAVG    = ckp.tr_lhAVG;
	  tr->lhDEC    = ckp.tr_lhDEC;   	 
	  tr->itCount = ckp.tr_itCount;	 

	  adef->useCheckpoint = FALSE;
	}
      else
	/* otherwise we restore the currently best tree and load it from bestT into our tree data 
	   structuire tr */
	recallBestTree(bestT, 1, tr);

      /* now, we write a checkpoint */
      /* Andre I believe that the code below, except for
	 writeCheckpoint cann still only be executed by process 0
	 => Andre: see above */
	{              
	  ckp.state = SLOW_SPRS;  
	  ckp.optimizeRateCategoryInvocations = optimizeRateCategoryInvocations;              
	  
	  
	  ckp.impr = impr;
	  ckp.Thorough = Thorough;
	  ckp.bestTrav = bestTrav;
	  ckp.treeVectorLength = treeVectorLength;
	  ckp.rearrangementsMax = rearrangementsMax;
	  ckp.rearrangementsMin = rearrangementsMin;
	  ckp.thoroughIterations = thoroughIterations;
	  ckp.fastIterations = fastIterations;
	  
	  
	  ckp.lh = lh;
	  ckp.previousLh = previousLh;
	  ckp.difference = difference;
	  ckp.epsilon    = epsilon; 
	  
	  
	  ckp.bestTrav = bestTrav;       
	  ckp.impr = impr;
	  
	  ckp.tr_startLH  = tr->startLH;
	  ckp.tr_endLH    = tr->endLH;
	  ckp.tr_likelihood = tr->likelihood;
	  ckp.tr_bestOfNode = tr->bestOfNode;
	  
	  ckp.tr_lhCutoff = tr->lhCutoff;
	  ckp.tr_lhAVG    = tr->lhAVG;
	  ckp.tr_lhDEC    = tr->lhDEC;     
	  ckp.tr_itCount  = tr->itCount;	
	  
	  /* write binary checkpoint to file */
	  
	  writeCheckpoint(tr, adef); 
	}
    
      if(impr)
	{	    
	  /* if the logl has improved write out some stuff and adapt the rearrangement radii */
	  printResult(tr, adef, FALSE);
	  /* minimum rearrangement radius */
	  rearrangementsMin = 1;
	  /* max radius, this is probably something I need to explain at the whiteboard */
	  rearrangementsMax = adef->stepwidth;	
	 
	  /* once again the convergence criterion */

	  if(tr->searchConvergenceCriterion && processID == 0)
	    {
	      int 
		bCounter = 0;	   	
	      
	      char 
		*buffer = (char*)calloc(tr->treeStringLength, sizeof(char));   

	      if(thoroughIterations > 1)
		cleanupHashTable(tr->h, (thoroughIterations % 2));		
		
	      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, thoroughIterations % 2, BIPARTITIONS_RF, (branchInfo *)NULL,
				      &bCounter, 1, FALSE, FALSE);	    
	      	     	    	   
	      	      	
#ifdef _DEBUG_CHECKPOINTING		
	      printf("Storing tree in slot %d\n", thoroughIterations % 2);
#endif
	      
	      Tree2String(buffer, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);
	      
	      if(thoroughIterations % 2 == 0)	      
		memcpy(tr->tree0, buffer, tr->treeStringLength * sizeof(char));
	      else
		memcpy(tr->tree1, buffer, tr->treeStringLength * sizeof(char));	    
	      
	      free(buffer);	      

	      assert(bCounter == tr->mxtips - 3);

	      if(thoroughIterations > 0)
		{
		  double 
		    rrf = convergenceCriterion(tr->h, tr->mxtips);
		  
		  MPI_Bcast(&rrf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		  
		  if(rrf <= 0.01) /* 1% cutoff */
		    {
		      printBothOpen("ML search converged at thorough SPR cycle %d with stopping criterion\n", thoroughIterations);
		      printBothOpen("Relative Robinson-Foulds (RF) distance between respective best trees after one succseful SPR cycle: %f%s\n", rrf, "%");
		      goto cleanup;
		    }
		  else		    
		    printBothOpen("ML search convergence criterion thorough cycle %d->%d Relative Robinson-Foulds %f\n", thoroughIterations - 1, thoroughIterations, rrf);
		}
	    }
	  
	  if(tr->searchConvergenceCriterion && processID != 0 && thoroughIterations > 0)
	    {
	      double 
		rrf;
	      
	      MPI_Bcast(&rrf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	      
	      if(rrf <= 0.01) /* 1% cutoff */		   
		goto cleanup;	      
	    }

	 
	   	  
	  thoroughIterations++;	  
	}			  			
      else
	{
	  /* if the lnl has not imrpved by the current SPR cycle adapt the min and max rearrangemnt radii and try again */
		       	   
	  rearrangementsMax += adef->stepwidth;
	  rearrangementsMin += adef->stepwidth; 	        	      

	  /* if we have already tried them then abandon this loop, the search has converged */
	  if(rearrangementsMax > adef->max_rearrange)	     	     	 
	    goto cleanup; 	   
	}
      
      /* optimize branch lengths of best tree */

      treeEvaluate(tr, 1.0);
     
      /* do some bokkeeping and printouts again */
      previousLh = lh = tr->likelihood;	      
      saveBestTree(bestT, tr, TRUE);     
      printLog(tr);

      /* do a cycle of thorough SPR moves with the minimum and maximum rearrangement radii */

      treeOptimizeRapid(tr, rearrangementsMin, rearrangementsMax, adef, bt, bestML);
	
      impr = 0;			      		            

      /* once again get the best 20 trees produced by the SPR cycle, load them from the bt tree list into tr
	 optimize their branch lengths and figure out if the LnL of the tree has improved */

      for(i = 1; i <= bt->nvalid; i++)
	{		 
	  recallBestTree(bt, i, tr);	 	    	    	
	  
	  treeEvaluate(tr, 0.25);	    	 
	    
	  difference = ((tr->likelihood > previousLh)? 
			tr->likelihood - previousLh: 
			previousLh - tr->likelihood); 	    
	  if(tr->likelihood > lh && difference > epsilon)
	    {
	      impr = 1;	       
	      lh = tr->likelihood;	  	     
	      saveBestTree(bestT, tr, TRUE);
	    }	   	   
	}  

#ifdef _DEBUG_CHECKPOINTING
      printBothOpen("SLOW LH: %f\n", lh);              
#endif
    }

 cleanup: 
  
  /* do a final full tree traversal, not sure if this is required here */
  
  evaluateGeneric(tr, tr->start, TRUE);
    
#ifdef _DEBUG_CHECKPOINTING
  printBothOpen("After SLOW SPRs Final %f\n", tr->likelihood);   
#endif
   
  printBothOpen("\nLikelihood of best tree: %f\n", tr->likelihood);
  /* print the absolut best tree */

  printLog(tr);
  printResult(tr, adef, TRUE);

  /* print other good trees encountered during the search */

  if(tr->saveBestTrees > 0)
    { 
      char 
	fileName[2048] = "",
	buf[64] = "";
     
      printBothOpen("\n\nEvaluating %d other good ML trees\n\n", bestML->nvalid);
      
      for(i = 1; i <= bestML->nvalid; i++)
	{		 
	  recallBestTree(bestML, i, tr);	 	    	    		  
	  /*treeEvaluate(tr, 0.25);*/
	  printBothOpen("tree %d likelihood %1.80f\n", i, tr->likelihood);

	  if(processID == 0)
	    { 	      		
	      FILE 
		*treeFile;

	      strcpy(fileName,       workdir);
	      strcat(fileName, "RAxML_");
	      sprintf(buf, "%d", bestML->nvalid);
	      strcat(fileName, buf);
	      strcat(fileName, "_goodTrees.");
	      strcat(fileName, run_id);

	      treeFile = myfopen(fileName, "a");
	     
	      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, SUMMARIZE_LH, FALSE, FALSE);
	       
	      fprintf(treeFile, "%s", tr->tree_string);
	      fclose(treeFile);	      
	    }

	  
	}      
	
      printBothOpen("\n\nOther good trees written to file %s\n", fileName);			
    }


  /* free data structures */

  if(tr->searchConvergenceCriterion && processID == 0)
    {
      freeBitVectors(tr->bitVectors, 2 * tr->mxtips);
      free(tr->bitVectors);
      freeHashTable(tr->h);
      free(tr->h);
    }
  
  freeBestTree(bestT);
  free(bestT);
  freeBestTree(bt);
  free(bt);
  freeInfoList();  
  

  /* and we are done, return to main() in axml.c  */

}



boolean treeEvaluate (tree *tr, double smoothFactor)       /* Evaluate a user tree */
{
  boolean result;

 
  result = smoothTree(tr, (int)((double)smoothings * smoothFactor));
  
  assert(result); 

  //make sure that all vectors are oriented correctly !

  evaluateGeneric(tr, tr->start, TRUE);   
    

  return TRUE;
}


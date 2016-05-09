/*  RAxML-HPC, a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright March 2006 by Alexandros Stamatakis
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
 *  stamatak@ics.forth.gr
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *  
 *  Alexandros Stamatakis: "An Efficient Program for phylogenetic Inference Using Simulated Annealing". 
 *  Proceedings of IPDPS2005,  Denver, Colorado, April 2005.
 *  
 *  AND
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>
#include <sys/types.h>
#include "axml.h"

extern double masterTime;
extern char workdir[1024];
extern char run_id[128];
extern char quartetGroupingFileName[1024];
extern char quartetFileName[1024];
extern checkPointState ckp;
extern int processID;

/* a parser error function */

static void parseError(int c)
{
  printBothOpen("Quartet grouping parser expecting symbol: %c\n", c);
  assert(0);
}

/* parser for the taxon grouping format, one has to specify 4 groups in a newick-like 
   format from which quartets (a substantially smaller number compared to ungrouped quartets) 
   will be drawn */

static void groupingParser(char *quartetGroupFileName, int *groups[4], int groupSize[4], tree *tr)
{
  FILE 
    *f = myfopen(quartetGroupFileName, "r");
  
  int 
    taxonCounter = 0,
    n,
    state = 0,
    groupCounter = 0,
    ch,
    i;

  printBothOpen("%s\n", quartetGroupFileName);

  for(i = 0; i < 4; i++)
    {
      groups[i] = (int*)malloc(sizeof(int) * (tr->mxtips + 1));
      groupSize[i] = 0;
    }
  
  while((ch = getc(f)) != EOF)
    {
      if(!whitechar(ch))
	{
	  switch(state)
	    {
	    case 0:
	      if(ch != '(')
		parseError('(');
	      state = 1;
	      break;
	    case 1:
	      ungetc(ch, f);
	      n = treeFindTipName(f, tr, FALSE);  
	      if(n <= 0 || n > tr->mxtips)		
		printBothOpen("parsing error, raxml is expecting to read a taxon name, found \"%c\" instead\n", ch);		
	      assert(n > 0 && n <= tr->mxtips);	     
	      taxonCounter++;
	      groups[groupCounter][groupSize[groupCounter]] = n;
	      groupSize[groupCounter] = groupSize[groupCounter] + 1;	    
	      state = 2;
	      break;
	    case 2:
	      if(ch == ',')
		state = 1;
	      else
		{
		  if(ch == ')')
		    {
		      groupCounter++;
		      state = 3;
		    }
		  else
		    parseError('?');
		}
	      break;
	    case 3:
	      if(groupCounter == 4)
		{
		  if(ch == ';')
		    state = 4;
		  else
		    parseError(';');
		}
	      else
		{
		  if(ch != ',')
		    parseError(',');
		  state = 0;
		}
	      break; 
	    case 4:
	      printBothOpen("Error: extra char after ; %c\n", ch);
	      assert(0);
	    default:
	      assert(0);
	    }
	}
    }

  assert(state == 4);
  assert(groupCounter == 4); 
  assert(taxonCounter == tr->mxtips);

  printBothOpen("Successfully parsed quartet groups\n\n");

  /* print out the taxa that have been assigned to the 4 groups */

  for(i = 0; i < 4; i++)
    {
      int 
	j;
      
      printBothOpen("group %d has %d members\n", i, groupSize[i]);

      for(j = 0; j < groupSize[i]; j++)
	printBothOpen("%s\n", tr->nameList[groups[i][j]]);

      printBothOpen("\n");
    }

  fclose(f);
}

/*****************************/

static void nniSmooth(tree *tr, nodeptr p, int maxtimes)
{
  int
    i;

  for(i = 0; i < tr->numBranches; i++)	
    tr->partitionConverged[i] = FALSE;	
 
  while (--maxtimes >= 0) 
    {           
      for(i = 0; i < tr->numBranches; i++)	
	tr->partitionSmoothed[i] = TRUE;
            
      assert(!isTip(p->number, tr->mxtips)); 	
      assert(!isTip(p->back->number, tr->mxtips));  
      
      update(tr, p);
     
      update(tr, p->next);
     
      update(tr, p->next->next);
      
      update(tr, p->back->next);
      
      update(tr, p->back->next->next);           
     
      if (allSmoothed(tr)) 
	break;      
    }

  for(i = 0; i < tr->numBranches; i++)
    {
      tr->partitionSmoothed[i] = FALSE; 
      tr->partitionConverged[i] = FALSE;
    }  
}





static double quartetLikelihood(tree *tr, nodeptr p1, nodeptr p2, nodeptr p3, nodeptr p4, nodeptr q1, nodeptr q2, analdef *adef, boolean firstQuartet)
{
  /* 
     build a quartet tree, where q1 and q2 are the inner nodes and p1, p2, p3, p4
     are the tips of the quartet where the sequence data is located.

     initially set all branch lengths to the default value.
  */

  /* 
     for the tree and node data structure used, please see one of the last chapter's of Joe 
     Felsensteins book. 
  */

  hookupDefault(q1, q2, tr->numBranches);
  
  hookupDefault(q1->next,       p1, tr->numBranches);
  hookupDefault(q1->next->next, p2, tr->numBranches);
  
  hookupDefault(q2->next,       p3, tr->numBranches);
  hookupDefault(q2->next->next, p4, tr->numBranches);
  
  /* now compute the likelihood vectors at the two inner nodes of the tree,
     here the virtual root is located between the two inner nodes q1 and q2.
  */

  newviewGeneric(tr, q1, FALSE);
  newviewGeneric(tr, q2, FALSE);
  
  /* call a function that is also used for NNIs that iteratively optimizes all 
     5 branch lengths in the tree.

     Note that 16 is an important tuning parameter, this integer value determines 
     how many times we visit all branches until we give up further optimizing the branch length 
     configuration.
  */

  nniSmooth(tr, q1, 16);

  /* now compute the log likelihood of the tree for the virtual root located between inner nodes q1 and q2 */
  
  /* debugging code 
     {
    double l;
  */
  
  evaluateGeneric(tr, q1->back->next->next, FALSE);
  
  /* debugging code 
     
     l = tr->likelihood;

     newviewGeneric(tr, q1);
     newviewGeneric(tr, q2);
     evaluateGeneric(tr, q1);
     
   
     assert(ABS(l - tr->likelihood) < 0.00001);
     }
  */

  return (tr->likelihood);
}



static void computeAllThreeQuartets(tree *tr, nodeptr q1, nodeptr q2, int t1, int t2, int t3, int t4, FILE *f, analdef *adef)
{
  /* set the tip nodes to different sequences 
     with the tip indices t1, t2, t3, t4 */
	       
  nodeptr 
    p1 = tr->nodep[t1],
    p2 = tr->nodep[t2],
    p3 = tr->nodep[t3], 
    p4 = tr->nodep[t4];
  
  double 
    l;
  
  /* first quartet */	    
  
  /* compute the likelihood of tree ((p1, p2), (p3, p4)) */
  
  l = quartetLikelihood(tr, p1, p2, p3, p4, q1, q2, adef, TRUE);
 
  if(processID == 0)
    fprintf(f, "%d %d | %d %d: %f\n", p1->number, p2->number, p3->number, p4->number, l);

  /* second quartet */	    
  
  /* compute the likelihood of tree ((p1, p3), (p2, p4)) */
  
  l = quartetLikelihood(tr, p1, p3, p2, p4, q1, q2, adef, FALSE);

  if(processID == 0)
    fprintf(f, "%d %d | %d %d: %f\n", p1->number, p3->number, p2->number, p4->number, l);

  /* third quartet */	    
  
  /* compute the likelihood of tree ((p1, p4), (p2, p3)) */
  
  l = quartetLikelihood(tr, p1, p4, p2, p3, q1, q2, adef, FALSE);

  if(processID == 0)
    fprintf(f, "%d %d | %d %d: %f\n", p1->number, p4->number, p2->number, p3->number, l);	    	   
}

/* the three quartet options: all quartets, randomly sub-sample a certain number n of quartets, 
   subsample all quartets from 4 pre-defined groups of quartets */


static void writeQuartetCheckpoint(uint64_t quartetCounter, FILE *f, tree *tr, analdef *adef)
{
  if(quartetCounter % adef->quartetCkpInterval == 0)
    {     
      ckp.quartetCounter = quartetCounter;
      if(processID == 0)
        { 
          fflush(f);
          ckp.filePosition = ftell(f);
        }
      printBothOpen("\nPrinting checkpoint after %f seconds of run-time\n", gettime() - masterTime);      
      writeCheckpoint(tr, adef);      
    }
}


#define ALL_QUARTETS 0
#define RANDOM_QUARTETS 1
#define GROUPED_QUARTETS 2

void computeQuartets(tree *tr, analdef *adef)
{
  /* some indices for generating quartets in an arbitrary way */

  int
    flavor = ALL_QUARTETS, //type of quartet calculation 
    i, 
    t1, 
    t2, 
    t3, 
    t4, 
    *groups[4],
    groupSize[4];    

  double
    fraction = 0.0;

  uint64_t
    randomQuartets = (uint64_t)(adef->numberRandomQuartets), //number of random quartets to compute 
    quartetCounter = 0, 
    //total number of possible quartets, note that we count the following ((A,B),(C,D)), ((A,C),(B,D)), ((A,D),(B,C)) as one quartet here 
    numberOfQuartets = ((uint64_t)tr->mxtips * ((uint64_t)tr->mxtips - 1) * ((uint64_t)tr->mxtips - 2) * ((uint64_t)tr->mxtips - 3)) / 24; 
  
  /* use two inner tree nodes for building quartet trees */

  nodeptr 	
    q1 = tr->nodep[tr->mxtips + 1],
    q2 = tr->nodep[tr->mxtips + 2];

  FILE 
    *f;

  long
    seed = (long)(tr->randomSeed);
       
  /***********************************/  

  /* get a starting tree on which we optimize the likelihood model parameters: either reads in a tree or computes a randomized stepwise addition parsimony tree */
  if(adef->useCheckpoint)
    { 
      /* read checkpoint file */
      restart(tr, adef);
      strncpy(quartetFileName, ckp.quartetFileName, 1024);
      printBothOpen("Time for reading checkpoint file: %f\n\n", gettime() - masterTime); 

      seed = ckp.seed;
   
      if(processID == 0)
        {   
           f = myfopen(quartetFileName, "r+");
        
           fseek(f, ckp.filePosition, SEEK_SET);  
           if(ftruncate(fileno(f),  ckp.filePosition) != 0)
  	    assert(0);
        }
    }
  else
    {
      getStartingTree(tr);
      evaluateGeneric(tr, tr->start, TRUE);
      treeEvaluate(tr, 1);

      /* optimize model parameters on that comprehensive tree that can subsequently be used for evaluation of quartet likelihoods */

      modOpt(tr, adef->likelihoodEpsilon, adef, 0);

      printBothOpen("Time for parsing input tree or building parsimony tree and optimizing model parameters: %f\n\n", gettime() - masterTime); 
      printBothOpen("Tree likelihood: %f\n\n", tr->likelihood);  

      if(processID == 0)
        f = myfopen(quartetFileName, "w");
    }

  /* figure out which flavor of quartets we want to compute */

  if(adef->useQuartetGrouping)
    {
      //quartet grouping evaluates all possible quartets from four disjoint 
      //sets of user-specified taxon names 

      flavor = GROUPED_QUARTETS;
      
      //parse the four disjoint sets of taxon names specified by the user from file      
      groupingParser(quartetGroupingFileName, groups, groupSize, tr);
    }
  else
    {
      //if the user specified more random quartets to sample than there actually 
      //exist for the number of taxa, then fix this.

      if(randomQuartets == 0 || randomQuartets >= numberOfQuartets)
	//TODO add usre warning? if second case above true? 
	flavor = ALL_QUARTETS;
      else
	{     	  
	  //compute the fraction of random quartets to sample 
	  //there may be an issue here with the unit64_t <-> double cast
	  fraction = (double)randomQuartets / (double)numberOfQuartets;      
	  assert(fraction < 1.0);
	  flavor = RANDOM_QUARTETS;
	}
    }

  ckp.state = QUARTETS;
  ckp.seed = seed; 
  strncpy(ckp.quartetFileName, quartetFileName, 1024);

  /* print some output on what we are doing*/

  switch(flavor)
    {
    case ALL_QUARTETS:
      printBothOpen("There are %" PRIu64 " quartet sets for which RAxML will evaluate all %" PRIu64 " quartet trees\n", numberOfQuartets, numberOfQuartets * 3);
      break;
    case RANDOM_QUARTETS:
      printBothOpen("There are %" PRIu64 " quartet sets for which RAxML will randomly sub-sambple %" PRIu64 " sets (%f per cent), i.e., compute %" PRIu64 " quartet trees\n", 
		    numberOfQuartets, randomQuartets, 100 * fraction, randomQuartets * 3);
      break;
    case GROUPED_QUARTETS:  
      printBothOpen("There are 4 quartet groups from which RAxML will evaluate all %u quartet trees\n", 
		    (unsigned int)groupSize[0] * (unsigned int)groupSize[1] * (unsigned int)groupSize[2] * (unsigned int)groupSize[3] * 3);
      break;
    default:
      assert(0);
    }

  /* print taxon name to taxon number correspondance table to output file */

  if(!adef->useCheckpoint)
    {
      if(processID == 0)
	fprintf(f, "Taxon names and indices:\n\n");
      
      for(i = 1; i <= tr->mxtips; i++)
	{
	  if(processID == 0)
	    fprintf(f, "%s %d\n", tr->nameList[i], i);
	  assert(tr->nodep[i]->number == i);
	}
      
      if(processID == 0)
        {
	  fprintf(f, "\n\n");
          fflush(f);
        }
    }
  
  
  /* do a loop to generate some quartets to test.
     note that tip nodes/sequences in RAxML are indexed from 1,...,n
     and not from 0,...,n-1 as one might expect 
     
     tr->mxtips is the maximum number of tips in the alignment/tree
  */


  //now do the respective quartet evaluations by switching over the three distinct flavors 
     
  switch(flavor)
    {
    case ALL_QUARTETS:
      {		    
	/* compute all possible quartets */	   	   
	
	for(t1 = 1; t1 <= tr->mxtips; t1++)
	  for(t2 = t1 + 1; t2 <= tr->mxtips; t2++)
	    for(t3 = t2 + 1; t3 <= tr->mxtips; t3++)
	      for(t4 = t3 + 1; t4 <= tr->mxtips; t4++)
		{		      
		  if((adef->useCheckpoint && quartetCounter >= ckp.quartetCounter) || !adef->useCheckpoint)
		    {
		      writeQuartetCheckpoint(quartetCounter, f, tr, adef);		      
		      
		      computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, f, adef);
		    }
		  quartetCounter++;
		}
	
	assert(quartetCounter == numberOfQuartets);
      }
      break;
    case RANDOM_QUARTETS:
      {	 
	
	//endless loop ta make sure we randomly sub-sample exactly as many quartets as the user specified
	
	//This is not very elegant, but it works, note however, that especially when the number of 
	//random quartets to be sampled is large, that is, close to the total number of quartets 
	//some quartets may be sampled twice by pure chance. To randomly sample unique quartets 
	//using hashes or bitmaps to store which quartets have already been sampled is not memory efficient.
	//Insetad, we need to use a random number generator that can generate a unique series of random numbers 
	//and then have a function f() that maps those random numbers to the corresponding index quartet (t1, t2, t3, t4).
	
	do
	  {	      
	    //loop over all quartets 
	    for(t1 = 1; t1 <= tr->mxtips; t1++)
	      for(t2 = t1 + 1; t2 <= tr->mxtips; t2++)
		for(t3 = t2 + 1; t3 <= tr->mxtips; t3++)
		  for(t4 = t3 + 1; t4 <= tr->mxtips; t4++)
		    {
		      //chose a random number
		      double
			r = randum(&seed);
			  			  
		      //if the random number is smaller than the fraction of quartets to subsample
		      //evaluate the likelihood of the current quartet
		      if(r < fraction)
			{
			  if((adef->useCheckpoint && quartetCounter >= ckp.quartetCounter) || !adef->useCheckpoint)
			    {
			      writeQuartetCheckpoint(quartetCounter, f, tr, adef);			     			      
			      
			      //function that computes the likelihood for all three possible unrooted trees 
			      //defined by the given quartet of taxa 
			      computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, f, adef);
			    }
			  //increment quartet counter that counts how many quartets we have evaluated
			  quartetCounter++;
			}
		      
		      //exit endless loop if we have randomly sub-sampled as many quartets as the user specified
		      if(quartetCounter == randomQuartets)
			goto DONE;
		    }
	  }
	while(1);
	
      DONE:
	assert(quartetCounter == randomQuartets);	  
      }
      break;
    case GROUPED_QUARTETS:
      {
	/* compute all quartets that can be built out of the four pre-defined groups */
	
	for(t1 = 0; t1 < groupSize[0]; t1++)
	  for(t2 = 0; t2 < groupSize[1]; t2++)
	    for(t3 = 0; t3 < groupSize[2]; t3++)
	      for(t4 = 0; t4 < groupSize[3]; t4++)
		{
		  int
		    i1 = groups[0][t1],
		    i2 = groups[1][t2],
		    i3 = groups[2][t3],
		    i4 = groups[3][t4];
		  
		  if((adef->useCheckpoint && quartetCounter >= ckp.quartetCounter) || !adef->useCheckpoint)
		    {
		      writeQuartetCheckpoint(quartetCounter, f, tr, adef);
		      		     
		      computeAllThreeQuartets(tr, q1, q2, i1, i2, i3, i4, f, adef);
		    }
		  quartetCounter++;
		}
	
	printBothOpen("\nComputed all %" PRIu64 " possible grouped quartets\n", quartetCounter); 	    
      }
      break;
    default:
      assert(0);
    }
  
  fclose(f);
}

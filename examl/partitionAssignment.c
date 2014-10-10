#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "partitionAssignment.h" 

extern int processID; 


void initializePartitionAssignment( PartitionAssignment **pAssPtr, pInfo **partitions, int numPart, int numProc)
{
  int 
    i; 
  
  PartitionAssignment 
    *pAss;

  *pAssPtr = (PartitionAssignment*)calloc(1, sizeof(PartitionAssignment)); 
  
  pAss = *pAssPtr; 

  pAss->numProc = numProc; 
  pAss->numPartitions = numPart; 

  pAss->partitions = (Partition *)calloc((size_t)pAss->numPartitions, sizeof(Partition)); 
  
  for(i = 0; i < numPart; ++i)
    {
      Partition 
	*p = pAss->partitions + i; 
      p->id = i; 
      p->width = partitions[i]->upper - partitions[i]->lower;
      p->type = partitions[i]->states;
    }
  
  pAss->assignPerProc = (Assignment **)calloc((size_t)pAss->numProc , sizeof(Assignment*)); 
  pAss->numAssignPerProc = (int *)calloc((size_t)pAss->numProc, sizeof(int)); 
}


void deletePartitionAssignment(PartitionAssignment *pAss)
{
  int 
    i; 

  free(pAss->partitions); 
  for(i = 0; i < pAss->numProc; ++i)
    free(pAss->assignPerProc[i]);
  free(pAss->assignPerProc); 
  free(pAss);
}


static int partSort(const void *a, const void *b ) 
{
  return  ((const Partition*)a)->width - ((const Partition*) b)->width ; 
}


/** 
    helper function that executes the assignment of a partial assignment (only numElem character are assigned)
 */ 
static void assignPartitionPartial(PartitionAssignment *pa, Partition* p, int procId, int *numAssigned, size_t *sizeAssigned, size_t offset, size_t numElem)
{
  Assignment 
    *a;
  
  int 
    newArrayLen; 
 
  ++numAssigned[procId]; 
  ++pa->numAssignPerProc[procId]; 

  newArrayLen = pa->numAssignPerProc[procId]; 

  pa->assignPerProc[procId] = (Assignment*)realloc(pa->assignPerProc[procId], newArrayLen * sizeof(Assignment)); 
  
  a = pa->assignPerProc[procId] + (newArrayLen-1); 

  a->offset = offset; 
  a->partId = p->id; 
  a->width = numElem; 
  sizeAssigned[procId] += numElem; 
}


/** 
    helper function that executes the assignment of a full partition 
 */ 
static void assignPartitionFull(PartitionAssignment* pa, Partition* p, int procId, int *numAssigned, size_t *sizeAssigned)
{
  assignPartitionPartial(pa, p, procId, numAssigned, sizeAssigned, 0, p->width);
} 


/** 
    Request a process to which a part of the partition should
    be assigned to. 

    At this stage there are processes that have one more partition
    assignment than others. These have been categorized into two
    stacks (high and low).  For each stack, we have in iter variable
    (of type int**) that gets decremented, if an element is removed
    from the stack and a the start of the array that allows us to
    determine when the stack is empty.

    popAndYield tries to satisfy the request for a process with more
    or less assignments (see wantLow), but will return any process, if
    the request cannot be fulfilled.

    If both queues are empty, it will return -1.
 */ 
static int popAndYield(int **procsHighIter, int *procsHighStart, int **procsLowIter, int *procsLowStart, boolean wantLow)
{
  boolean 
    fromHigh = FALSE, 
    fromLow = FALSE;

  int 
    result = -1; 
  
  if(wantLow)
    {
      if(*procsLowIter - procsLowStart > 0) 
	fromLow = TRUE; 
      else 
	if(*procsHighIter - procsHighStart > 0)
	  fromHigh = TRUE; 
    }
  else 
    {
      if(*procsHighIter - procsHighStart > 0 )
	fromHigh = TRUE; 
      else 
	if(*procsLowIter - procsLowStart > 0)
	  fromLow = TRUE; 
    }

  if(fromHigh) 
    {
      result = **procsHighIter;
      --(*procsHighIter); 
    }
  else 
    if(fromLow)
      {
	result = **procsLowIter; 
	--(*procsLowIter); 
      }

  return result; 
}



static void assignThesePartitions(PartitionAssignment* pa, Partition *partitions, int numCur)
{
  int
    proc, 
    remainder,			/* number of processes that receive 1 character less than other s */
    i,
    numFull = 0,		/* number of processes that cannot take any more   */
    numLow = 0,			/*  */
    *numAssigned = (int *)NULL,	/* number of characters assigned to a process */
    *procsHighIter = (int *)NULL, /* stack of processes that have one assignment more than others  */
    *procsLowIter = (int *)NULL,  /* stack of processes that have one assignment less than others  */
    *procsLowStart = (int *)NULL, /* start of stack */
    *procsHighStart = (int *)NULL, /* start of stack */
    highProc,			   /* id of a process that potentially has more partitions assigned than others   */
    lowProc;			   /* id of a process that potentially has less partitions assigned than others   */
  
  size_t
    totalElems = 0,
    *sizeAssigned = (size_t *)NULL,
    toAdd,
    cap; 			/* defines a cap: once we have
				   assigned this many characters, the
				   remaining processes will get one
				   character less */
  
  boolean 
    iterate = TRUE; 
  
  Partition 
    *partIter = partitions,
    *partEnd = partitions +numCur; 

  /* The following implements Kassian's algorithm. Originally, his
     algorithm consists of 5 phases */

  /*
    Sorts partitions according to their size. According to Kassians
    algorithm, this step is NOT obligatory. However, if we do it, then
    phase 3 (called "top-up" phase) is not necessary (this has been
    clarified with Kassian).
  */ 
  qsort(partitions, numCur, sizeof(Partition), partSort);

  for(i = 0; i < numCur; ++i)
    totalElems += partitions[i].width; 

  cap = ceil( (double)totalElems / (double)pa->numProc ); 
  
  remainder = cap * pa->numProc - totalElems; 
  
  assert(remainder >= 0 ); 
  
  numAssigned = (int *)calloc((size_t)pa->numProc, sizeof(int));
  
  sizeAssigned = (size_t *)calloc((size_t)pa->numProc, sizeof(size_t)); 
  
  /* phase 2: initial distribution of full partitions to procesess. We
     distribute full partitions until for the first time, we cannot
     assign an entire partition any more, because this would exceed
     the number of characters we want to assign to this process */
  while(iterate)
    {           
      for(proc = 0; proc < pa->numProc;++proc)
	{
	  if(partIter < partEnd && sizeAssigned[proc] + partIter->width <= cap)
	    {
	      assignPartitionFull(pa, partIter, proc, numAssigned, sizeAssigned); 
	      
	      if(sizeAssigned[proc] == cap)
		{
		  ++numFull;
		  if(numFull == pa->numProc - remainder)
		    --cap; 
		}
	      
	      ++partIter; 
	    }
	  else 
	    {
	      numLow = numAssigned[proc]; 
	      iterate = FALSE; 
	      break; 
	    }
	}
    }

  /* phase 3: top-up => not necessary because of previous sorting */

  
  /* 
     phase 4: stick breaking

     Here we partially assign the remaining partitions to processes
     until every process has as many characters as it can take.
  */

  
  /* first categorize processes into two stacks, dependent on whether
     they have gotton one more partition than others */
  procsHighIter =  (int*)calloc((size_t)pa->numProc + 1, sizeof(int)); 
  procsLowIter = (int *)calloc((size_t)pa->numProc + 1, sizeof(int)); 
  procsLowStart = procsLowIter; 
  procsHighStart = procsHighIter;
    

  numFull = 0; 
  
  for(proc = 0; proc < pa->numProc; ++proc)
    {
      if(sizeAssigned[proc] < cap)
 	{
 	  if(numAssigned[proc] == numLow)
 	    {
 	      ++procsLowIter; 
 	      *procsLowIter = proc; 
 	    }
 	  else  
 	    {
 	      ++procsHighIter;
 	      *procsHighIter = proc; 
 	    }
 	}
      else 	
	++numFull; 
    }
  
  assert((procsHighIter - procsHighStart) + (procsLowIter - procsLowStart) + numFull == pa->numProc); 

  toAdd = (partIter < partEnd) ? partIter->width : 0  ; 
  highProc = popAndYield(&procsHighIter, procsHighStart, &procsLowIter, procsLowStart, FALSE); 
  lowProc = popAndYield(&procsHighIter, procsHighStart, &procsLowIter, procsLowStart, TRUE); 


  /* 
     now assign as long as there is something to assign. Once both
     stacks are empty, popAndYield yields -1. This then breaks the
     loop condition here.
   */ 
  while(  ! (highProc == -1 && lowProc == -1 
	     && (procsHighIter - procsHighStart <= 0 ) && (procsLowIter - procsLowStart <= 0 )) )
    {
      /* try to finish a assignments for a process that has many partitions  */
      if(highProc != -1 && sizeAssigned[highProc] + toAdd >= cap)
	{ 
	  size_t
	    toTransfer = cap - sizeAssigned[highProc],
	    offset = partIter->width - toAdd; 
	  
	  assignPartitionPartial( pa, partIter, highProc, numAssigned, sizeAssigned, offset, toTransfer);
	  
	  toAdd -= toTransfer; 

	  if(toAdd == 0 && partIter < partEnd ) 
	    {
	      ++partIter; 
	      toAdd = partIter < partEnd ?  partIter->width : 0; 
	    }
	  ++numFull; 
	  
	  if(numFull == pa->numProc - remainder)
	    --cap; 

	  highProc = popAndYield(&procsHighIter, procsHighStart, &procsLowIter, procsLowStart, FALSE); 
	}
      else 
	if(lowProc != -1)
	  {
	    /* assign the enitre remaining portion to a process that
	       still has fewer partitions */
	    if(sizeAssigned[lowProc] + toAdd < cap)
	      {
		size_t
		  offset = partIter->width - toAdd; 
		
		assignPartitionPartial(pa, partIter, lowProc, numAssigned, sizeAssigned, offset, toAdd); 
		
		if(highProc != -1 )
		  {
		    ++procsHighIter; 
		    *procsHighIter = highProc; 
		  }
		
		highProc = lowProc; 
		
		toAdd = 0; 
		
		if( partIter != partEnd )
		  {
		    ++partIter;
		    toAdd = partIter < partEnd ? partIter->width : 0  ;
		  }
		
		lowProc = popAndYield(&procsHighIter, procsHighStart, &procsLowIter, procsLowStart, TRUE); 
	      }
	    else 
	      { 
		/* assign as much as possible to a process with less
		   partitions (the rest probably needs to be assigned
		   to the next process) */

		size_t
		  toTransfer = cap - sizeAssigned[lowProc],
		  offset = partIter->width - toAdd; 
		
		assignPartitionPartial(pa, partIter, lowProc, numAssigned, sizeAssigned, offset, toTransfer); 
		
		toAdd -= toTransfer;
		
		if(toAdd == 0 && partIter < partEnd)
		  {
		    ++partIter; 
		    toAdd = partIter < partEnd ? partIter->width : 0; 
		  }
		
		++numFull;
		if(numFull == pa->numProc - remainder)
		  --cap ;
		
		lowProc = popAndYield(&procsHighIter, procsHighStart, &procsLowIter, procsLowStart, FALSE); 
	      }
	  }
	else 
	  {
	    /* should not occurr, but I am not entirely happly with
	       this assert. */
	    assert(0);
	  }
    }

  assert(toAdd == 0 ); 
  assert(partIter == partEnd); 
  
  free(numAssigned); 
  free(sizeAssigned); 
}

/** 
    Assigns all partitions. Notice that for each data type (currently
    only AA and DNA), we execute the algorithm separately. Thus, in
    the worst case imbalances for AA and DNA could hit the same
    processes (but this is probably not worth bothering). 
 */ 
void assign(PartitionAssignment *pa)
{
  int 
    partitionsHandled = 0,
    curType = -1,
    i; 

  while(partitionsHandled < pa->numPartitions)
    {
      size_t 
	cnt; 
      
      Partition 
	*curPartitions = (Partition *)NULL; 

      /* search next type */
      for( i = 0; i < pa->numPartitions; ++i)
	{
	  if(curType < pa->partitions[i].type )
	    {
	      curType = pa->partitions[i].type; 
	      break; 
	    }
	}
      
      /* count number of type */
      cnt = 0; 
      for(i = 0; i < pa->numPartitions; ++i)      
	{
	  if(pa->partitions[i].type == curType)
	    ++cnt; 
	}

      curPartitions = (Partition*)calloc((size_t)cnt, sizeof(Partition)); 
      
      cnt = 0; 
      for(i = 0; i< pa->numPartitions; ++i)
	{
	  if(pa->partitions[i].type == curType)
	    curPartitions[cnt++] = pa->partitions[i]; 
	}

      assignThesePartitions(pa, curPartitions, cnt);
      free(curPartitions); 
      
      partitionsHandled += cnt; 
    }
}




void printAssignment(Assignment a, int procid)
{
  printf("p: %d\t(%lu,%lu) -> proc %d\n", a.partId, a.offset, a.width , procid); 
} 


void printAssignments(PartitionAssignment *pa)
{
  int i,j; 
  printf("proc\toffset\tlength\tpart\n"); 
  for(i = 0; i < pa->numProc; ++i)
    {
      for(j = 0; j < pa->numAssignPerProc[i] ; ++j)
	{
	  Assignment a = pa->assignPerProc[i][j]; 
	  printf("%d\t%lu\t%lu\t%d\n", i,a.offset, a.width, a.partId); 
	}
    }
}


void printLoad(PartitionAssignment *pa)
{
  int 
    i,
    j, 
    *numsPerProc = (int *)calloc((size_t)pa->numProc, sizeof(int)); 
  
  size_t
    *sitesPerProc = (size_t *)calloc((size_t)pa->numProc, sizeof(size_t)); 

  for(i = 0; i< pa->numProc; ++i)
    {
      for(j = 0; j < pa->numAssignPerProc[i]; ++j)
	{
	  Assignment a = pa->assignPerProc[i][j]; 
	  sitesPerProc[i] += a.width; 
	  ++numsPerProc[i]; 
	}
    }

  printf("#proc\t#part\t#sites\n"); 
  for( i = 0; i < pa->numProc ; ++i)
    printf("%d\t%d\t%lu\n", i, numsPerProc[i], sitesPerProc[i]); 
}



/**
   allocates global arrays for CAT and sets the pointers in each pInfo
   instance.
 */ 
static void setupBasePointersInTree(tree *tr)
{
  size_t 
    len = 0;
  int
    i; 
  
  for(i = 0; i < tr->NumberOfModels; ++i)
    len += tr->partitionData[i].width ; 

  tr->patrat_basePtr = (double*) calloc((size_t)len, sizeof(double));
  tr->rateCategory_basePtr = (int*) calloc((size_t)len, sizeof(int)); 
  tr->lhs_basePtr = (double*) calloc((size_t)len, sizeof(double)); 

  len = 0; 
  for(i = 0; i < tr->NumberOfModels; ++i)
    {
      if(tr->partitionData[i].width > 0)
	{
	  tr->partitionData[i].rateCategory = tr->rateCategory_basePtr + len; 
	  tr->partitionData[i].patrat = tr->patrat_basePtr + len; 
	  tr->partitionData[i].lhs = tr->lhs_basePtr + len ; 
	}
      else 
	{
	  tr->partitionData[i].rateCategory = (int *)NULL; 
	  tr->partitionData[i].patrat = (double*)NULL; 
	  tr->partitionData[i].lhs = (double*)NULL;
	}

      len += tr->partitionData[i].width; 
    }
}



static int sortById(const void *a, const void *b)
{
  return ((Assign*) a)->partitionId -  ((Assign*) b)->partitionId ;
}



void copyAssignmentInfoToTree(PartitionAssignment *pa, tree *tr)
{
  int 
    i,
    numAssign = 0; 
  
  Assign 
    *assIter;
  
  for(i = 0; i < pa->numProc; ++i)
    numAssign += pa->numAssignPerProc[i]; 
  
  /* copy the partition assignment to the tree structure */

  tr->numAssignments  = numAssign; 
  tr->partAssigns = (Assign *)calloc((size_t)numAssign, sizeof(Assign)); 

  assIter = tr->partAssigns; 
  
  for(i = 0; i < pa->numProc; ++i)
    {
      int j; 
      for( j = 0 ; j < pa->numAssignPerProc[i]; ++j)
	{
	  Assignment *ass = pa->assignPerProc[i] + j; 
	  assIter->procId = i; 
	  assIter->offset = ass->offset; 
	  assIter->width = ass->width; 
	  assIter->partitionId = ass->partId; 
	  ++assIter; 
	}
    }

  /* 
     the sorting makes it easier to deal with the assignments later in
     case of a gather/scatter at the master. Thus, we do not need to
     jump around in the array that be obtained or send to a process,
     because we are sure that the data is ordered the same we as we
     obtained.
   */ 
  qsort(tr->partAssigns, tr->numAssignments, sizeof(Assign), sortById);
  
  if(tr->rateHetModel == CAT)
    setupBasePointersInTree( tr);
}

#ifdef _USE_OMP
void copyThreadAssignmentInfoToTree(PartitionAssignment *pa, tree *tr)
{
  int
    i, j;

  /* we want to know max number of partitions assigned to a single thread -> mainly for memory allocation */

  int
    *numsPerProc = (int *)calloc((size_t)pa->numProc, sizeof(int)),
    *numsPerPart = (int *)calloc((size_t)pa->numPartitions, sizeof(int));

  for(i = 0; i< pa->numProc; ++i)
    {
      for(j = 0; j < pa->numAssignPerProc[i]; ++j)
	{
	  Assignment *a = &pa->assignPerProc[i][j];
	  ++numsPerProc[i];
	  ++numsPerPart[a->partId];
	}
    }

  int
    pmax = 0;

  for(i = 1; i< pa->numProc; ++i)
    {
      if (numsPerProc[i] > numsPerProc[pmax])
	pmax = i;
    }

  /* save max partition count to the tree structure */
  tr->maxModelsPerThread = numsPerProc[pmax];

  assert(tr->maxModelsPerThread > 0 && tr->maxModelsPerThread <= pa->numPartitions);

  pmax = 0;
  for(i = 1; i< pa->numPartitions; ++i)
    {
      if (numsPerPart[i] > numsPerPart[pmax])
	pmax = i;
    }

  /* save max threads count to the tree structure */
  tr->maxThreadsPerModel = numsPerPart[pmax];

  assert(tr->maxThreadsPerModel > 0 && tr->maxThreadsPerModel <= pa->numProc);

  free(numsPerProc);
  free(numsPerPart);

  printf("\n maxModelsPerThread: %d,   maxThreadsPerModel: %d\n", tr->maxModelsPerThread, tr->maxThreadsPerModel);

  /* copy the partition assignment to the tree structure */

  int
    threadPartSize = pa->numProc * tr->maxModelsPerThread,
    partThreadSize = pa->numPartitions * tr->maxThreadsPerModel;

  tr->threadPartAssigns = (Assign **)calloc((size_t)threadPartSize, sizeof(Assign*));
  tr->partThreadAssigns = (Assign **)calloc((size_t)partThreadSize, sizeof(Assign*));

  for(i = 0; i < pa->numProc; ++i)
    {
      int
	partCount = 0;

      for( j = 0 ; j < pa->numAssignPerProc[i]; ++j)
        {
	  Assignment *ass = pa->assignPerProc[i] + j;
	  Assign* pTreeAss = (Assign *)calloc(1, sizeof(Assign));
	  pTreeAss->procId = i;
	  pTreeAss->offset = ass->offset;
	  pTreeAss->width = ass->width;
	  pTreeAss->partitionId = ass->partId;

	  size_t
	    ind = i * tr->maxModelsPerThread + partCount;

	  assert(ind < (i+1) * tr->maxModelsPerThread);

	  tr->threadPartAssigns[ind] = pTreeAss;
	  ++partCount;

	  ind = ass->partId * tr->maxThreadsPerModel;
	  while (tr->partThreadAssigns[ind])
	    ++ind;

	  assert( ind < (ass->partId+1) * tr->maxThreadsPerModel);

	  tr->partThreadAssigns[ind] = pTreeAss;
        }
    }
}
#endif

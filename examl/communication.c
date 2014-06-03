#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <mpi.h>

#include "axml.h"


extern int processes; 
extern int processID; 



/** 
    computes the count and displacement for gatherv/scatterv, assuming
    that the new partition assignment algorithm was used
*/ 
void calculateLengthAndDisplPerProcess(tree *tr, int **length_result, int **disp_result)
{
  int i; 

  *length_result = (int*) calloc((size_t) processes , sizeof(int)); 
  *disp_result = (int*) calloc((size_t) processes, sizeof(int)); 

  int* numPerProc = *length_result; 
  int* displPerProc= *disp_result;
  
  for(i = 0; i < tr->numAssignments; ++i)
    {
      Assign* ass = &(tr->partAssigns[i]); 
      numPerProc[ass->procId] += ass->width; 
    }

  displPerProc[0] = 0; 
  for(i = 1; i < processes  ; ++i)
    displPerProc[i] = displPerProc[i-1] + numPerProc[i-1];   
}


static size_t mapMpiTypeToSize(MPI_Datatype type)
{
  if(type == MPI_INT)
    return sizeof(int); 
  else if(type == MPI_DOUBLE)
    return sizeof(double); 
  else 
    {
      assert(0); 
      return 0; 
    }
}


/** 
    scatters a distributed array (e.g., what used to be
    tr->rateCategory) to partition-specfic arrays (e.g.,
    tr->partitionData[i].rateCategory). 

    This works, because tr->partitionData[i].rateCategory is a
    non-owning pointer to a position in the global resource array
    (e.g., tr->rateCategory_basePtr).
*/ 
void scatterDistrbutedArray(tree *tr, void *src, void *destination, MPI_Datatype type, int *countPerProc, int *displPerProc)
{
  int 
    i; 
  
  size_t 
    typeLen = mapMpiTypeToSize(type); 
  
  char 
    *srcReordered = (char *)NULL; 

  /* master must reorder the data   */
  if(processID == 0)
    {
      srcReordered = (char *)malloc(tr->originalCrunchedLength * typeLen); 
      int *seenPerProcesses = (int *)calloc((size_t) processes, sizeof(int)); 
      
      Assign *aIter = tr->partAssigns; 
      Assign *aEnd = &(tr->partAssigns[ tr->numAssignments ] ); 

      while(aIter != aEnd)
	{
	  pInfo *partition = &(tr->partitionData[ aIter->partitionId ]) ; 
	  memcpy( srcReordered +  ( (size_t) displPerProc[aIter->procId] + (size_t) seenPerProcesses[aIter->procId] )  * typeLen , 
		  ((char*) src) + (partition->lower + aIter->offset) * typeLen, 
		  aIter->width * typeLen); 
	  seenPerProcesses[aIter->procId] += aIter->width; 
	  ++aIter; 
	}

      for(i = 0; i < processes; ++i)
	assert(seenPerProcesses[i] == countPerProc[i]) ; 

      free(seenPerProcesses); 
    }
  
  MPI_Scatterv(srcReordered, countPerProc, displPerProc, type, destination, countPerProc[processID], type, 0, MPI_COMM_WORLD); 
 
  /* after this scatter, every process already has the data correctly
     ordered at its repective base pointer */
 
  if(processID == 0)
    free(srcReordered); 
}


/** 
    gathers a distributed array (e.g., what used to be
    tr->rateCategory) to partition-specfic arrays (e.g.,
    tr->partitionData[i].rateCategory).

    This works, because tr->partitionData[i].rateCategory is a
    non-owning pointer to a position in the global resource array
    (e.g., tr->rateCategory_basePtr).
*/ 
void gatherDistributedArray(tree *tr, void **destinationPtr, void *src, MPI_Datatype type, int* countPerProc, int *displPerProc)
{
  /* this is the raw array that the master will obtain from his
     peers. Data in this arrays are ordered per process */
  char 
    *destinationUnordered = (char*)NULL; 
  
  char
    *destination = (char*)NULL; 
  
  size_t
    typeLen = mapMpiTypeToSize(type); 
  
  if(processID == 0)
    {
      //TODO one pointer is of type void the other of type char, not really nice
      *destinationPtr = (void *)malloc( tr->originalCrunchedLength *  typeLen); 
      destinationUnordered = (char *)malloc( tr->originalCrunchedLength * typeLen); 
      destination = *destinationPtr; 
    }
  
  MPI_Gatherv(src, countPerProc[processID], type, destinationUnordered, countPerProc, displPerProc, type,0 , MPI_COMM_WORLD ); 
  
  /*
    here the master reorders the array it has obtained. Afterwards,
    destinationPtr is a pointer to the array that contains the global
    array that can be indexed by alignment position (i.e., if we have
    gathered tr->partitionData[i].lhs, then *destinationPtr
    corresponds to what previously was tr->lhs). This strongly couples
    the respective distributed array to tr->partAssigns.
   */ 
  if( processID == 0 )
    {
      int
	i, 
	*seenPerProcesses = (int*) calloc(processes, sizeof(int)); 

      Assign
	*aIter = tr->partAssigns; 

      Assign
	*aEnd = tr->partAssigns + tr->numAssignments; 

      while(aIter != aEnd)
	{
	  pInfo
	    *partition  = &(tr->partitionData[aIter->partitionId]); 

	  memcpy(destination + (size_t) (partition->lower + aIter->offset) * typeLen,
		 destinationUnordered +  (size_t) (displPerProc[aIter->procId] + seenPerProcesses[aIter->procId]) * typeLen ,
		 typeLen * aIter->width);
	  seenPerProcesses[aIter->procId] += aIter->width; 
	  ++aIter  ; 
	}
      
      /* check, if everything has been reordered */
      for(i = 0; i < processes; ++i)
      	assert(seenPerProcesses[i] == countPerProc[i]);
      free(seenPerProcesses); 
    }
}

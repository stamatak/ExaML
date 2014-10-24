#ifndef _PARTITIION_ASSIGNMENT
#define _PARTITIION_ASSIGNMENT

#include "axml.h"

#define not !  


typedef struct 
{
  int partId;
  size_t width; 
  size_t offset; 
} Assignment; 


typedef struct
{
  int id; 
  size_t width; 
  int type;
}  Partition; 


typedef struct
{
  int numProc; 
  int numPartitions; 
  Partition *partitions;
  Assignment **assignPerProc;  	/* procid -> array of assignments  */
  int *numAssignPerProc; 	/* procid -> size of above array */
} PartitionAssignment; 

/*
  constructor
*/ 
void initializePartitionAssignment( PartitionAssignment **pAssPtr, pInfo **partitions, int numPart, int numProc);
/* 
   deletor 
 */ 
void deletePartitionAssignment(PartitionAssignment *pAss); 
/*
  assign partitions to all proceses  
 */ 
void assign(PartitionAssignment *pa); 
/* 
   prints a single assignment 
 */ 
void printAssignment(Assignment a, int procid); 
/* 
   calculates and prints the load (number of partitions and number of
   sites) for each process
 */ 
void printLoad(PartitionAssignment *pa); 

void copyAssignmentInfoToTree(PartitionAssignment *pa, tree *tr); 

#ifdef _USE_OMP
void copyThreadAssignmentInfoToTree(PartitionAssignment *pa, tree *tr);
#endif

void printAssignments(PartitionAssignment *pa); 

#endif

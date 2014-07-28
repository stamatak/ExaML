#ifndef _BYTE_FILE
#define _BYTE_FILE

#include "axml.h"

#include "partitionAssignment.h"

#define ALN_HEAD 1 
#define ALN_WEIGHTS 2 
#define ALN_TAXA 4 
#define ALN_PARTITIONS 8  
#define ALN_ALIGNMENT 16 



typedef struct 
{
  int numTax; 
  size_t numPattern; 
  int numPartitions; 
  double gappyness; 
  pInfo **partitions;
  char **taxaNames; 
  FILE *fh; 
  char hasRead ; 
} ByteFile; 

/* 
   constructor  
*/ 
void initializeByteFile(ByteFile **bf, char *name); 
/* 
   deletor 
*/ 
void deleteByteFile(ByteFile *bf) ; 
/* 
   reads the header of a byte file   
*/  
void readHeader(ByteFile* bf); 
/* 
   reads partition information in a byte file 
*/ 
void readPartitions(ByteFile *bf); 
/* 
   reads the taxon names in a byte file   
*/
void readTaxa(ByteFile *bf); 
/* 
   reads weights and alignment characters in a byte file 
*/ 
void readMyData(ByteFile *bf, PartitionAssignment *pa, int procId); 
/*
  initializes a tree from a byte file.  

  @notice Since shallow copies are involved, you cannot copy the
  information from a byte file into multiple tree instances.
 */
void initializeTreeFromByteFile(ByteFile *bf, tree *tr); 

#endif

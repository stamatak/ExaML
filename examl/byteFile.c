#include <string.h> 
#include <malloc.h>

#include "byteFile.h"


#define READ_VAR(file,var) assert( fread(&var, sizeof(var),1, file ) == 1  )
#define READ_ARRAY(file, arrPtr, numElem, size)  assert( fread(arrPtr, size, numElem, file) ==  (unsigned int) numElem)



/** 
    seekPos finds the position in the byte file where a certain type
    of information is stored. See byteFile.h for possible values of
    "pos" .

    Notice, that this is a "fall-through" switch statement: if -- for
    instance -- we want to get to the position of the taxa, we have to
    skip everything that comes prior to the taxa in the file (but
    naturally not the taxa themselves).
 */ 
static void seekPos(ByteFile *bf, int pos)
{
  exa_off_t
    toSkip = 0;
  
  int 
    i; 

  switch(pos)
    {
    case ALN_ALIGNMENT: 	/* skips partitions */
      {
	assert(bf->hasRead & ALN_PARTITIONS); 
	pInfo p ; 

	toSkip +=  bf->numPartitions * ( sizeof(p.states) + sizeof(p.maxTipStates) + sizeof(p.lower) 
					 + sizeof(p.upper) + sizeof(p.width) + sizeof(p.dataType) + sizeof(p.protModels) 
					 + sizeof(p.protFreqs) + sizeof(p.nonGTR) + sizeof(p.optimizeBaseFrequencies)); 
	
	/* skip the names and their lengths */
	for( i = 0 ; i < bf->numPartitions; ++i)
	  {
	    pInfo *p = bf->partitions[i]; 
	    toSkip += (strlen(p->partitionName)+1 ) * sizeof(char) + sizeof(int);
	    toSkip += sizeof(double) * p->states; /* also skip frequncies */
	  }
      }
    case ALN_PARTITIONS: 	/* skips taxa  */
      {
	assert(bf->hasRead & ALN_TAXA); 
	for(i = 0; i < bf->numTax; ++i)
	  toSkip += (strlen(bf->taxaNames[i]) + 1)  * sizeof(char) + sizeof(int); 
      }
    case ALN_TAXA: 		/* skips weights */
      {
	assert(bf->hasRead & ALN_HEAD); 
	toSkip += bf->numPattern * sizeof(int); 
      }
    case ALN_WEIGHTS: 		/* skips header */
      {
	toSkip += 
	  sizeof(bf->numTax) + sizeof(bf->numPattern) 
	  + sizeof(bf->numPartitions) + sizeof(bf->gappyness); 
      }
    case ALN_HEAD : 
      toSkip += sizeof(int); 		/* skips the initial int that tells us how many bytes a size_t has */
      break; 
    default : 
      assert(0); 
    }
  
  exa_fseek(bf->fh, toSkip, SEEK_SET); 
} 

/** 
    initializes ByteFile **bf 
 */ 
void initializeByteFile(ByteFile **bf, char *name)
{
  *bf = (ByteFile *)calloc(1,sizeof(ByteFile)); 
  ByteFile *result = *bf; 
  result->fh  = myfopen(name, "rb"); 

  int sizeOfSizeT = 0; 
  READ_VAR(result->fh, sizeOfSizeT); 

  if(sizeOfSizeT != sizeof(size_t))
    {
      printf("Error: the address data type has a size of %d bits on the current system while on the system on which you created the binary alignment file using the parser the address size is %d bits!\n", 
	     8 * (int)sizeof(size_t), 8 * sizeOfSizeT);
      printf("Usually this indicates that the parser was executed on a 32-bit system while you are trying to run ExaML on a 64-bit system.\n");
      printf("Please parse the binary alignment file on the same hardware on which you intend to run ExaML.\n"); 
      MPI_Finalize();
    }
} 


/** 
    a shallow cleanup of ByteFile *bf. Notice, that various data may
    have been copied (by pointer value) to our tree instance and
    therefore should not be clean up.
 */ 
void deleteByteFile(ByteFile *bf) 
{
  /* only a shallow free! pointers inside the pInfo must persist */
  int i; 

  if(bf->partitions)
    {
      for( i = 0; i < bf->numPartitions; ++i)
      	free(bf->partitions[i]);
      free(bf->partitions); 
    }

  if(bf->fh)
    fclose(bf->fh); 

  if(bf->taxaNames )
    {
      for(i = 0; i < bf->numTax; ++i)
	free(bf->taxaNames[i] ); 
    }
  free(bf->taxaNames); 
  free(bf);
} 




/** 
    only reads initial header information 
 */ 
void readHeader(ByteFile* bf)
{
  seekPos(bf, ALN_HEAD); 
  READ_VAR(bf->fh, bf->numTax); 
  READ_VAR(bf->fh, bf->numPattern); 
  READ_VAR(bf->fh, bf->numPartitions); 
  READ_VAR(bf->fh, bf->gappyness) ;
  bf->hasRead |= ALN_HEAD;

}


/** 
    reads partition information from the byte file.
 */
void readPartitions(ByteFile *bf)
{
  int i ; 
  
  seekPos(bf, ALN_PARTITIONS); 
  
  assert(bf->partitions == (pInfo **)NULL); 
  bf->partitions = (pInfo **)calloc(bf->numPartitions, sizeof(pInfo*) ); 
  for(i = 0; i < bf->numPartitions; ++i)
    {
      bf->partitions[i] = (pInfo*)calloc(1,sizeof(pInfo)); 
      pInfo* p = bf->partitions[i]; 

      READ_VAR(bf->fh, p->states);
      READ_VAR(bf->fh, p->maxTipStates);
      READ_VAR(bf->fh, p->lower);
      READ_VAR(bf->fh, p->upper);

      /* DONT use this value! */
      READ_VAR(bf->fh, p->width);
      p->width = 0; 

      READ_VAR(bf->fh, p->dataType);
      READ_VAR(bf->fh, p->protModels);
      //READ_VAR(bf->fh, p->autoProtModels);
      READ_VAR(bf->fh, p->protFreqs);
      READ_VAR(bf->fh, p->nonGTR);
      READ_VAR(bf->fh, p->optimizeBaseFrequencies);
      //      READ_VAR(bf->fh, p->numberOfCategories);

      /* read string */
      unsigned int len = 0; 
      READ_VAR(bf->fh, len); 
      p->partitionName = (char*)calloc(len,sizeof(char));
      READ_ARRAY(bf->fh, p->partitionName, len, sizeof(char)); 

      p->frequencies = (double*)calloc(p->states, sizeof(double)); 
      READ_ARRAY(bf->fh, p->frequencies, p->states , sizeof(double)); 
    }
  
  bf->hasRead |= ALN_PARTITIONS; 
}


/** 
    reads the taxon names from the byte file  
 */ 
void readTaxa(ByteFile *bf)
{
  int i; 

  assert(bf->taxaNames == (char **)NULL);
  seekPos(bf,  ALN_TAXA); 

  bf->taxaNames = (char **)calloc(bf->numTax, sizeof(char*));
  for(i = 0; i < bf->numTax; ++i)
    {
      int len = 0; 
      READ_VAR(bf->fh, len ); 
      bf->taxaNames[i] = (char*)calloc(len, sizeof(char)); 
      READ_ARRAY(bf->fh, bf->taxaNames[i], len, sizeof(char)); 
    }

  bf->hasRead |= ALN_TAXA; 
}


/** 
    uses the information in the PartitionAssignment to only extract
    data relevant to this process (weights and alignment characters).
 */ 
void readMyData(ByteFile *bf, PartitionAssignment *pa, int procId)
{
  seekPos(bf, ALN_ALIGNMENT); 

  exa_off_t
    alnPos = exa_ftell(bf->fh); 

  size_t 
    len; 

  int numAssign = pa->numAssignPerProc[procId];
  Assignment *myAssigns = pa->assignPerProc[procId];

  /* first read aln characters   */
  int i,j ; 
  for(i = 0; i < numAssign; ++i )
    {
      Assignment a = myAssigns[i]; 
      /* printf("reading for: ") ;  */
      /* printAssignment(a, procId);  */

      pInfo *partition = bf->partitions[a.partId]; 
      partition->width = a.width; 
      partition->offset = a.offset; 
      len = bf->numTax * a.width; 
      partition->yResource = (unsigned char*)malloc_aligned( len * sizeof(unsigned char)); 
      memset(partition->yResource,0,len * sizeof(unsigned char)); 
      partition->yVector = (unsigned char**) calloc(bf->numTax + 1 , sizeof(unsigned char*)); 
      for(j = 1; j <= bf->numTax; ++j)
	partition->yVector[j] = partition->yResource + (j-1) * a.width; 

      for(j = 1; j <= bf->numTax; ++j )
	{
	  exa_off_t pos = alnPos + (  bf->numPattern * (j-1)    +  partition->lower + a.offset ) * sizeof(unsigned char); 
	  assert(alnPos <= pos); 
	  exa_fseek(bf->fh, pos, SEEK_SET); 
	  READ_ARRAY(bf->fh, partition->yVector[j], a.width, sizeof(unsigned char));
	}
    }

  
  /* now read weights  */
  seekPos(bf, ALN_WEIGHTS); 

  exa_off_t
    wgtPos = exa_ftell(bf->fh); 
  assert( ! (wgtPos <  0) );

  for(i = 0; i < numAssign; ++i)
    {
      Assignment a = myAssigns[i]; 
      pInfo *partition = bf->partitions[a.partId]; 
      len = a.width * sizeof(int); 
      partition->wgt = (int*)malloc_aligned( len); 
      memset(partition->wgt, 0, len); 

      exa_off_t pos = wgtPos +  (partition->lower  + a.offset) * sizeof(int); 
      assert(wgtPos <= pos );
      
      exa_fseek(bf->fh, pos, SEEK_SET); 
      READ_ARRAY(bf->fh, partition->wgt, a.width, sizeof(int)); 

    }

  bf->hasRead |= ALN_ALIGNMENT; 
  bf->hasRead |= ALN_WEIGHTS; 
} 


/** 
    copies all relevant information from our byte file to the tree
    instance.
 */ 
void initializeTreeFromByteFile(ByteFile *bf, tree *tr)
{
  assert( ( bf->hasRead & ALN_HEAD )
	  && (bf->hasRead & ALN_WEIGHTS)
	  && (bf->hasRead & ALN_TAXA) 
	  && (bf->hasRead & ALN_PARTITIONS)
	  && (bf->hasRead & ALN_ALIGNMENT ) ); 
 
  /* some additional stuff we read */
  tr->mxtips = bf->numTax;
  tr->originalCrunchedLength = bf->numPattern; 
  tr->NumberOfModels = bf->numPartitions; 
  tr->gapyness = bf->gappyness; 
  
  /* deep copy of taxa */
  int i ; 
  tr->nameList = (char **)calloc((size_t)(tr->mxtips + 1), sizeof(char *)  );
  
  tr->nameList[0] = (char *)NULL;

  for(i = 1; i <= bf->numTax; ++i)
    {
      tr->nameList[i] = (char*)calloc(strlen(bf->taxaNames[i-1]) + 1, sizeof(char)); 
      strcpy(tr->nameList[i], bf->taxaNames[i-1]);      
    }

  /* 
   * shallow copy of partitions 
   * 
   * partition contains only shallow copies of a few data arrays that
   * needed to be initialized at this point
   */
  int myLength = 0; 
  tr->partitionData = (pInfo*)calloc( tr->NumberOfModels, sizeof(pInfo));
  for(i = 0; i < tr->NumberOfModels; ++i)
    {
      tr->partitionData[i] = *(bf->partitions[i]);
      myLength += tr->partitionData[i].width; 
      assert( bf->partitions[i]->wgt != (int*)NULL || bf->partitions[i]->width == 0); 
      assert( ( tr->partitionData[i].wgt != (int*)NULL)  || ( tr->partitionData[i].width == 0 ) ); 
    }
} 



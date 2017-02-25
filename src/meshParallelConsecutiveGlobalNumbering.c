#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include <mpi.h>

#include "mesh.h"

typedef struct {

  iint localId;
  iint globalId;
  iint recvId;
  iint newGlobalId;
  iint originalRank;
  iint ownerRank;
  
}parallelNode_t;

// compare on global indices 
int parallelCompareGlobalIndices(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(fa->globalId < fb->globalId) return -1;
  if(fa->globalId > fb->globalId) return +1;

  return 0;  
}

// compare on global indices 
int parallelCompareSourceIndices(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;
  
  if(fa->originalRank < fb->originalRank) return -1;
  if(fa->originalRank > fb->originalRank) return +1;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  return 0;

}

// compare on global indices 
int parallelCompareOwners(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(fa->ownerRank < fb->ownerRank) return -1;
  if(fa->ownerRank > fb->ownerRank) return +1;

  return 0;  
}



// squeeze gaps out of a globalNumbering of local nodes (arranged in NpNum blocks
void meshParallelConsecutiveGlobalNumbering(iint Nnum,
					    iint *globalNumbering, iint *globalOwners, 
              iint *globalStarts, iint *sendSortId,
              iint *globalSortId, iint **compressId,
              iint *sendCounts, iint *sendOffsets,
              iint *recvCounts, iint *recvOffsets){

  // need to handle globalNumbering = 0
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // build GS for this numbering
  void *gsh = gsParallelGatherScatterSetup(Nnum, globalNumbering);

  iint *ranks = (iint*) calloc(Nnum, sizeof(iint));
  for(iint n=0;n<Nnum;++n)
    ranks[n] = rank;
  
  // find lowest rank process that contains each node (decides ownership)
  gsParallelGatherScatter(gsh, ranks, "int", "min"); // should use iint

  // clean up
  gsParallelGatherScatterDestroy(gsh);
  
  // count how many nodes to send to each process
  
  iint *allCounts   = (iint*) calloc(size, sizeof(iint));
  iint *allOffsets  = (iint*) calloc(size+1, sizeof(iint));

  for(iint n=0;n<Nnum;++n)
    sendCounts[ranks[n]] += sizeof(parallelNode_t);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_IINT, recvCounts, 1, MPI_IINT, MPI_COMM_WORLD);
  
  // find send and recv offsets for gather
  iint recvNtotal = 0;
  for(iint r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r]/sizeof(parallelNode_t);
  }

  // populate parallel nodes to send
  parallelNode_t *sendNodes = (parallelNode_t*) calloc(Nnum, sizeof(parallelNode_t));
  for(iint n=0;n<Nnum;++n){
    sendNodes[n].localId = n;
    sendNodes[n].globalId = globalNumbering[n];
    sendNodes[n].newGlobalId = -1;
    sendNodes[n].originalRank = rank;
    sendNodes[n].ownerRank = ranks[n];
  }

  // sort by global index
  qsort(sendNodes, Nnum, sizeof(parallelNode_t), parallelCompareOwners);

  //record send indexing
  for (iint n=0;n<Nnum;n++) sendSortId[n] = sendNodes[n].localId;
  
  parallelNode_t *recvNodes = (parallelNode_t*) calloc(recvNtotal, sizeof(parallelNode_t));
  
  // load up node data to send (NEED TO SCALE sendCounts, sendOffsets etc by sizeof(parallelNode_t)
  MPI_Alltoallv(sendNodes, sendCounts, sendOffsets, MPI_CHAR,
		recvNodes, recvCounts, recvOffsets, MPI_CHAR,
		MPI_COMM_WORLD);

  for (iint n = 0; n<recvNtotal;n++) recvNodes[n].recvId = n;

  // sort by global index
  qsort(recvNodes, recvNtotal, sizeof(parallelNode_t), parallelCompareGlobalIndices);

  //record global indexing
  for (iint n=0;n<recvNtotal;n++) globalSortId[n] = recvNodes[n].recvId;

  // renumber unique nodes starting from 0 (need to be careful about zeros)
  iint cnt = 0;
  recvNodes[0].newGlobalId = cnt;
  for(iint n=1;n<recvNtotal;++n){
    if(recvNodes[n].globalId!=recvNodes[n-1].globalId){ // new node
      ++cnt;
    }
    recvNodes[n].newGlobalId = cnt;
  }
  ++cnt; // increment to actual number of unique nodes on this rank

  *compressId = (iint*) calloc(cnt+1,sizeof(iint));
  cnt = 0;
  *compressId[0] = 0;
  for(iint n=1;n<recvNtotal;++n){
    if(recvNodes[n].globalId!=recvNodes[n-1].globalId){ // new node
      (*compressId)[++cnt] = n;
    }
  }
  (*compressId)[++cnt] = recvNtotal;

  //  printf("cnt %d on rank %d\n", cnt, rank);

  // collect unique node counts from all processes
  MPI_Allgather(&cnt, 1, MPI_IINT, allCounts, 1, MPI_IINT, MPI_COMM_WORLD);

  // cumulative sum of unique node counts => starting node index for each process
  for(iint r=0;r<size;++r)
    allOffsets[r+1] = allOffsets[r] + allCounts[r];

  memcpy(globalStarts, allOffsets, (size+1)*sizeof(iint));
  
  // shift numbering
  for(iint n=0;n<recvNtotal;++n)
    recvNodes[n].newGlobalId += allOffsets[rank];

#if 0
  for(iint r=0;r<size;++r) fflush(stdout);
  for(iint r=0;r<size;++r){
    MPI_Barrier(MPI_COMM_WORLD);
    if(r==rank){
      printf("allOFfsets = %d\n", allOffsets[r]);
      for(iint n=0;n<recvNtotal;++n)
	printf("rank %d old id %d new id %d\n", rank, recvNodes[n].globalId, recvNodes[n].newGlobalId);
      fflush(stdout);
    }
  }
#endif
  
  // sort by rank, local index
  qsort(recvNodes, recvNtotal, sizeof(parallelNode_t), parallelCompareSourceIndices);

  // reverse all to all to reclaim nodes
  MPI_Alltoallv(recvNodes, recvCounts, recvOffsets, MPI_CHAR,
		sendNodes, sendCounts, sendOffsets, MPI_CHAR,
		MPI_COMM_WORLD);

  // extract new global indices and push back to original numbering array
  for(iint n=0;n<Nnum;++n){
    // shuffle incoming nodes based on local id
    iint id = sendNodes[n].localId;
    globalNumbering[id] = sendNodes[n].newGlobalId;
    globalOwners[id] = sendNodes[n].ownerRank;
  }

  for (iint i=0;i<size;i++) {
    sendCounts[i] /= sizeof(parallelNode_t);
    recvCounts[i] /= sizeof(parallelNode_t);
  }
  for (iint i=0;i<size+1;i++) {
    sendOffsets[i] /= sizeof(parallelNode_t);
    recvOffsets[i] /= sizeof(parallelNode_t);
  }

  free(ranks);

  free(allCounts);
  free(allOffsets);
  free(sendNodes);
  free(recvNodes);
}

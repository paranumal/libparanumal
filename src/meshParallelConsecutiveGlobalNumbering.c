#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"

typedef struct {

  iint localId;
  iint globalId;
  iint recvId;
  iint newGlobalId;
  iint originalRank;
  iint ownerRank;
  
}parallelNode2_t;

// compare on global indices 
int parallelCompareGlobalIndices(const void *a, const void *b){

  parallelNode2_t *fa = (parallelNode2_t*) a;
  parallelNode2_t *fb = (parallelNode2_t*) b;

  if(fa->globalId < fb->globalId) return -1;
  if(fa->globalId > fb->globalId) return +1;

  return 0;  
}

// compare on global indices 
int parallelCompareSourceIndices(const void *a, const void *b){

  parallelNode2_t *fa = (parallelNode2_t*) a;
  parallelNode2_t *fb = (parallelNode2_t*) b;
  
  if(fa->originalRank < fb->originalRank) return -1;
  if(fa->originalRank > fb->originalRank) return +1;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  return 0;

}

// compare on global indices 
int parallelCompareOwners2(const void *a, const void *b){

  parallelNode2_t *fa = (parallelNode2_t*) a;
  parallelNode2_t *fb = (parallelNode2_t*) b;

  if(fa->ownerRank < fb->ownerRank) return -1;
  if(fa->ownerRank > fb->ownerRank) return +1;

  return 0;  
}

// squeeze gaps out of a globalNumbering of local nodes (arranged in NpNum blocks
void meshParallelConsecutiveGlobalNumbering(mesh_t *mesh,
                                            iint Nnum,
                                            iint *globalNumbering, 
                                            iint *globalOwners, 
                                            iint *globalStarts){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  // count how many nodes to send to each process
  iint *allCounts   = (iint*) calloc(size, sizeof(iint));

  iint *sendCounts = (iint *) calloc(size,sizeof(iint));
  iint *recvCounts = (iint *) calloc(size,sizeof(iint));
  iint *sendOffsets = (iint *) calloc(size+1,sizeof(iint));
  iint *recvOffsets = (iint *) calloc(size+1,sizeof(iint));
  
  iint cnt = 0;
  for(iint n=0;n<Nnum;++n) {
    if (globalNumbering[n] < 0) continue; //skip negative ids
    sendCounts[globalOwners[n]] += sizeof(parallelNode2_t);
    cnt++;
  }

  iint Nlocal = cnt; //number of unmasked nodes

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_IINT, recvCounts, 1, MPI_IINT, MPI_COMM_WORLD);
  
  // find send and recv offsets for gather
  iint recvNtotal = 0;
  for(iint r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r]/sizeof(parallelNode2_t);
  }

  // populate parallel nodes to send
  parallelNode2_t *sendNodes;
  if (Nlocal)
    sendNodes = (parallelNode2_t*) calloc(Nlocal, sizeof(parallelNode2_t));
  cnt = 0;
  for(iint n=0;n<Nnum;++n){
    if (globalNumbering[n] < 0) continue; //skip negative ids
    sendNodes[cnt].localId = n;
    sendNodes[cnt].globalId = globalNumbering[n];
    sendNodes[cnt].newGlobalId = -1;
    sendNodes[cnt].originalRank = rank;
    sendNodes[cnt].ownerRank = globalOwners[n];
    cnt++;
  }

  // sort by global index
  qsort(sendNodes, Nlocal, sizeof(parallelNode2_t), parallelCompareOwners2);
  
  parallelNode2_t *recvNodes;
  if (recvNtotal)
    recvNodes = (parallelNode2_t*) calloc(recvNtotal, sizeof(parallelNode2_t));
  
  // load up node data to send (NEED TO SCALE sendCounts, sendOffsets etc by sizeof(parallelNode2_t)
  MPI_Alltoallv(sendNodes, sendCounts, sendOffsets, MPI_CHAR,
    recvNodes, recvCounts, recvOffsets, MPI_CHAR,
    MPI_COMM_WORLD);

  for (iint n = 0; n<recvNtotal;n++) recvNodes[n].recvId = n;

  // sort by global index
  qsort(recvNodes, recvNtotal, sizeof(parallelNode2_t), parallelCompareGlobalIndices);

  // renumber unique nodes starting from 0 (need to be careful about zeros)
  cnt = 0;
  if (recvNtotal) recvNodes[0].newGlobalId = cnt;
  for(iint n=1;n<recvNtotal;++n){
    if(recvNodes[n].globalId!=recvNodes[n-1].globalId){ // new node
      ++cnt;
    }
    recvNodes[n].newGlobalId = cnt;
  }
  if (recvNtotal) ++cnt; // increment to actual number of unique nodes on this rank

  // collect unique node counts from all processes
  MPI_Allgather(&cnt, 1, MPI_IINT, allCounts, 1, MPI_IINT, MPI_COMM_WORLD);

  // cumulative sum of unique node counts => starting node index for each process
  for(iint r=0;r<size;++r)
    globalStarts[r+1] = globalStarts[r] + allCounts[r];
  
  // shift numbering
  for(iint n=0;n<recvNtotal;++n)
    recvNodes[n].newGlobalId += globalStarts[rank];
  
  // sort by rank, local index
  qsort(recvNodes, recvNtotal, sizeof(parallelNode2_t), parallelCompareSourceIndices);

  // reverse all to all to reclaim nodes
  MPI_Alltoallv(recvNodes, recvCounts, recvOffsets, MPI_CHAR,
    sendNodes, sendCounts, sendOffsets, MPI_CHAR,
    MPI_COMM_WORLD);

  // extract new global indices and push back to original numbering array
  for(iint n=0;n<Nlocal;++n){
    // shuffle incoming nodes based on local id
    iint id = sendNodes[n].localId;
    globalNumbering[id] = sendNodes[n].newGlobalId;
  }

  free(sendCounts);
  free(recvCounts);
  free(sendOffsets);
  free(recvOffsets);
  free(allCounts);
  free(sendNodes);
  free(recvNodes);
}
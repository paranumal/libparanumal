#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"

typedef struct{

  int element; // local element id
  int node;    // local node id
  int rank;    // rank of original node
  int id;      // original id
  int haloFlag;

  // info on base node (lowest rank node)
  int baseElement; 
  int baseNode;
  int baseRank;
  int baseId;

  int newGlobalId;

}parallelNode_t;

// compare on base rank then by globalId
int parallelCompareOwners(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(fa->baseRank < fb->baseRank) return -1;
  if(fa->baseRank > fb->baseRank) return +1;

  return 0;
}

// compare on base rank then by globalId
int parallelCompareSourceRank(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(fa->rank < fb->rank) return -1;
  if(fa->rank > fb->rank) return +1;

  return 0;
}

// compare on base rank for sorting suitable for destination
int parallelCompareBaseNodes(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if ((fa->baseRank==rank)&&(fb->baseRank!=rank)) return -1; //move locally-owned nodes to the beginning
  if ((fa->baseRank!=rank)&&(fb->baseRank==rank)) return  1;

  if(fa->baseRank < fb->baseRank) return -1;
  if(fa->baseRank > fb->baseRank) return +1;

  if(fa->haloFlag < fb->haloFlag) return -1;
  if(fa->haloFlag > fb->haloFlag) return +1;

  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;

  return 0;
}

// iteratively find a gather numbering for all local element nodes
void meshParallelConnectNodes(mesh_t *mesh){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int localNodeCount = mesh->Np*mesh->Nelements;
  int *allLocalNodeCounts = (int*) calloc(size, sizeof(int));

  MPI_Allgather(&localNodeCount, 1, MPI_INT,
                allLocalNodeCounts, 1, MPI_INT,
                MPI_COMM_WORLD);
  
  int gatherNodeStart = 0;
  for(int r=0;r<rank;++r)
    gatherNodeStart += allLocalNodeCounts[r];
  
  free(allLocalNodeCounts);

  // form continuous node numbering (local=>virtual gather)
  parallelNode_t *sendNodes =
    (parallelNode_t*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,
                             sizeof(parallelNode_t));

  // use local numbering
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      int id = e*mesh->Np+n;
      sendNodes[id].element = e;
      sendNodes[id].node = n;
      sendNodes[id].rank = rank;
      sendNodes[id].id = 1 + id + mesh->Nnodes + gatherNodeStart;

      sendNodes[id].baseElement = e;
      sendNodes[id].baseNode = n;
      sendNodes[id].baseRank = rank;
      sendNodes[id].baseId = 1 + id + mesh->Nnodes + gatherNodeStart;

    }

    // use vertex ids for vertex nodes to reduce iterations
    for(int v=0;v<mesh->Nverts;++v){
      int id = e*mesh->Np + mesh->vertexNodes[v];
      int gid = mesh->EToV[e*mesh->Nverts+v] + 1;
      sendNodes[id].id = gid; 
      sendNodes[id].baseId = gid; 
    }

    // label halo flags
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1){
        for(int n=0;n<mesh->Nfp;++n){
          int id = e*mesh->Np+mesh->faceNodes[f*mesh->Nfp+n];
          sendNodes[id].haloFlag = 1;
        }
      }
    }
  }

  int localChange = 0, gatherChange = 1;

  parallelNode_t *sendBuffer =
    (parallelNode_t*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(parallelNode_t));

  // keep comparing numbers on positive and negative traces until convergence
  while(gatherChange>0){

    // reset change counter
    localChange = 0;

    // send halo data and recv into extension of buffer
    meshHaloExchange(mesh, mesh->Np*sizeof(parallelNode_t),
                     sendNodes, sendBuffer, sendNodes+localNodeCount);

    // compare trace nodes
    for(int e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
        int id  = e*mesh->Nfp*mesh->Nfaces + n;
        int idM = mesh->vmapM[id];
        int idP = mesh->vmapP[id];
        int gidM = sendNodes[idM].baseId;
        int gidP = sendNodes[idP].baseId;

        int baseRankM = sendNodes[idM].baseRank;
        int baseRankP = sendNodes[idP].baseRank;

        // use minimum of trace variables
        int haloM = sendNodes[idM].haloFlag;
        int haloP = sendNodes[idP].haloFlag;
        if(haloM!=haloP) {
          ++localChange;
          sendNodes[idM].haloFlag = mymax(haloM, haloP);
          sendNodes[idP].haloFlag = mymax(haloM, haloP);  
        }
        
        if(gidM<gidP || (gidP==gidM && baseRankM<baseRankP)){
          ++localChange;
          sendNodes[idP].baseElement = sendNodes[idM].baseElement;
          sendNodes[idP].baseNode    = sendNodes[idM].baseNode;
          sendNodes[idP].baseRank    = sendNodes[idM].baseRank;
          sendNodes[idP].baseId      = sendNodes[idM].baseId;
        }
        
        if(gidP<gidM || (gidP==gidM && baseRankP<baseRankM)){
          ++localChange;
          sendNodes[idM].baseElement = sendNodes[idP].baseElement;
          sendNodes[idM].baseNode    = sendNodes[idP].baseNode;
          sendNodes[idM].baseRank    = sendNodes[idP].baseRank;
          sendNodes[idM].baseId      = sendNodes[idP].baseId;
        }
      }
    }

    // sum up changes
    MPI_Allreduce(&localChange, &gatherChange, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }

  // sort based on base nodes (rank,element,node at base)
  qsort(sendNodes, localNodeCount, sizeof(parallelNode_t), parallelCompareOwners);
 
  // count how many nodes to send to each process
  int *sendCounts = (int *) calloc(size,sizeof(int));
  int *recvCounts = (int *) calloc(size,sizeof(int));
  int *sendOffsets = (int *) calloc(size+1,sizeof(int));
  int *recvOffsets = (int *) calloc(size+1,sizeof(int));

  for(int n=0;n<localNodeCount;++n) 
    sendCounts[sendNodes[n].baseRank] += sizeof(parallelNode_t);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, MPI_COMM_WORLD);
  
  // find send and recv offsets for gather
  int recvNtotal = 0;
  for(int r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r]/sizeof(parallelNode_t);
  }

  parallelNode_t *recvNodes;
  if (recvNtotal) recvNodes = (parallelNode_t*) calloc(recvNtotal, sizeof(parallelNode_t));
  
  // load up node data to send 
  MPI_Alltoallv(sendNodes, sendCounts, sendOffsets, MPI_CHAR,
                recvNodes, recvCounts, recvOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  // sort by global index shifting halo nodes to the end
  qsort(recvNodes, recvNtotal, sizeof(parallelNode_t), parallelCompareBaseNodes);

  // renumber unique nodes starting from 0 (need to be careful about zeros)
  int Ngather = 0;
  if (recvNtotal) recvNodes[0].newGlobalId = Ngather;
  for(int n=1;n<recvNtotal;++n){
    if(recvNodes[n].baseId!=recvNodes[n-1].baseId){ // new node
      ++Ngather;
    }
    recvNodes[n].newGlobalId = Ngather;
  }
  if (recvNtotal) ++Ngather; // increment to actual number of unique nodes on this rank

  // collect unique node counts from all processes
  int *allGather   = (int*) calloc(size, sizeof(int));
  MPI_Allgather(&Ngather, 1, MPI_INT, allGather, 1, MPI_INT, MPI_COMM_WORLD);

  // cumulative sum of unique node counts => starting node index for each process
  mesh->gatherGlobalStarts = (int*) calloc(size+1, sizeof(int));
  for(int r=0;r<size;++r)
    mesh->gatherGlobalStarts[r+1] = mesh->gatherGlobalStarts[r] + allGather[r];

  // shift numbering
  for(int n=0;n<recvNtotal;++n)
    recvNodes[n].newGlobalId += mesh->gatherGlobalStarts[rank];
  
  // sort by rank, local index
  qsort(recvNodes, recvNtotal, sizeof(parallelNode_t), parallelCompareSourceRank);

  // reverse all to all to reclaim nodes
  MPI_Alltoallv(recvNodes, recvCounts, recvOffsets, MPI_CHAR,
                sendNodes, sendCounts, sendOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  // sort by rank, local index
  qsort(sendNodes, localNodeCount, sizeof(parallelNode_t), parallelCompareBaseNodes);

  // extract base index of each node (i.e. gather numbering)
  mesh->gatherLocalIds  = (int*) calloc(localNodeCount, sizeof(int));
  mesh->gatherBaseIds   = (int*) calloc(localNodeCount, sizeof(int));
  mesh->gatherBaseRanks = (int*) calloc(localNodeCount, sizeof(int));
  mesh->gatherHaloFlags = (int*) calloc(localNodeCount, sizeof(int));

  for(int id=0;id<localNodeCount;++id){
    mesh->gatherLocalIds[id]  = sendNodes[id].element*mesh->Np+sendNodes[id].node;
    mesh->gatherBaseIds[id]   = sendNodes[id].newGlobalId+1;
    mesh->gatherBaseRanks[id] = sendNodes[id].baseRank;
    mesh->gatherHaloFlags[id] = sendNodes[id].haloFlag;
  }
  
  free(sendBuffer);
  free(sendNodes);

  //make a locally-ordered version
  mesh->globalIds      = (int*) calloc(localNodeCount, sizeof(int));
  mesh->globalOwners   = (int*) calloc(localNodeCount, sizeof(int));
  mesh->globalHaloFlags= (int*) calloc(localNodeCount, sizeof(int));

  for(int id=0;id<localNodeCount;++id){
    int localId = mesh->gatherLocalIds[id];
    mesh->globalIds[localId]      = mesh->gatherBaseIds[id];    
    mesh->globalOwners[localId]   = mesh->gatherBaseRanks[id];
    mesh->globalHaloFlags[localId]= mesh->gatherHaloFlags[id];
  }
}

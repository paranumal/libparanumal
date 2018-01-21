#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"

typedef struct{

  iint element; // local element id
  iint node;    // local node id
  iint rank;    // rank of original node
  iint id;      // original id
  int haloFlag;

  // info on base node (lowest rank node)
  iint baseElement; 
  iint baseNode;
  iint baseRank;
  iint baseId;

  iint newGlobalId;

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

  iint rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if ((fa->baseRank==rank)&&(fb->baseRank!=rank)) return -1; //move locally-owned nodes to the beginning
  if ((fa->baseRank!=rank)&&(fb->baseRank==rank)) return  1;

  if(fa->baseRank < fb->baseRank) return -1;
  if(fa->baseRank > fb->baseRank) return +1;

  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;

  return 0;
}

// iteratively find a gather numbering for all local element nodes
void meshParallelConnectNodes(mesh_t *mesh){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  iint localNodeCount = mesh->Np*mesh->Nelements;
  iint *allLocalNodeCounts = (iint*) calloc(size, sizeof(iint));

  MPI_Allgather(&localNodeCount, 1, MPI_IINT,
                allLocalNodeCounts, 1, MPI_IINT,
                MPI_COMM_WORLD);
  
  iint gatherNodeStart = 0;
  for(iint r=0;r<rank;++r)
    gatherNodeStart += allLocalNodeCounts[r];
  
  free(allLocalNodeCounts);

  // form continuous node numbering (local=>virtual gather)
  parallelNode_t *sendNodes =
    (parallelNode_t*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,
                             sizeof(parallelNode_t));

  // use local numbering
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint id = e*mesh->Np+n;
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
    for(iint v=0;v<mesh->Nverts;++v){
      iint id = e*mesh->Np + mesh->vertexNodes[v];
      iint gid = mesh->EToV[e*mesh->Nverts+v] + 1;
      sendNodes[id].id = gid; 
      sendNodes[id].baseId = gid; 
    }

    // label halo flags
    for(iint f=0;f<mesh->Nfaces;++f){
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1){
        for(iint n=0;n<mesh->Nfp;++n){
          iint id = e*mesh->Np+mesh->faceNodes[f*mesh->Nfp+n];
          sendNodes[id].haloFlag = 1;
        }
      }
    }
  }

  iint localChange = 0, gatherChange = 1;

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
    for(iint e=0;e<mesh->Nelements;++e){
      for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){
        iint id  = e*mesh->Nfp*mesh->Nfaces + n;
        iint idM = mesh->vmapM[id];
        iint idP = mesh->vmapP[id];
        iint gidM = sendNodes[idM].baseId;
        iint gidP = sendNodes[idP].baseId;

        iint baseRankM = sendNodes[idM].baseRank;
        iint baseRankP = sendNodes[idP].baseRank;

        // use minimum of trace variables
        iint haloM = sendNodes[idM].haloFlag;
        iint haloP = sendNodes[idP].haloFlag;
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
    MPI_Allreduce(&localChange, &gatherChange, 1, MPI_IINT, MPI_SUM, MPI_COMM_WORLD);
  }

  // sort based on base nodes (rank,element,node at base)
  qsort(sendNodes, localNodeCount, sizeof(parallelNode_t), parallelCompareOwners);

  // now squeeze gaps out of a globalNumbering of local nodes
 
  // count how many nodes to send to each process
  iint *sendCounts = (iint *) calloc(size,sizeof(iint));
  iint *recvCounts = (iint *) calloc(size,sizeof(iint));
  iint *sendOffsets = (iint *) calloc(size+1,sizeof(iint));
  iint *recvOffsets = (iint *) calloc(size+1,sizeof(iint));

  for(iint n=0;n<localNodeCount;++n) 
    sendCounts[sendNodes[n].baseRank] += sizeof(parallelNode_t);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_IINT, recvCounts, 1, MPI_IINT, MPI_COMM_WORLD);
  
  // find send and recv offsets for gather
  iint recvNtotal = 0;
  for(iint r=0;r<size;++r){
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

  // sort by global index
  qsort(recvNodes, recvNtotal, sizeof(parallelNode_t), parallelCompareBaseNodes);

  // renumber unique nodes starting from 0 (need to be careful about zeros)
  iint Ngather = 0;
  if (recvNtotal) recvNodes[0].newGlobalId = Ngather;
  for(iint n=1;n<recvNtotal;++n){
    if(recvNodes[n].baseId!=recvNodes[n-1].baseId){ // new node
      ++Ngather;
    }
    recvNodes[n].newGlobalId = Ngather;
  }
  if (recvNtotal) ++Ngather; // increment to actual number of unique nodes on this rank

  // collect unique node counts from all processes
  iint *allGather   = (iint*) calloc(size, sizeof(iint));
  MPI_Allgather(&Ngather, 1, MPI_IINT, allGather, 1, MPI_IINT, MPI_COMM_WORLD);

  // cumulative sum of unique node counts => starting node index for each process
  mesh->gatherGlobalStarts = (iint*) calloc(size+1, sizeof(iint));
  for(iint r=0;r<size;++r)
    mesh->gatherGlobalStarts[r+1] = mesh->gatherGlobalStarts[r] + allGather[r];

  // shift numbering
  for(iint n=0;n<recvNtotal;++n)
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
  mesh->gatherLocalIds  = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->gatherBaseIds   = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->gatherBaseRanks = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->gatherHaloFlags = (iint*) calloc(localNodeCount, sizeof(iint));

  for(iint id=0;id<localNodeCount;++id){
    mesh->gatherLocalIds[id]  = sendNodes[id].element*mesh->Np+sendNodes[id].node;
    mesh->gatherBaseIds[id]   = sendNodes[id].newGlobalId;
    mesh->gatherBaseRanks[id] = sendNodes[id].baseRank;
    mesh->gatherHaloFlags[id] = sendNodes[id].haloFlag;
  }
  
  free(sendBuffer);
  free(sendNodes);

  //make a locally-ordered version
  mesh->globalIds      = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->globalOwners   = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->globalHaloFlags= (int*) calloc(localNodeCount, sizeof(int));

  for(iint id=0;id<localNodeCount;++id){
    iint localId = mesh->gatherLocalIds[id];
    mesh->globalIds[localId]      = mesh->gatherBaseIds[id];    
    mesh->globalOwners[localId]   = mesh->gatherBaseRanks[id];
    mesh->globalHaloFlags[localId]= mesh->gatherHaloFlags[id];
  }
}

/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"


// int rank for this process (not host thread safe)
int localRank = -1;

typedef struct{

  dlong element; // local element id
  int node;    // local node id
  int rank;    // rank of original node
  hlong id;      // original id
  int haloFlag;

  // info on base node (lowest rank node)
  dlong baseElement; 
  int baseNode;
  int baseRank;
  hlong baseId;

  hlong newGlobalId;

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

  int rank = localRank;

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
  rank = mesh->rank; 
  size = mesh->size; 
  localRank = rank;

  dlong localNodeCount = mesh->Np*mesh->Nelements;
  dlong *allLocalNodeCounts = (dlong*) calloc(size, sizeof(dlong));

  MPI_Allgather(&localNodeCount,    1, MPI_DLONG,
                allLocalNodeCounts, 1, MPI_DLONG,
                mesh->comm);
  
  hlong gatherNodeStart = 0;
  for(int r=0;r<rank;++r)
    gatherNodeStart += allLocalNodeCounts[r];
  
  free(allLocalNodeCounts);

  // form continuous node numbering (local=>virtual gather)
  parallelNode_t *sendNodes =
    (parallelNode_t*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,
                             sizeof(parallelNode_t));

  // use local numbering
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dlong id = e*mesh->Np+n;
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
      dlong id = e*mesh->Np + mesh->vertexNodes[v];
      hlong gid = mesh->EToV[e*mesh->Nverts+v] + 1;
      sendNodes[id].id = gid; 
      sendNodes[id].baseId = gid; 
    }

    // label halo flags
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1){
        for(int n=0;n<mesh->Nfp;++n){
          dlong id = e*mesh->Np+mesh->faceNodes[f*mesh->Nfp+n];
          sendNodes[id].haloFlag = 1;
        }
      }
    }
  }

  dlong localChange = 0, gatherChange = 1;

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
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
        dlong id  = e*mesh->Nfp*mesh->Nfaces + n;
        dlong idM = mesh->vmapM[id];
        dlong idP = mesh->vmapP[id];
        hlong gidM = sendNodes[idM].baseId;
        hlong gidP = sendNodes[idP].baseId;

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
    MPI_Allreduce(&localChange, &gatherChange, 1, MPI_DLONG, MPI_SUM, mesh->comm);
  }

  // sort based on base nodes (rank,element,node at base)
  qsort(sendNodes, localNodeCount, sizeof(parallelNode_t), parallelCompareOwners);

  // Make the MPI_PARALLELNODE_T data type
  MPI_Datatype MPI_PARALLELNODE_T;
  MPI_Datatype dtype[10] = {MPI_DLONG, MPI_INT, MPI_INT, MPI_HLONG, MPI_INT,
                            MPI_DLONG, MPI_INT, MPI_INT, MPI_HLONG, MPI_HLONG};
  int blength[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Aint addr[10], displ[10];
  MPI_Get_address ( &(sendNodes[0]            ), addr+0);
  MPI_Get_address ( &(sendNodes[0].node       ), addr+1);
  MPI_Get_address ( &(sendNodes[0].rank       ), addr+2);
  MPI_Get_address ( &(sendNodes[0].id         ), addr+3);
  MPI_Get_address ( &(sendNodes[0].haloFlag   ), addr+4);
  MPI_Get_address ( &(sendNodes[0].baseElement), addr+5);
  MPI_Get_address ( &(sendNodes[0].baseNode   ), addr+6);
  MPI_Get_address ( &(sendNodes[0].baseRank   ), addr+7);
  MPI_Get_address ( &(sendNodes[0].baseId     ), addr+8);
  MPI_Get_address ( &(sendNodes[0].newGlobalId), addr+9);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  displ[4] = addr[4] - addr[0];
  displ[5] = addr[5] - addr[0];
  displ[6] = addr[6] - addr[0];
  displ[7] = addr[7] - addr[0];
  displ[8] = addr[8] - addr[0];
  displ[9] = addr[9] - addr[0];
  MPI_Type_create_struct (10, blength, displ, dtype, &MPI_PARALLELNODE_T);
  MPI_Type_commit (&MPI_PARALLELNODE_T);

  // count how many nodes to send to each process
  int *sendCounts = (int *) calloc(size,sizeof(int));
  int *recvCounts = (int *) calloc(size,sizeof(int));
  int *sendOffsets = (int *) calloc(size+1,sizeof(int));
  int *recvOffsets = (int *) calloc(size+1,sizeof(int));

  for(dlong n=0;n<localNodeCount;++n) 
    sendCounts[sendNodes[n].baseRank]++;

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, mesh->comm);
  
  // find send and recv offsets for gather
  dlong recvNtotal = 0;
  for(int r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r];
  }

  parallelNode_t *recvNodes;
  if (recvNtotal) recvNodes = (parallelNode_t*) calloc(recvNtotal, sizeof(parallelNode_t));
  
  // load up node data to send 
  MPI_Alltoallv(sendNodes, sendCounts, sendOffsets, MPI_PARALLELNODE_T,
                recvNodes, recvCounts, recvOffsets, MPI_PARALLELNODE_T,
                mesh->comm);

  // sort by global index shifting halo nodes to the end
  qsort(recvNodes, recvNtotal, sizeof(parallelNode_t), parallelCompareBaseNodes);

  // renumber unique nodes starting from 0 (need to be careful about zeros)
  dlong Ngather = 0;
  if (recvNtotal) recvNodes[0].newGlobalId = Ngather;
  for(dlong n=1;n<recvNtotal;++n){
    if(recvNodes[n].baseId!=recvNodes[n-1].baseId){ // new node
      ++Ngather;
    }
    recvNodes[n].newGlobalId = Ngather;
  }
  if (recvNtotal) ++Ngather; // increment to actual number of unique nodes on this rank

  // collect unique node counts from all processes
  dlong *allGather   = (dlong*) calloc(size, sizeof(dlong));
  MPI_Allgather(&Ngather, 1, MPI_DLONG, allGather, 1, MPI_DLONG, mesh->comm);

  // cumulative sum of unique node counts => starting node index for each process
  mesh->gatherGlobalStarts = (hlong*) calloc(size+1, sizeof(hlong));
  for(int r=0;r<size;++r)
    mesh->gatherGlobalStarts[r+1] = mesh->gatherGlobalStarts[r] + allGather[r];

  // shift numbering
  for(dlong n=0;n<recvNtotal;++n)
    recvNodes[n].newGlobalId += mesh->gatherGlobalStarts[rank];
  
  // sort by rank, local index
  qsort(recvNodes, recvNtotal, sizeof(parallelNode_t), parallelCompareSourceRank);

  // reverse all to all to reclaim nodes
  MPI_Alltoallv(recvNodes, recvCounts, recvOffsets, MPI_PARALLELNODE_T,
                sendNodes, sendCounts, sendOffsets, MPI_PARALLELNODE_T,
                mesh->comm);

  // sort by rank, local index
  qsort(sendNodes, localNodeCount, sizeof(parallelNode_t), parallelCompareBaseNodes);

  // extract base index of each node (i.e. gather numbering)
  mesh->gatherLocalIds  = (dlong*) calloc(localNodeCount, sizeof(dlong));
  mesh->gatherBaseIds   = (hlong*) calloc(localNodeCount, sizeof(hlong));
  mesh->gatherBaseRanks = (int*) calloc(localNodeCount, sizeof(int));
  mesh->gatherHaloFlags = (int*) calloc(localNodeCount, sizeof(int));

  for(dlong id=0;id<localNodeCount;++id){
    mesh->gatherLocalIds[id]  = sendNodes[id].element*mesh->Np+sendNodes[id].node;
    mesh->gatherBaseIds[id]   = sendNodes[id].newGlobalId+1;
    mesh->gatherBaseRanks[id] = sendNodes[id].baseRank;
    mesh->gatherHaloFlags[id] = sendNodes[id].haloFlag;
  }
  
  MPI_Barrier(mesh->comm);
  MPI_Type_free(&MPI_PARALLELNODE_T);
  free(sendBuffer);
  free(sendNodes);

  //make a locally-ordered version
  mesh->globalIds      = (hlong*) calloc(localNodeCount, sizeof(hlong));
  mesh->globalOwners   = (int*) calloc(localNodeCount, sizeof(int));
  mesh->globalHaloFlags= (int*) calloc(localNodeCount, sizeof(int));

  for(dlong id=0;id<localNodeCount;++id){
    dlong localId = mesh->gatherLocalIds[id];
    mesh->globalIds[localId]      = mesh->gatherBaseIds[id];    
    mesh->globalOwners[localId]   = mesh->gatherBaseRanks[id];
    mesh->globalHaloFlags[localId]= mesh->gatherHaloFlags[id];
  }
}

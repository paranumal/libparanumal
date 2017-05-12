#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include <mpi.h>

#include "mesh.h"

typedef struct {

  iint localId;
  iint globalId;
  iint ownerRank;
  
}parallelNode_t;

// compare on global owners 
int parallelCompareOwnersAndGlobalId(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  iint rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if ((fa->ownerRank==rank)&&(fb->ownerRank!=rank)) return -1;
  if ((fa->ownerRank!=rank)&&(fb->ownerRank==rank)) return  1;

  if(fa->ownerRank < fb->ownerRank) return -1;
  if(fa->ownerRank > fb->ownerRank) return +1;

  if(fa->globalId < fb->globalId) return -1;
  if(fa->globalId > fb->globalId) return +1;

  return 0;  
}

// compare on global indices 
int parallelCompareGlobalId(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(fa->globalId < fb->globalId) return -1;
  if(fa->globalId > fb->globalId) return +1;

  return 0;  
}

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
hgs_t *meshParallelGatherSetup(mesh_t *mesh,    // provides DEVICE
				      iint Nlocal,     // number of local nodes
				      iint *globalNumbering,  // global index of nodes
				      iint *globalOwners){   // node owners

  iint rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  hgs_t *hgs = (hgs_t*) calloc(1, sizeof(hgs_t));

  hgs->Nlocal = Nlocal;

  parallelNode_t *nodes = (parallelNode_t *) calloc(Nlocal,sizeof(parallelNode_t));

  //populate node list
  for (iint n=0;n<Nlocal;n++) {
    nodes[n].localId = n;
    nodes[n].globalId = globalNumbering[n];
    nodes[n].ownerRank = globalOwners[n];
  }

  //sort by owning rank
  qsort(nodes, Nlocal, sizeof(parallelNode_t), parallelCompareOwnersAndGlobalId);

  //record this ordering
  hgs->sortIds = (iint *) calloc(Nlocal,sizeof(iint));
  for (iint n=0;n<Nlocal;n++)
    hgs->sortIds[n] = nodes[n].localId;

  //count how many nodes we're sending
  hgs->Nsend = (iint *) calloc(size,sizeof(iint));
  for (iint n=0;n<Nlocal;n++)    
    hgs->Nsend[nodes[n].ownerRank]++;

  //share how many nodes we're sending
  hgs->Nrecv = (iint *) calloc(size,sizeof(iint));
  MPI_Alltoall(hgs->Nsend, 1, MPI_IINT, hgs->Nrecv, 1, MPI_IINT, MPI_COMM_WORLD);

  hgs->NsendTotal=0;
  hgs->NrecvTotal=0;
  iint recvMessage=0,sendMessage=0;
  for (iint r=0;r<size;r++) {
    if (r!=rank) {
      hgs->NsendTotal += hgs->Nsend[r];
      hgs->NrecvTotal += hgs->Nrecv[r];
      sendMessage += (hgs->Nsend[r] > 0);
      recvMessage += (hgs->Nrecv[r] > 0);
    }
  }

  //these nodes are at the beginning of the list and don't need to be shared
  iint Nowned = Nlocal - hgs->NsendTotal;
  hgs->haloOffset = Nowned;

  hgs->haloSendRequests = calloc(sendMessage,sizeof(MPI_Request));
  hgs->haloRecvRequests = calloc(recvMessage,sizeof(MPI_Request));

  //this is how many nodes we'll end up with
  iint NownedTotal = Nowned + hgs->NrecvTotal;

  parallelNode_t *newNodes = (parallelNode_t *) calloc(NownedTotal,sizeof(parallelNode_t));
  memcpy(newNodes,nodes,Nowned*sizeof(parallelNode_t));

  char *sendBuffer = (char *)(nodes+Nowned);
  char *recvBuffer = (char *)(newNodes+Nowned);

  // initiate immediate send  and receives to each other process as needed
  int tag =999;
  recvMessage=0;
  iint recvOffset=0;
  for(iint r=0;r<size;++r){
    if(r!=rank){
      if(hgs->Nrecv[r]){
        MPI_Irecv(recvBuffer+recvOffset, hgs->Nrecv[r]*sizeof(parallelNode_t), MPI_CHAR, r, tag,
            MPI_COMM_WORLD, (MPI_Request*)hgs->haloRecvRequests+recvMessage);
        recvOffset += hgs->Nrecv[r]*sizeof(parallelNode_t);
        ++recvMessage;
      }
    }
  }
  sendMessage=0;
  iint sendOffset=0;
  for(iint r=0;r<size;++r){
    if(r!=rank){
      if(hgs->Nsend[r]){
        MPI_Isend(sendBuffer+sendOffset, hgs->Nsend[r]*sizeof(parallelNode_t), MPI_CHAR, r, tag,
            MPI_COMM_WORLD, (MPI_Request*)hgs->haloSendRequests+sendMessage);
        sendOffset += hgs->Nsend[r]*sizeof(parallelNode_t);
        ++sendMessage;
      }
    }
  }

  // Wait for all sent messages to have left and received messages to have arrived
  MPI_Status *sendStatus = (MPI_Status*) calloc(sendMessage, sizeof(MPI_Status));
  MPI_Status *recvStatus = (MPI_Status*) calloc(recvMessage, sizeof(MPI_Status));
  
  MPI_Waitall(recvMessage, (MPI_Request*)hgs->haloRecvRequests, recvStatus);
  MPI_Waitall(sendMessage, (MPI_Request*)hgs->haloSendRequests, sendStatus);
  
  free(recvStatus);
  free(sendStatus);   

  free(nodes);

  //newNodes now conatains all the nodes which are owned by this rank. We just need to define the local gatherIds
  for (iint n=0;n<NownedTotal;n++)
    newNodes[n].localId = n;

  //sort by global index
  qsort(newNodes, NownedTotal, sizeof(parallelNode_t), parallelCompareGlobalId);

  //count the number of unique globalIds
  hgs->Ngather=0;
  if (NownedTotal) hgs->Ngather++;
  for (iint n=1;n<NownedTotal;n++){
    if (newNodes[n].globalId!=newNodes[n-1].globalId) hgs->Ngather++;
  }

  //count the number of globalIds to gather
  hgs->gatherOffsets = (iint *) calloc(hgs->Ngather+1,sizeof(iint));
  iint cnt =1;
  if (NownedTotal) hgs->gatherOffsets[1]++;
  for (iint n=1;n<NownedTotal;n++) {
    if (newNodes[n].globalId!=newNodes[n-1].globalId) {
      hgs->gatherOffsets[++cnt]++;
    } else {
      hgs->gatherOffsets[cnt]++;
    }
  }

  //create the degree vector using these counts
  hgs->invDegree = (dfloat *) calloc(hgs->Ngather,sizeof(dfloat));
  for (iint n=0;n<hgs->Ngather;n++)
    hgs->invDegree[n] = 1.0/hgs->gatherOffsets[n+1];

  //make offsets from running total of counts
  for (iint n=1;n<hgs->Ngather+1;n++)
    hgs->gatherOffsets[n] += hgs->gatherOffsets[n-1];
  
  //finally, record the new local index of each node to gather
  hgs->gatherLocalIds = (iint *) calloc(NownedTotal,sizeof(iint));
  for (iint n=0;n<NownedTotal;n++) 
    hgs->gatherLocalIds[n] = newNodes[n].localId;
  
  free(newNodes);

  //allocate device storage buffers
  hgs->o_sortIds = mesh->device.malloc(Nlocal*sizeof(iint),hgs->sortIds);
  hgs->o_gatherOffsets = mesh->device.malloc((hgs->Ngather+1)*sizeof(iint),hgs->gatherOffsets);
  hgs->o_gatherLocalIds = mesh->device.malloc(NownedTotal*sizeof(iint),hgs->gatherLocalIds);
  hgs->o_invDegree = mesh->device.malloc(hgs->Ngather*sizeof(dfloat),hgs->invDegree);



  iint tmpSize = mymax(Nlocal,NownedTotal);
  void *tmpBuffer = calloc(tmpSize,sizeof(dfloat));
  hgs->o_gatherTmp = mesh->device.malloc(tmpSize*sizeof(dfloat),tmpBuffer);
  free(tmpBuffer);

  hgs->sendBuffer = (dfloat *) calloc(hgs->NsendTotal,sizeof(dfloat));
  hgs->recvBuffer = (dfloat *) calloc(hgs->NrecvTotal,sizeof(dfloat));

  return hgs;
}

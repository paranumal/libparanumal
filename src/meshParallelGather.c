#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"
#include "mpi.h"


void meshParallelGather(mesh_t *mesh,
			       hgs_t *hgs, 
			       occa::memory &o_v,
			       occa::memory &o_gv){

  // MPI info
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  occaTimerTic(mesh->device,"hgsGatherKernel");
  
  // arrange nodes for exchange + gather
  mesh->getKernel(hgs->Nlocal, o_v, hgs->o_sortIds, hgs->o_gatherTmp); 
  
  if(hgs->NsendTotal||hgs->NrecvTotal){  
    // copy partially gathered halo data from DEVICE to HOST
    if (hgs->NsendTotal)
      hgs->o_gatherTmp.copyTo(hgs->sendBuffer, (uintptr_t) hgs->NsendTotal*sizeof(dfloat), hgs->haloOffset*sizeof(dfloat));

    iint tag = 999;

    // initiate immediate send  and receives to each other process as needed
    iint recvOffset = 0, recvMessage=0;
    for(iint r=0;r<size;++r){
      if(r!=rank){
        if(hgs->Nrecv[r]){
          MPI_Irecv(hgs->recvBuffer+recvOffset, hgs->Nrecv[r], MPI_DFLOAT, r, tag,
              MPI_COMM_WORLD, (MPI_Request*)hgs->haloRecvRequests+recvMessage);
          recvOffset += hgs->Nrecv[r];
          ++recvMessage;
        }
      }
    }
    iint sendOffset = 0, sendMessage=0;
    for(iint r=0;r<size;++r){
      if(r!=rank){
        if(hgs->Nsend[r]){
          MPI_Isend(hgs->sendBuffer+sendOffset, hgs->Nsend[r], MPI_DFLOAT, r, tag,
              MPI_COMM_WORLD, (MPI_Request*)hgs->haloSendRequests+sendMessage);
          sendOffset += hgs->Nsend[r];
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

    // copy halo data back from HOST to DEVICE
    if (hgs->NrecvTotal)
      hgs->o_gatherTmp.copyFrom(hgs->recvBuffer, hgs->NrecvTotal*sizeof(dfloat), hgs->haloOffset*sizeof(dfloat));
  }

  // gather on DEVICE
  mesh->gatherKernel(hgs->Ngather, hgs->o_gatherOffsets, hgs->o_gatherLocalIds, hgs->o_gatherTmp, o_gv);

  occaTimerToc(mesh->device,"hgsGatherKernel");
}

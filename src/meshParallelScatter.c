#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"

void meshParallelScatter(mesh_t *mesh,
			       hgs_t *hgs, 
			       occa::memory &o_v,
			       occa::memory &o_sv){

  // MPI info
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  occaTimerTic(mesh->device,"hgsScatterKernel");

  // scatter on DEVICE
  mesh->scatterKernel(hgs->Ngather, hgs->o_gatherOffsets, hgs->o_gatherLocalIds, o_v, hgs->o_gatherTmp);
  
  if(hgs->NsendTotal||hgs->NrecvTotal){  
    // copy halo data to HOST to DEVICE
    if (hgs->NrecvTotal)
      hgs->o_gatherTmp.copyTo(hgs->recvBuffer, hgs->NrecvTotal*sizeof(dfloat), hgs->haloOffset*sizeof(dfloat));

    int tag = 999;

    // initiate immediate send  and receives to each other process as needed
    int sendOffset = 0,sendMessage=0;
    for(int r=0;r<size;++r){
      if(r!=rank){
        if(hgs->Nsend[r]){
          MPI_Irecv(hgs->sendBuffer+sendOffset, hgs->Nsend[r], MPI_DFLOAT, r, tag,
              MPI_COMM_WORLD, (MPI_Request*)hgs->haloSendRequests+sendMessage);
          sendOffset += hgs->Nsend[r];
          ++sendMessage;
        }
      }
    }
    int recvOffset = 0,recvMessage=0;
    for(int r=0;r<size;++r){
      if(r!=rank){
        if(hgs->Nrecv[r]){
          MPI_Isend(hgs->recvBuffer+recvOffset, hgs->Nrecv[r], MPI_DFLOAT, r, tag,
              MPI_COMM_WORLD, (MPI_Request*)hgs->haloRecvRequests+recvMessage);
          recvOffset += hgs->Nrecv[r];
          ++recvMessage;
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

    // copy scattered halo data from DEVICE to HOST
    if (hgs->NsendTotal)
      hgs->o_gatherTmp.copyFrom(hgs->sendBuffer, hgs->NsendTotal*sizeof(dfloat), hgs->haloOffset*sizeof(dfloat));
  }

  // extract scattered nodes
  mesh->putKernel(hgs->Nlocal, hgs->o_gatherTmp, hgs->o_sortIds, o_sv); // subv = v[ids]

  occaTimerToc(mesh->device,"hgsScatterKernel");
}

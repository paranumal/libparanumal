#include "parAlmond.h"


dcsr *newDCSR(almond_t *almond, csr *B){

  dcsr *A = (dcsr *) calloc(1,sizeof(dcsr));

  A->Nrows	= B->Nrows;
  A->Ncols	= B->Ncols;
  A->nnz		= B->nnz;

  //copy to device
  if(B->nnz){
    A->o_rowStarts = almond->device.malloc((A->Nrows+1)*sizeof(iint), B->rowStarts);
    A->o_cols      = almond->device.malloc(A->nnz*sizeof(iint), B->cols);
    A->o_coefs     = almond->device.malloc(A->nnz*sizeof(dfloat), B->coefs);
  }

  A->NsendTotal = B->NsendTotal;
  A->NrecvTotal = B->NrecvTotal;

  A->numLocalIds = B->numLocalIds;
  A->NHalo = B->NHalo;
  A->colMap = B->colMap;
  A->NrecvTotal = B->NrecvTotal;
  A->NsendTotal = B->NsendTotal;
  A->haloElementList = B->haloElementList;
  if (A->NsendTotal) 
    A->o_haloElementList = almond->device.malloc(A->NsendTotal*sizeof(iint),A->haloElementList);
  A->NsendPairs = B->NsendPairs;
  A->NrecvPairs = B->NrecvPairs;
  A->NsendMessages = B->NsendMessages;
  A->NrecvMessages = B->NrecvMessages;
  A->sendBuffer = B->sendBuffer;
  if (A->NrecvTotal) 
    A->recvBuffer = (dfloat*) malloc(A->NrecvTotal*sizeof(dfloat));
  if (A->NsendTotal) 
    A->o_haloBuffer = almond->device.malloc(A->NsendTotal*sizeof(dfloat),A->sendBuffer);

  A->haloSendRequests = B->haloSendRequests;
  A->haloRecvRequests = B->haloRecvRequests;

  return A;
}

void dcsrHaloExchangeStart(dcsr *A, size_t Nbytes, void *sendBuffer, void *recvBuffer) {
  // MPI info
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // count outgoing and incoming meshes
  iint tag = 999;

  // initiate immediate send  and receives to each other process as needed
  iint recvOffset = 0;
  iint sendOffset = 0;
  iint sendMessage = 0, recvMessage = 0;
  for(iint r=0;r<size;++r){
    if (A->NrecvTotal) {
      if(A->NrecvPairs[r]) {
        MPI_Irecv(((char*)A->recvBuffer)+recvOffset, A->NrecvPairs[r]*Nbytes, MPI_CHAR, r, tag,
            MPI_COMM_WORLD, (MPI_Request*)A->haloRecvRequests+recvMessage);
        recvOffset += A->NrecvPairs[r]*Nbytes;
        ++recvMessage;
      }
    }
    if (A->NsendTotal) {
      if(A->NsendPairs[r]){
        MPI_Isend(((char*)A->sendBuffer)+sendOffset, A->NsendPairs[r]*Nbytes, MPI_CHAR, r, tag,
            MPI_COMM_WORLD, (MPI_Request*)A->haloSendRequests+sendMessage);
        sendOffset += A->NsendPairs[r]*Nbytes;
        ++sendMessage;
      }
    }
  }
}

void dcsrHaloExchangeFinish(dcsr *A) {
  // Wait for all sent messages to have left and received messages to have arrived
  if (A->NsendTotal) {
    MPI_Status *sendStatus = (MPI_Status*) calloc(A->NsendMessages, sizeof(MPI_Status));
    MPI_Waitall(A->NsendMessages, (MPI_Request*)A->haloSendRequests, sendStatus);
    free(sendStatus);
  }
  if (A->NrecvTotal) {
    MPI_Status *recvStatus = (MPI_Status*) calloc(A->NrecvMessages, sizeof(MPI_Status));
    MPI_Waitall(A->NrecvMessages, (MPI_Request*)A->haloRecvRequests, recvStatus);
    free(recvStatus);
  }
}



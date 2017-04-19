#include "parAlmond.h"


hyb * newHYB(almond_t *almond, csr *csrA) {
  
  hyb *A = (hyb *) calloc(1,sizeof(hyb));

  A->Nrows	= csrA->Nrows;
  A->Ncols	= csrA->Ncols;


  A->diagInv = (dfloat *) calloc(A->Nrows, sizeof(dfloat));

  for(iint i=0; i<A->Nrows; i++)
  	A->diagInv[i] = 1.0/csrA->coefs[csrA->rowStarts[i]];

  dfloat *rowCounters = (dfloat*) calloc(Nrows, sizeof(dfloat));

  iint maxNnzPerRow = 0;
  iint minNnzPerRow = A->Ncols;
  for(iint i=0; i<Nrows; i++) {
    iint rowNnz = csrA->rowStarts[i+1] - csrA->rowStarts[i];
    rowCounters[i] = rowNnz;

    maxNnzPerRow = (rowNnz > maxNnzPerRow) ? rowNnz : maxNnzPerRow;
    minNnzPerRow = (rowNnz < minNnzPerRow) ? rowNnz : minNnzPerRow;
  }

  // create bins
  iint numBins = maxNnzPerRow - minNnzPerRow + 1;

  //zero row check
  if (numBins<0) numBins =0;

  iint *bins = (iint *) calloc(numBins, sizeof(iint));
  for(iint i=0; i<A->Nrows; i++){
    bins[rowCounters[i]-minNnzPerRow]++;
  }

  iint nnzPerRow = 0;
  iint nnz = 0;
  for(iint i=0; i<numBins; i++){
    nnz += bins[i] * (i+minNnzPerRow);
    if(nnz > 2.0*csrA->nnz/3.0){
    	nnzPerRow = i+minNnzPerRow;
    	break;
    }
  }

  A->E->Nrows	= A->Nrows;
  A->E->Ncols	= A->Ncols;
  A->E->nnzPerRow	= nnzPerRow;
  A->E->strideLength = A->Nrows;

  if(nnzPerRow){
    A->E->cols  = (iint *) calloc(A->Nrows*nnzPerRow, sizeof(iint));
    A->E->coefs = (dfloat *) calloc(A->Nrows*nnzPerRow, sizeof(dfloat));
  }

  iint nnzC = 0;

  // count the number of nonzeros to be stored in coo format
  for(iint i=0; i< A->Nrows; i++)
    if(rowCounters[i] > nnzPerRow) nnzC += (rowCounters[i] - nnzPerRow);

  A->E->actualNNZ  = csrA->nnz - nnzC;

  A->C->Nrows	= A->Nrows;
  A->C->Ncols	= A->Ncols;
  A->C->nnz		= nnzC;

  A->C->rows  = (iint *) calloc(nnzC, sizeof(iint));
  A->C->cols  = (iint *) calloc(nnzC, sizeof(iint));
  A->C->coefs = (dfloat *) calloc(nnzC, sizeof(dfloat));

  nnzC = 0;
  for(iint i=0; i<A->Nrows; i++){
    iint Jstart	= csrA->rowStarts[i];
    iint Jend	  = csrA->rowStarts[i+1];
    iint rowNnz	= Jend - Jstart;

    // store only min of nnzPerRow and rowNnz
    iint maxNnz = (nnzPerRow >= rowNnz) ? rowNnz : nnzPerRow;

    for(iint c=0; c<maxNnz; c++){
    	A->E->cols [i+c*E.strideLength()]	= csrA->cols[Jstart+c];
    	A->E->coefs[i+c*E.strideLength()]	= csrA->coefs[Jstart+c];
    }

    // store the remaining in coo format
    if(rowNnz > nnzPerRow){
    	for(iint c=nnzPerRow; c<rowNnz; c++){
    	  A->C->rows[nnzC]	= i;
    	  A->C->cols[nnzC]	= csrA->cols[Jstart+c];
    	  A->C->coefs[nnzC]	= csrA->coefs[Jstart+c];
    	  nnzC++;
    	}
    }
  }

  // copy the data to device memory
  if(A->Nrows)
    o_diagInv = almond->device.malloc(A->Nrows*sizeof(dfloat), A->diagInv);

  if(A->E->nnzPerRow){
    A->E->o_cols  = almond->device.malloc(A->Nrows*nnzPerRow*sizeof(iint), A->E->cols);
    A->E->o_coefs = almond->device.malloc(A->Nrows*nnzPerRow*sizeof(dfloat), A->E->coefs);
  }
  if(A->C->nnz){
    A->C->o_rows  = almond->device.malloc(A->C->nnz*sizeof(iint), A->C->rows);
    A->C->o_cols  = almond->device.malloc(A->C->nnz*sizeof(iint), A->C->cols);
    A->C->o_coefs = almond->device.malloc(A->C->nnz*sizeof(dfloat), A->C->coefs);

    const iint numBlocks = (A->C->nnz + almond->agmgBdim - 1)/almond->agmgBdim;

    iint *dummyRows = (iint *) calloc(numBlocks, sizeof(iint));
    dfloat *dummyAx = (dfloat *) calloc(numBlocks, sizeof(dfloat));

    A->C->temp_rows = almond->device.malloc(numBlocks*sizeof(iint), dummyRows);
    A->C->temp_Ax   = almond->device.malloc(numBlocks*sizeof(dfloat), dummyAx);
  }

  A->NsendTotal = csrA->NsendTotal;
  A->NrecvTotal = csrA->NrecvTotal;

  A->numLocalIds = csrA->numLocalIds;
  A->NHalo = csrA->NHalo;
  A->colMap = csrA->colMap;
  A->NrecvTotal = csrA->NrecvTotal;
  A->NsendTotal = csrA->NsendTotal;
  A->haloElementList = csrA->haloElementList;
  if (A->NsendTotal) A->o_haloElementList = almond->device.malloc(A->NsendTotal*sizeof(iint),A->haloElementList);
  A->NsendPairs = csrA->NsendPairs;
  A->NrecvPairs = csrA->NrecvPairs;
  A->NsendMessages = csrA->NsendMessages;
  A->NrecvMessages = csrA->NrecvMessages;
  A->sendBuffer = csrA->sendBuffer;
  if (A->NrecvTotal) A->recvBuffer = (dfloat *) malloc(A->NrecvTotal*sizeof(dfloat));
  A->haloSendRequests = csrA->haloSendRequests;
  A->haloRecvRequests = csrA->haloRecvRequests;

  if (A->NsendTotal) A->o_haloBuffer = almond->device.malloc(A->NsendTotal*sizeof(dfloat),A->sendBuffer);
}


void axpy(almond_t *almond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y) {

  if (A->NsendTotal) {
    haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);
  
    //copy from device
    A->o_haloBuffer.copyTo(A->sendBuffer);
  }

  if (A->NsendTotal+A->NrecvTotal) {
    // start halo exchange
    hybHaloExchangeStart(A, sizeof(dfloat),A->sendBuffer, A->recvBuffer);

    // immediately end the exchange TODO: make the exchange async using the A = E + C hybrid form
    hybHaloExchangeFinish(A);
  }

  if (A->NrecvTotal) {
    //copy back to device
    o_x.copyFrom(A->recvBuffer,A->NrecvTotal*sizeof(dfloat),A->numLocalIds*sizeof(dfloat));
  }

  // y <-- alpha*E*x+beta*y
  axpy(almond, A->E, alpha, o_x, beta, o_y);

  // y <-- alpha*C*x + y
  if(A->C->nnz)
    ax(almond, A->C, alpha, o_x, o_y);
}


void zeqaxpy(almond_t *almond, hyb *A, dfloat alpha, occa::memory o_x,
	  dfloat beta,  occa::memory o_y, occa::memory o_z) {

  if (A->NsendTotal) {
    haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);
  
    //copy from device
    A->o_haloBuffer.copyTo(A->sendBuffer);
  }

  if (A->NsendTotal+A->NrecvTotal) {
    // start halo exchange
    hybHaloExchangeStart(A, sizeof(dfloat),A->sendBuffer, A->recvBuffer);

    // immediately end the exchange TODO: make the exchange async using the A = E + C hybrid form
    hybHaloExchangeFinish(A);
  }

  if (A->NrecvTotal) {
    //copy back to device
    o_x.copyFrom(A->recvBuffer,A->NrecvTotal*sizeof(dfloat),A->numLocalIds*sizeof(dfloat));
  }

  // z <-- alpha*E*x+ beta*y
  zeqaxpy(almond, A->E, alpha, o_x, beta, o_y, o_z);

  // z <-- alpha*C*x + z
  if(A->C->nnz)
    ax(almond, A->C, alpha, o_x, o_z);
}


void smoothJacobi(almond_t *almond, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero) {

  if(x_is_zero){
    dfloat alpha = 1.0;
    dfloat beta = 0.0;
    dotStar(almond, A->Nrows, alpha, A->o_diagInv, o_r, beta, o_x);
    return;
  }


  if (A->NsendTotal) {
    haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);
  
    //copy from device
    A->o_haloBuffer.copyTo(A->sendBuffer);
  }

  if (A->NsendTotal+A->NrecvTotal) {
    // start halo exchange
    hybHaloExchangeStart(A, sizeof(dfloat),A->sendBuffer, A->recvBuffer);

    // immediately end the exchange TODO: make the exchange async using the A = E + C hybrid form
    hybHaloExchangeFinish(A);
  }

  if (A->NrecvTotal) {
    //copy back to device
    o_x.copyFrom(A->recvBuffer,A->NrecvTotal*sizeof(dfloat),A->numLocalIds*sizeof(dfloat));
  }


  if(A->o_temp1.bytes()== 0){
    dfloat *dummy = calloc(A->Nrows, sizeof(dfloat));
    A->o_temp1 = almond->device.malloc(A->Nrows*sizeof(dfloat), dummy);
  }


  const iint numBlocks = (A->Nrows+almond->agmgBdim-1)/almond->agmgBdim;

  //occa::tic("ellJacobi1");
  almond->ellJacobi1Kernel(numBlocks, almond->agmgBdim, A->Nrows, A->E->nnzPerRow, A->E->strideLength, 
                  A->E->o_cols, A->E->o_coefs, o_x, o_r, A->o_temp1);

  // temp1 += -C*x
  if(A->C->nnz)
    ax(almond, A->C, -1.0, o_x, A->o_temp1);

  // x = invD*temp1
  dotStar(almond, A->Nrows, 1.0, A->o_diagInv, A->o_temp1, 0., o_x);
}


void smoothDampedJacobi(almond_t *almond,
         hyb *A,
         occa::memory o_r,
		     occa::memory o_x,
		     dfloat alpha,
		     bool x_is_zero){

  if(x_is_zero){
    dfloat beta = 0.0;
    dotStar(almond, A->Nrows, alpha, A->o_diagInv, o_r, beta, o_x);
    return;
  }

  if (A->NsendTotal) {
    haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);
  
    //copy from device
    A->o_haloBuffer.copyTo(A->sendBuffer);
  }

  if (A->NsendTotal+A->NrecvTotal) {
    // start halo exchange
    hybHaloExchangeStart(A, sizeof(dfloat),A->sendBuffer, A->recvBuffer);

    // immediately end the exchange TODO: make the exchange async using the A = E + C hybrid form
    hybHaloExchangeFinish(A);
  }

  if (A->NrecvTotal) {
    //copy back to device
    o_x.copyFrom(A->recvBuffer,A->NrecvTotal*sizeof(dfloat),A->numLocalIds*sizeof(dfloat));
  }


  if(A->o_temp1.bytes()== 0){
    dfloat *dummy = calloc(A->Nrows, sizeof(dfloat));
    A->o_temp1 = almond->device.malloc(A->Nrows*sizeof(dfloat), dummy);
  }


  const iint numBlocks = (A->Nrows+almond->agmgBdim-1)/almond->agmgBdim;

  ellJacobi1Kernel(numBlocks, almond->agmgBdim, A->Nrows, A->E->nnzPerRow, A->E->strideLength, 
                    A->E->o_cols, A->E->o_coefs, o_x, o_r, A->o_temp1);

  // temp1 += -C*x
  if(A->C->nnz)
    ax(almond, A->C, -1.0, o_x, A->o_temp1);

  // x = alpha*invD*temp1 + (1-alpha)*x
  const dfloat beta = 1.0 - alpha;
  dotStar(almond, A->Nrows, alpha, A->o_diagInv, A->o_temp1, beta, o_x);
}


void hybHaloExchangeStart(hyb *A, size_t Nbytes, void *sendBuffer, void *recvBuffer) {
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
        MPI_Irecv(((char*)recvBuffer)+recvOffset, A->NrecvPairs[r]*Nbytes, MPI_CHAR, r, tag,
            MPI_COMM_WORLD, (MPI_Request*)A->haloRecvRequests+recvMessage);
        recvOffset += A->NrecvPairs[r]*Nbytes;
        ++recvMessage;
      }
    }
    if (A->NsendTotal) {
      if(A->NsendPairs[r]){
        MPI_Isend(((char*)sendBuffer)+sendOffset, A->NsendPairs[r]*Nbytes, MPI_CHAR, r, tag,
            MPI_COMM_WORLD, (MPI_Request*)A->haloSendRequests+sendMessage);
        sendOffset += A->NsendPairs[r]*Nbytes;
        ++sendMessage;
      }
    }
  }
}

void hybHaloExchangeFinish(hyb *A) {
  // Wait for all sent messages to have left and received messages to have arrived
  if (A->NsendTotal) {
    MPI_Status *sendStatus = (MPI_Status*) calloc(A->NsendMessages, sizeof(MPI_Status));
    MPI_Waitall(NsendMessages, (MPI_Request*)A->haloSendRequests, sendStatus);
    free(sendStatus);
  }
  if (A->NrecvTotal) {
    MPI_Status *recvStatus = (MPI_Status*) calloc(A->NrecvMessages, sizeof(MPI_Status));
    MPI_Waitall(NrecvMessages, (MPI_Request*)A->haloRecvRequests, recvStatus);
    free(recvStatus);
  }
}


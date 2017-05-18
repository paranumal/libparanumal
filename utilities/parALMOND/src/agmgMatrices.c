#include "parAlmond.h"


csr * newCSR(iint Nrows, iint Ncols, iint nnz,
      iint *rowStarts, iint *cols, dfloat *coefs){

  csr *A = (csr *) calloc(1,sizeof(csr));

  A->Nrows = Nrows;
  A->Ncols = Ncols;
  A->nnz   = nnz;

  A->rowStarts = (iint *) calloc(A->Nrows+1,sizeof(iint));
  A->cols      = (iint *) calloc(A->nnz, sizeof(iint));
  A->coefs     = (dfloat *) calloc(A->nnz, sizeof(dfloat));


  for (iint i=0; i<Nrows; i++) {
    iint start = rowStarts[i];
    A->rowStarts[i] = start; 
    iint cnt = 1;
    for (iint j=rowStarts[i]; j<rowStarts[i+1]; j++) {
      if (cols[j] == i) { //diagonal entry
        A->cols[start] = cols[j];
        A->coefs[start] = coefs[j];
      } else {
        A->cols[start+cnt] = cols[j];
        A->coefs[start+cnt] = coefs[j];
        cnt++;
      }
    }
  }
  A->rowStarts[Nrows] = rowStarts[Nrows];

  A->NsendTotal = 0;
  A->NrecvTotal = 0;
  A->NHalo = 0;

  return A;
}

void freeCSR(csr *A) {
  free(A->rowStarts);
  free(A->cols);
  free(A->coefs);

  free(A);
}

dcsr *newDCSR(parAlmond_t *parAlmond, csr *B){

  dcsr *A = (dcsr *) calloc(1,sizeof(dcsr));

  A->Nrows  = B->Nrows;
  A->Ncols  = B->Ncols;
  A->nnz    = B->nnz;

  //copy to device
  if(B->nnz){
    A->o_rowStarts = parAlmond->device.malloc((A->Nrows+1)*sizeof(iint), B->rowStarts);
    A->o_cols      = parAlmond->device.malloc(A->nnz*sizeof(iint), B->cols);
    A->o_coefs     = parAlmond->device.malloc(A->nnz*sizeof(dfloat), B->coefs);
  }

  if(B->Nrows) {
    dfloat *diagInv = (dfloat *) calloc(A->Nrows, sizeof(dfloat));

    for(iint i=0; i<B->Nrows; i++)
      diagInv[i] = 1.0/B->coefs[B->rowStarts[i]];
    A->o_diagInv = parAlmond->device.malloc(A->Nrows*sizeof(dfloat), diagInv);
    A->o_temp1 = parAlmond->device.malloc(A->Nrows*sizeof(dfloat), diagInv);
    free(diagInv); 
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
    A->o_haloElementList = parAlmond->device.malloc(A->NsendTotal*sizeof(iint),A->haloElementList);
  A->NsendPairs = B->NsendPairs;
  A->NrecvPairs = B->NrecvPairs;
  A->NsendMessages = B->NsendMessages;
  A->NrecvMessages = B->NrecvMessages;
  A->sendBuffer = B->sendBuffer;
  if (A->NrecvTotal) 
    A->recvBuffer = (dfloat*) malloc(A->NrecvTotal*sizeof(dfloat));
  if (A->NsendTotal) 
    A->o_haloBuffer = parAlmond->device.malloc(A->NsendTotal*sizeof(dfloat),A->sendBuffer);

  A->haloSendRequests = B->haloSendRequests;
  A->haloRecvRequests = B->haloRecvRequests;

  return A;
}

hyb * newHYB(parAlmond_t *parAlmond, csr *csrA) {

  hyb *A = (hyb *) calloc(1,sizeof(hyb));

  A->Nrows  = csrA->Nrows;
  A->Ncols  = csrA->Ncols;

  dfloat *diagInv;
  iint *rowCounters;
  if (csrA->Nrows) {
    diagInv = (dfloat *) calloc(csrA->Nrows, sizeof(dfloat));
    rowCounters = (iint*) calloc(csrA->Nrows, sizeof(iint));
  }

  for(iint i=0; i<csrA->Nrows; i++)
    diagInv[i] = 1.0/csrA->coefs[csrA->rowStarts[i]];

  iint maxNnzPerRow = 0;
  iint minNnzPerRow = csrA->Ncols;
  for(iint i=0; i<csrA->Nrows; i++) {
    iint rowNnz = csrA->rowStarts[i+1] - csrA->rowStarts[i];
    rowCounters[i] = rowNnz;

    maxNnzPerRow = (rowNnz > maxNnzPerRow) ? rowNnz : maxNnzPerRow;
    minNnzPerRow = (rowNnz < minNnzPerRow) ? rowNnz : minNnzPerRow;
  }

  // create bins
  iint numBins = maxNnzPerRow - minNnzPerRow + 1;

  //zero row check
  if (numBins<0) numBins =0;

  iint *bins;
  if (numBins)
    bins = (iint *) calloc(numBins, sizeof(iint));

  for(iint i=0; i<csrA->Nrows; i++){
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

  A->E = (ell *) calloc(1, sizeof(ell));

  A->E->Nrows = csrA->Nrows;
  A->E->Ncols = csrA->Ncols;
  A->E->nnzPerRow = nnzPerRow;
  A->E->strideLength = csrA->Nrows;

  iint *Ecols;
  dfloat *Ecoefs;
  if(nnzPerRow){
    Ecols  = (iint *) calloc(csrA->Nrows*nnzPerRow, sizeof(iint));
    Ecoefs = (dfloat *) calloc(csrA->Nrows*nnzPerRow, sizeof(dfloat));
  }

  iint nnzC = 0;

  // count the number of nonzeros to be stored in coo format
  for(iint i=0; i< csrA->Nrows; i++)
    if(rowCounters[i] > nnzPerRow) nnzC += (rowCounters[i] - nnzPerRow);

  A->E->actualNNZ  = csrA->nnz - nnzC;

  A->C = (coo *) calloc(1, sizeof(coo));

  A->C->Nrows = csrA->Nrows;
  A->C->Ncols = csrA->Ncols;
  A->C->nnz   = nnzC;

  iint *Coffsets;
  iint *Ccols;
  dfloat *Ccoefs;

  Coffsets = (iint *) calloc(csrA->Nrows+1, sizeof(iint));
  if (nnzC) {
    Ccols    = (iint *) calloc(nnzC, sizeof(iint));
    Ccoefs   = (dfloat *) calloc(nnzC, sizeof(dfloat));
  }

  nnzC = 0;
  for(iint i=0; i<csrA->Nrows; i++){
    iint Jstart = csrA->rowStarts[i];
    iint Jend   = csrA->rowStarts[i+1];
    iint rowNnz = Jend - Jstart;

    // store only min of nnzPerRow and rowNnz
    iint maxNnz = (nnzPerRow >= rowNnz) ? rowNnz : nnzPerRow;

    for(iint c=0; c<maxNnz; c++){
      Ecols [i+c*A->E->strideLength]  = csrA->cols[Jstart+c];
      Ecoefs[i+c*A->E->strideLength]  = csrA->coefs[Jstart+c];
    }

    // store the remaining in coo format
    if(rowNnz > nnzPerRow){
      for(iint c=nnzPerRow; c<rowNnz; c++){
        Coffsets[i+1]++;
        Ccols[nnzC]   = csrA->cols[Jstart+c];
        Ccoefs[nnzC]  = csrA->coefs[Jstart+c];
        nnzC++;
      }
    }
  }

  //use counts to create offsets
  for (iint i=0;i<csrA->Nrows;i++)
    Coffsets[i+1] += Coffsets[i];

  // copy the data to device memory
  if(csrA->Nrows) {
    A->o_diagInv = parAlmond->device.malloc(csrA->Nrows*sizeof(dfloat), diagInv);
    A->o_temp1 = parAlmond->device.malloc(csrA->Nrows*sizeof(dfloat), diagInv);
    free(diagInv); free(rowCounters); free(bins);
  }

  if(A->E->nnzPerRow){
    A->E->o_cols  = parAlmond->device.malloc(csrA->Nrows*nnzPerRow*sizeof(iint), Ecols);
    A->E->o_coefs = parAlmond->device.malloc(csrA->Nrows*nnzPerRow*sizeof(dfloat), Ecoefs);
    free(Ecols); free(Ecoefs);
  }
  
  if(A->C->nnz){
    A->C->o_offsets = parAlmond->device.malloc((csrA->Nrows+1)*sizeof(iint), Coffsets);
    A->C->o_cols    = parAlmond->device.malloc(A->C->nnz*sizeof(iint), Ccols);
    A->C->o_coefs   = parAlmond->device.malloc(A->C->nnz*sizeof(dfloat), Ccoefs);

    free(Ccols); free(Ccoefs);
  }

  free(Coffsets);

  A->NsendTotal = csrA->NsendTotal;
  A->NrecvTotal = csrA->NrecvTotal;

  A->numLocalIds = csrA->numLocalIds;
  A->NHalo = csrA->NHalo;
  A->colMap = csrA->colMap;
  A->NrecvTotal = csrA->NrecvTotal;
  A->NsendTotal = csrA->NsendTotal;
  A->haloElementList = csrA->haloElementList;
  if (A->NsendTotal) A->o_haloElementList = parAlmond->device.malloc(A->NsendTotal*sizeof(iint),A->haloElementList);
  A->NsendPairs = csrA->NsendPairs;
  A->NrecvPairs = csrA->NrecvPairs;
  A->NsendMessages = csrA->NsendMessages;
  A->NrecvMessages = csrA->NrecvMessages;
  A->sendBuffer = csrA->sendBuffer;
  if (A->NrecvTotal) A->recvBuffer = (dfloat *) malloc(A->NrecvTotal*sizeof(dfloat));
  A->haloSendRequests = csrA->haloSendRequests;
  A->haloRecvRequests = csrA->haloRecvRequests;

  if (A->NsendTotal) A->o_haloBuffer = parAlmond->device.malloc(A->NsendTotal*sizeof(dfloat),A->sendBuffer);

  return A;
}


void axpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y) {

  occa::tic("csr axpy");
  csrHaloExchange(A, sizeof(dfloat), x, A->sendBuffer, x+A->numLocalIds);

  // y[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  #pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){
    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    dfloat result = 0.0;
    for(iint jj=Jstart; jj<Jend; jj++)
      result += (A->coefs[jj]*x[A->cols[jj]]);
    
    y[i] = alpha*result + beta*y[i];
  }
  occa::toc("csr axpy");
}

void zeqaxpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y, dfloat *z) {

  occa::tic("csr zeqaxpy");
  csrHaloExchange(A, sizeof(dfloat), x, A->sendBuffer, x+A->numLocalIds);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  #pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){
    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    dfloat result = 0.0;

    for(iint jj=Jstart; jj<Jend; jj++)
      result += A->coefs[jj]*x[A->cols[jj]];

    z[i] = alpha*result + beta*y[i];
  }
  occa::toc("csr zeqaxpy");
}

void axpy(parAlmond_t *parAlmond, dcsr *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y) {

  occaTimerTic(parAlmond->device,"dcsr axpy");
  if (A->NsendTotal) {
    parAlmond->haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);
    
    //copy from device
    A->o_haloBuffer.copyTo(A->sendBuffer);
  }

  if (A->NsendTotal + A->NrecvTotal) {
    // start halo exchange
    dcsrHaloExchangeStart(A, sizeof(dfloat), A->sendBuffer, A->recvBuffer);

    // immediately end the exchange TODO: make the exchange async using the A = E + C hybrid form
    dcsrHaloExchangeFinish(A);
  }

  if(A->NrecvTotal) {
    //copy back to device
    o_x.copyFrom(A->recvBuffer,A->NrecvTotal*sizeof(dfloat),
                  A->numLocalIds*sizeof(dfloat));
  }

  parAlmond->agg_interpolateKernel(A->Nrows, A->o_cols, A->o_coefs, o_x, o_y);
  //parAlmond->dcsrAXPYKernel(A->Nrows, alpha,beta, A->o_rowStarts, 
  //                                A->o_cols, A->o_coefs, o_x, o_y);
  occaTimerToc(parAlmond->device,"dcsr axpy");
}

void axpy(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y) {

  occaTimerTic(parAlmond->device,"hyb axpy");
  if (A->NsendTotal) {
    parAlmond->haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);
  
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
  axpy(parAlmond, A->E, alpha, o_x, beta, o_y);

  // y <-- alpha*C*x + y
  if(A->C->nnz)
    ax(parAlmond, A->C, alpha, o_x, o_y);

  occaTimerToc(parAlmond->device,"hyb axpy");
}


void zeqaxpy(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x,
    dfloat beta,  occa::memory o_y, occa::memory o_z) {

  occaTimerTic(parAlmond->device,"hyb zeqaxpy");
  if (A->NsendTotal) {
    parAlmond->haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);
  
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
  zeqaxpy(parAlmond, A->E, alpha, o_x, beta, o_y, o_z);

  // z <-- alpha*C*x + z
  if(A->C->nnz)
    ax(parAlmond, A->C, alpha, o_x, o_z);

  occaTimerToc(parAlmond->device,"hyb zeqaxpy");
}

void axpy(parAlmond_t *parAlmond, ell *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y) {

  if(A->Nrows){
    occaTimerTic(parAlmond->device,"ell axpy");
    parAlmond->ellAXPYKernel(A->Nrows, A->nnzPerRow, A->strideLength, 
                          alpha, beta, A->o_cols, A->o_coefs, o_x, o_y);
    occaTimerToc(parAlmond->device,"ell axpy");
  }
}

void zeqaxpy(parAlmond_t *parAlmond, ell *A, dfloat alpha, occa::memory o_x, 
            dfloat beta, occa::memory o_y,  occa::memory o_z) {

  if(A->Nrows){
    occaTimerTic(parAlmond->device,"ell zeqaxpy");
    parAlmond->ellZeqAXPYKernel(A->Nrows, A->nnzPerRow, A->strideLength, 
                          alpha, beta, A->o_cols, A->o_coefs, o_x, o_y, o_z);
    occaTimerToc(parAlmond->device,"ell zeqaxpy");
  }
}

void ax(parAlmond_t *parAlmond, coo *C, dfloat alpha, occa::memory o_x, occa::memory o_y) {

  // do block-wise product
  if(C->nnz){
    occaTimerTic(parAlmond->device,"coo ax");
    parAlmond->cooAXKernel(C->Nrows, alpha, C->o_offsets, C->o_cols, C->o_coefs,o_x, o_y);
    occaTimerToc(parAlmond->device,"coo ax");
  }
}

void matFreeZeqAXPY(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, 
              occa::memory o_y, occa::memory o_z) {

  occaTimerTic(parAlmond->device,"matfree zeqaxpy");
  
  //matrix free A action, temp1 = Ax
  parAlmondMatrixFreeAX(parAlmond, o_x, A->o_temp1);

  //z = alpha*Ax + beta*y
  vectorAdd(parAlmond, A->Nrows, alpha, A->o_temp1, beta, o_y, o_z);

  occaTimerToc(parAlmond->device,"matfree zeqaxpy");
}

void smoothJacobi(csr *A, dfloat *r, dfloat *x, bool x_is_zero) {

  occa::tic("csr smoothJacobi");
  // x = inv(D)*(b-R*x)  where R = A-D
  if(x_is_zero){
    #pragma omp parallel for
    for(iint i=0; i<A->Nrows; i++){
      dfloat invD = 1.0/A->coefs[A->rowStarts[i]];
      x[i] = invD*r[i];
    }
    occa::toc("csr smoothJacobi");
    return;
  }

  csrHaloExchange(A, sizeof(dfloat), x, A->sendBuffer, x+A->numLocalIds);

  dfloat y[A->Nrows]; //TODO get rid of this stack alloc

#pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){
    dfloat result = r[i];

    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    iint jj = Jstart;

    const dfloat invD = 1./A->coefs[jj++];

    for(; jj<Jend; jj++)
      result -= A->coefs[jj]*x[A->cols[jj]];

    y[i] = invD*result;
  }

  // copy the buffer vector to x
#pragma omp parallel for
  for (iint i=0;i<A->Nrows;i++)
    x[i] = y[i];

  occa::toc("csr smoothJacobi");
}


void smoothDampedJacobi(csr *A, dfloat *r, dfloat *x, dfloat alpha, bool x_is_zero) {
  
  occa::tic("csr smoothDampedJacobi");
  if(x_is_zero){
#pragma omp parallel for
    for(iint i=0; i<A->Nrows; i++){
      dfloat invD = 1.0/A->coefs[A->rowStarts[i]];
      x[i] = alpha*invD*r[i];
    }
    occa::toc("csr smoothDampedJacobi");
    return;
  }

  csrHaloExchange(A, sizeof(dfloat), x, A->sendBuffer, x+A->numLocalIds);

  // x = (1-alpha)*x + alpha*inv(D) * (b-R*x) where R = A-D
  dfloat y[A->Nrows];

  const dfloat oneMalpha = 1. - alpha;

#pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){
    dfloat result = r[i];

    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    iint jj = Jstart;

    const dfloat invD = 1./A->coefs[jj++];

    for(; jj<Jend; jj++)
      result -= A->coefs[jj]*x[A->cols[jj]];

    y[i] = oneMalpha*x[i] + alpha*invD*result;
  }

  // copy the buffer vector to x
#pragma omp parallel for
  for (iint i=0;i<A->Nrows;i++)
    x[i] = y[i];

  occa::toc("csr smoothDampedJacobi");
}

void smoothJacobi(parAlmond_t *parAlmond, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero) {

  occaTimerTic(parAlmond->device,"hyb smoothJacobi");
  if(x_is_zero){
    if (A->Nrows)
      dotStar(parAlmond, A->Nrows, 1.0, A->o_diagInv, o_r, 0.0, o_x);
    occaTimerToc(parAlmond->device,"hyb smoothJacobi");
    return;
  }

  if (A->NsendTotal) {
    parAlmond->haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);
  
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

  occaTimerTic(parAlmond->device,"ellJacobi1");
  if (A->Nrows)
    parAlmond->ellJacobiKernel(A->Nrows, A->E->nnzPerRow, A->E->strideLength, 
                            A->E->o_cols, A->E->o_coefs, o_x, o_r, A->o_temp1);
  occaTimerToc(parAlmond->device,"ellJacobi1");

  // temp1 += -C*x
  if(A->C->nnz)
    ax(parAlmond, A->C, -1.0, o_x, A->o_temp1);

  // x = invD*temp1
  if (A->Nrows)
    dotStar(parAlmond, A->Nrows, 1.0, A->o_diagInv, A->o_temp1, 0., o_x);
  occaTimerToc(parAlmond->device,"hyb smoothJacobi");
}

void smoothDampedJacobi(parAlmond_t *parAlmond,
         hyb *A,
         occa::memory o_r,
         occa::memory o_x,
         dfloat alpha,
         bool x_is_zero){

  occaTimerTic(parAlmond->device,"hyb smoothDampedJacobi");
  if(x_is_zero){
    if (A->Nrows)
      dotStar(parAlmond, A->Nrows, alpha, A->o_diagInv, o_r, 0.0, o_x);
    occaTimerToc(parAlmond->device,"hyb smoothDampedJacobi");
    return;
  }

  if (A->NsendTotal) {
    parAlmond->haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);
  
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

  occaTimerTic(parAlmond->device,"ellJacobi1");
  if (A->Nrows)
    parAlmond->ellJacobiKernel(A->Nrows, A->E->nnzPerRow, A->E->strideLength, 
                              A->E->o_cols, A->E->o_coefs, o_x, o_r, A->o_temp1);
  occaTimerToc(parAlmond->device,"ellJacobi1");

  // temp1 += -C*x
  if(A->C->nnz)
    ax(parAlmond, A->C, -1.0, o_x, A->o_temp1);

  // x = alpha*invD*temp1 + (1-alpha)*x
  const dfloat beta = 1.0 - alpha;
  if (A->Nrows)
    dotStar(parAlmond, A->Nrows, alpha, A->o_diagInv, A->o_temp1, beta, o_x);

  occaTimerToc(parAlmond->device,"hyb smoothDampedJacobi");
}

void matFreeSmoothJacobi(parAlmond_t *parAlmond, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero) {

  occaTimerTic(parAlmond->device,"matfree smoothJacobi");
  if(x_is_zero){
    dotStar(parAlmond, A->Nrows, 1., A->o_diagInv, o_r, 0., o_x);
    occaTimerToc(parAlmond->device,"matfree smoothJacobi");
    return;
  }

  //matrix free A action, temp1 = Ax
  parAlmondMatrixFreeAX(parAlmond, o_x, A->o_temp1);

  //temp1 = r-Ax
  vectorAdd(parAlmond, A->Nrows, 1.0, o_r, -1.0, A->o_temp1);

  // x = x + invD*temp1
  dotStar(parAlmond, A->Nrows, 1.0, A->o_diagInv, A->o_temp1, 1., o_x);
  occaTimerToc(parAlmond->device,"matfree smoothJacobi");
}


void matFreeSmoothDampedJacobi(parAlmond_t *parAlmond, hyb* A, occa::memory o_r,
                            occa::memory o_x, dfloat alpha, bool x_is_zero){

  occaTimerTic(parAlmond->device,"matfree smoothDampedJacobi");
  if(x_is_zero){
    dotStar(parAlmond, A->Nrows, alpha, A->o_diagInv, o_r, 0., o_x);
    occaTimerToc(parAlmond->device,"matfree smoothDampedJacobi");
    return;
  }

  //matrix free A action, temp1 = Ax
  parAlmondMatrixFreeAX(parAlmond, o_x, A->o_temp1);

  //temp1 = r-Ax
  vectorAdd(parAlmond, A->Nrows, 1.0, o_r, -1.0, A->o_temp1);

  // x = x + alpha*invD*temp1
  dotStar(parAlmond, A->Nrows, alpha, A->o_diagInv, A->o_temp1, 1., o_x);
  occaTimerToc(parAlmond->device,"matfree smoothDampedJacobi");
}

// set up halo infomation for inter-processor MPI 
// exchange of trace nodes
void csrHaloSetup(csr *A, iint *globalColStarts){

  // MPI info
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // non-blocking MPI isend/irecv requests (used in meshHaloExchange)
  A->haloSendRequests = calloc(size, sizeof(MPI_Request));
  A->haloRecvRequests = calloc(size, sizeof(MPI_Request));
  
  // count number of halo element nodes to swap
  A->NrecvTotal = 0;
  A->NsendPairs = (iint*) calloc(size, sizeof(int));
  A->NrecvPairs = (iint*) calloc(size, sizeof(int));
  for(iint n=A->numLocalIds;n<A->Ncols;++n){ //for just the halo
    iint id = A->colMap[n]; // global index
    for (iint r=0;r<size;r++) { //find owner's rank
      if (globalColStarts[r]-1<id && id < globalColStarts[r+1]) {
        A->NrecvTotal++;
        A->NrecvPairs[r]++;  
      }
    }
  }

  MPI_Alltoall(A->NrecvPairs, 1, MPI_IINT, A->NsendPairs, 1, MPI_IINT, MPI_COMM_WORLD);

  A->NsendTotal = 0;
  for (iint r=0;r<size;r++)
    A->NsendTotal += A->NsendPairs[r];

  A->haloElementList = (iint *) calloc(A->NsendTotal,sizeof(iint));

  // count number of MPI messages in halo exchange
  A->NsendMessages = 0;
  A->NrecvMessages = 0;
  for(iint r=0;r<size;++r) {
    if(A->NsendPairs[r])
      A->NsendMessages++;
    if(A->NrecvPairs[r])
      A->NrecvMessages++;
  }

  //exchange the needed ids
  iint tag = 999;
  iint recvOffset = A->numLocalIds;
  iint sendOffset = 0;
  iint sendMessage = 0, recvMessage = 0;
  for(iint r=0;r<size;++r){
     if(A->NsendPairs[r]) {
      MPI_Irecv(A->haloElementList+sendOffset, A->NsendPairs[r], MPI_IINT, r, tag,
          MPI_COMM_WORLD, (MPI_Request*)A->haloSendRequests+sendMessage);
      sendOffset += A->NsendPairs[r];
      ++sendMessage;
    }
    if(A->NrecvPairs[r]){
      MPI_Isend(A->colMap+recvOffset, A->NrecvPairs[r], MPI_IINT, r, tag,
          MPI_COMM_WORLD, (MPI_Request*)A->haloRecvRequests+recvMessage);
      recvOffset += A->NrecvPairs[r];
      ++recvMessage;
    }
  }

  // Wait for all sent messages to have left and received messages to have arrived
  MPI_Status *sendStatus = (MPI_Status*) calloc(A->NsendMessages, sizeof(MPI_Status));
  MPI_Status *recvStatus = (MPI_Status*) calloc(A->NrecvMessages, sizeof(MPI_Status));
  
  MPI_Waitall(A->NrecvMessages, (MPI_Request*)A->haloRecvRequests, recvStatus);
  MPI_Waitall(A->NsendMessages, (MPI_Request*)A->haloSendRequests, sendStatus);
  
  free(recvStatus);
  free(sendStatus);

  //shift to local ids
  for (iint n=0;n<A->NsendTotal;n++)
    A->haloElementList[n] -= globalColStarts[rank];

  A->sendBuffer = (dfloat *) calloc(A->NsendTotal,sizeof(dfloat));    
}

void csrHaloExchange(csr *A, 
                    size_t Nbytes,         // message size per element
                    void *sourceBuffer,  
                    void *sendBuffer,    // temporary buffer
                    void *recvBuffer) {
  // MPI info
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // count outgoing and incoming meshes
  iint tag = 999;

  // copy data from outgoing elements into temporary send buffer
  for(iint i=0;i<A->NsendTotal;++i){
    // outgoing element
    iint id = A->haloElementList[i];

    memcpy(((char*)sendBuffer)+i*Nbytes, ((char*)sourceBuffer)+id*Nbytes, Nbytes);
  }

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
    MPI_Waitall(A->NsendMessages, (MPI_Request*)A->haloSendRequests, sendStatus);
    free(sendStatus);
  }
  if (A->NrecvTotal) {
    MPI_Status *recvStatus = (MPI_Status*) calloc(A->NrecvMessages, sizeof(MPI_Status));
    MPI_Waitall(A->NrecvMessages, (MPI_Request*)A->haloRecvRequests, recvStatus);
    free(recvStatus);
  }
}


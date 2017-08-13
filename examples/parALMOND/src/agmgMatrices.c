#include "agmg.h"

csr * newCSRfromCOO(iint N, iint* globalRowStarts,
            iint nnz, iint *Ai, iint *Aj, dfloat *Avals){

  iint size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  csr *A = (csr *) calloc(1,sizeof(csr));

  A->Nrows = N;
  A->Ncols = N;

  A->NlocalCols = N;

  iint globalOffset = globalRowStarts[rank];

  //first, count number of local, and non-local non-zeros
  iint diagNNZ=0;
  iint offdNNZ=0;
  for (iint n=0;n<nnz;n++) {
    if ((Aj[n] < globalOffset) || (Aj[n]>globalOffset+N-1)) offdNNZ++;
    else diagNNZ++;
  }

  iint   *diagAi, *diagAj;
  iint   *offdAi, *offdAj;
  dfloat *diagAvals, *offdAvals;

  if (diagNNZ) {
    diagAi        = (iint *)   calloc(diagNNZ, sizeof(iint));
    diagAj        = (iint *)   calloc(diagNNZ, sizeof(iint));
    diagAvals     = (dfloat *) calloc(diagNNZ, sizeof(dfloat));
  }
  if (offdNNZ) {
    offdAi        = (iint *)   calloc(offdNNZ, sizeof(iint));
    offdAj        = (iint *)   calloc(offdNNZ, sizeof(iint));
    offdAvals     = (dfloat *) calloc(offdNNZ, sizeof(dfloat));
  }

  //split into local and non-local COO matrices
  diagNNZ =0;
  offdNNZ =0;
  for (iint n=0;n<nnz;n++) {
    if ((Aj[n] < globalOffset) || (Aj[n]>globalOffset+N-1)) {
      offdAi[offdNNZ] = Ai[n] - globalOffset; //local index
      offdAj[offdNNZ] = Aj[n];                //global index
      offdAvals[offdNNZ] = Avals[n];
      offdNNZ++;
    } else {
      diagAi[diagNNZ] = Ai[n] - globalOffset; //local index
      diagAj[diagNNZ] = Aj[n] - globalOffset; //local index
      diagAvals[diagNNZ] = Avals[n];
      diagNNZ++;
    }
  }

  A->diagNNZ   = diagNNZ;
  A->offdNNZ   = offdNNZ;

  if (N) {
    A->diagRowStarts = (iint *)   calloc(N+1,sizeof(iint));
    A->offdRowStarts = (iint *)   calloc(N+1,sizeof(iint));
  }
  if (diagNNZ) {
    A->diagCols  = (iint *)   calloc(diagNNZ, sizeof(iint));
    A->diagCoefs = (dfloat *) calloc(diagNNZ, sizeof(dfloat));
  }
  if (offdNNZ) {
    A->offdCols  = (iint *)   calloc(offdNNZ,sizeof(iint));
    A->offdCoefs = (dfloat *) calloc(offdNNZ, sizeof(dfloat));
  }

  // Convert to csr storage, assumes orginal matrix was presorted by rows
  for(iint n=0;n<diagNNZ;++n) {
    iint row = diagAi[n];
    A->diagRowStarts[row+1]++;
  }
  for(iint n=0;n<offdNNZ;++n) {
    iint row = offdAi[n];
    A->offdRowStarts[row+1]++;
  }
  //cumulative sum
  for (iint i=0;i<A->Nrows;i++) {
    A->diagRowStarts[i+1] += A->diagRowStarts[i];
    A->offdRowStarts[i+1] += A->offdRowStarts[i];
  }

  //copy input data into struct
  if (diagNNZ) {
    for (iint i=0; i<N; i++) {
      iint start = A->diagRowStarts[i];
      iint cnt = 1;
      for (iint j=A->diagRowStarts[i]; j<A->diagRowStarts[i+1]; j++) {
        if (diagAj[j] == i) { //move diagonal to first entry
          A->diagCols[start]  = diagAj[j];
          A->diagCoefs[start] = diagAvals[j];
        } else {
          A->diagCols[start+cnt]  = diagAj[j];
          A->diagCoefs[start+cnt] = diagAvals[j];
          cnt++;
        }
      }
    }
  }

  //record global indexing of columns
  A->colMap = (iint *)   calloc(A->Ncols, sizeof(iint));
  for (iint i=0;i<A->Ncols;i++)
    A->colMap[i] = i + globalOffset;

  if (offdNNZ) {
    for (iint i=0; i<N; i++) {
      iint start = A->offdRowStarts[i];
      iint cnt = 0;
      for (iint j=A->offdRowStarts[i]; j<A->offdRowStarts[i+1]; j++) {
        A->offdCols[start+cnt]  = offdAj[j];
        A->offdCoefs[start+cnt] = offdAvals[j];
        cnt++;
      }
    }

    //we now need to reorder the x vector for the halo, and shift the column indices
    iint *col = (iint *) calloc(A->offdNNZ,sizeof(iint));
    for (iint n=0;n<offdNNZ;n++)
      col[n] = A->offdCols[n]; //copy non-local column global ids

    //sort by global index
    std::sort(col,col+offdNNZ);

    //count unique non-local column ids
    A->NHalo = 0;
    for (iint n=1;n<offdNNZ;n++)
      if (col[n]!=col[n-1])  col[++A->NHalo] = col[n];
    A->NHalo++; //number of unique columns

    A->Ncols += A->NHalo;

    //save global column ids in colMap
    A->colMap    = (iint *) realloc(A->colMap, A->Ncols*sizeof(iint));
    for (iint n=0; n<A->NHalo; n++)
      A->colMap[n+A->NlocalCols] = col[n];
    free(col);

    //shift the column indices to local indexing
    for (iint n=0;n<offdNNZ;n++) {
      iint gcol = A->offdCols[n];
      for (iint m=A->NlocalCols;m<A->Ncols;m++) {
        if (gcol == A->colMap[m])
          A->offdCols[n] = m;
      }
    }
  }

  if (diagNNZ) {
    free(diagAi);
    free(diagAj);
    free(diagAvals);
  }
  if (offdNNZ) {
    free(offdAi);
    free(offdAj);
    free(offdAvals);
  }

  csrHaloSetup(A,globalRowStarts);

  return A;
}

void freeCSR(csr *A) {
  if (A->diagNNZ) {
    free(A->diagRowStarts);
    free(A->diagCols);
    free(A->diagCoefs);
  }
  if (A->offdNNZ) {
    free(A->offdRowStarts);
    free(A->offdCols);
    free(A->offdCoefs);
  }
  if (A->Ncols) {
    free(A->colMap);
  }
  free(A->haloSendRequests);
  free(A->haloRecvRequests);
  free(A->NsendPairs);
  free(A->NrecvPairs);
  if (A->NsendTotal) {
    free(A->sendBuffer);
    free(A->haloElementList);
  }

  free(A);
}

//create a device version of a coo matrix
dcoo *newDCOO(parAlmond_t *parAlmond, csr *B){

  dcoo *A = (dcoo *) calloc(1,sizeof(dcoo));

  A->Nrows  = B->Nrows;
  A->Ncols  = B->Ncols;

  A->NHalo = B->NHalo;
  A->NlocalCols = B->NlocalCols;

  A->diagNNZ = B->diagNNZ;
  A->offdNNZ = B->offdNNZ;

  iint *diagRows;
  iint *offdRows;
  if (B->diagNNZ)
    diagRows = (iint *) calloc(B->diagNNZ,sizeof(iint));
  if (B->offdNNZ)
    offdRows = (iint *) calloc(B->offdNNZ,sizeof(iint));

  iint diagCnt =0;
  iint offdCnt =0;
  for (iint i=0;i<B->Nrows;i++) {
    for (iint j=B->diagRowStarts[i];j<B->diagRowStarts[i+1];j++)
      diagRows[diagCnt++] = i;

    for (iint j=B->offdRowStarts[i];j<B->offdRowStarts[i+1];j++)
      offdRows[offdCnt++] = i;
  }

  //copy to device
  if(B->diagNNZ){
    A->o_diagRows  = parAlmond->device.malloc(A->diagNNZ*sizeof(iint),   diagRows);
    A->o_diagCols  = parAlmond->device.malloc(A->diagNNZ*sizeof(iint),   B->diagCols);
    A->o_diagCoefs = parAlmond->device.malloc(A->diagNNZ*sizeof(dfloat), B->diagCoefs);
  }
  if(B->offdNNZ){
    A->o_offdRows  = parAlmond->device.malloc(A->offdNNZ*sizeof(iint), offdRows);
    A->o_offdCols  = parAlmond->device.malloc(A->offdNNZ*sizeof(iint),   B->offdCols);
    A->o_offdCoefs = parAlmond->device.malloc(A->offdNNZ*sizeof(dfloat), B->offdCoefs);
  }

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

  A->NlocalCols = csrA->NlocalCols;
  A->NHalo = csrA->NHalo;


  dfloat *diagInv;
  iint *rowCounters;
  if (csrA->Nrows) {
    diagInv = (dfloat *) calloc(csrA->Nrows, sizeof(dfloat));
    rowCounters = (iint*) calloc(csrA->Nrows, sizeof(iint));
  }

  for(iint i=0; i<csrA->Nrows; i++)
    diagInv[i] = 1.0/csrA->diagCoefs[csrA->diagRowStarts[i]];

  iint maxNnzPerRow = 0;
  iint minNnzPerRow = csrA->Ncols;
  for(iint i=0; i<csrA->Nrows; i++) {
    iint rowNnz = csrA->diagRowStarts[i+1] - csrA->diagRowStarts[i];
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

  dfloat threshold = 2.0/3.0;
  iint totalNNZ = csrA->diagNNZ+csrA->offdNNZ;
  iint nnzPerRow = 0;
  iint nnz = 0;

  //increase the nnz per row in E until it holds threshold*totalnnz nonzeros
  for(iint i=0; i<numBins; i++){
    nnz += bins[i] * (i+minNnzPerRow);
    if((nnz > threshold*totalNNZ)||(i==numBins-1)){
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
  for(iint i=0; i<csrA->Nrows; i++) {
    //excess from row in diag
    if(rowCounters[i] > nnzPerRow) nnzC += (rowCounters[i] - nnzPerRow);

    //all of offd
    iint offdRowNnz = csrA->offdRowStarts[i+1]-csrA->offdRowStarts[i];

    nnzC += offdRowNnz;
  }

  A->E->actualNNZ  = totalNNZ - nnzC;

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
    iint Jstart = csrA->diagRowStarts[i];
    iint Jend   = csrA->diagRowStarts[i+1];
    iint rowNnz = Jend - Jstart;

    // store only min of nnzPerRow and rowNnz
    iint maxNnz = (nnzPerRow >= rowNnz) ? rowNnz : nnzPerRow;

    for(iint c=0; c<maxNnz; c++){
      Ecols [i+c*A->E->strideLength]  = csrA->diagCols[Jstart+c];
      Ecoefs[i+c*A->E->strideLength]  = csrA->diagCoefs[Jstart+c];
    }

    // store the remaining in coo format
    if(rowNnz > nnzPerRow){
      for(iint c=nnzPerRow; c<rowNnz; c++){
        Coffsets[i+1]++;
        Ccols[nnzC]   = csrA->diagCols[Jstart+c];
        Ccoefs[nnzC]  = csrA->diagCoefs[Jstart+c];
        nnzC++;
      }
    }

    //add the offd non-zeros
    for (iint j=csrA->offdRowStarts[i];j<csrA->offdRowStarts[i+1];j++) {
      Coffsets[i+1]++;
      Ccols[nnzC]   = csrA->offdCols[j];
      Ccoefs[nnzC]  = csrA->offdCoefs[j];
      nnzC++;
    }
  }

  //use counts to create offsets
  for (iint i=0;i<csrA->Nrows;i++)
    Coffsets[i+1] += Coffsets[i];

  // copy the data to device memory
  if(csrA->Nrows) {
    A->o_diagInv = parAlmond->device.malloc(csrA->Nrows*sizeof(dfloat), diagInv);
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

  if (A->NsendTotal + A->NrecvTotal)
    csrHaloExchangeStart(A, sizeof(dfloat), x, A->sendBuffer, x+A->NlocalCols);

  // y[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  #pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){ //local
    dfloat result = 0.0;
    for(iint jj=A->diagRowStarts[i]; jj<A->diagRowStarts[i+1]; jj++)
      result += (A->diagCoefs[jj]*x[A->diagCols[jj]]);

    y[i] = alpha*result + beta*y[i];
  }

  if (A->NsendTotal + A->NrecvTotal)
    csrHaloExchangeFinish(A);

  #pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){ //nonlocal
    dfloat result = 0.0;
    for(iint jj=A->offdRowStarts[i]; jj<A->offdRowStarts[i+1]; jj++)
      result += (A->offdCoefs[jj]*x[A->offdCols[jj]]);

    y[i] += alpha*result;
  }
}

void zeqaxpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y, dfloat *z) {

  if (A->NsendTotal + A->NrecvTotal)
    csrHaloExchangeStart(A, sizeof(dfloat), x, A->sendBuffer, x+A->NlocalCols);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  #pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){
    dfloat result = 0.0;
    for(iint jj=A->diagRowStarts[i]; jj<A->diagRowStarts[i+1]; jj++)
      result += (A->diagCoefs[jj]*x[A->diagCols[jj]]);

    z[i] = alpha*result + beta*y[i];
  }

  if (A->NsendTotal + A->NrecvTotal)
    csrHaloExchangeFinish(A);

  #pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){
    dfloat result = 0.0;
    for(iint jj=A->offdRowStarts[i]; jj<A->offdRowStarts[i+1]; jj++)
      result += (A->offdCoefs[jj]*x[A->offdCols[jj]]);

    z[i] += alpha*result;
  }
}

void axpy(parAlmond_t *parAlmond, dcoo *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y) {

  occaTimerTic(parAlmond->device,"dcoo axpy");
  if (A->NsendTotal) {
    parAlmond->haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);

    //copy from device
    A->o_haloBuffer.copyTo(A->sendBuffer);
  }

  if (A->NsendTotal + A->NrecvTotal)
    dcooHaloExchangeStart(A, sizeof(dfloat), A->sendBuffer, A->recvBuffer);

  if (A->diagNNZ)
    parAlmond->agg_interpolateKernel(A->diagNNZ, A->o_diagRows, A->o_diagCols, A->o_diagCoefs, o_x, o_y);

  if (A->NsendTotal + A->NrecvTotal)
    dcooHaloExchangeFinish(A);

  //copy back to device
  if(A->NrecvTotal)
    o_x.copyFrom(A->recvBuffer,A->NrecvTotal*sizeof(dfloat),
                  A->NlocalCols*sizeof(dfloat));

  if (A->offdNNZ)
    parAlmond->agg_interpolateKernel(A->offdNNZ, A->o_offdRows, A->o_offdCols, A->o_offdCoefs, o_x, o_y);

  occaTimerToc(parAlmond->device,"dcoo axpy");
}

void axpy(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y) {

  occaTimerTic(parAlmond->device,"hyb axpy");
  if (A->NsendTotal) {
    parAlmond->haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);

    //copy from device
    A->o_haloBuffer.copyTo(A->sendBuffer);
  }

  if (A->NsendTotal+A->NrecvTotal)
    hybHaloExchangeStart(A, sizeof(dfloat),A->sendBuffer, A->recvBuffer);

  // y <-- alpha*E*x+beta*y
  axpy(parAlmond, A->E, alpha, o_x, beta, o_y);

  if (A->NsendTotal+A->NrecvTotal)
    hybHaloExchangeFinish(A);

  //copy back to device
  if (A->NrecvTotal)
    o_x.copyFrom(A->recvBuffer,A->NrecvTotal*sizeof(dfloat),A->NlocalCols*sizeof(dfloat));

  // y <-- alpha*C*x + y
  if (A->C->nnz)
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

  if (A->NsendTotal+A->NrecvTotal)
    hybHaloExchangeStart(A, sizeof(dfloat),A->sendBuffer, A->recvBuffer);

  // z <-- alpha*E*x+ beta*y
  zeqaxpy(parAlmond, A->E, alpha, o_x, beta, o_y, o_z);

  if (A->NsendTotal+A->NrecvTotal)
    hybHaloExchangeFinish(A);

  if (A->NrecvTotal)
    o_x.copyFrom(A->recvBuffer,A->NrecvTotal*sizeof(dfloat),A->NlocalCols*sizeof(dfloat));

  // z <-- alpha*C*x + z
  if (A->C->nnz)
    ax(parAlmond, A->C, alpha, o_x, o_z);

  occaTimerToc(parAlmond->device,"hyb zeqaxpy");
}

void axpy(parAlmond_t *parAlmond, ell *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y) {

  if(A->actualNNZ){
    occaTimerTic(parAlmond->device,"ell axpy");
    parAlmond->ellAXPYKernel(A->Nrows, A->nnzPerRow, A->strideLength,
                          alpha, beta, A->o_cols, A->o_coefs, o_x, o_y);
    occaTimerToc(parAlmond->device,"ell axpy");
  }
}

void zeqaxpy(parAlmond_t *parAlmond, ell *A, dfloat alpha, occa::memory o_x,
            dfloat beta, occa::memory o_y,  occa::memory o_z) {

  if(A->actualNNZ){
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

void smoothJacobi(agmgLevel *level, csr *A, dfloat *r, dfloat *x, bool x_is_zero) {

  // x = inv(D)*(b-R*x)  where R = A-D
  if(x_is_zero){
    #pragma omp parallel for
    for(iint i=0; i<A->Nrows; i++){
      dfloat invD = 1.0/A->diagCoefs[A->diagRowStarts[i]];
      x[i] = invD*r[i];
    }
    return;
  }

  if (A->NsendTotal + A->NrecvTotal)
    csrHaloExchangeStart(A, sizeof(dfloat), x, A->sendBuffer, x+A->NlocalCols);

  dfloat *y = level->smootherResidual;
  #pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){
    dfloat result = r[i];
    for(iint jj = A->diagRowStarts[i]; jj<A->diagRowStarts[i+1]; jj++)
      result -= A->diagCoefs[jj]*x[A->diagCols[jj]];

    y[i] = result;
  }

  if (A->NsendTotal + A->NrecvTotal)
    csrHaloExchangeFinish(A);

  #pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){
    dfloat result = y[i];
    const dfloat invD = 1./A->diagCoefs[A->diagRowStarts[i]];
    for(iint jj = A->offdRowStarts[i]; jj<A->offdRowStarts[i+1]; jj++)
      result -= A->offdCoefs[jj]*x[A->offdCols[jj]];

    y[i] = invD*result;
  }

  // copy the buffer vector to x
  #pragma omp parallel for
  for (iint i=0;i<A->Nrows;i++)
    x[i] = y[i];

}


void smoothDampedJacobi(agmgLevel *level, csr *A, dfloat *r, dfloat *x, bool x_is_zero) {

  dfloat alpha = level->smoother_params[0];

  if(x_is_zero){
  #pragma omp parallel for
    for(iint i=0; i<A->Nrows; i++){
      dfloat invD = 1.0/A->diagCoefs[A->diagRowStarts[i]];
      x[i] = alpha*invD*r[i];
    }
    return;
  }

  if (A->NsendTotal + A->NrecvTotal)
    csrHaloExchangeStart(A, sizeof(dfloat), x, A->sendBuffer, x+A->NlocalCols);

  // x = (1-alpha)*x + alpha*inv(D) * (b-R*x) where R = A-D
  dfloat *y = level->smootherResidual;

  const dfloat oneMalpha = 1. - alpha;
  #pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){
    dfloat result = r[i];

    for(iint jj =A->diagRowStarts[i]; jj<A->diagRowStarts[i+1]; jj++)
      result -= A->diagCoefs[jj]*x[A->diagCols[jj]];

    y[i] = result;
  }

  if (A->NsendTotal + A->NrecvTotal)
    csrHaloExchangeFinish(A);

  #pragma omp parallel for
  for(iint i=0; i<A->Nrows; i++){
    dfloat result = y[i];
    const dfloat invD = 1./A->diagCoefs[A->diagRowStarts[i]];
    for(iint jj =A->offdRowStarts[i]; jj<A->offdRowStarts[i+1]; jj++)
      result -= A->offdCoefs[jj]*x[A->offdCols[jj]];

    y[i] = oneMalpha*x[i] + alpha*invD*result;
  }

  // copy the buffer vector to x
  #pragma omp parallel for
  for (iint i=0;i<A->Nrows;i++)
    x[i] = y[i];
}

void smoothChebyshev(agmgLevel *level, csr *A, dfloat *r, dfloat *x, bool x_is_zero) {

  dfloat lambdaN = level->smoother_params[0];
  dfloat lambda1 = level->smoother_params[1];

  dfloat theta = 0.5*(lambdaN+lambda1);
  dfloat delta = 0.5*(lambdaN-lambda1);
  dfloat invTheta = 1.0/theta;
  dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  dfloat *res = level->smootherResidual;
  dfloat *Ad  = level->smootherResidual2;
  dfloat *d   = level->smootherUpdate;

  if(x_is_zero){ //skip the Ax if x is zero
    #pragma omp parallel for
    for(iint i=0; i<A->Nrows; i++){
      dfloat invD = 1.0/A->diagCoefs[A->diagRowStarts[i]];
      res[i] = invD*r[i];
      x[i] = 0.;
      d[i] = invTheta*res[i];
    }
  } else {

    level->Ax(level->AxArgs,x,res);

    #pragma omp parallel for
    for(iint i=0; i<A->Nrows; i++){
      dfloat invD = 1.0/A->diagCoefs[A->diagRowStarts[i]];
      res[i] = invD*(r[i]-res[i]);
      d[i]   = invTheta*res[i];
    }
  }

  for (int k=0;k<level->ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    vectorAdd(A->Nrows, 1.0, d, 1.0, x);

    //r_k+1 = r_k - D^{-1}Ad_k
    level->Ax(level->AxArgs,d,Ad);
    #pragma omp parallel for
    for(iint i=0; i<A->Nrows; i++) {
      dfloat invD = 1.0/A->diagCoefs[A->diagRowStarts[i]];
      res[i] = res[i] - invD*Ad[i];
    }

    rho_np1 = 1.0/(2.*sigma-rho_n);

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    vectorAdd(A->Nrows, 2.0*rho_np1/delta, res, rho_np1*rho_n, d);
    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  vectorAdd(A->Nrows, 1.0, d, 1.0, x);
}

void smoothJacobi(parAlmond_t *parAlmond, agmgLevel *level, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero) {

  occaTimerTic(parAlmond->device,"hyb smoothJacobi");
  if(x_is_zero){
    if (A->Nrows)
      dotStar(parAlmond, A->Nrows, 1.0, A->o_diagInv, o_r, 0.0, o_x);
    occaTimerToc(parAlmond->device,"hyb smoothJacobi");
    return;
  }

  occa::memory o_res = level->o_smootherResidual;

  if (A->NsendTotal) {
    parAlmond->haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);

    //copy from device
    A->o_haloBuffer.copyTo(A->sendBuffer);
  }

  if (A->NsendTotal+A->NrecvTotal)
    hybHaloExchangeStart(A, sizeof(dfloat),A->sendBuffer, A->recvBuffer);

  occaTimerTic(parAlmond->device,"ellJacobi1");
  if (A->Nrows)
    parAlmond->ellJacobiKernel(A->Nrows, A->E->nnzPerRow, A->E->strideLength,
                            A->E->o_cols, A->E->o_coefs, o_x, o_r, o_res);
  occaTimerToc(parAlmond->device,"ellJacobi1");

  if (A->NsendTotal+A->NrecvTotal)
    hybHaloExchangeFinish(A);

  //copy back to device
  if (A->NrecvTotal)
    o_x.copyFrom(A->recvBuffer,A->NrecvTotal*sizeof(dfloat),A->NlocalCols*sizeof(dfloat));

  // temp1 += -C*x
  ax(parAlmond, A->C, -1.0, o_x, o_res);

  // x = invD*temp1
  if (A->Nrows)
    dotStar(parAlmond, A->Nrows, 1.0, A->o_diagInv, o_res, 0., o_x);
  occaTimerToc(parAlmond->device,"hyb smoothJacobi");
}

void smoothDampedJacobi(parAlmond_t *parAlmond, agmgLevel *level, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero){

  dfloat alpha = level->smoother_params[0];

  occaTimerTic(parAlmond->device,"hyb smoothDampedJacobi");
  if(x_is_zero){
    if (A->Nrows)
      dotStar(parAlmond, A->Nrows, alpha, A->o_diagInv, o_r, 0.0, o_x);
    occaTimerToc(parAlmond->device,"hyb smoothDampedJacobi");
    return;
  }

  occa::memory o_res = level->o_smootherResidual;

  if (A->NsendTotal) {
    parAlmond->haloExtract(A->NsendTotal, 1, A->o_haloElementList, o_x, A->o_haloBuffer);

    //copy from device
    A->o_haloBuffer.copyTo(A->sendBuffer);
  }

  if (A->NsendTotal+A->NrecvTotal)
    hybHaloExchangeStart(A, sizeof(dfloat),A->sendBuffer, A->recvBuffer);

  occaTimerTic(parAlmond->device,"ellJacobi1");
  if (A->Nrows)
    parAlmond->ellJacobiKernel(A->Nrows, A->E->nnzPerRow, A->E->strideLength,
                              A->E->o_cols, A->E->o_coefs, o_x, o_r, o_res);
  occaTimerToc(parAlmond->device,"ellJacobi1");

  if (A->NsendTotal+A->NrecvTotal)
    hybHaloExchangeFinish(A);

  //copy back to device
  if (A->NrecvTotal)
    o_x.copyFrom(A->recvBuffer,A->NrecvTotal*sizeof(dfloat),A->NlocalCols*sizeof(dfloat));


  // temp1 += -C*x
  ax(parAlmond, A->C, -1.0, o_x, o_res);

  // x = alpha*invD*temp1 + (1-alpha)*x
  const dfloat beta = 1.0 - alpha;
  dotStar(parAlmond, A->Nrows, alpha, A->o_diagInv, o_res, beta, o_x);

  occaTimerToc(parAlmond->device,"hyb smoothDampedJacobi");
}

void smoothChebyshev(parAlmond_t *parAlmond, agmgLevel *level, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero) {

  dfloat lambdaN = level->smoother_params[0];
  dfloat lambda1 = level->smoother_params[1];

  dfloat theta = 0.5*(lambdaN+lambda1);
  dfloat delta = 0.5*(lambdaN-lambda1);
  dfloat invTheta = 1.0/theta;
  dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  occa::memory o_res = level->o_smootherResidual;
  occa::memory o_Ad  = level->o_smootherResidual2;
  occa::memory o_d   = level->o_smootherUpdate;

  occaTimerToc(parAlmond->device,"hyb smoothChebyshev");

  if(x_is_zero){ //skip the Ax if x is zero
    //res = D^{-1}r
    dotStar(parAlmond, A->Nrows, 1.0, A->o_diagInv, o_r, 0.0, o_res);
    setVector(parAlmond, A->Nrows, o_x, 0.0);
    //d = invTheta*res
    vectorAdd(parAlmond, A->Nrows, invTheta, o_res, 0.0, o_d);

  } else {

    //res = D^{-1}(r-Ax)
    level->device_Ax(level->AxArgs,o_x,o_res);
    vectorAdd(parAlmond, A->Nrows, 1.0, o_r, -1.0, o_res);
    dotStar(parAlmond, A->Nrows, A->o_diagInv, o_res);

    //d = invTheta*res
    vectorAdd(parAlmond, A->Nrows, invTheta, o_res, 0.0, o_d);
  }

  for (int k=0;k<level->ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    vectorAdd(parAlmond, A->Nrows, 1.0, o_d, 1.0, o_x);

    //r_k+1 = r_k - D^{-1}Ad_k
    level->device_Ax(level->AxArgs,o_d,o_Ad);
    dotStar(parAlmond, A->Nrows, -1.0, A->o_diagInv, o_Ad, 1.0, o_res);

    rho_np1 = 1.0/(2.*sigma-rho_n);

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    vectorAdd(parAlmond, A->Nrows, 2.0*rho_np1/delta, o_res, rho_np1*rho_n, o_d);
    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  vectorAdd(parAlmond, A->Nrows, 1.0, o_d, 1.0, o_x);

  occaTimerToc(parAlmond->device,"hyb smoothChebyshev");
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
  A->NsendPairs = (iint*) calloc(size, sizeof(iint));
  A->NrecvPairs = (iint*) calloc(size, sizeof(iint));
  for(iint n=A->NlocalCols;n<A->Ncols;++n){ //for just the halo
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

  if (A->NsendTotal)
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
  iint recvOffset = A->NlocalCols;
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

  if (A->NsendTotal)
    A->sendBuffer = (dfloat *) calloc(A->NsendTotal,sizeof(dfloat));

  A->totalHaloPairs = A->NsendTotal+A->NrecvTotal;
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

void csrHaloExchangeStart(csr *A,
                    size_t Nbytes,         // message size per element
                    void *sourceBuffer,
                    void *sendBuffer,    // temporary buffer
                    void *recvBuffer) {
  // MPI info
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
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
}

void csrHaloExchangeFinish(csr *A) {
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

void dcooHaloExchangeStart(dcoo *A, size_t Nbytes, void *sendBuffer, void *recvBuffer) {
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

void dcooHaloExchangeFinish(dcoo *A) {
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


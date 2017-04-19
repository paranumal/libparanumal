#include "parAlmond.h"


csr * newCSR(iint Nrows, iint Ncolumns, iint nnz,
      iint *rowStarts, iint *cols, dfloat *vals){


  csr *A = (csr *) calloc(1,sizeof(csr));

  A->Nrows = Nrows;
  A->Ncolumns = Ncolumns;
  A->nnz = nnz;

  A->rowStarts = (iint *) calloc(A->Nrows+1,sizeof(iint));
  A->cols      = (iint *) calloc(A->nnz, sizeof(iint));
  A->coefs     = (iint *) calloc(A->nnz, sizeof(dfloat));


  for (iint i=0; i<Nrows; i++) {
    start = rowStarts[i];
    A->rowStarts[i] = start; 
    iint cnt = 1;
    for (iint j=rowStarts[i]; j<rowStarts[i+1]; j++) {
      if (cols[j] == i) { //diagonal entry
        A->cols[start] = cols[j];
        A->vals[start] = vals[j];
      } else {
        A->cols[start+cnt] = cols[j];
        A->vals[start+cnt] = vals[j];
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


void axpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y) {

  HaloExchange(sizeof(dfloat), x, A->sendBuffer, x+A->numLocalIds);

  // y[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  for(iint i=0; i<A->Nrows; i++){
    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    dfloat result = 0.0;
    for(iint jj=Jstart; jj<Jend; jj++)
      result += (A->coefs[jj]*x[A->cols[jj]]);
    
    y[i] = alpha*result + beta*y[i];
  }
}

void zeqaxpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y, dfloat *z) {

  HaloExchange(sizeof(dfloat), x, A->sendBuffer, x+A->numLocalIds);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  for(iint i=0; i<A->Nrows; i++){
    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    dfloat result = (T) 0.0;

    for(iint jj=Jstart; jj<Jend; jj++)
      result += A->coefs[jj]*x[A->cols[jj]];

    z[i] = alpha*result + beta*y[i];
  }
}


void smoothJacobi(csr *A, dfloat *r, dfloat *x, const bool x_is_zero) {

  // x = inv(D)*(b-R*x)  where R = A-D
  if(x_is_zero){
    for(iint i=0; i<A->Nrows; i++){
      dfloat invD = 1.0/A->coefs[A->rowStarts[i]];
      x[i] = invD*r[i];
    }
    return;
  }

  HaloExchange(sizeof(dfloat), x, A->sendBuffer, x+A->numLocalIds);

  dfloat y[Nrows];

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
  for (iint i=0;i<A->Nrows;i++)
    x[i] = y[i];
}


void smoothDampedJacobi(csr *A, dfloat *r, dfloat *x, dfloat alpha, const bool x_is_zero) {
  
  if(x_is_zero){
    for(iint i=0; i<A->Nrows; i++){
      dfloat invD = 1.0/A->coefs[A->rowStarts[i]];
      x[i] = alpha*invD*r[i];
    }
    return;
  }

  HaloExchange(sizeof(dfloat), x, A->sendBuffer, x+A->numLocalIds);

  // x = (1-alpha)*x + alpha*inv(D) * (b-R*x) where R = A-D
  dfloat y[Nrows];

  const dfloat oneMalpha = 1. - alpha;

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
  for (iint i=0;i<A->Nrows;i++)
    x[i] = y[i];
}


csr * transpose(csr *A){

  csr *A = (csr *) calloc(1,sizeof(csr));

  At->Nrows    = A->Ncolumns;
  At->Ncolumns = A->Nrows;
  At->nnz      = A->NNZ;

  At->rowStarts = (iint *) calloc(At->Nrows+1, sizeof(iint));
  At->cols      = (iint *) calloc(At->nnz, sizeof(iint));
  At->coefs     = (dfloat *) calloc(At->nnz, sizeof(dfloat));

  // count the no of nonzeros per row for transpose
  for(iint i=0; i<A->nnz; i++){
    iint row = A->cols[i];
    At->rowStarts[row+1]++;
  }

  // cumulative sum for rows
  for(iint i=1; i<=At->Nrows; i++)
    At->rowStarts[i] += At->rowStarts[i-1];

  iint counter[At->Nrows+1];
  for (iint i=0; i<At->Nrows+1; i++)
    counter[i] = At->rowStarts[i];

  if(A->Nrows != A->Ncolumns){
    for(iint i=0; i<A->Nrows; i++){
      const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];
      
      for(iint jj=Jstart; jj<Jend; jj++){
        iint row = A->cols[jj];
        At.cols[counter[row]] = i;
        At.coefs[counter[row]] = A->coefs[jj];

        counter[row]++;
      }
    }
  }
  else{
    // fill in diagonals first
    for(iint i=0; i<A->Nrows; i++){
      At->cols[counter[i]] = i;

      At->coefs[counter[i]] = A->coefs[A->rowStarts[i]];
      counter[i]++;
    }

    // fill the remaining ones
    for(iint i=0; i<Nrows; i++){
      const iint Jstart = A->rowStarts[i]+1, Jend = A->rowStarts[i+1];

      for(iint jj=Jstart; jj<Jend; jj++){
        iint row = A->cols[jj];

        At->cols[counter[row]] = i;
        At->coefs[counter[row]] = A->coefs[jj];

        counter[row]++;
      }
    }
  }

  return At;
}

struct key_value_pair{
  iint key;
  dfloat value;
};


struct compare_key{
  bool operator()(key_value_pair a, key_value_pair b){
    return a.key < b.key;
  }
};


// TODO : this has to be implemented in a better way
csr * spmm(csr *A, csr *B){

  csr *C = (csr*) calloc(1,sizeof(csr));

  C->Nrows = A->Nrows;
  C->Ncols = B->Ncols;

  C->rowStarts = (iint *) calloc(C->Nrows,sizeof(iint));

  iint   **C_cols =    (iint **) calloc(C->Nrows,sizeof(iint*));
  dfloat **C_coefs = (dfloat **) calloc(C->Nrows,sizeof(dfloat*));

  C->nnz = 0;
  iint Rcounter[C->Nrows];

  struct key_value_pair *col_val_pair;

  for(iint i=0; i<C->Nrows; i++){
    // count all possible nonzeros for each row
    iint max_rowsize = 0;
    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    for(iint jj=Jstart; jj<Jend; jj++){

      const iint k = A->cols[jj];

      max_rowsize += (B->rowStarts[k+1] - B->rowStarts[k]);
    }

    // do the product
    col_val_pair = (struct key_value_pair *) calloc(max_rowsize,sizeof(struct key_value_pair));

    iint counter = 0;
    for(iint jj=Jstart; jj<Jend; jj++){

      const iint k = A->cols[jj];
      const dfloat Aik = A->coefs[jj];

      const iint JstartB = B->rowStarts[k], JendB = B->rowStarts[k+1];
      for (iint jjB=JstartB; jjB<JendB; jjB++){

        // if it corresponds to diagonal keep it at the beginning of row
        // setting the key to be -ve ensures that it appears first after sorting
        iint j = B->cols[jjB];
        if(j == i) j = -1;

        // j
        col_val_pair[counter].key = j;

        // Aik*Bkj
        col_val_pair[counter].value = Aik*B->coefs[jjB];

        counter++;
      }
    }

    std::sort(col_val_pair, col_val_pair+max_rowsize, compare_key());

    // find unique entries
    iint rowsize = 1;

    for(iint r=1; r<max_rowsize; r++){
      if(col_val_pair[r].key != col_val_pair[r-1].key)
        rowsize++;
    }


    C->rowStarts[i+1] = rowsize;
    C->nnz += rowsize;

    C_cols[i] = (iint *) calloc(rowsize, sizeof(iint));
    C_coefs[i] = (dfloat *) calloc(rowsize, sizeof(dfloat));

    Rcounter[i] = 0;
    iint col = col_val_pair[0].key;
    dfloat coef = col_val_pair[0].value;
    if(col == -1) col = i;

    C_cols[i][0]  = col;
    C_coefs[i][0] = coef;


    for(iint r=1; r<max_rowsize; r++){
      if(col_val_pair[r].key == col_val_pair[r-1].key){
        coef = col_val_pair[r].value;
        // if it already exists, add the value to the previous one
        C_coefs[i][Rcounter[i]] += coef;
      } else {
        col = col_val_pair[r].key;
        if(col == -1) col = i;

        coef = col_val_pair[r].value;

        Rcounter[i]++;
        C_cols[i][Rcounter[i]] = col;
        C_coefs[i][Rcounter[i]] = coef;
      }
    }
  }

  // copy it to global
  C->cols  = (iint *) calloc(C->nnz,sizeof(iint));
  C->coefs = (dfloat *) calloc(C->nnz,sizeof(dfloat));

  for(iint i=0; i<C->Nrows; i++){
    for (iint j=0;j<Rcounter[i];j++) {
      C->cols[C->rowStarts[i]+j] = C_cols[i][j];
      C->coefs[C->rowStarts[i]+j] = C_coefs[i][j];
    }

    // cumulative sum
    C->rowStarts[i+1] += C->rowStarts[i];
  }

  for(iint i=0; i<C->Nrows; i++){
    free(C_cols[i]); 
    free(C_coefs[i]);
  }

  return C;
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
  for(iint n=A->numLocalIds;n<A->Ncolumns;++n){ //for just the halo
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
  MPI_Status *sendStatus = (MPI_Status*) calloc(NsendMessages, sizeof(MPI_Status));
  MPI_Status *recvStatus = (MPI_Status*) calloc(NrecvMessages, sizeof(MPI_Status));
  
  MPI_Waitall(NrecvMessages, (MPI_Request*)haloRecvRequests, recvStatus);
  MPI_Waitall(NsendMessages, (MPI_Request*)haloSendRequests, sendStatus);
  
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
    if (NsendTotal) {
      if(NsendPairs[r]){
        MPI_Isend(((char*)sendBuffer)+sendOffset, NsendPairs[r]*Nbytes, MPI_CHAR, r, tag,
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


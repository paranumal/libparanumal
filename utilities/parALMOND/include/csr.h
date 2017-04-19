

csr * newCSR(iint Nrows, iint Ncolumns, iint nnz,
      iint *rowStarts, iint *cols, dfloat *vals);

void axpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

void zeqaxpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y, dfloat *z);

void smoothJacobi(csr *A, dfloat *r, dfloat *x, bool x_is_zero);

void smoothDampedJacobi(csr *A, dfloat *r, dfloat *x, dfloat alpha, bool x_is_zero);

csr * transpose(csr *A);

csr *galerkinProd(csr *A, csr *P, iint ag);

csr * spmm(csr *A, csr *B);

void csrHaloSetup(csr *A, iint *globalColStarts);

void csrHaloExchange(csr *A, 
                    size_t Nbytes,         // message size per element
                    void *sourceBuffer,  
                    void *sendBuffer,    // temporary buffer
                    void *recvBuffer);
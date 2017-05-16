
//creators
csr * newCSR(iint Nrows, iint Ncolumns, iint nnz,
      iint *rowStarts, iint *cols, dfloat *vals);
void freeCSR(csr *A);
dcsr *newDCSR(parAlmond_t *parAlmond, csr *B);
hyb * newHYB(parAlmond_t *parAlmond, csr *csrA);


void axpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

void zeqaxpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y, dfloat *z);

void axpy(parAlmond_t *parAlmond, dcsr *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void axpy(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void zeqaxpy(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x,
      dfloat beta,  occa::memory o_y, occa::memory o_z);

void axpy(parAlmond_t *parAlmond, ell *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void zeqaxpy(parAlmond_t *parAlmond, ell *A, dfloat alpha, occa::memory o_x, 
            dfloat beta, occa::memory o_y,  occa::memory o_z);

void zeqaxpy(parAlmond_t *parAlmond, dcsr *A, dfloat alpha, occa::memory o_x, dfloat beta, 
              occa::memory o_y, occa::memory o_z);

void ax(parAlmond_t *parAlmond, coo *C, dfloat alpha, occa::memory o_x, occa::memory o_y);

void matFreeZeqAXPY(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, 
              occa::memory o_y, occa::memory o_z);

//smoothing
void smoothJacobi(csr *A, dfloat *r, dfloat *x, bool x_is_zero);
void smoothDampedJacobi(csr *A, dfloat *r, dfloat *x, dfloat alpha, bool x_is_zero);
void smoothJacobi(parAlmond_t *parAlmond, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);
void smoothDampedJacobi(parAlmond_t *parAlmond, hyb *A, occa::memory o_r, occa::memory o_x,
             				dfloat alpha, bool x_is_zero);
void smoothJacobi(parAlmond_t *parAlmond, dcsr *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);
void smoothDampedJacobi(parAlmond_t *parAlmond, dcsr *A, occa::memory o_r, occa::memory o_x,
         					dfloat alpha, bool x_is_zero);

void matFreeSmoothJacobi(parAlmond_t *parAlmond, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);
void matFreeSmoothDampedJacobi(parAlmond_t *parAlmond, hyb* A, occa::memory o_r,
                            occa::memory o_x, dfloat alpha, bool x_is_zero);

//halo exchange
void csrHaloSetup(csr *A, iint *globalColStarts);
void csrHaloExchange(csr *A, size_t Nbytes, void *sourceBuffer, void *sendBuffer, void *recvBuffer);
void dcsrHaloExchangeStart(dcsr *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);
void dcsrHaloExchangeFinish(dcsr *A);
void hybHaloExchangeStart(hyb *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);
void hybHaloExchangeFinish(hyb *A);
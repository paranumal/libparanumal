
//creators
csr * newCSR(iint Nrows, iint Ncolumns, iint nnz,
      iint *rowStarts, iint *cols, dfloat *vals);
dcsr *newDCSR(almond_t *almond, csr *B);
hyb * newHYB(almond_t *almond, csr *csrA);


void axpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

void zeqaxpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y, dfloat *z);

void axpy(almond_t *almond, dcsr *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void axpy(almond_t *almond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void zeqaxpy(almond_t *almond, hyb *A, dfloat alpha, occa::memory o_x,
      dfloat beta,  occa::memory o_y, occa::memory o_z);

void axpy(almond_t *almond, ell *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void zeqaxpy(almond_t *almond, ell *A, dfloat alpha, occa::memory o_x, 
            dfloat beta, occa::memory o_y,  occa::memory o_z);

void ax(almond_t *almond, coo *C, dfloat alpha, occa::memory o_x, occa::memory o_y);

//smoothing
void smoothJacobi(csr *A, dfloat *r, dfloat *x, bool x_is_zero);
void smoothDampedJacobi(csr *A, dfloat *r, dfloat *x, dfloat alpha, bool x_is_zero);
void smoothJacobi(almond_t *almond, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);
void smoothDampedJacobi(almond_t *almond, hyb *A, occa::memory o_r, occa::memory o_x,
             				dfloat alpha, bool x_is_zero);

//halo exchange
void csrHaloSetup(csr *A, iint *globalColStarts);
void csrHaloExchange(csr *A, size_t Nbytes, void *sourceBuffer, void *sendBuffer, void *recvBuffer);
void dcsrHaloExchangeStart(dcsr *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);
void dcsrHaloExchangeFinish(dcsr *A);
void hybHaloExchangeStart(hyb *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);
void hybHaloExchangeFinish(hyb *A);
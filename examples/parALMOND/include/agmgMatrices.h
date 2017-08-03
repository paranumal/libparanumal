
//creators
csr * newCSRfromCOO(iint N, iint* globalRowStarts,
            iint NNZ,   iint *Ai, iint *Aj, dfloat *Avals);
void freeCSR(csr *A);
dcoo *newDCOO(parAlmond_t *parAlmond, csr *B);
hyb * newHYB(parAlmond_t *parAlmond, csr *csrA);


void axpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

void zeqaxpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y, dfloat *z);

void axpy(parAlmond_t *parAlmond, dcoo *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void axpy(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void zeqaxpy(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x,
      dfloat beta,  occa::memory o_y, occa::memory o_z);

void axpy(parAlmond_t *parAlmond, ell *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void zeqaxpy(parAlmond_t *parAlmond, ell *A, dfloat alpha, occa::memory o_x,
            dfloat beta, occa::memory o_y,  occa::memory o_z);

void ax(parAlmond_t *parAlmond, coo *C, dfloat alpha, occa::memory o_x, occa::memory o_y);

void matFreeZeqAXPY(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta,
              occa::memory o_y, occa::memory o_z);

//smoothing
void smoothJacobi(csr *A, dfloat *r, dfloat *x, bool x_is_zero);
void smoothDampedJacobi(csr *A, dfloat *r, dfloat *x, dfloat alpha, bool x_is_zero);
void smoothJacobi(parAlmond_t *parAlmond, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);
void smoothDampedJacobi(parAlmond_t *parAlmond, hyb *A, occa::memory o_r, occa::memory o_x,
             				dfloat alpha, bool x_is_zero);

void matFreeSmoothJacobi(parAlmond_t *parAlmond, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);
void matFreeSmoothDampedJacobi(parAlmond_t *parAlmond, hyb* A, occa::memory o_r,
                            occa::memory o_x, dfloat alpha, bool x_is_zero);

//halo exchange
void csrHaloSetup(csr *A, iint *globalColStarts);
void csrHaloExchange(csr *A, size_t Nbytes, void *sourceBuffer, void *sendBuffer, void *recvBuffer);
void csrHaloExchangeStart(csr *A, size_t Nbytes, void *sourceBuffer, void *sendBuffer, void *recvBuffer);
void csrHaloExchangeFinish(csr *A);
void dcooHaloExchangeStart(dcoo *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);
void dcooHaloExchangeFinish(dcoo *A);
void hybHaloExchangeStart(hyb *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);
void hybHaloExchangeFinish(hyb *A);
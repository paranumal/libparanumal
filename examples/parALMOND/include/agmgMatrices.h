
//creators
csr * newCSRfromCOO(dlong N, hlong* globalRowStarts,
            		dlong NNZ,   hlong *Ai, hlong *Aj, dfloat *Avals);
void freeCSR(csr *A);
dcoo *newDCOO(parAlmond_t *parAlmond, csr *B);
hyb * newHYB(parAlmond_t *parAlmond, csr *csrA);


void axpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y, bool nullSpace, dfloat nullSpacePenalty);

void axpy(parAlmond_t *parAlmond, dcoo *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void axpy(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y, bool nullSpace, dfloat nullSpacePenalty);

void axpy(parAlmond_t *parAlmond, ell *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void ax(parAlmond_t *parAlmond, coo *C, dfloat alpha, occa::memory o_x, occa::memory o_y);


//smoothing
void smoothJacobi      (parAlmond_t *parAlmond, agmgLevel *level, csr *A, dfloat *r, dfloat *x, bool x_is_zero);
void smoothDampedJacobi(parAlmond_t *parAlmond, agmgLevel *level, csr *A, dfloat *r, dfloat *x, bool x_is_zero);
void smoothChebyshev   (parAlmond_t *parAlmond, agmgLevel *level, csr *A, dfloat *r, dfloat *x, bool x_is_zero);
void smoothJacobi      (parAlmond_t *parAlmond, agmgLevel *level, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);
void smoothDampedJacobi(parAlmond_t *parAlmond, agmgLevel *level, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);
void smoothChebyshev   (parAlmond_t *parAlmond, agmgLevel *level, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);

//halo exchange
void csrHaloSetup(csr *A, hlong *globalColStarts);
void csrHaloExchange(csr *A, size_t Nbytes, void *sourceBuffer, void *sendBuffer, void *recvBuffer);
void csrHaloExchangeStart(csr *A, size_t Nbytes, void *sourceBuffer, void *sendBuffer, void *recvBuffer);
void csrHaloExchangeFinish(csr *A);
void dcooHaloExchangeStart(dcoo *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);
void dcooHaloExchangeFinish(dcoo *A);
void hybHaloExchangeStart(hyb *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);
void hybHaloExchangeFinish(hyb *A);
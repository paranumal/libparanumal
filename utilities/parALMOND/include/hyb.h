

hyb * newHYB(almond_t *almond, csr *csrA);

void axpy(almond_t *almond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void zeqaxpy(almond_t *almond, hyb *A, dfloat alpha, occa::memory o_x,
      dfloat beta,  occa::memory o_y, occa::memory o_z);

void smoothJacobi(almond_t *almond, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);

void smoothDampedJacobi(almond_t *almond, hyb *A, occa::memory o_r, occa::memory o_x,
             dfloat alpha, bool x_is_zero);

void hybHaloExchangeStart(hyb *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);

void hybHaloExchangeFinish(hyb *A);
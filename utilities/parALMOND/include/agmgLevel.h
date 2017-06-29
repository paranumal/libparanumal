
void restrict(agmgLevel *level, dfloat *r, dfloat *Rr);

void restrict(parAlmond_t *parAlmond, agmgLevel *level, occa::memory o_r, occa::memory o_Rr);

void interpolate(agmgLevel *level, dfloat *x, dfloat *Px);

void interpolate(parAlmond_t *parAlmond, agmgLevel *level, occa::memory o_x, occa::memory o_Px);

void setupSmoother(agmgLevel *level, SmoothType s);

dfloat rhoDinvA(csr *A, dfloat *invD);

void smooth(agmgLevel *level, dfloat *rhs, dfloat *x, bool x_is_zero);

void smooth(parAlmond_t *parAlmond, agmgLevel *level, occa::memory o_rhs, occa::memory o_x, bool x_is_zero);

void matFreeSmooth(parAlmond_t *parAlmond, agmgLevel *level, occa::memory &o_r, occa::memory &o_x, bool x_is_zero);

void setupExactSolve(parAlmond_t *parAlmond, agmgLevel *level);
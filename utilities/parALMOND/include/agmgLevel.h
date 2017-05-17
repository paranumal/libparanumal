
void restrict(agmgLevel *level, dfloat *r, dfloat *Rr);

void restrict(parAlmond_t *parAlmond, agmgLevel *level, occa::memory o_r, occa::memory o_Rr);

void interpolate(agmgLevel *level, dfloat *x, dfloat *Px);

void interpolate(parAlmond_t *parAlmond, agmgLevel *level, occa::memory o_x, occa::memory o_Px);

void allocate(agmgLevel *level);

csr* distribute(csr *A, iint *globalRowStarts, iint *globalColStarts);

void setup_smoother(agmgLevel *level, SmoothType s);

dfloat rhoDinvA(csr *A, dfloat *invD);

void smooth(agmgLevel *level, dfloat *rhs, dfloat *x, bool x_is_zero);

void smooth(parAlmond_t *parAlmond, agmgLevel *level, occa::memory o_rhs, occa::memory o_x, bool x_is_zero);

void matFreeSmooth(parAlmond_t *parAlmond_t, agmgLevel *level, occa::memory &o_r, occa::memory &o_x, bool x_is_zero);

csr *strong_graph(csr *A, dfloat threshold);

bool customLess(iint smax, dfloat rmax, iint imax, iint s, dfloat r, iint i);

iint *form_aggregates(agmgLevel *level, csr *C);

void construct_interpolator(agmgLevel *level, iint *FineToCoarse, dfloat **nullCoarseA);

csr *galerkinProd(agmgLevel *level);

csr *transpose(csr *A);

void coarsen(agmgLevel *level, csr **coarseA, dfloat **nullCoarseA);



dfloat norm(iint n, dfloat *a);

dfloat innerProd(iint n, dfloat *a, dfloat *b);

void vectorAdd(iint n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

void dotStar(iint m, dfloat *a, dfloat *b);

void scaleVector(iint m, dfloat *a, dfloat alpha);

void randomize(iint m, dfloat *a);

dfloat maxEntry(iint n, dfloat *a);

void copyVector(parAlmond_t *parAlmond, iint N, occa::memory o_a, occa::memory o_b);

void scaleVector(parAlmond_t *parAlmond, iint N, occa::memory o_a, dfloat alpha);

void dotStar(parAlmond_t *parAlmond, iint N, occa::memory o_a, occa::memory o_b);

void dotStar(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_a,
	           occa::memory o_b, dfloat beta, occa::memory o_c);

dfloat innerProd(parAlmond_t *parAlmond, iint N, occa::memory o_a, occa::memory o_b);

void vectorAdd(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void vectorAdd(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_x,
	 dfloat beta, occa::memory o_y, occa::memory o_z);

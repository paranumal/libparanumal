

dfloat norm(iint n, dfloat *a);

dfloat innerProd(iint n, dfloat *a, dfloat *b);

void vectorAdd(iint n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

void dotStar(iint m, dfloat *a, dfloat *b);

void scaleVector(iint m, dfloat *a, dfloat alpha);

void randomize(iint m, dfloat *a);

dfloat maxEntry(iint n, dfloat *a);

void copyVector(almond_t *almond, iint N, occa::memory o_a, occa::memory o_b);

void scaleVector(almond_t *almond, iint N, occa::memory o_a, dfloat alpha);

void dotStar(almond_t *almond, iint N, occa::memory o_a, occa::memory o_b);

void dotStar(almond_t *almond, iint N, dfloat alpha, occa::memory o_a,
	           occa::memory o_b, dfloat beta, occa::memory o_c);

dfloat innerProd(almond_t *almond, iint N, occa::memory o_a, occa::memory o_b);

void vectorAdd(almond_t *almond, iint N, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void vectorAdd(almond_t *almond, iint N, dfloat alpha, occa::memory o_x,
	 dfloat beta, occa::memory o_y, occa::memory o_z);



dfloat innerProd(iint n, dfloat *a, dfloat *b);

void doubleInnerProd(iint n, dfloat *aDotbc, dfloat *a, dfloat *b, dfloat *c);

void kcycleCombinedOp1(iint n, dfloat *aDotbc, dfloat *a, dfloat *b, dfloat *c);

void kcycleCombinedOp2(iint n, dfloat *aDotbcd, dfloat *a, dfloat *b, dfloat *c, dfloat* d) ;

void vectorAdd(iint n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

dfloat vectorAddInnerProd(iint n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

void dotStar(iint m, dfloat *a, dfloat *b);

void scaleVector(iint m, dfloat *a, dfloat alpha);

void setVector(iint m, dfloat *a, dfloat alpha);

void randomize(iint m, dfloat *a);

dfloat maxEntry(iint n, dfloat *a);

void scaleVector(parAlmond_t *parAlmond, iint N, occa::memory o_a, dfloat alpha);

void setVector(parAlmond_t *parAlmond, iint N, occa::memory o_a, dfloat alpha);

void dotStar(parAlmond_t *parAlmond, iint N, occa::memory o_a, occa::memory o_b);

void dotStar(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_a,
	           occa::memory o_b, dfloat beta, occa::memory o_c);

dfloat innerProd(parAlmond_t *parAlmond, iint N, occa::memory o_x, occa::memory o_y);

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(parAlmond_t *parAlmond, iint n, dfloat *aDotbc, occa::memory o_a,
                                        occa::memory o_b, occa::memory o_c);

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(parAlmond_t *parAlmond, iint n, dfloat *aDotbcd, occa::memory o_a,
                                              occa::memory o_b, occa::memory o_c, occa::memory o_d);

// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(parAlmond_t *parAlmond, iint n, dfloat alpha, occa::memory o_x,
                                                          dfloat beta, occa::memory o_y);

void vectorAdd(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void vectorAdd(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_x,
	 dfloat beta, occa::memory o_y, occa::memory o_z);

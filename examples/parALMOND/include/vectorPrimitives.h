

dfloat innerProd(int n, dfloat *a, dfloat *b);

void doubleInnerProd(int n, dfloat *aDotbc, dfloat *a, dfloat *b, dfloat *c);

void kcycleCombinedOp1(int n, dfloat *aDotbc, dfloat *a, dfloat *b, dfloat *c);

void kcycleCombinedOp2(int n, dfloat *aDotbcd, dfloat *a, dfloat *b, dfloat *c, dfloat* d) ;

void vectorAdd(int n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

dfloat vectorAddInnerProd(int n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

void dotStar(int m, dfloat *a, dfloat *b);

void scaleVector(int m, dfloat *a, dfloat alpha);

void setVector(int m, dfloat *a, dfloat alpha);

dfloat sumVector(int m, dfloat *a);

void addScalar(int m, dfloat alpha, dfloat *a);

void randomize(int m, dfloat *a);

dfloat maxEntry(int n, dfloat *a);

void scaleVector(parAlmond_t *parAlmond, int N, occa::memory o_a, dfloat alpha);

void setVector(parAlmond_t *parAlmond, int N, occa::memory o_a, dfloat alpha);

dfloat sumVector(parAlmond_t *parAlmond, int N, occa::memory o_a);

void addScalar(parAlmond_t *parAlmond, int N, dfloat alpha, occa::memory o_a);

void dotStar(parAlmond_t *parAlmond, int N, occa::memory o_a, occa::memory o_b);

void dotStar(parAlmond_t *parAlmond, int N, dfloat alpha, occa::memory o_a,
	           occa::memory o_b, dfloat beta, occa::memory o_c);

dfloat innerProd(parAlmond_t *parAlmond, int N, occa::memory o_x, occa::memory o_y);

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(parAlmond_t *parAlmond, int n, dfloat *aDotbc, occa::memory o_a,
                                        occa::memory o_b, occa::memory o_c);

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(parAlmond_t *parAlmond, int n, dfloat *aDotbcd, occa::memory o_a,
                                              occa::memory o_b, occa::memory o_c, occa::memory o_d);

// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(parAlmond_t *parAlmond, int n, dfloat alpha, occa::memory o_x,
                                                          dfloat beta, occa::memory o_y);

void vectorAdd(parAlmond_t *parAlmond, int N, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void vectorAdd(parAlmond_t *parAlmond, int N, dfloat alpha, occa::memory o_x,
	 dfloat beta, occa::memory o_y, occa::memory o_z);

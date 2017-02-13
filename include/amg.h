
// needs mesh.h for iint and dfloat

typedef struct{

  iint row, col, own;
  dfloat val;

}nonZero_t;

void *amgSetup(dfloat tol, iint *starts, nonZero_t *nonZeros);

iint amgSolve(iint maxIterations, void *amg, dfloat *x, dfloat *rhs);

//iint amgSolve(iint maxIterations, void *amg, occa::memory &x, occa::memory &rhs);

void amgDestroy(void *amg);

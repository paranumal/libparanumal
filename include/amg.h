
// needs mesh.h for int and dfloat

typedef struct{

  int row, col, own;
  dfloat val;

}nonZero_t;

void *amgSetup(dfloat tol, int *starts, nonZero_t *nonZeros);

int amgSolve(int maxIterations, void *amg, dfloat *x, dfloat *rhs);

//int amgSolve(int maxIterations, void *amg, occa::memory &x, occa::memory &rhs);

void amgDestroy(void *amg);

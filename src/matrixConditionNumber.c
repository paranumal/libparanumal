#include "mesh.h"

dfloat matrixConditionNumber(int N, dfloat *A) {

  int lwork = 4*N;
  int info;

  char norm = '1';

  double Acond;
  double Anorm;

  double *tmpLU = (double*) calloc(N*N, sizeof(double));

  int *ipiv = (int*) calloc(N, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));
  int  *iwork = (int*) calloc(N, sizeof(int));

  for(int n=0;n<N*N;++n){
    tmpLU[n] = (double) A[n];
  }

  //get the matrix norm of A
  Anorm = dlange_(&norm, &N, &N, tmpLU, &N, work);

  //compute LU factorization
  dgetrf_ (&N, &N, tmpLU, &N, ipiv, &info);

  //compute inverse condition number
  dgecon_(&norm, &N, tmpLU, &N, &Anorm, &Acond, work, iwork, &info);

  if(info)
    printf("inv: dgetrf/dgecon reports info = %d when computing condition number\n", info);

  free(work);
  free(iwork);
  free(ipiv);
  free(tmpLU);

  return (dfloat) 1.0/Acond;
}
#include "mesh.h"

void matrixInverse(int N, dfloat *A){
  int lwork = N*N;
  int info;

  // compute inverse mass matrix
  double *tmpInvA = (double*) calloc(N*N, sizeof(double));

  int *ipiv = (int*) calloc(N, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));

  for(int n=0;n<N*N;++n){
    tmpInvA[n] = A[n];
  }

  dgetrf_ (&N, &N, tmpInvA, &N, ipiv, &info);
  dgetri_ (&N, tmpInvA, &N, ipiv, work, &lwork, &info);

  if(info)
    printf("inv: dgetrf/dgetri reports info = %d when inverting matrix\n", info);

  for(int n=0;n<N*N;++n)
    A[n] = tmpInvA[n];

  free(work);
  free(ipiv);
  free(tmpInvA);
}

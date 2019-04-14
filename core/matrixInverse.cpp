/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

extern "C" {
  void dgesv_ ( int     *N, int     *NRHS, double  *A,
                int     *LDA,
                int     *IPIV,
                double  *B,
                int     *LDB,
                int     *INFO );

  void sgesv_(int *N, int *NRHS,float  *A, int *LDA, int *IPIV, float  *B, int *LDB,int *INFO);

  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

  void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);
  void sgetri_(int* N, float* A, int* lda, int* IPIV, float* WORK, int* lwork, int* INFO);

  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
              double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );

  double dlange_(char *NORM, int *M, int *N, double *A, int *LDA, double *WORK);
  void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM,
                double *RCOND, double *WORK, int *IWORK, int *INFO );
}

void matrixInverse(int N, double *A){
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

void matrixInverse(int N, float *A){
  int lwork = N*N;
  int info;

  // compute inverse mass matrix
  float *tmpInvA = (float*) calloc(N*N, sizeof(float));

  int *ipiv = (int*) calloc(N, sizeof(int));
  float *work = (float*) calloc(lwork, sizeof(float));

  for(int n=0;n<N*N;++n){
    tmpInvA[n] = A[n];
  }

  sgetrf_ (&N, &N, tmpInvA, &N, ipiv, &info);
  sgetri_ (&N, tmpInvA, &N, ipiv, work, &lwork, &info);

  if(info)
    printf("inv: sgetrf/sgetri reports info = %d when inverting matrix\n", info);

  for(int n=0;n<N*N;++n)
    A[n] = tmpInvA[n];

  free(work);
  free(ipiv);
  free(tmpInvA);
}
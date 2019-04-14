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

  void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);

  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
              double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );

  double dlange_(char *NORM, int *M, int *N, double *A, int *LDA, double *WORK);
  void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM,
                double *RCOND, double *WORK, int *IWORK, int *INFO );

  double slange_(char *NORM, int *M, int *N, float *A, int *LDA, float *WORK);
  void sgecon_(char *NORM, int *N, float *A, int *LDA, float *ANORM,
                float *RCOND, float *WORK, int *IWORK, int *INFO );
}

double matrixConditionNumber(int N, double *A) {

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

  return (double) 1.0/Acond;
}

float matrixConditionNumber(int N, float *A) {

  int lwork = 4*N;
  int info;

  char norm = '1';

  float Acond;
  float Anorm;

  float *tmpLU = (float*) calloc(N*N, sizeof(float));

  int *ipiv = (int*) calloc(N, sizeof(int));
  float *work = (float*) calloc(lwork, sizeof(float));
  int  *iwork = (int*) calloc(N, sizeof(int));

  for(int n=0;n<N*N;++n){
    tmpLU[n] = (float) A[n];
  }

  //get the matrix norm of A
  Anorm = slange_(&norm, &N, &N, tmpLU, &N, work);

  //compute LU factorization
  sgetrf_ (&N, &N, tmpLU, &N, ipiv, &info);

  //compute inverse condition number
  sgecon_(&norm, &N, tmpLU, &N, &Anorm, &Acond, work, iwork, &info);

  if(info)
    printf("inv: sgetrf/sgecon reports info = %d when computing condition number\n", info);

  free(work);
  free(iwork);
  free(ipiv);
  free(tmpLU);

  return (float) 1.0/Acond;
}
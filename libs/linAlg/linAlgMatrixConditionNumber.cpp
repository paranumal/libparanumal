/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "linAlg.hpp"

extern "C" {
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);

  double dlange_(char *NORM, int *M, int *N, double *A, int *LDA, double *WORK);
  double slange_(char *NORM, int *M, int *N, float *A, int *LDA, float *WORK);

  void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM,
                double *RCOND, double *WORK, int *IWORK, int *INFO );
  void sgecon_(char *NORM, int *N, float *A, int *LDA, float *ANORM,
                float *RCOND, float *WORK, int *IWORK, int *INFO );
}

namespace libp {

double linAlg_t::matrixConditionNumber(const int N, const memory<double> A) {

  int n = N;
  int lwork = 4*N;
  int info;

  char norm = '1';

  double Acond;
  double Anorm;

  memory<double> tmpLU(N*N);

  memory<int> ipiv(N);
  memory<double> work(lwork);
  memory<int> iwork(N);

  tmpLU.copyFrom(A, N*N);

  //get the matrix norm of A
  Anorm = dlange_(&norm, &n, &n, tmpLU.ptr(), &n, work.ptr());

  //compute LU factorization
  dgetrf_ (&n, &n, tmpLU.ptr(), &n, ipiv.ptr(), &info);

  LIBP_ABORT("dgetrf reports info = " << info << " when computing condition number",
             info);

  //compute inverse condition number
  dgecon_(&norm, &n, tmpLU.ptr(), &n, &Anorm, &Acond, work.ptr(), iwork.ptr(), &info);

  LIBP_ABORT("dgecon reports info = " << info << " when computing condition number",
             info);

  return 1.0/Acond;
}

float linAlg_t::matrixConditionNumber(const int N, const memory<float> A) {

  int n = N;
  int lwork = 4*N;
  int info;

  char norm = '1';

  float Acond;
  float Anorm;

  memory<float> tmpLU(N*N);

  memory<int> ipiv(N);
  memory<float> work(lwork);
  memory<int> iwork(N);

  tmpLU.copyFrom(A, N*N);

  //get the matrix norm of A
  Anorm = slange_(&norm, &n, &n, tmpLU.ptr(), &n, work.ptr());

  //compute LU factorization
  sgetrf_ (&n, &n, tmpLU.ptr(), &n, ipiv.ptr(), &info);

  LIBP_ABORT("sgetrf reports info = " << info << " when computing condition number",
             info);

  //compute inverse condition number
  sgecon_(&norm, &n, tmpLU.ptr(), &n, &Anorm, &Acond, work.ptr(), iwork.ptr(), &info);

  LIBP_ABORT("sgecon reports info = " << info << " when computing condition number",
             info);

  return 1.0/Acond;
}

} //namespace libp

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

#include "core.hpp"

extern "C" {
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

  void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);
  void sgetri_(int* N, float* A, int* lda, int* IPIV, float* WORK, int* lwork, int* INFO);
}

void matrixInverse(int N, double *A){
  int lwork = N*N;
  int info;

  // compute inverse mass matrix
  int *ipiv = (int*) calloc(N, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));

  dgetrf_ (&N, &N, A, &N, ipiv, &info);

  if(info) {
    std::stringstream ss;
    ss << "dgetrf_ reports info = " << info;
    LIBP_ABORT(ss.str());
  }

  dgetri_ (&N, A, &N, ipiv, work, &lwork, &info);

  if(info) {
    std::stringstream ss;
    ss << "dgetri_ reports info = " << info;
    LIBP_ABORT(ss.str());
  }

  free(work);
  free(ipiv);
}

void matrixInverse(int N, float *A){
  int lwork = N*N;
  int info;

  // compute inverse mass matrix
  int *ipiv = (int*) calloc(N, sizeof(int));
  float *work = (float*) calloc(lwork, sizeof(float));

  sgetrf_ (&N, &N, A, &N, ipiv, &info);

  if(info) {
    std::stringstream ss;
    ss << "sgetrf_ reports info = " << info;
    LIBP_ABORT(ss.str());
  }

  sgetri_ (&N, A, &N, ipiv, work, &lwork, &info);

  if(info) {
    std::stringstream ss;
    ss << "sgetri_ reports info = " << info;
    LIBP_ABORT(ss.str());
  }

  free(work);
  free(ipiv);
}

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
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

  void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);
  void sgetri_(int* N, float* A, int* lda, int* IPIV, float* WORK, int* lwork, int* INFO);
}

namespace libp {

void linAlg_t::matrixInverse(const int N, memory<double> A){
  int n = N;
  int lwork = N*N;
  int info;

  // compute inverse matrix in-place
  memory<int> ipiv(N);
  memory<double> work(lwork);

  dgetrf_ (&n, &n, A.ptr(), &n, ipiv.ptr(), &info);

  LIBP_ABORT("dgetrf_ reports info = " << info, info);

  dgetri_ (&n, A.ptr(), &n, ipiv.ptr(), work.ptr(), &lwork, &info);

  LIBP_ABORT("dgetri_ reports info = " << info, info);
}

void linAlg_t::matrixInverse(const int N, memory<float> A){
  int n = N;
  int lwork = N*N;
  int info;

  // compute inverse matrix in-place
  memory<int> ipiv(N);
  memory<float> work(lwork);

  sgetrf_ (&n, &n, A.ptr(), &n, ipiv.ptr(), &info);

  LIBP_ABORT("sgetrf_ reports info = " << info, info);

  sgetri_ (&n, A.ptr(), &n, ipiv.ptr(), work.ptr(), &lwork, &info);

  LIBP_ABORT("sgetri_ reports info = " << info, info);
}

// Pseudo Invverse
void linAlg_t::matrixPseudoInverse(const int M, const int N, memory<double> A){

  memory<double> tmpA(M*N);
  memory<double> V(M*N);

  // Copy the matrix
  for(int m=0; m<M*N; m++)
    tmpA[m] = A[m]; 

  // First compute  inv(trans(A)*A)
  for(int n=0; n<N; n++){
    for(int m=0; m<N; m++){
      double tmp = 0; 
      for(int i = 0; i<M; i++){
       tmp += tmpA[i*N +n]*tmpA[i*N +m];
     }
     V[n*N + m] = tmp; 
   }
 }

 matrixInverse(N, V); 

// (A^T*A)^-1 * A^T
 for(int n=0; n<N; n++){
  for(int m=0; m<M; m++){
    double tmp = 0; 
    for(int i=0; i<N; i++){
      tmp += V[n*N + i]*tmpA[m*N + i];
    }
    A[n*M + m] = tmp; 
  }
}

}


// Pseudo Inverse 
void linAlg_t::matrixPseudoInverse(const int M, const int N, memory<float> A){

  memory<float> tmpA(M*N);
  memory<float> V(M*N);

  // Copy the matrix
  for(int m=0; m<M*N; m++)
    tmpA[m] = A[m]; 

  // First compute  inv(trans(A)*A)
  for(int n=0; n<N; n++){
    for(int m=0; m<N; m++){
      float tmp = 0; 
      for(int i = 0; i<M; i++){
       tmp += tmpA[i*N +n]*tmpA[i*N +m];
     }
     V[n*N + m] = tmp; 
   }
 }
 matrixInverse(N, V); 

// (A^T*A)^-1 * A^T
 for(int n=0; n<N; n++){
  for(int m=0; m<M; m++){
    float tmp = 0; 
    for(int i=0; i<N; i++){
      tmp += V[n*N + i]*tmpA[m*N + i];
    }
    A[n*M + m] = tmp; 
  }
}
}


} //namespace libp

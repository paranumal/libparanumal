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
  void dgemm_(char *transA, char *transB, int* M, int *N, int *K, double *alpha, double* A, int* lda, double * B, int *ldb, double *beta, double *C, int *ldc);
  void sgemm_(char *transA, char *transB, int* M, int *N, int *K, float *alpha, float* A, int* lda, float * B, int *ldb, float *beta, float *C, int *ldc);
}

namespace libp {

void linAlg_t::matrixMultiply(const int NrowsA, const int NcolsA, memory<double> A,
                              const int NrowsB, const int NcolsB, memory<double> B,
                              memory<double> C){

  if(NcolsA!=NrowsB){
    std::cout << "In linAlg_t::matrixMultiply incompatible column size for A and row size for B" << std::endl;
  }

  char transA = 'N';
  char transB = 'N';
  int M = NcolsB;
  int N = NrowsA;
  int K = NrowsB;
  int lda = NcolsB, ldb = NcolsA, ldc = NcolsB;
  double alpha = 1;
  double beta = 0;


  
  dgemm_ (&transA, &transB, &M, &N, &K, &alpha, B.ptr(), &lda, A.ptr(), &ldb, &beta, C.ptr(), &ldc);
}

void linAlg_t::matrixMultiply(const int NrowsA, const int NcolsA, memory<float> A,
                              const int NrowsB, const int NcolsB, memory<float> B,
                              memory<float> C){

  if(NcolsA!=NrowsB){
    std::cout << "In linAlg_t::matrixMultiply incompatible column size for A and row size for B" << std::endl;
  }

  
  char transA = 'N';
  char transB = 'N';
  int M = NcolsB;
  int N = NrowsA;
  int K = NrowsB;
  int lda = NcolsB, ldb = NcolsA, ldc = NcolsA;

  float alpha = 1;
  float beta = 0;
  
  sgemm_ (&transA, &transB, &M, &N, &K, &alpha, B.ptr(), &lda, A.ptr(), &ldb, &beta, C.ptr(), &ldc);
}

} //namespace libp

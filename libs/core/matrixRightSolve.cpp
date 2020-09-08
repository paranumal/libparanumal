/*

The MIT License (MIT)

Copyright (c) 2020 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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
  void dgesv_ ( int     *N, int     *NRHS, double  *A,
                int     *LDA,
                int     *IPIV,
                double  *B,
                int     *LDB,
                int     *INFO );
  void sgesv_ ( int     *N, int     *NRHS, float  *A,
                int     *LDA,
                int     *IPIV,
                float  *B,
                int     *LDB,
                int     *INFO );
}

// C = A/B  = trans(trans(B)\trans(A))
// assume row major
void matrixRightSolve(int NrowsA, int NcolsA, double *A, int NrowsB, int NcolsB, double *B, double *C){

  int info;

  int NrowsX = NcolsB;
  int NcolsX = NrowsB;

  int NrowsY = NcolsA;
  int NcolsY = NrowsA;

  int lwork = NrowsX*NcolsX;

  // compute inverse mass matrix
  double *tmpX = (double*) calloc(NrowsX*NcolsX, sizeof(double));
  double *tmpY = (double*) calloc(NrowsY*NcolsY, sizeof(double));

  int    *ipiv = (int*) calloc(NrowsX, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));

  for(int n=0;n<NrowsX*NcolsX;++n){
    tmpX[n] = B[n];
  }

  for(int n=0;n<NrowsY*NcolsY;++n){
    tmpY[n] =A[n];
  }

  dgesv_(&NrowsX, &NcolsY, tmpX, &NrowsX, ipiv, tmpY, &NrowsY, &info); // ?

  if(info) {
    std::stringstream ss;
    ss << "dgesv_ reports info = " << info;
    LIBP_ABORT(ss.str());
  }

  for(int n=0;n<NrowsY*NcolsY;++n){
    C[n] = tmpY[n];
  }

  free(work);
  free(ipiv);
  free(tmpX);
  free(tmpY);
}

// C = A/B  = trans(trans(B)\trans(A))
// assume row major
void matrixRightSolve(int NrowsA, int NcolsA, float *A, int NrowsB, int NcolsB, float *B, float *C){

  int info;

  int NrowsX = NcolsB;
  int NcolsX = NrowsB;

  int NrowsY = NcolsA;
  int NcolsY = NrowsA;

  int lwork = NrowsX*NcolsX;

  // compute inverse mass matrix
  float *tmpX = (float*) calloc(NrowsX*NcolsX, sizeof(float));
  float *tmpY = (float*) calloc(NrowsY*NcolsY, sizeof(float));

  int    *ipiv = (int*) calloc(NrowsX, sizeof(int));
  float *work = (float*) calloc(lwork, sizeof(float));

  for(int n=0;n<NrowsX*NcolsX;++n){
    tmpX[n] = B[n];
  }

  for(int n=0;n<NrowsY*NcolsY;++n){
    tmpY[n] =A[n];
  }

  sgesv_(&NrowsX, &NcolsY, tmpX, &NrowsX, ipiv, tmpY, &NrowsY, &info); // ?

  if(info) {
    std::stringstream ss;
    ss << "sgesv_ reports info = " << info;
    LIBP_ABORT(ss.str());
  }

  for(int n=0;n<NrowsY*NcolsY;++n){
    C[n] = tmpY[n];
  }

  free(work);
  free(ipiv);
  free(tmpX);
  free(tmpY);
}
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
  void sgeev_(char *JOBVL, char *JOBVR, int *N, float *A, int *LDA, float *WR, float *WI,
              float *VL, int *LDVL, float *VR, int *LDVR, float *WORK, int *LWORK, int *INFO );
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
              double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
}

// compute right eigenvectors
void matrixEigenVectors(int N, double *A, double *VR, double *WR, double *WI){

  char JOBVL = 'N';
  char JOBVR = 'V';
  int LDA = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 8*N;

  double *VL = NULL;
  double *WORK  = (double*) calloc(LWORK,sizeof(double));

  double *tmpA  = (double*) calloc(N*N,sizeof(double));
  double *tmpVR = (double*) calloc(N*N,sizeof(double));

  for(int n=0;n<N;++n){
    for(int m=0;m<N;++m){
      tmpA[n+m*N] = A[n*N+m];
    }
  }

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &N, tmpA, &LDA, WR, WI,
          VL, &LDVL, tmpVR, &LDVR, WORK, &LWORK, &INFO);

  if(INFO) {
    std::stringstream ss;
    ss << "dgeev_ reports info = " << INFO;
    LIBP_ABORT(ss.str());
  }

  for(int n=0;n<N;++n){
    for(int m=0;m<N;++m){
      VR[n+m*N] = tmpVR[n*N+m];
    }
  }

  free(tmpVR);
  free(tmpA);
  free(WORK);
}

// compute right eigenvectors
void matrixEigenVectors(int N, float *A, float *VR, float *WR, float *WI){

  char JOBVL = 'N';
  char JOBVR = 'V';
  int LDA = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 8*N;

  float *VL = NULL;
  float *WORK  = (float*) calloc(LWORK,sizeof(float));

  float *tmpA  = (float*) calloc(N*N,sizeof(float));
  float *tmpVR = (float*) calloc(N*N,sizeof(float));

  for(int n=0;n<N;++n){
    for(int m=0;m<N;++m){
      tmpA[n+m*N] = A[n*N+m];
    }
  }

  int INFO = -999;

  sgeev_ (&JOBVL, &JOBVR, &N, tmpA, &LDA, WR, WI,
          VL, &LDVL, tmpVR, &LDVR, WORK, &LWORK, &INFO);

  if(INFO) {
    std::stringstream ss;
    ss << "sgeev_ reports info = " << INFO;
    LIBP_ABORT(ss.str());
  }

  for(int n=0;n<N;++n){
    for(int m=0;m<N;++m){
      VR[n+m*N] = tmpVR[n*N+m];
    }
  }

  free(tmpVR);
  free(tmpA);
  free(WORK);
}

// compute eigenvalues
void matrixEigenValues(int N, double *A, double *WR, double *WI){

  char JOBVL  = 'N';
  char JOBVR  = 'N';
  int LDA = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 8*N;

  double *VR = nullptr;
  double *VL = nullptr;
  double *WORK  = (double*) calloc(LWORK,sizeof(double));

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
          VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);

  if(INFO) {
    std::stringstream ss;
    ss << "dgeev_ reports info = " << INFO;
    LIBP_ABORT(ss.str());
  }

  free(WORK);
}

// compute eigenvalues
void matrixEigenValues(int N, float *A, float *WR, float *WI){

  char JOBVL  = 'N';
  char JOBVR  = 'N';
  int LDA = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 8*N;

  float *VR = nullptr;
  float *VL = nullptr;
  float *WORK  = (float*) calloc(LWORK,sizeof(float));

  int INFO = -999;

  sgeev_ (&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
          VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);

  if(INFO) {
    std::stringstream ss;
    ss << "sgeev_ reports info = " << INFO;
    LIBP_ABORT(ss.str());
  }

  free(WORK);
}
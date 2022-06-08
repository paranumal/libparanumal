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
  void sgeev_(char *JOBVL, char *JOBVR, int *N, float *A, int *LDA, float *WR, float *WI,
              float *VL, int *LDVL, float *VR, int *LDVR, float *WORK, int *LWORK, int *INFO );
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
              double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
}

namespace libp {

// compute right eigenvectors
void linAlg_t::matrixEigenVectors(const int N, const memory<double> A,
                                  memory<double> VR,
                                  memory<double> WR,
                                  memory<double> WI){

  int n = N;
  char JOBVL = 'N';
  char JOBVR = 'V';
  int LDA = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 8*N;

  memory<double> WORK(LWORK);
  memory<double> tmpA(N*LDA);
  memory<double> tmpVR(N*LDVR);

  //tmpA = A^T (row major to column-major)
  linAlg_t::matrixTranspose(N, N, A, LDA, tmpA, LDA);

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &n, tmpA.ptr(), &LDA, WR.ptr(), WI.ptr(),
          nullptr, &LDVL, tmpVR.ptr(), &LDVR, WORK.ptr(), &LWORK, &INFO);

  LIBP_ABORT("dgeev_ reports info = " << INFO, INFO);

  //VR = tmpVR^T (column major to row major)
  linAlg_t::matrixTranspose(N, N, tmpVR, LDVR, VR, LDVR);
}

// compute right eigenvectors
void linAlg_t::matrixEigenVectors(const int N, const memory<float> A,
                                  memory<float> VR,
                                  memory<float> WR,
                                  memory<float> WI){

  int n = N;
  char JOBVL = 'N';
  char JOBVR = 'V';
  int LDA = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 8*N;

  memory<float> WORK(LWORK);
  memory<float> tmpA(N*LDA);
  memory<float> tmpVR(N*LDVR);

  //tmpA = A^T (row major to column-major)
  linAlg_t::matrixTranspose(N, N, A, LDA, tmpA, LDA);

  int INFO = -999;

  sgeev_ (&JOBVL, &JOBVR, &n, tmpA.ptr(), &LDA, WR.ptr(), WI.ptr(),
          nullptr, &LDVL, tmpVR.ptr(), &LDVR, WORK.ptr(), &LWORK, &INFO);

  LIBP_ABORT("sgeev_ reports info = " << INFO, INFO);

  //VR = tmpVR^T (column major to row major)
  linAlg_t::matrixTranspose(N, N, tmpVR, LDVR, VR, LDVR);
}

// compute eigenvalues
void linAlg_t::matrixEigenValues(const int N, const memory<double> A,
                                 memory<double> WR,
                                 memory<double> WI){

  int n = N;
  char JOBVL  = 'N';
  char JOBVR  = 'N';
  int LDA = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 8*N;

  double* Aptr = const_cast<double*>(A.ptr());

  memory<double> WORK(LWORK);

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &n, Aptr, &LDA, WR.ptr(), WI.ptr(),
          nullptr, &LDVL, nullptr, &LDVR, WORK.ptr(), &LWORK, &INFO);

  LIBP_ABORT("dgeev_ reports info = " << INFO, INFO);
}

// compute eigenvalues
void linAlg_t::matrixEigenValues(const int N, const memory<float> A,
                                 memory<float> WR,
                                 memory<float> WI){

  int n = N;
  char JOBVL  = 'N';
  char JOBVR  = 'N';
  int LDA = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 8*N;

  float* Aptr = const_cast<float*>(A.ptr());

  memory<float> WORK(LWORK);

  int INFO = -999;

  sgeev_ (&JOBVL, &JOBVR, &n, Aptr, &LDA, WR.ptr(), WI.ptr(),
          nullptr, &LDVL, nullptr, &LDVR, WORK.ptr(), &LWORK, &INFO);

  LIBP_ABORT("sgeev_ reports info = " << INFO, INFO);
}

} //namespace libp

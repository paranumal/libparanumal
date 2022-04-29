/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim WarburtonTim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

  void dgels_ ( char   *TRANS,
                int    *M,
                int    *N,
                int    *NRHS,
                double *A,
                int    *LDA,
                double *B,
                int    *LDB,
                double *WORK,
                int    *LWORK,
                int    *INFO);

  void sgels_ ( char   *TRANS,
                int    *M,
                int    *N,
                int    *NRHS,
                float  *A,
                int    *LDA,
                float  *B,
                int    *LDB,
                float  *WORK,
                int    *LWORK,
                int    *INFO);

  void dgeqp3_( int    *M,
                int    *N,
                double *A,
                int    *LDA,
                int    *JPVT,
                double *TAU,
                double *WORK,
                int    *LWORK,
                int    *INFO);

  void sgeqp3_( int    *M,
                int    *N,
                float  *A,
                int    *LDA,
                int    *JPVT,
                float  *TAU,
                float  *WORK,
                int    *LWORK,
                int    *INFO);

  void dormqr_( char   *SIDE,
                char   *TRANS,
                int    *M,
                int    *N,
                int    *K,
                double *A,
                int    *LDA,
                double *TAU,
                double *C,
                int    *LDC,
                double *WORK,
                int    *LWORK,
                int    *INFO);

  void sormqr_( char   *SIDE,
                char   *TRANS,
                int    *M,
                int    *N,
                int    *K,
                float  *A,
                int    *LDA,
                float  *TAU,
                float  *C,
                int    *LDC,
                float  *WORK,
                int    *LWORK,
                int    *INFO);

  void dtrsm_ ( char   *SIDE,
                char   *UPLO,
                char   *TRANSA,
                char   *DIAG,
                int    *M,
                int    *N,
                double *ALPHA,
                double *A,
                int    *LDA,
                double *B,
                int    *LDB);

  void strsm_ ( char   *SIDE,
                char   *UPLO,
                char   *TRANSA,
                char   *DIAG,
                int    *M,
                int    *N,
                float  *ALPHA,
                float  *A,
                int    *LDA,
                float  *B,
                int    *LDB);
}

namespace libp {

// C = A/B  = trans(trans(B)\trans(A))
// assume row major
void linAlg_t::matrixRightSolve(const int NrowsA, const int NcolsA, const memory<double> A,
                                const int NrowsB, const int NcolsB, const memory<double> B,
                                memory<double> C){

  int info;

  int NrowsX = NcolsB;
  int NcolsX = NrowsB;

  int NrowsY = NcolsA;
  int NcolsY = NrowsA;

  int lwork = NrowsX*NcolsX;

  // compute inverse mass matrix
  memory<double> tmpX(NrowsX*NcolsX);
  memory<int>    ipiv(NrowsX);
  memory<double> work(lwork);

  tmpX.copyFrom(B, NrowsX*NcolsX);
  C.copyFrom(A, NrowsY*NcolsY);

  dgesv_(&NrowsX, &NcolsY, tmpX.ptr(), &NrowsX, ipiv.ptr(), C.ptr(), &NrowsY, &info);

  LIBP_ABORT("dgesv_ reports info = " << info, info);
}

// C = A/B  = trans(trans(B)\trans(A))
// assume row major
void linAlg_t::matrixRightSolve(const int NrowsA, const int NcolsA, const memory<float> A,
                                const int NrowsB, const int NcolsB, const memory<float> B,
                                memory<float> C){

  int info;

  int NrowsX = NcolsB;
  int NcolsX = NrowsB;

  int NrowsY = NcolsA;
  int NcolsY = NrowsA;

  int lwork = NrowsX*NcolsX;

  // compute inverse mass matrix
  memory<float> tmpX(NrowsX*NcolsX);
  memory<int>   ipiv(NrowsX);
  memory<float> work(lwork);

  tmpX.copyFrom(B, NrowsX*NcolsX);
  C.copyFrom(A, NrowsY*NcolsY);

  sgesv_(&NrowsX, &NcolsY, tmpX.ptr(), &NrowsX, ipiv.ptr(), C.ptr(), &NrowsY, &info); // ?

  LIBP_ABORT("sgesv_ reports info = " << info, info);
}

// Find minimum-norm solution to xA = b with NrowsA > NcolsA (underdetermined).
//
// NB:  A must be stored ROW MAJOR.
void linAlg_t::matrixUnderdeterminedRightSolveMinNorm(const int NrowsA, const int NcolsA,
                                                      const memory<double> A, const memory<double> b,
                                                      memory<double> x) {
  // Solve A^T x^T = b^T.  Note TRANS = 'N', since A is row major.
  int  INFO  = 0;
  char TRANS = 'N';
  int  NRHS  = 1;
  int  LWORK = 2*NrowsA*NcolsA;
  int  Nrows = NrowsA;
  int  Ncols = NcolsA;

  memory<double> WORK(LWORK);
  memory<double> tmpA(NrowsA*NcolsA);
  memory<double> tmpb(NrowsA);

  tmpA.copyFrom(A, NrowsA*NcolsA);
  tmpb.copyFrom(b, NcolsA);

  dgels_(&TRANS, &Ncols, &Nrows, &NRHS, tmpA.ptr(), &Ncols, tmpb.ptr(), &Nrows, WORK.ptr(), &LWORK, &INFO);

  LIBP_ABORT("dgels_ returned INFO = " << INFO, INFO);

  // Copy to output.
  x.copyFrom(tmpb, NrowsA);
}

// Find minimum-norm solution to xA = b with NrowsA > NcolsA (underdetermined).
//
// NB:  A must be stored ROW MAJOR.
void linAlg_t::matrixUnderdeterminedRightSolveMinNorm(const int NrowsA, const int NcolsA,
                                                      const memory<float> A, const memory<float> b,
                                                      memory<float> x) {
  // Solve A^T x^T = b^T.  Note TRANS = 'N', since A is row major.
  int  INFO  = 0;
  char TRANS = 'N';
  int  NRHS  = 1;
  int  LWORK = 2*NrowsA*NcolsA;
  int  Nrows = NrowsA;
  int  Ncols = NcolsA;

  memory<float> WORK(LWORK);
  memory<float> tmpA(NrowsA*NcolsA);
  memory<float> tmpb(NrowsA);

  tmpA.copyFrom(A, NrowsA*NcolsA);
  tmpb.copyFrom(b, NcolsA);

  sgels_(&TRANS, &Ncols, &Nrows, &NRHS, tmpA.ptr(), &Ncols, tmpb.ptr(), &Nrows, WORK.ptr(), &LWORK, &INFO);

  LIBP_ABORT("dgels_ returned INFO = " << INFO, INFO);

  // Copy to output.
  x.copyFrom(tmpb, NrowsA);
}

// Solve xA = b with NrowsA > NcolsA (underdetermined) using column-pivoted QR.
//
// Done by solving A^T x^T = b^T in 4 steps:
//   1.  Decompose A^T * P = Q * R.  -->  Q * R * P^T x^T = b^T
//   2.  Multiply by Q^T.            -->  R * P^T x^T = Q^T b^T
//   3.  Backsolve with R1.          -->  P^T * x^T = R1^{-1} Q^T b^T
//       where R1 = leading NcolsA * NcolsA submatrix of R.
//   4.  Apply permutation.          -->  x^T = P R1^{-1} Q^T b^T
//
// NB:  A must be stored ROW MAJOR.
void linAlg_t::matrixUnderdeterminedRightSolveCPQR(const int NrowsA, const int NcolsA,
                                                   const memory<double> A, const memory<double> b,
                                                   memory<double> x) {
  int INFO  = 0;
  int LWORK = 3*NrowsA + 1;
  int Nrows = NrowsA;
  int Ncols = NcolsA;

  memory<int>    JPVT(NrowsA, 0);
  memory<double> TAU(std::min(NrowsA, NcolsA));

  memory<double> WORK;
  memory<double> tmpA(NrowsA*NcolsA);
  memory<double> tmpb(NrowsA, 0.0);

  WORK.malloc(LWORK);
  tmpA.copyFrom(A, NrowsA*NcolsA);
  tmpb.copyFrom(b, NcolsA);

  // Compute A^T * P = Q * R.
  dgeqp3_(&Ncols, &Nrows, tmpA.ptr(), &Ncols, JPVT.ptr(), TAU.ptr(), WORK.ptr(), &LWORK, &INFO);

  LIBP_ABORT("dgeqp3_ returned INFO = " << INFO, INFO);

  // Compute Q^T * b^T.
  char SIDE = 'L';
  char TRANS = 'T';
  int  NRHS = 1;
  int  NREFLS = NcolsA;

  LWORK = 1;
  WORK.malloc(LWORK);
  dormqr_(&SIDE, &TRANS, &Ncols, &NRHS, &NREFLS, tmpA.ptr(), &Ncols, TAU.ptr(), tmpb.ptr(), &Ncols, WORK.ptr(), &LWORK, &INFO);

  LIBP_ABORT("dormqr_ returned INFO = " << INFO, INFO);

  // Compute R1^{-1} * Q^T * b^T
  SIDE = 'L';
  char UPLO = 'U';
  char TRANSA = 'N';
  char DIAG = 'N';
  NRHS = 1;
  double ALPHA = 1.0;

  dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &Ncols, &NRHS, &ALPHA, tmpA.ptr(), &Ncols, tmpb.ptr(), &Ncols);

  // Apply the permutation.
  for (int i = 0; i < NrowsA; i++)
    x[JPVT[i] - 1] = tmpb[i];
}

// Solve xA = b with NrowsA > NcolsA (underdetermined) using column-pivoted QR.
//
// Done by solving A^T x^T = b^T in 4 steps:
//   1.  Decompose A^T * P = Q * R.  -->  Q * R * P^T x^T = b^T
//   2.  Multiply by Q^T.            -->  R * P^T x^T = Q^T b^T
//   3.  Backsolve with R1.          -->  P^T * x^T = R1^{-1} Q^T b^T
//       where R1 = leading NcolsA * NcolsA submatrix of R.
//   4.  Apply permutation.          -->  x^T = P R1^{-1} Q^T b^T
//
// NB:  A must be stored ROW MAJOR.
void linAlg_t::matrixUnderdeterminedRightSolveCPQR(const int NrowsA, const int NcolsA,
                                                   const memory<float> A, const memory<float> b,
                                                   memory<float> x) {
  int INFO  = 0;
  int LWORK = 3*NrowsA + 1;
  int Nrows = NrowsA;
  int Ncols = NcolsA;

  memory<int>    JPVT(NrowsA, 0);
  memory<float> TAU(std::min(NrowsA, NcolsA));

  memory<float> WORK;
  memory<float> tmpA(NrowsA*NcolsA);
  memory<float> tmpb(NrowsA, 0.0);

  WORK.malloc(LWORK);
  tmpA.copyFrom(A, NrowsA*NcolsA);
  tmpb.copyFrom(b, NcolsA);

  // Compute A^T * P = Q * R.
  sgeqp3_(&Ncols, &Nrows, tmpA.ptr(), &Ncols, JPVT.ptr(), TAU.ptr(), WORK.ptr(), &LWORK, &INFO);

  LIBP_ABORT("dgeqp3_ returned INFO = " << INFO, INFO);

  // Compute Q^T * b^T.
  char SIDE = 'L';
  char TRANS = 'T';
  int  NRHS = 1;
  int  NREFLS = NcolsA;

  LWORK = 1;
  WORK.malloc(LWORK);
  sormqr_(&SIDE, &TRANS, &Ncols, &NRHS, &NREFLS, tmpA.ptr(), &Ncols, TAU.ptr(), tmpb.ptr(), &Ncols, WORK.ptr(), &LWORK, &INFO);

  LIBP_ABORT("dormqr_ returned INFO = " << INFO, INFO);

  // Compute R1^{-1} * Q^T * b^T
  SIDE = 'L';
  char UPLO = 'U';
  char TRANSA = 'N';
  char DIAG = 'N';
  NRHS = 1;
  float ALPHA = 1.0;

  strsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &Ncols, &NRHS, &ALPHA, tmpA.ptr(), &Ncols, tmpb.ptr(), &Ncols);

  // Apply the permutation.
  for (int i = 0; i < NrowsA; i++)
    x[JPVT[i] - 1] = tmpb[i];
}

} //namespace libp

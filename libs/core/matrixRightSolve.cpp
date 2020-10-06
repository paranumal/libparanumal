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

  void dgeqp3_( int    *M,
                int    *N,
                double *A,
                int    *LDA,
                int    *JPVT,
                double *TAU,
                double *WORK,
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

// Find minimum-norm solution to xA = b with NrowsA > NcolsA (underdetermined).
//
// NB:  A must be stored ROW MAJOR.
void matrixUnderdeterminedRightSolveMinNorm(int NrowsA, int NcolsA, dfloat *A, dfloat *b, dfloat *x)
{
  int     LWORK, INFO = 0;
  dfloat* WORK;

  dfloat* tmpA = new dfloat[NrowsA*NcolsA]();
  for (int i = 0; i < NrowsA*NcolsA; i++)
    tmpA[i] = A[i];

  dfloat* tmpb = new dfloat[NrowsA]();
  for (int i = 0; i < NcolsA; i++)
    tmpb[i] = b[i];

  // Solve A^T x^T = b^T.  Note TRANS = 'N', since A is row major.
  char TRANS = 'N';
  int  NRHS = 1;

  LWORK = 2*NrowsA*NcolsA;
  WORK = new dfloat[LWORK]();
  dgels_(&TRANS, &NcolsA, &NrowsA, &NRHS, tmpA, &NcolsA, tmpb, &NrowsA, WORK, &LWORK, &INFO);

  if (INFO != 0) {
    std::stringstream ss;
    ss << "dgels_ returned INFO = " << INFO;
    LIBP_ABORT(ss.str());
  }

  // Copy to output.
  for (int i = 0; i < NrowsA; i++)
    x[i] = tmpb[i];

  delete[] WORK;
  delete[] tmpA;
  delete[] tmpb;
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
void matrixUnderdeterminedRightSolveCPQR(int NrowsA, int NcolsA, dfloat *A, dfloat *b, dfloat *x)
{
  int     LWORK, INFO = 0;
  dfloat* WORK;

  dfloat* tmpA = new dfloat[NrowsA*NcolsA]();
  for (int i = 0; i < NrowsA*NcolsA; i++)
    tmpA[i] = A[i];

  dfloat* tmpb = new dfloat[NrowsA]();
  for (int i = 0; i < NcolsA; i++)
    tmpb[i] = b[i];

  // Compute A^T * P = Q * R.
  int*    JPVT = new int[NrowsA]();
  dfloat* TAU = new dfloat[mymin(NrowsA, NcolsA)]();

  LWORK = 3*NrowsA + 1;
  WORK  = new dfloat[LWORK]();
  dgeqp3_(&NcolsA, &NrowsA, tmpA, &NcolsA, JPVT, TAU, WORK, &LWORK, &INFO);

  if (INFO != 0) {
    std::stringstream ss;
    ss << "dgeqp3_ returned INFO = " << INFO;
    LIBP_ABORT(ss.str());
  }

  delete[] WORK;

  // Compute Q^T * b^T.
  char SIDE = 'L';
  char TRANS = 'T';
  int  NRHS = 1;
  int  NREFLS = NcolsA;

  LWORK = 1;
  WORK  = new dfloat[LWORK]();
  dormqr_(&SIDE, &TRANS, &NcolsA, &NRHS, &NREFLS, tmpA, &NcolsA, TAU, tmpb, &NcolsA, WORK, &LWORK, &INFO);

  if (INFO != 0) {
    std::stringstream ss;
    ss << "dormqr_ returned INFO = " << INFO;
    LIBP_ABORT(ss.str());
  }

  delete[] WORK;

  // Compute R1^{-1} * Q^T * b^T
  SIDE = 'L';
  char UPLO = 'U';
  char TRANSA = 'N';
  char DIAG = 'N';
  NRHS = 1;
  dfloat ALPHA = 1.0;

  dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &NcolsA, &NRHS, &ALPHA, tmpA, &NcolsA, tmpb, &NcolsA);

  // Apply the permutation.
  for (int i = 0; i < NrowsA; i++)
    x[JPVT[i] - 1] = tmpb[i];

  delete[] JPVT;
  delete[] TAU;
  delete[] tmpA;
  delete[] tmpb;
}

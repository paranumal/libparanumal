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

#include "c99.h"
#include "name.h"
#include "types.h"

#if !defined(USE_CBLAS)

#define tensor_dot  PREFIXED_NAME(tensor_dot )
#define tensor_mtxm PREFIXED_NAME(tensor_mtxm)

/* Matrices are always column-major (FORTRAN style) */

double tensor_dot(const double *a, const double *b, uint n)
{
  double sum = 0;
  for(;n;--n) sum += *a++ * *b++;
  return sum;
}

#  if defined(USE_NAIVE_BLAS)
#    define tensor_mxv  PREFIXED_NAME(tensor_mxv )
#    define tensor_mtxv PREFIXED_NAME(tensor_mtxv)
#    define tensor_mxm  PREFIXED_NAME(tensor_mxm )

/* y = A x */
void tensor_mxv(
  double *restrict y, uint ny,
  const double *restrict A, const double *restrict x, uint nx)
{
  uint i;
  for(i=0;i<ny;++i) y[i]=0;
  for(;nx;--nx) {
    const double xk = *x++;
    for(i=0;i<ny;++i) y[i] += (*A++)*xk;
  }
}

/* y = A^T x */
void tensor_mtxv(
  double *restrict y, uint ny,
  const double *restrict A, const double *restrict x, uint nx)
{
  for(;ny;--ny) {
    const double *restrict xp = x;
    uint n = nx;
    double sum = *A++ * *xp++;
    for(--n;n;--n) sum += *A++ * *xp++;
    *y++ = sum;
  }
}

/* C = A * B */
void tensor_mxm(
  double *restrict C, uint nc,
  const double *restrict A, uint na, const double *restrict B, uint nb)
{
  uint i,j,k;
  for(i=0;i<nc*nb;++i) C[i]=0;
  for(j=0;j<nb;++j,C+=nc) {
    const double *restrict A_ = A;
    for(k=0;k<na;++k) {
      const double b = *B++;
      for(i=0;i<nc;++i) C[i] += (*A_++) * b;
    }
  }
}

#  endif

/* C = A^T * B */
void tensor_mtxm(
  double *restrict C, uint nc,
  const double *restrict A, uint na, const double *restrict B, uint nb)
{
  uint i,j;
  for(j=0;j<nb;++j,B+=na) {
    const double *restrict A_ = A;
    for(i=0;i<nc;++i,A_+=na) *C++ = tensor_dot(A_,B,na);
  }
}

#endif


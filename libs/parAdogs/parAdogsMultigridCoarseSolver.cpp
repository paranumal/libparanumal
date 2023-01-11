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

#include "parAdogs.hpp"
#include "parAdogs/parAdogsGraph.hpp"

namespace libp {

namespace paradogs {

void coarseSolver_t::Solve(memory<dfloat>& rhs, memory<dfloat>& x) {

  //gather the global rhs
  comm.Allgatherv(rhs, N,
                  grhs, coarseCounts, coarseOffsets);

  #pragma omp parallel for
  for (int n=0;n<N;++n) {
    dfloat xn=0.0;
    for (int m=0;m<coarseTotal;++m) {
      xn += invA[n*coarseTotal + m]*grhs[m];
    }
    x[n] = xn;
  }
}

void coarseSolver_t::Setup(parCSR& A, memory<dfloat>& null) {

  comm = A.comm;
  int size = comm.size();

  N = static_cast<int>(A.Nrows);
  Nrows = A.Nrows;
  Ncols = A.Ncols;

  coarseCounts.malloc(size);
  coarseOffsets.malloc(size);

  //collect partitioning info
  comm.Allgather(N, coarseCounts);

  coarseTotal=0;
  for (int r=0;r<size;++r) {
    coarseTotal+=coarseCounts[r];
  }
  coarseOffsets[0]=0;
  for (int r=1;r<size;++r) {
    coarseOffsets[r]= coarseOffsets[r-1] + coarseCounts[r-1];
  }

  //gather global null vector
  memory<dfloat> gnull(coarseTotal);

  comm.Allgatherv( null, N,
                  gnull, coarseCounts, coarseOffsets);

  //populate local dense matrix
  memory<dfloat> localA(N*coarseTotal);

  /*Fill the matrix with the null boost*/
  #pragma omp parallel for
  for (int n=0;n<N;n++) {
    for (int m=0;m<coarseTotal;m++) {
      localA[n*coarseTotal+m] = null[n]*gnull[m];
    }
  }
  gnull.free();

  /*Add sparse entries*/
  #pragma omp parallel for
  for (int n=0;n<N;n++) {
    const int start = static_cast<int>(A.diag.rowStarts[n]);
    const int end   = static_cast<int>(A.diag.rowStarts[n+1]);
    for (int m=start;m<end;m++) {
      const int col = static_cast<int>(A.diag.cols[m] + A.colOffsetL);
      localA[n*coarseTotal+col] += A.diag.vals[m];
    }
  }
  #pragma omp parallel for
  for (int n=0;n<A.offd.nzRows;n++) {
    const int row   = static_cast<int>(A.offd.rows[n]);
    const int start = static_cast<int>(A.offd.mRowStarts[n]);
    const int end   = static_cast<int>(A.offd.mRowStarts[n+1]);
    for (int m=start;m<end;m++) {
      const int col = static_cast<int>(A.colMap[A.offd.cols[m]]);
      localA[row*coarseTotal+col] += A.offd.vals[m];
    }
  }

  //assemble the full matrix
  memory<dfloat> gA(coarseTotal*coarseTotal);

  for (int r=0;r<size;++r) {
    coarseCounts[r] *= coarseTotal;
    coarseOffsets[r] *= coarseTotal;
  }

  comm.Allgatherv(localA, N*coarseTotal,
                      gA, coarseCounts, coarseOffsets);
  localA.free();

  for (int r=0;r<size;++r) {
    coarseCounts[r]  /= coarseTotal;
    coarseOffsets[r] /= coarseTotal;
  }

  linAlg_t::matrixInverse(coarseTotal, gA);

  //diag piece of invA
  invA.malloc(N*coarseTotal);

  #pragma omp parallel for
  for (int n=0;n<N;n++) {
    for (int m=0;m<coarseTotal;m++) {
      invA[n*coarseTotal+m] = gA[(n+A.rowOffsetL)*coarseTotal+m];
    }
  }

  /*Space for global rhs*/
  grhs.malloc(coarseTotal);
}

} //namespace paradogs

} //namespace libp

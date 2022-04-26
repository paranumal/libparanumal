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
#include <limits>

extern "C" {
  void dsyev_ (char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO);
}

namespace libp {

namespace paradogs {

/*Compute Fiedler vector of graph via multilevel heirarchy*/
memory<dfloat>& graph_t::FiedlerVector() {

  /*Fiedler vector on coarsest level*/
  L[Nlevels-1].FiedlerVector();

  /*Project and improve the Fiedler vector to the fine level*/
  for (int l=Nlevels-2;l>=0;--l) {
    /*Prolongate Fiedler vector to fine graph*/
    L[l].P.SpMV(1.0, L[l+1].Fiedler, 0.0, L[l].Fiedler);

    /*Refine the Fiedler vector*/
    Refine(l);
  }

  return L[0].Fiedler;
}



/*Compute Fiedler vector of graph Laplacian*/
void mgLevel_t::FiedlerVector() {

  const int N = static_cast<int>(A.Nrows);

  int size = A.comm.size();
  memory<int> counts(size);
  memory<int> offsets(size);

  //collect partitioning info
  A.comm.Allgather(N, counts);

  int Ntotal=0;
  for (int r=0;r<size;++r) {
    Ntotal+=counts[r];
  }
  offsets[0]=0;
  for (int r=1;r<size;++r) {
    offsets[r]= offsets[r-1] + counts[r-1];
  }

  //populate local dense matrix
  memory<double> localA(N*Ntotal, 0.0);

  /*Add sparse entries*/
  #pragma omp parallel for
  for (int n=0;n<N;n++) {
    const int start = static_cast<int>(A.diag.rowStarts[n]);
    const int end   = static_cast<int>(A.diag.rowStarts[n+1]);
    for (int m=start;m<end;m++) {
      const int col = static_cast<int>(A.diag.cols[m] + A.colOffsetL);
      localA[n*Ntotal+col] += A.diag.vals[m];
    }
  }
  #pragma omp parallel for
  for (int n=0;n<A.offd.nzRows;n++) {
    const int row   = static_cast<int>(A.offd.rows[n]);
    const int start = static_cast<int>(A.offd.mRowStarts[n]);
    const int end   = static_cast<int>(A.offd.mRowStarts[n+1]);
    for (int m=start;m<end;m++) {
      const int col = static_cast<int>(A.colMap[A.offd.cols[m]]);
      localA[row*Ntotal+col] += A.offd.vals[m];
    }
  }

  //assemble the full matrix
  memory<double> M(Ntotal*Ntotal);

  for (int r=0;r<size;++r) {
    counts[r] *= Ntotal;
    offsets[r] *= Ntotal;
  }

  A.comm.Allgatherv(localA, N*Ntotal,
                    M, counts, offsets);

  localA.free();
  counts.free();
  offsets.free();

  /*Call LaPack to find eigen pairs*/
  int INFO = -999;
  char JOBZ='V';
  char UPLO='L';
  int LWORK = -1;
  int LDA = Ntotal;
  double WORKSIZE=0.0;
  memory<double> W(Ntotal);
  dsyev_(&JOBZ, &UPLO, &Ntotal, M.ptr(), &LDA, W.ptr(), &WORKSIZE, &LWORK, &INFO); //Size query

  LWORK = int(WORKSIZE);
  double *WORK= new double[LWORK];
  dsyev_(&JOBZ, &UPLO, &Ntotal, M.ptr(), &LDA, W.ptr(), WORK, &LWORK, &INFO);
  delete[] WORK;

  LIBP_ABORT("Paradogs: dsyev_ reports info = " << INFO << " in FiedlerVector",
             INFO);

  /*Find the second smallest eigenvalue (the smallest is 0)*/
  double min0 = std::numeric_limits<double>::max();
  double min1 = std::numeric_limits<double>::max();
  int minloc0 = -1;
  int minloc1 = -1;
  for (int i=0;i<Ntotal;++i) {
    // printf("Eig[%d] = %f\n", i, W[i]);

    if (W[i]<min0) {
      min1 = min0;
      min0 = W[i];
      minloc1 = minloc0;
      minloc0 = i;
    } else if (W[i]<min1) {
      min1 = W[i];
      minloc1 = i;
    }
  }

  // printf("min1 = %f, minloc1 = %d \n", min1, minloc1);

  memory<double> minV = M + minloc1*Ntotal;
  for (int i=0;i<N;++i) {
    Fiedler[i] = minV[i+A.rowOffsetL];
  }

  // Fiedler vector is already orthogonal to null

  /* Fiedler vector is probably already normalized, but just in case */
  dfloat norm = 0.0;
  for (dlong n=0;n<N;++n) norm += Fiedler[n]*Fiedler[n];
  A.comm.Allreduce(norm);
  norm = sqrt(norm);

  for (dlong n=0;n<N;++n) Fiedler[n] /= norm;
}

} //namespace paradogs

} //namespace libp

/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"
#include "parAlmond/parAlmondAMGSetup.hpp"

namespace libp {

namespace parAlmond {

static strongGraph_t RugeStubenStrength(parCSR& A, pfloat theta);
static strongGraph_t SymmetricStrength(parCSR& A, pfloat theta);

strongGraph_t strongGraph(parCSR& A, StrengthType type, pfloat theta){

  if (type==RUGESTUBEN) {
    return RugeStubenStrength(A, theta);
  } else { // (type==SYMMETRIC)
    return SymmetricStrength(A, theta);
  }

}

static strongGraph_t RugeStubenStrength(parCSR& A, pfloat theta) {

  const dlong N = A.Nrows;
  const dlong M = A.Ncols;

  strongGraph_t C(N, M, A.platform, A.comm);

  C.rowStarts.malloc(N+1);

  memory<pfloat> maxOD(N,0.0);

  memory<pfloat> diagA = A.diagA;

  //find maxOD
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    const int sign = (diagA[i] >= 0) ? 1:-1;

    //local entries
    dlong Jstart = A.diag.rowStarts[i];
    dlong Jend   = A.diag.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = A.diag.cols[jj];
      if (col==i) continue;
      const pfloat OD = -sign*A.diag.vals[jj];
      if(OD > maxOD[i]) maxOD[i] = OD;
    }
    //non-local entries
    Jstart = A.offd.rowStarts[i];
    Jend   = A.offd.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      pfloat OD = -sign*A.offd.vals[jj];
      if(OD > maxOD[i]) maxOD[i] = OD;
    }

    int strong_per_row = 1; // diagonal entry

    //local entries
    Jstart = A.diag.rowStarts[i];
    Jend   = A.diag.rowStarts[i+1];
    for(dlong jj = Jstart; jj<Jend; jj++){
      const dlong col = A.diag.cols[jj];
      if (col==i) continue;
      const pfloat OD = -sign*A.diag.vals[jj];
      if(OD > theta*maxOD[i]) strong_per_row++;
    }
    //non-local entries
    Jstart = A.offd.rowStarts[i];
    Jend   = A.offd.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const pfloat OD = -sign*A.offd.vals[jj];
      if(OD > theta*maxOD[i]) strong_per_row++;
    }
    C.rowStarts[i+1] = strong_per_row;
  }

  // cumulative sum
  C.rowStarts[0] = 0.0;
  for(dlong i=1; i<N+1 ; i++) {
    C.rowStarts[i] += C.rowStarts[i-1];
  }
  C.nnz = C.rowStarts[N];
  C.cols.malloc(C.nnz);


  // fill in the columns for strong connections
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    const int sign = (diagA[i] >= 0) ? 1:-1;

    dlong counter = C.rowStarts[i];

    //local entries
    dlong Jstart = A.diag.rowStarts[i];
    dlong Jend   = A.diag.rowStarts[i+1];
    for(dlong jj = Jstart; jj<Jend; jj++){
      const dlong col = A.diag.cols[jj];
      if (col==i) {
        C.cols[counter++] = col;// diag entry
        continue;
      }

      const pfloat OD = -sign*A.diag.vals[jj];
      if(OD > theta*maxOD[i])
        C.cols[counter++] = col;
    }
    //nonlocal entries
    Jstart = A.offd.rowStarts[i];
    Jend = A.offd.rowStarts[i+1];
    for(dlong jj = Jstart; jj<Jend; jj++){
      const dlong col = A.offd.cols[jj];
      const pfloat OD = -sign*A.offd.vals[jj];
      if(OD > theta*maxOD[i])
        C.cols[counter++] = col;
    }
  }

  return C;
}

static strongGraph_t SymmetricStrength(parCSR& A, pfloat theta) {

  const dlong N = A.Nrows;
  const dlong M = A.Ncols;

  strongGraph_t C(N, M, A.platform, A.comm);

  C.rowStarts.malloc(N+1);

  memory<pfloat> diagA = A.diagA;

  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int strong_per_row = 1; // diagonal entry

    const pfloat Aii = std::abs(diagA[i]);

    //local entries
    dlong Jstart = A.diag.rowStarts[i];
    dlong Jend   = A.diag.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = A.diag.cols[jj];
      if (col==i) continue;

      const pfloat Ajj = std::abs(diagA[col]);

      if(std::abs(A.diag.vals[jj]) > theta*(sqrt(Aii*Ajj)))
        strong_per_row++;
    }
    //non-local entries
    Jstart = A.offd.rowStarts[i];
    Jend   = A.offd.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = A.offd.cols[jj];
      const pfloat Ajj = std::abs(diagA[col]);

      if(std::abs(A.offd.vals[jj]) > theta*(sqrt(Aii*Ajj)))
        strong_per_row++;
    }

    C.rowStarts[i+1] = strong_per_row;
  }

  // cumulative sum
  C.rowStarts[0] = 0;
  for(dlong i=1; i<N+1 ; i++) {
    C.rowStarts[i] += C.rowStarts[i-1];
  }
  C.nnz = C.rowStarts[N];
  C.cols.malloc(C.nnz);


  // fill in the columns for strong connections
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    const pfloat Aii = std::abs(diagA[i]);

    dlong counter = C.rowStarts[i];

    //local entries
    dlong Jstart = A.diag.rowStarts[i];
    dlong Jend   = A.diag.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = A.diag.cols[jj];
      if (col==i) {
        C.cols[counter++] = col;// diag entry
        continue;
      }

      const pfloat Ajj = std::abs(diagA[col]);

      if(std::abs(A.diag.vals[jj]) > theta*(sqrt(Aii*Ajj)))
        C.cols[counter++] = col;
    }
    //non-local entries
    Jstart = A.offd.rowStarts[i];
    Jend   = A.offd.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = A.offd.cols[jj];

      const pfloat Ajj = std::abs(diagA[col]);

      if(std::abs(A.offd.vals[jj]) > theta*(sqrt(Aii*Ajj)))
        C.cols[counter++] = col;
    }
  }

  return C;
}

} //namespace parAlmond

} //namespace libp

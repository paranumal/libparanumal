/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

namespace parAlmond {

parCSR *constructProlongation(parCSR *A, hlong *FineToCoarse,
                            hlong *globalAggStarts, dfloat *null){

  int rank;
  MPI_Comm_rank(A->comm, &rank);

  const dlong N = A->Nrows;

  const hlong globalAggOffset = globalAggStarts[rank];
  const dlong NCoarse = (dlong) (globalAggStarts[rank+1]-globalAggStarts[rank]); //local num agg

  parCSR* P = new parCSR(N, NCoarse, A->platform, A->comm);

  P->globalRowStarts = A->globalRowStarts;
  P->globalColStarts = globalAggStarts;

  P->diag.rowStarts = (dlong *) calloc(N+1, sizeof(dlong));

  dlong* offdRowCounts = (dlong *) calloc(N+1, sizeof(dlong));

  // each row has exactly one nonzero
  for(dlong i=0; i<N; i++) {
    const hlong col = FineToCoarse[i];
    if ((col>globalAggOffset-1)&&(col<globalAggOffset+NCoarse)) {
      P->diag.rowStarts[i+1]++;
    } else {
      offdRowCounts[i+1]++;
    }
  }

  // count how many rows are shared
  P->offd.nzRows=0;
  for(dlong i=0; i<N; i++)
    if (offdRowCounts[i+1]>0) P->offd.nzRows++;

  P->offd.rows      = (dlong *) calloc(P->offd.nzRows, sizeof(dlong));
  P->offd.rowStarts = (dlong *) calloc(P->offd.nzRows+1, sizeof(dlong));

  // cumulative sum
  dlong cnt=0;
  for(dlong i=0; i<N; i++) {

    P->diag.rowStarts[i+1] += P->diag.rowStarts[i];

    if (offdRowCounts[i+1]>0) {
      P->offd.rows[cnt] = i; //record row id
      P->offd.rowStarts[cnt+1] = P->offd.rowStarts[cnt] + offdRowCounts[i+1];
      cnt++;
    }
  }
  P->diag.nnz = P->diag.rowStarts[N];
  P->offd.nnz = P->offd.rowStarts[P->offd.nzRows];

  free(offdRowCounts);

  // Halo setup
  cnt=0;
  hlong *colIds = (hlong *) malloc(P->offd.nnz*sizeof(hlong));
  for (dlong i=0;i<N;i++) {
    hlong col = FineToCoarse[i];
    if ((col<globalAggOffset)||(col>globalAggOffset+NCoarse-1))
      colIds[cnt++] = col;
  }
  P->haloSetup(colIds); //setup halo, and transform colIds to a local indexing

  //fill entries of P with null vector
  P->diag.cols = (dlong *)  calloc(P->diag.nnz, sizeof(dlong));
  P->diag.vals = (dfloat *) calloc(P->diag.nnz, sizeof(dfloat));
  P->offd.cols = (dlong *)  calloc(P->offd.nnz, sizeof(dlong));
  P->offd.vals = (dfloat *) calloc(P->offd.nnz, sizeof(dfloat));

  dlong diagCnt = 0;
  dlong offdCnt = 0;
  for(dlong i=0; i<N; i++) {
    const hlong col = FineToCoarse[i];
    if ((col>globalAggStarts[rank]-1)&&(col<globalAggStarts[rank+1])) {
      P->diag.cols[diagCnt  ] = (dlong) (col - globalAggOffset); //local index
      P->diag.vals[diagCnt++] = null[i];
    } else {
      P->offd.cols[offdCnt  ] = colIds[offdCnt];
      P->offd.vals[offdCnt++] = null[i];
    }
  }

  // normalize the columns of P

  //check size. If this ever triggers, we'll have to implement a re-alloc of null
  if (P->Ncols > A->Ncols)
    LIBP_ABORT(string("Size of Coarse nullvector is too large, need to re-alloc"))

  //set coarse null to 0
  for(dlong i=0; i<P->Ncols; i++) null[i] = 0.0;

  //add local nonzeros
  for(dlong i=0; i<P->diag.nnz; i++)
    null[P->diag.cols[i]] += P->diag.vals[i] * P->diag.vals[i];

  //add nonlocal nonzeros
  for(dlong i=0; i<P->offd.nnz; i++)
    null[P->offd.cols[i]] += P->offd.vals[i] * P->offd.vals[i];

  //add the halo values to their origins
  P->halo->Combine(null, 1, ogs_dfloat);

  for(dlong i=0; i<NCoarse; i++)
    null[i] = sqrt(null[i]);

  //share the results
  P->halo->Exchange(null, 1, ogs_dfloat);

  for(dlong i=0; i<P->diag.nnz; i++)
    P->diag.vals[i] /= null[P->diag.cols[i]];
  for(dlong i=0; i<P->offd.nnz; i++)
    P->offd.vals[i] /= null[P->offd.cols[i]];

  return P;
}

} //namespace parAlmond
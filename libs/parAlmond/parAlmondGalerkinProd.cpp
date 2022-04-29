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

parCSR galerkinProd(parCSR& A, parCSR& P){

  // MPI info
  int rank = A.comm.rank();
  int size = A.comm.size();

  memory<hlong> globalAggStarts = P.globalColStarts;
  hlong globalAggOffset = globalAggStarts[rank];

  //The galerkin product can be computed as
  // (P^T A P)_IJ = sum_{i in Agg_I} sum_{j in Agg_J} P_iI A_ij P_jJ
  // Since each row of P has only one entry, we can share the necessary
  // P entries, form the products, and send them to their destination rank

  const dlong N = A.Nrows;
  const dlong M = A.Ncols;

  //printf("Level has %d rows, and is making %d aggregates\n", N, globalAggStarts[rank+1]-globalAggStarts[rank]);

  // Exploit the fact that we know P has one non-zero per row to
  // compress the global Ids of the columns and nonzero values to
  // single vectors
  memory<hlong>  Pcols(M);
  memory<pfloat> Pvals(M);

  //record the entries of P that this rank has
  for (dlong i=0;i<N;i++) {
    for (dlong j=P.diag.rowStarts[i];j<P.diag.rowStarts[i+1];j++) {
      Pcols[i] = P.diag.cols[j] + globalAggOffset; //global ID
      Pvals[i] = P.diag.vals[j];
    }
  }
  for (dlong i=0;i<P.offd.nzRows;i++) {
    const dlong row = P.offd.rows[i];
    for (dlong j=P.offd.mRowStarts[i];j<P.offd.mRowStarts[i+1];j++) {
      Pcols[row] = P.colMap[P.offd.cols[j]]; //global ID
      Pvals[row] = P.offd.vals[j];
    }
  }

  //fill the halo region
  A.halo.Exchange(Pcols, 1);
  A.halo.Exchange(Pvals, 1);

  dlong sendNtotal = A.diag.nnz+A.offd.nnz;
  memory<parCOO::nonZero_t> sendPTAP(sendNtotal);

  //form the fine PTAP products
  dlong cnt =0;
  for (dlong i=0;i<N;i++) {
    const dlong start = A.diag.rowStarts[i];
    const dlong end   = A.diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong  col = A.diag.cols[j];
      const dfloat val = A.diag.vals[j];

      sendPTAP[cnt].row = Pcols[i];
      sendPTAP[cnt].col = Pcols[col];
      sendPTAP[cnt].val = val*Pvals[i]*Pvals[col];
      cnt++;
    }
  }
  for (dlong i=0;i<A.offd.nzRows;i++) {
    const dlong row   = A.offd.rows[i];
    const dlong start = A.offd.mRowStarts[i];
    const dlong end   = A.offd.mRowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong  col = A.offd.cols[j];
      const dfloat val = A.offd.vals[j];

      sendPTAP[cnt].row = Pcols[row];
      sendPTAP[cnt].col = Pcols[col];
      sendPTAP[cnt].val = val*Pvals[row]*Pvals[col];
      cnt++;
    }
  }

  //sort entries by the coarse row and col
  std::sort(sendPTAP.ptr(), sendPTAP.ptr()+sendNtotal,
            [](const parCOO::nonZero_t& a, const parCOO::nonZero_t& b) {
              if (a.row < b.row) return true;
              if (a.row > b.row) return false;

              return a.col < b.col;
            });

  //count number of non-zeros we're sending
  memory<int> sendCounts(size,0);
  memory<int> recvCounts(size);
  memory<int> sendOffsets(size+1);
  memory<int> recvOffsets(size+1);

  int r=0;
  for(dlong i=0;i<sendNtotal;++i) {
    hlong id = sendPTAP[i].row;
    while(id>=globalAggStarts[r+1]) r++;
    sendCounts[r]++;
  }

  // find how many nodes to expect (should use sparse version)
  A.comm.Alltoall(sendCounts, recvCounts);

  // find send and recv offsets for gather
  sendOffsets[0] = 0;
  recvOffsets[0] = 0;
  for(int rr=0;rr<size;++rr){
    sendOffsets[rr+1] = sendOffsets[rr] + sendCounts[rr];
    recvOffsets[rr+1] = recvOffsets[rr] + recvCounts[rr];
  }
  dlong recvNtotal = recvOffsets[size];

  memory<parCOO::nonZero_t> recvPTAP(recvNtotal);

  A.comm.Alltoallv(sendPTAP, sendCounts, sendOffsets,
                   recvPTAP, recvCounts, recvOffsets);

  //sort entries by the coarse row and col
  std::sort(recvPTAP.ptr(), recvPTAP.ptr()+recvNtotal,
            [](const parCOO::nonZero_t& a, const parCOO::nonZero_t& b) {
              if (a.row < b.row) return true;
              if (a.row > b.row) return false;

              return a.col < b.col;
            });

  //count total number of nonzeros;
  dlong nnz =0;
  if (recvNtotal) nnz++;
  for (dlong i=1;i<recvNtotal;i++)
    if ((recvPTAP[i].row!=recvPTAP[i-1].row)||
        (recvPTAP[i].col!=recvPTAP[i-1].col)) nnz++;


  parCOO PTAP(A.platform, A.comm);

  //copy global partition
  PTAP.globalRowStarts = globalAggStarts;
  PTAP.globalColStarts = globalAggStarts;

  PTAP.nnz = nnz;
  PTAP.entries.malloc(PTAP.nnz);

  //compress nonzeros
  nnz = 0;
  if (recvNtotal) PTAP.entries[nnz++] = recvPTAP[0];
  for (dlong i=1;i<recvNtotal;i++) {
    if ((recvPTAP[i].row!=recvPTAP[i-1].row)||
        (recvPTAP[i].col!=recvPTAP[i-1].col)) {
      PTAP.entries[nnz++] = recvPTAP[i];
    } else {
      PTAP.entries[nnz-1].val += recvPTAP[i].val;
    }
  }

  //build Ac from coo matrix
  return parCSR(PTAP);
}

} //namespace parAlmond

} //namespace libp

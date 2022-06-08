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

parCSR transpose(parCSR& A){

  // MPI info
  int rank = A.comm.rank();
  int size = A.comm.size();

  // copy data from nonlocal entries into send buffer
  memory<parCOO::nonZero_t> sendNonZeros(A.offd.nnz);
  for(dlong i=0;i<A.offd.nzRows;++i){
    const hlong row = A.offd.rows[i] + A.globalRowStarts[rank]; //global ids
    for (dlong j=A.offd.mRowStarts[i];j<A.offd.mRowStarts[i+1];j++) {
      const hlong col =  A.colMap[A.offd.cols[j]]; //global ids
      sendNonZeros[j].row = col;
      sendNonZeros[j].col = row;
      sendNonZeros[j].val = A.offd.vals[j];
    }
  }

  //sort by destination row
  std::sort(sendNonZeros.ptr(), sendNonZeros.ptr()+A.offd.nnz,
            [](const parCOO::nonZero_t& a, const parCOO::nonZero_t& b) {
              if (a.row < b.row) return true;
              if (a.row > b.row) return false;

              return a.col < b.col;
            });

  //count number of non-zeros we're sending
  memory<int> sendCounts(size, 0);
  memory<int> recvCounts(size);
  memory<int> sendOffsets(size+1);
  memory<int> recvOffsets(size+1);

  int r=0;
  for (dlong n=0;n<A.offd.nnz;n++) {
    dlong row = sendNonZeros[n].row;
    while(row>=A.globalColStarts[r+1]) r++;
    sendCounts[r]++;
  }

  A.comm.Alltoall(sendCounts, recvCounts);

  sendOffsets[0] = 0;
  recvOffsets[0] = 0;
  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }
  dlong offdnnz = recvOffsets[size]; //total offd nonzeros


  parCOO cooAt(A.platform, A.comm);

  //copy global partition
  cooAt.globalRowStarts = A.globalColStarts;
  cooAt.globalColStarts = A.globalRowStarts;

  cooAt.nnz = A.diag.nnz+offdnnz;
  cooAt.entries.malloc(cooAt.nnz);

  //fill local nonzeros
  for(dlong i=0; i<A.Nrows; i++){
    const dlong Jstart = A.diag.rowStarts[i];
    const dlong Jend   = A.diag.rowStarts[i+1];

    for(dlong jj=Jstart; jj<Jend; jj++){
      cooAt.entries[jj].row = A.diag.cols[jj] + A.globalColStarts[rank];
      cooAt.entries[jj].col = i + A.globalRowStarts[rank];
      cooAt.entries[jj].val = A.diag.vals[jj];
    }
  }

  // receive non-local nonzeros
  A.comm.Alltoallv(sendNonZeros,             sendCounts, sendOffsets,
                   cooAt.entries+A.diag.nnz, recvCounts, recvOffsets);

  //sort by row
  std::sort(cooAt.entries.ptr(), cooAt.entries.ptr()+cooAt.nnz,
            [](const parCOO::nonZero_t& a, const parCOO::nonZero_t& b) {
              if (a.row < b.row) return true;
              if (a.row > b.row) return false;

              return a.col < b.col;
            });

  return parCSR(cooAt);
}

} //namespace parAlmond

} //namespace libp

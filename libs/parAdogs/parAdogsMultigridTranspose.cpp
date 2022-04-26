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

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

namespace libp {

namespace paradogs {

parCSR Transpose(const parCSR& A) {

  // MPI info
  int size = A.comm.size();

  // copy data from nonlocal entries into send buffer
  memory<nonZero_t> sendNonZeros(A.offd.nnz);
  for(dlong i=0;i<A.offd.nzRows;++i){
    const hlong row = A.offd.rows[i] + A.rowOffsetL; //global ids
    for (dlong j=A.offd.mRowStarts[i];j<A.offd.mRowStarts[i+1];j++) {
      const hlong col =  A.colMap[A.offd.cols[j]]; //global ids
      sendNonZeros[j].row = col;
      sendNonZeros[j].col = row;
      sendNonZeros[j].val = A.offd.vals[j];
    }
  }

  //sort by destination row
  std::sort(sendNonZeros.ptr(), sendNonZeros.ptr()+A.offd.nnz,
            [](const nonZero_t& a, const nonZero_t& b) {
              if (a.row < b.row) return true;
              if (a.row > b.row) return false;

              return a.col < b.col;
            });

  // //count number of non-zeros we're sending
  memory<int> sendCounts(size, 0);
  memory<int> recvCounts(size);
  memory<int> sendOffsets(size+1);
  memory<int> recvOffsets(size+1);

  memory<hlong> globalColStarts(size+1);
  globalColStarts[0]=0;
  A.comm.Allgather(A.colOffsetU, globalColStarts+1);

  int r=0;
  for (dlong n=0;n<A.offd.nnz;n++) {
    dlong row = sendNonZeros[n].row;
    while(row>=globalColStarts[r+1]) r++;
    sendCounts[r]++;
  }
  globalColStarts.free();

  A.comm.Alltoall(sendCounts, recvCounts);

  sendOffsets[0]=0;
  recvOffsets[0]=0;
  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }
  dlong offdnnz = recvOffsets[size]; //total offd nonzeros

  memory<nonZero_t> offdNonZeros(offdnnz);

  // receive non-local nonzeros
  A.comm.Alltoallv(sendNonZeros, sendCounts, sendOffsets,
                   offdNonZeros, recvCounts, recvOffsets);

  //clean up
  sendNonZeros.free();
  sendCounts.free();
  recvCounts.free();
  sendOffsets.free();
  recvOffsets.free();

  dlong NNZ = A.diag.nnz+offdnnz;

  memory<nonZero_t> entries(NNZ);

  memory<dlong> rowStarts(A.NlocalCols+1, 0);
  memory<dlong> rowCounts(A.NlocalCols, 0);

  /*Count entries per row*/
  for(dlong i=0; i<A.Nrows; i++) {
    const dlong Jstart = A.diag.rowStarts[i];
    const dlong Jend   = A.diag.rowStarts[i+1];

    for(dlong jj=Jstart; jj<Jend; jj++){
      rowStarts[A.diag.cols[jj]+1]++;
    }
  }
  for(dlong i=0; i<offdnnz; i++) {
    const dlong row = static_cast<dlong>(offdNonZeros[i].row-A.colOffsetL);
    rowStarts[row+1]++;
  }

  /*Cumulative sum*/
  for(dlong i=1; i<A.NlocalCols+1; i++) {
    rowStarts[i] += rowStarts[i-1];
  }

  //fill local nonzeros
  // #pragma omp parallel for
  for(dlong i=0; i<A.Nrows; i++){
    const dlong Jstart = A.diag.rowStarts[i];
    const dlong Jend   = A.diag.rowStarts[i+1];

    for(dlong jj=Jstart; jj<Jend; jj++){
      const dlong row = A.diag.cols[jj];
      const dlong c = rowStarts[row] + rowCounts[row];

      entries[c].row = row + A.colOffsetL;
      entries[c].col = i + A.rowOffsetL;
      entries[c].val = A.diag.vals[jj];
      rowCounts[row]++;
    }
  }
  for(dlong i=0; i<offdnnz; i++) {
    const dlong row = static_cast<dlong>(offdNonZeros[i].row-A.colOffsetL);
    const dlong c = rowStarts[row] + rowCounts[row];
    entries[c] = offdNonZeros[i];
    rowCounts[row]++;
  }

  offdNonZeros.free();

  //sort each row by column id
  #pragma omp parallel for
  for(dlong i=0; i<A.NlocalCols; i++){
    const dlong Nentries = rowStarts[i+1]-rowStarts[i];
    const dlong c = rowStarts[i];
    std::sort(entries.ptr()+c, entries.ptr()+c+Nentries,
          [](const nonZero_t& a, const nonZero_t& b) {
            return a.col < b.col;
          });
  }

  return parCSR(A.NlocalCols, A.Nrows,
                NNZ, entries,
                A.platform, A.comm);
}

} //namespace paradogs

} //namespace libp

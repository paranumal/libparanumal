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
#include "parAdogs/parAdogsMatrix.hpp"
#include "parAdogs/parAdogsPartition.hpp"

namespace libp {

namespace paradogs {

parCSR SpMM(const parCSR& A, const parCSR& B){

  // MPI info
  int size = A.comm.size();

  // To compute C = A*B we need all the rows B(j,:) for which
  // j is a column index for the nonzeros of A on this rank.
  // For all local column indices in A.diag, we will already
  // have the row of B on this rank, so we just need to gather
  // the offd colIds

  memory<hlong> recvRows(A.Ncols-A.NlocalCols);
  memory<int> sendCounts(size);
  memory<int> recvCounts(size, 0);
  memory<int> sendOffsets(size+1);
  memory<int> recvOffsets(size+1);

  memory<hlong> globalRowStarts(size+1);
  globalRowStarts[0]=0;
  B.comm.Allgather(B.rowOffsetU, globalRowStarts+1);

  //use the colMap of A to list the needed rows of B
  int r=0;
  for (dlong n=A.NlocalCols;n<A.Ncols;n++) {
    const hlong id = A.colMap[n];
    while (id>=globalRowStarts[r+1]) r++; //assumes the halo is sorted
    recvCounts[r]++;
    recvRows[n-A.NlocalCols] = id; //record the row to recv
  }
  globalRowStarts.free();

  //share the counts
  A.comm.Alltoall(recvCounts, sendCounts);

  sendOffsets[0]=0;
  recvOffsets[0]=0;
  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }

  int sendTotal = sendOffsets[size];
  memory<hlong> sendRows(sendTotal);

  //share the rowIds
  B.comm.Alltoallv(recvRows, recvCounts, recvOffsets,
                   sendRows, sendCounts, sendOffsets);

  //we now have a list of rows to send, count the nnz to send
  dlong NNZ=0;
  for (r=0;r<size;r++) {
    sendCounts[r] =0; //reset
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      const dlong i = static_cast<dlong>(sendRows[n]-B.rowOffsetL); //local row id
      sendCounts[r]+= B.diag.rowStarts[i+1]-B.diag.rowStarts[i]; //count entries in this row
      sendCounts[r]+= B.offd.rowStarts[i+1]-B.offd.rowStarts[i]; //count entries in this row
    }
    NNZ += sendCounts[r]; //tally the total
  }

  memory<nonZero_t> sendNonZeros(NNZ);

  NNZ=0; //reset
  for (r=0;r<size;r++) {
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      const dlong i = static_cast<dlong>(sendRows[n] - B.rowOffsetL); //local row id
      for (dlong jj=B.diag.rowStarts[i]; jj<B.diag.rowStarts[i+1];jj++){
        sendNonZeros[NNZ].row = sendRows[n];
        sendNonZeros[NNZ].col = B.diag.cols[jj] + B.colOffsetL;
        sendNonZeros[NNZ].val = B.diag.vals[jj];
        NNZ++;
      }
      for (dlong jj=B.offd.rowStarts[i]; jj<B.offd.rowStarts[i+1];jj++){
        sendNonZeros[NNZ].row = sendRows[n];
        sendNonZeros[NNZ].col = B.colMap[B.offd.cols[jj]];
        sendNonZeros[NNZ].val = B.offd.vals[jj];
        NNZ++;
      }
    }
  }

  A.comm.Alltoall(sendCounts, recvCounts);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }


  dlong Boffdnnz = recvOffsets[size]; //total nonzeros
  memory<nonZero_t> BoffdRows(Boffdnnz);

  B.comm.Alltoallv(sendNonZeros, sendCounts, sendOffsets,
                      BoffdRows, recvCounts, recvOffsets);

  //clean up
  sendNonZeros.free();
  sendRows.free();
  recvRows.free();
  sendCounts.free();
  recvCounts.free();
  sendOffsets.free();
  recvOffsets.free();

  //we now have all the needed nonlocal rows (should also be sorted by row then col)

  //make an array of row offsets so we know how large each row is
  memory<dlong> BoffdRowOffsets(A.Ncols-A.NlocalCols+1, 0);

  dlong id=0;
  for (dlong n=0;n<Boffdnnz;n++) {
    hlong row = BoffdRows[n].row;

    while(A.colMap[id+A.NlocalCols]!=row) id++;

    BoffdRowOffsets[id+1]++; //count entry in row
  }

  //cumulative sum
  for (dlong n=0;n<A.Ncols-A.NlocalCols;n++)
    BoffdRowOffsets[n+1] += BoffdRowOffsets[n];


  // The next step to compute C = A*B is to multiply each entry A(i,j) by the
  // row B(j,:), store the all the results, sort them by row+col, and compress
  // the entries

  // Find how big the intermediate form is
  memory<dlong> rowStarts(A.Nrows+1, 0);
  memory<dlong> rowCounts(A.Nrows, 0);

  /*Count entries per row*/
  #pragma omp parallel for
  for (dlong i=0;i<A.Nrows;i++) {
    //local entries
    dlong start = A.diag.rowStarts[i];
    dlong end   = A.diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A.diag.cols[j];
      rowStarts[i+1] +=  B.diag.rowStarts[col+1]-B.diag.rowStarts[col]
                        +B.offd.rowStarts[col+1]-B.offd.rowStarts[col];
    }
    //non-local entries
    start = A.offd.rowStarts[i];
    end   = A.offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A.offd.cols[j]-A.NlocalCols;
      rowStarts[i+1] += BoffdRowOffsets[col+1] - BoffdRowOffsets[col];
    }
  }

  /*Cumulative sum*/
  for(dlong i=1; i<A.Nrows+1; i++) {
    rowStarts[i] += rowStarts[i-1];
  }

  NNZ = rowStarts[A.Nrows];

  memory<nonZero_t> Ctmp(NNZ);

  //count total number of nonzeros;
  dlong nnz =0;

  // Fill the intermediate form of C
  // #pragma omp parallel for
  for (dlong i=0;i<A.Nrows;i++) {
    const dlong cStart = rowStarts[i];
    dlong& c = rowCounts[i];

    //local A entries
    dlong start = A.diag.rowStarts[i];
    dlong end   = A.diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A.diag.cols[j];
      const dfloat Aval = A.diag.vals[j];

      //local B entries
      dlong Bstart = B.diag.rowStarts[col];
      dlong Bend   = B.diag.rowStarts[col+1];
      for (dlong jj=Bstart;jj<Bend;jj++) {
        Ctmp[cStart+c].row = i + A.rowOffsetL;
        Ctmp[cStart+c].col = B.diag.cols[jj] + B.colOffsetL; //global id
        Ctmp[cStart+c].val = Aval*B.diag.vals[jj];
        c++;
      }
      //non-local B entries
      Bstart = B.offd.rowStarts[col];
      Bend   = B.offd.rowStarts[col+1];
      for (dlong jj=Bstart;jj<Bend;jj++) {
        Ctmp[cStart+c].row = i + A.rowOffsetL;
        Ctmp[cStart+c].col = B.colMap[B.offd.cols[jj]]; //global id
        Ctmp[cStart+c].val = Aval*B.offd.vals[jj];
        c++;
      }
    }
    //non-local A entries
    start = A.offd.rowStarts[i];
    end   = A.offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A.offd.cols[j]-A.NlocalCols;
      const dfloat Aval = A.offd.vals[j];

      // entries from recived rows of B
      dlong Bstart = BoffdRowOffsets[col];
      dlong Bend   = BoffdRowOffsets[col+1];
      for (dlong jj=Bstart;jj<Bend;jj++) {
        Ctmp[cStart+c].row = i + A.rowOffsetL;
        Ctmp[cStart+c].col = BoffdRows[jj].col; //global id
        Ctmp[cStart+c].val = Aval*BoffdRows[jj].val;
        c++;
      }
    }

    //sort entries in this row by col id
    std::sort(Ctmp.ptr()+cStart, Ctmp.ptr()+cStart+c,
              [](const nonZero_t& a, const nonZero_t& b) {
                return a.col < b.col;
              });

    /*Count how many actual nonzeros will be in this row*/
    dlong nnzRow=0;
    if (c>0) nnzRow++;
    for (dlong j=1;j<c;j++) {
      if ((Ctmp[cStart+j].col!=Ctmp[cStart+j-1].col)) nnzRow++;
    }

    nnz+=nnzRow; //Add to total
  }
  BoffdRowOffsets.free();
  BoffdRows.free();

  rowStarts.free();
  rowCounts.free();

  // cooC.nnz = nnz;
  memory<nonZero_t> entries(nnz);

  //compress nonzeros
  nnz = 0;
  if (NNZ) entries[nnz++] = Ctmp[0];
  for (dlong i=1;i<NNZ;i++) {
    if ((Ctmp[i].row!=Ctmp[i-1].row)||
        (Ctmp[i].col!=Ctmp[i-1].col)) {
      entries[nnz++] = Ctmp[i];
    } else {
      entries[nnz-1].val += Ctmp[i].val;
    }
  }
  //clean up
  Ctmp.free();

  //build C from coo matrix
  return parCSR(A.Nrows, B.NlocalCols,
                nnz, entries,
                A.platform, A.comm);
}

} //namespace paradogs

} //namespace libp

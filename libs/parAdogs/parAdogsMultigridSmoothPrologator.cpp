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

parCSR SmoothProlongator(const parCSR& A, const parCSR& T) {

  // MPI info
  int size = A.comm.size();

  // This function computes a smoothed prologation operator
  // via a single weighted Jacobi iteration on the tentative
  // prologator, i.e.,
  //
  //   P = (I - omega*D^{-1}*A)*T
  //
  // To compute D^{-1}*A*T we need all the rows T(j,:) for which
  // j is a column index for the nonzeros of A on this rank.
  // For all local column indices in A.diag, we will already
  // have the row of T on this rank, so we just need to gather
  // the offd colIds

  //Jacobi weight
  const dfloat omega = (4./3.)/A.rho;

  memory<hlong> recvRows(A.Ncols-A.NlocalCols);
  memory<int> sendCounts(size);
  memory<int> recvCounts(size, 0);
  memory<int> sendOffsets(size+1);
  memory<int> recvOffsets(size+1);

  memory<hlong> globalRowStarts(size+1);
  globalRowStarts[0]=0;
  T.comm.Allgather(T.rowOffsetU, globalRowStarts+1);

  //use the colMap of A to list the needed rows of T
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
  T.comm.Alltoallv(recvRows, recvCounts, recvOffsets,
                   sendRows, sendCounts, sendOffsets);

  //we now have a list of rows to send, count the nnz to send
  dlong nnzTotal=0;
  for (r=0;r<size;r++) {
    sendCounts[r] =0; //reset
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      const dlong i = static_cast<dlong>(sendRows[n]-T.rowOffsetL); //local row id
      sendCounts[r]+= T.diag.rowStarts[i+1]-T.diag.rowStarts[i]; //count entries in this row
      sendCounts[r]+= T.offd.rowStarts[i+1]-T.offd.rowStarts[i]; //count entries in this row
    }
    nnzTotal += sendCounts[r]; //tally the total
  }

  memory<nonZero_t> sendNonZeros(nnzTotal);

  nnzTotal=0; //reset
  for (r=0;r<size;r++) {
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      const dlong i = static_cast<dlong>(sendRows[n] - T.rowOffsetL); //local row id
      for (dlong jj=T.diag.rowStarts[i]; jj<T.diag.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = T.diag.cols[jj] + T.colOffsetL;
        sendNonZeros[nnzTotal].val = T.diag.vals[jj];
        nnzTotal++;
      }
      for (dlong jj=T.offd.rowStarts[i]; jj<T.offd.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = T.colMap[T.offd.cols[jj]];
        sendNonZeros[nnzTotal].val = T.offd.vals[jj];
        nnzTotal++;
      }
    }
  }

  A.comm.Alltoall(sendCounts, recvCounts);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }


  dlong Toffdnnz = recvOffsets[size]; //total nonzeros
  memory<nonZero_t> ToffdRows(Toffdnnz);

  T.comm.Alltoallv(sendNonZeros, sendCounts, sendOffsets,
                      ToffdRows, recvCounts, recvOffsets);

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
  memory<dlong> ToffdRowOffsets(A.Ncols-A.NlocalCols+1, 0);

  dlong id=0;
  for (dlong n=0;n<Toffdnnz;n++) {
    hlong row = ToffdRows[n].row;

    while(A.colMap[id+A.NlocalCols]!=row) id++;

    ToffdRowOffsets[id+1]++; //count entry in row
  }

  //cumulative sum
  for (dlong n=0;n<A.Ncols-A.NlocalCols;n++)
    ToffdRowOffsets[n+1] += ToffdRowOffsets[n];


  // The next step to compute D^{-1}*A*T is to multiply each entry A(i,j) by the
  // row T(j,:), store the all the results, sort them by row+col, and compress
  // the entries

  // Find how big the intermediate form is
  memory<dlong> rowStarts(A.Nrows+1, 0);
  memory<dlong> rowCounts(A.Nrows, 0);

  /*Count entries per row*/
  #pragma omp parallel for
  for(dlong i=0; i<A.Nrows; i++) {
    /*Start with entries for T*/
    rowStarts[i+1]+=T.diag.rowStarts[i+1]-T.diag.rowStarts[i] +
                    T.offd.rowStarts[i+1]-T.offd.rowStarts[i];

    /*Then add entries from A*T*/
    dlong Jstart = A.diag.rowStarts[i];
    dlong Jend   = A.diag.rowStarts[i+1];
    for(dlong jj=Jstart; jj<Jend; jj++){
      const dlong col = A.diag.cols[jj];
      rowStarts[i+1]+=T.diag.rowStarts[col+1]-T.diag.rowStarts[col] +
                      T.offd.rowStarts[col+1]-T.offd.rowStarts[col];
    }
    //non-local entries
    Jstart = A.offd.rowStarts[i];
    Jend   = A.offd.rowStarts[i+1];
    for (dlong jj=Jstart;jj<Jend;jj++) {
      const dlong col = A.offd.cols[jj]-A.NlocalCols;
      rowStarts[i+1]+= ToffdRowOffsets[col+1] - ToffdRowOffsets[col];
    }
  }

  /*Cumulative sum*/
  for(dlong i=1; i<A.Nrows+1; i++) {
    rowStarts[i] += rowStarts[i-1];
  }

  dlong NNZ = rowStarts[A.Nrows];

  memory<nonZero_t> Ptmp(NNZ);

  //count total number of nonzeros we find
  dlong nnz =0;

  // Fill the intermediate form of P
  // #pragma omp parallel for
  for (dlong i=0;i<A.Nrows;i++) {
    const dlong cStart = rowStarts[i];
    dlong& c = rowCounts[i];

    /*Start with P=T entries*/

    //local T entries
    dlong start = T.diag.rowStarts[i];
    dlong end   = T.diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      Ptmp[cStart+c].row = i + T.rowOffsetL;
      Ptmp[cStart+c].col = T.diag.cols[j] + T.colOffsetL; //global id
      Ptmp[cStart+c].val = T.diag.vals[j];
      c++;
    }
    //non-local T entries
    start = T.offd.rowStarts[i];
    end   = T.offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      Ptmp[cStart+c].row = i + T.rowOffsetL;
      Ptmp[cStart+c].col = T.colMap[T.offd.cols[j]];
      Ptmp[cStart+c].val = T.offd.vals[j];
      c++;
    }

    /*Then P -= omega*invD*A*T*/

    //local A entries
    start = A.diag.rowStarts[i];
    end   = A.diag.rowStarts[i+1];

    const dfloat invDi = 1.0/A.diagA[i];

    for (dlong j=start;j<end;j++) {
      const dlong col = A.diag.cols[j];
      const dfloat Aval = -omega*invDi*A.diag.vals[j];

      //local T entries
      dlong Tstart = T.diag.rowStarts[col];
      dlong Tend   = T.diag.rowStarts[col+1];
      for (dlong jj=Tstart;jj<Tend;jj++) {
        Ptmp[cStart+c].row = i + A.rowOffsetL;
        Ptmp[cStart+c].col = T.diag.cols[jj] + T.colOffsetL; //global id
        Ptmp[cStart+c].val = Aval*T.diag.vals[jj];
        c++;
      }
      //non-local T entries
      Tstart = T.offd.rowStarts[col];
      Tend   = T.offd.rowStarts[col+1];
      for (dlong jj=Tstart;jj<Tend;jj++) {
        Ptmp[cStart+c].row = i + A.rowOffsetL;
        Ptmp[cStart+c].col = T.colMap[T.offd.cols[jj]]; //global id
        Ptmp[cStart+c].val = Aval*T.offd.vals[jj];
        c++;
      }
    }
    //non-local A entries
    start = A.offd.rowStarts[i];
    end   = A.offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A.offd.cols[j]-A.NlocalCols;
      const dfloat Aval = -omega*invDi*A.offd.vals[j];

      // entries from recived rows of T
      dlong Tstart = ToffdRowOffsets[col];
      dlong Tend   = ToffdRowOffsets[col+1];
      for (dlong jj=Tstart;jj<Tend;jj++) {
        Ptmp[cStart+c].row = i + A.rowOffsetL;
        Ptmp[cStart+c].col = ToffdRows[jj].col; //global id
        Ptmp[cStart+c].val = Aval*ToffdRows[jj].val;
        c++;
      }
    }

    //sort entries in this row by col id
    std::sort(Ptmp.ptr()+cStart, Ptmp.ptr()+cStart+c,
              [](const nonZero_t& a, const nonZero_t& b) {
                return a.col < b.col;
              });

    /*Count how many actual nonzeros will be in this row*/
    dlong nnzRow=0;
    if (c>0) nnzRow++;
    for (dlong j=1;j<c;j++) {
      if ((Ptmp[cStart+j].col!=Ptmp[cStart+j-1].col)) nnzRow++;
    }

    nnz+=nnzRow; //Add to total
  }
  ToffdRowOffsets.free();
  ToffdRows.free();

  rowStarts.free();
  rowCounts.free();

  // cooP.nnz = nnz;
  memory<nonZero_t> entries(nnz);

  //compress nonzeros
  nnz = 0;
  if (NNZ) entries[nnz++] = Ptmp[0];
  for (dlong i=1;i<NNZ;i++) {
    if ((Ptmp[i].row!=Ptmp[i-1].row)||
        (Ptmp[i].col!=Ptmp[i-1].col)) {
      entries[nnz++] = Ptmp[i];
    } else {
      entries[nnz-1].val += Ptmp[i].val;
    }
  }
  //clean up
  Ptmp.free();

  //build P from coo matrix
  return parCSR(A.Nrows, T.NlocalCols,
                nnz, entries,
                A.platform, A.comm);
}

} //namespace paradogs

} //namespace libp

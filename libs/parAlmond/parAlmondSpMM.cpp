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

static int compareNonZeroByRow(const void *a, const void *b){
  parCOO::nonZero_t *pa = (parCOO::nonZero_t *) a;
  parCOO::nonZero_t *pb = (parCOO::nonZero_t *) b;

  if (pa->row < pb->row) return -1;
  if (pa->row > pb->row) return +1;

  if (pa->col < pb->col) return -1;
  if (pa->col > pb->col) return +1;

  return 0;
};

parCSR *SpMM(parCSR *A, parCSR *B){

  // MPI info
  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  // To compute C = A*B we need all the rows B(j,:) for which
  // j is a column index for the nonzeros of A on this rank.
  // For all local column indices in A->diag, we will already
  // have the row of B on this rank, so we just need to gather
  // the offd colIds

  hlong *recvRows = (hlong *) calloc(A->Ncols-A->NlocalCols, sizeof(hlong));
  int *sendCounts = (int*) calloc(size, sizeof(int));
  int *recvCounts = (int*) calloc(size, sizeof(int));
  int *sendOffsets = (int*) calloc(size+1, sizeof(int));
  int *recvOffsets = (int*) calloc(size+1, sizeof(int));

  //use the colMap of A to list the needed rows of B
  int r=0;
  for (dlong n=A->NlocalCols;n<A->Ncols;n++) {
    const hlong id = A->colMap[n];
    while (id>=B->globalRowStarts[r+1]) r++; //assumes the halo is sorted
    recvCounts[r]++;
    recvRows[n-A->NlocalCols] = id; //record the row to recv
  }

  //share the counts
  MPI_Alltoall(recvCounts, 1, MPI_INT,
               sendCounts, 1, MPI_INT, A->comm);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }

  int sendTotal = sendOffsets[size];
  hlong *sendRows = (hlong *) calloc(sendTotal, sizeof(hlong));

  //share the rowIds
  MPI_Alltoallv(recvRows, recvCounts, recvOffsets, MPI_HLONG,
                sendRows, sendCounts, sendOffsets, MPI_HLONG,
                B->comm);

  //we now have a list of rows to send, count the nnz to send
  dlong nnzTotal=0;
  for (r=0;r<size;r++) {
    sendCounts[r] =0; //reset
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      dlong i = (dlong) (sendRows[n]-B->globalRowStarts[rank]); //local row id
      sendCounts[r]+= B->diag.rowStarts[i+1]-B->diag.rowStarts[i]; //count entries in this row
      sendCounts[r]+= B->offd.rowStarts[i+1]-B->offd.rowStarts[i]; //count entries in this row
    }
    nnzTotal += sendCounts[r]; //tally the total
  }

  parCOO::nonZero_t *sendNonZeros = (parCOO::nonZero_t *) calloc(nnzTotal, sizeof(parCOO::nonZero_t));

  nnzTotal=0; //reset
  for (r=0;r<size;r++) {
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      dlong i = (dlong) (sendRows[n] - B->globalRowStarts[rank]); //local row id
      for (dlong jj=B->diag.rowStarts[i]; jj<B->diag.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = B->diag.cols[jj] + B->globalColStarts[rank];
        sendNonZeros[nnzTotal].val = B->diag.vals[jj];
        nnzTotal++;
      }
      for (dlong jj=B->offd.rowStarts[i]; jj<B->offd.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = B->colMap[B->offd.cols[jj]];
        sendNonZeros[nnzTotal].val = B->offd.vals[jj];
        nnzTotal++;
      }
    }
  }

  MPI_Alltoall(sendCounts, 1, MPI_INT,
               recvCounts, 1, MPI_INT, A->comm);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }


  dlong Boffdnnz = recvOffsets[size]; //total nonzeros
  parCOO::nonZero_t *BoffdRows = (parCOO::nonZero_t *)
                                 calloc(Boffdnnz, sizeof(parCOO::nonZero_t));

  MPI_Alltoallv(sendNonZeros, sendCounts, sendOffsets, MPI_NONZERO_T,
                BoffdRows, recvCounts, recvOffsets, MPI_NONZERO_T,
                B->comm);

  //clean up
  MPI_Barrier(B->comm);
  free(sendNonZeros);
  free(sendCounts);
  free(recvCounts);
  free(sendOffsets);
  free(recvOffsets);

  //we now have all the needed nonlocal rows (should also be sorted by row then col)

  //make an array of row offsets so we know how large each row is
  dlong *BoffdRowOffsets = (dlong *) calloc(A->Ncols-A->NlocalCols+1, sizeof(dlong));

  dlong id=0;
  for (dlong n=0;n<Boffdnnz;n++) {
    hlong row = BoffdRows[n].row;

    while(A->colMap[id+A->NlocalCols]!=row) id++;

    BoffdRowOffsets[id+1]++; //count entry in row
  }

  //cumulative sum
  for (dlong n=0;n<A->Ncols-A->NlocalCols;n++)
    BoffdRowOffsets[n+1] += BoffdRowOffsets[n];


  // The next step to compute C = A*B is to multiply each entry A(i,j) by the
  // row B(j,:), store the all the results, sort them by row+col, and compress
  // the entries

  // Find how big the intermediate form is
  nnzTotal = 0;
  for (dlong i=0;i<A->Nrows;i++) {
    //local entries
    dlong start = A->diag.rowStarts[i];
    dlong end   = A->diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->diag.cols[j];
      const int nnzBj =  B->diag.rowStarts[col+1]-B->diag.rowStarts[col]
                        +B->offd.rowStarts[col+1]-B->offd.rowStarts[col];
      nnzTotal += nnzBj;
    }
    //non-local entries
    start = A->offd.rowStarts[i];
    end   = A->offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->offd.cols[j]-A->NlocalCols;
      const int nnzBj = BoffdRowOffsets[col+1] - BoffdRowOffsets[col];
      nnzTotal += nnzBj;
    }
  }

  parCOO::nonZero_t *Ctmp = (parCOO::nonZero_t *)
                            calloc(nnzTotal, sizeof(parCOO::nonZero_t));

  // Fill the intermediate form of C
  dlong cnt = 0;
  for (dlong i=0;i<A->Nrows;i++) {
    //local A entries
    dlong start = A->diag.rowStarts[i];
    dlong end   = A->diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->diag.cols[j];
      const dfloat Aval = A->diag.vals[j];

      //local B entries
      dlong Bstart = B->diag.rowStarts[col];
      dlong Bend   = B->diag.rowStarts[col+1];
      for (dlong jj=Bstart;jj<Bend;jj++) {
        Ctmp[cnt].row = i + A->globalRowStarts[rank];
        Ctmp[cnt].col = B->diag.cols[jj]+B->globalColStarts[rank]; //global id
        Ctmp[cnt].val = Aval*B->diag.vals[jj];
        cnt++;
      }
      //non-local B entries
      Bstart = B->offd.rowStarts[col];
      Bend   = B->offd.rowStarts[col+1];
      for (dlong jj=Bstart;jj<Bend;jj++) {
        Ctmp[cnt].row = i + A->globalRowStarts[rank];
        Ctmp[cnt].col = B->colMap[B->offd.cols[jj]]; //global id
        Ctmp[cnt].val = Aval*B->offd.vals[jj];
        cnt++;
      }
    }
    //non-local A entries
    start = A->offd.rowStarts[i];
    end   = A->offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->offd.cols[j]-A->NlocalCols;
      const dfloat Aval = A->offd.vals[j];

      // entries from recived rows of B
      dlong Bstart = BoffdRowOffsets[col];
      dlong Bend   = BoffdRowOffsets[col+1];
      for (dlong jj=Bstart;jj<Bend;jj++) {
        Ctmp[cnt].row = i + A->globalRowStarts[rank];
        Ctmp[cnt].col = BoffdRows[jj].col; //global id
        Ctmp[cnt].val = Aval*BoffdRows[jj].val;
        cnt++;
      }
    }
  }
  free(BoffdRowOffsets);
  free(BoffdRows);


  //sort entries by the row and col
  qsort(Ctmp, nnzTotal, sizeof(parCOO::nonZero_t), compareNonZeroByRow);

  //count total number of nonzeros;
  dlong nnz =0;
  if (nnzTotal) nnz++;
  for (dlong i=1;i<nnzTotal;i++)
    if ((Ctmp[i].row!=Ctmp[i-1].row)||
        (Ctmp[i].col!=Ctmp[i-1].col)) nnz++;

  parCOO cooC(A->platform, A->comm);

  //copy global partition
  cooC.globalRowStarts = (hlong *) calloc(size+1,sizeof(hlong));
  cooC.globalColStarts = (hlong *) calloc(size+1,sizeof(hlong));
  memcpy(cooC.globalRowStarts, A->globalRowStarts, (size+1)*sizeof(hlong));
  memcpy(cooC.globalColStarts, B->globalColStarts, (size+1)*sizeof(hlong));

  cooC.nnz = nnz;
  cooC.entries = (parCOO::nonZero_t *) calloc(nnz,sizeof(parCOO::nonZero_t));

  //compress nonzeros
  nnz = 0;
  if (nnzTotal) cooC.entries[nnz++] = Ctmp[0];
  for (dlong i=1;i<nnzTotal;i++) {
    if ((Ctmp[i].row!=Ctmp[i-1].row)||
        (Ctmp[i].col!=Ctmp[i-1].col)) {
      cooC.entries[nnz++] = Ctmp[i];
    } else {
      cooC.entries[nnz-1].val += Ctmp[i].val;
    }
  }
  //clean up
  free(Ctmp);

  //build C from coo matrix
  return new parCSR(cooC);
}

} //namespace parAlmond
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

parCSR *smoothProlongator(parCSR *A, parCSR *T){

  // MPI info
  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  // This function computes a smoothed prologation operator
  // via a single weighted Jacobi iteration on the tentative
  // prologator, i.e.,
  //
  //   P = (I - omega*D^{-1}*A)*T
  //
  // To compute D^{-1}*A*T we need all the rows T(j,:) for which
  // j is a column index for the nonzeros of A on this rank.
  // For all local column indices in A->diag, we will already
  // have the row of T on this rank, so we just need to gather
  // the offd colIds

  //Jacobi weight
  const dfloat omega = (4./3.)/A->rho;

  hlong *recvRows = (hlong *) calloc(A->Ncols-A->NlocalCols, sizeof(hlong));
  int *sendCounts = (int*) calloc(size, sizeof(int));
  int *recvCounts = (int*) calloc(size, sizeof(int));
  int *sendOffsets = (int*) calloc(size+1, sizeof(int));
  int *recvOffsets = (int*) calloc(size+1, sizeof(int));

  //use the colMap of A to list the needed rows of T
  int r=0;
  for (dlong n=A->NlocalCols;n<A->Ncols;n++) {
    const hlong id = A->colMap[n];
    while (id>=T->globalRowStarts[r+1]) r++; //assumes the halo is sorted
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
                T->comm);

  //we now have a list of rows to send, count the nnz to send
  dlong nnzTotal=0;
  for (r=0;r<size;r++) {
    sendCounts[r] =0; //reset
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      dlong i = (dlong) (sendRows[n]-T->globalRowStarts[rank]); //local row id
      sendCounts[r]+= T->diag.rowStarts[i+1]-T->diag.rowStarts[i]; //count entries in this row
      sendCounts[r]+= T->offd.rowStarts[i+1]-T->offd.rowStarts[i]; //count entries in this row
    }
    nnzTotal += sendCounts[r]; //tally the total
  }

  parCOO::nonZero_t *sendNonZeros = (parCOO::nonZero_t *) calloc(nnzTotal, sizeof(parCOO::nonZero_t));

  nnzTotal=0; //reset
  for (r=0;r<size;r++) {
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      dlong i = (dlong) (sendRows[n] - T->globalRowStarts[rank]); //local row id
      for (dlong jj=T->diag.rowStarts[i]; jj<T->diag.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = T->diag.cols[jj] + T->globalColStarts[rank];
        sendNonZeros[nnzTotal].val = T->diag.vals[jj];
        nnzTotal++;
      }
      for (dlong jj=T->offd.rowStarts[i]; jj<T->offd.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = T->colMap[T->offd.cols[jj]];
        sendNonZeros[nnzTotal].val = T->offd.vals[jj];
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


  dlong Toffdnnz = recvOffsets[size]; //total nonzeros
  parCOO::nonZero_t *ToffdRows = (parCOO::nonZero_t *)
                                 calloc(Toffdnnz, sizeof(parCOO::nonZero_t));

  MPI_Alltoallv(sendNonZeros, sendCounts, sendOffsets, MPI_NONZERO_T,
                ToffdRows, recvCounts, recvOffsets, MPI_NONZERO_T,
                T->comm);

  //clean up
  MPI_Barrier(T->comm);
  free(sendNonZeros);
  free(sendCounts);
  free(recvCounts);
  free(sendOffsets);
  free(recvOffsets);

  //we now have all the needed nonlocal rows (should also be sorted by row then col)

  //make an array of row offsets so we know how large each row is
  dlong *ToffdRowOffsets = (dlong *) calloc(A->Ncols-A->NlocalCols+1, sizeof(dlong));

  dlong id=0;
  for (dlong n=0;n<Toffdnnz;n++) {
    hlong row = ToffdRows[n].row;

    while(A->colMap[id+A->NlocalCols]!=row) id++;

    ToffdRowOffsets[id+1]++; //count entry in row
  }

  //cumulative sum
  for (dlong n=0;n<A->Ncols-A->NlocalCols;n++)
    ToffdRowOffsets[n+1] += ToffdRowOffsets[n];


  // The next step to compute D^{-1}*A*T is to multiply each entry A(i,j) by the
  // row T(j,:), store the all the results, sort them by row+col, and compress
  // the entries

  // Find how big the intermediate form is
  nnzTotal = T->diag.nnz+T->offd.nnz; //start with T populated

  for (dlong i=0;i<A->Nrows;i++) {
    //local entries
    dlong start = A->diag.rowStarts[i];
    dlong end   = A->diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->diag.cols[j];
      const int nnzBj =  T->diag.rowStarts[col+1]-T->diag.rowStarts[col]
                        +T->offd.rowStarts[col+1]-T->offd.rowStarts[col];
      nnzTotal += nnzBj;
    }
    //non-local entries
    start = A->offd.rowStarts[i];
    end   = A->offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->offd.cols[j]-A->NlocalCols;
      const int nnzBj = ToffdRowOffsets[col+1] - ToffdRowOffsets[col];
      nnzTotal += nnzBj;
    }
  }

  parCOO::nonZero_t *Ptmp = (parCOO::nonZero_t *)
                            calloc(nnzTotal, sizeof(parCOO::nonZero_t));

  // Fill the intermediate form of P
  dlong cnt = 0;
  //First P = T
  for (dlong i=0;i<T->Nrows;i++) {
    //local T entries
    dlong start = T->diag.rowStarts[i];
    dlong end   = T->diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      Ptmp[cnt].row = i + T->globalRowStarts[rank];
      Ptmp[cnt].col = T->diag.cols[j]+T->globalColStarts[rank]; //global id
      Ptmp[cnt].val = T->diag.vals[j];
      cnt++;
    }
    //non-local T entries
    start = T->offd.rowStarts[i];
    end   = T->offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      Ptmp[cnt].row = i + T->globalRowStarts[rank];
      Ptmp[cnt].col = T->colMap[T->offd.cols[j]];
      Ptmp[cnt].val = T->offd.vals[j];
      cnt++;
    }
  }

  //Then P -= omega*invD*A*T
  for (dlong i=0;i<A->Nrows;i++) {
    //local A entries
    dlong start = A->diag.rowStarts[i];
    dlong end   = A->diag.rowStarts[i+1];

    const dfloat invDi = 1.0/A->diagA[i];

    for (dlong j=start;j<end;j++) {
      const dlong col = A->diag.cols[j];
      const dfloat Aval = -omega*invDi*A->diag.vals[j];

      //local T entries
      dlong Tstart = T->diag.rowStarts[col];
      dlong Tend   = T->diag.rowStarts[col+1];
      for (dlong jj=Tstart;jj<Tend;jj++) {
        Ptmp[cnt].row = i + A->globalRowStarts[rank];
        Ptmp[cnt].col = T->diag.cols[jj]+T->globalColStarts[rank]; //global id
        Ptmp[cnt].val = Aval*T->diag.vals[jj];
        cnt++;
      }
      //non-local T entries
      Tstart = T->offd.rowStarts[col];
      Tend   = T->offd.rowStarts[col+1];
      for (dlong jj=Tstart;jj<Tend;jj++) {
        Ptmp[cnt].row = i + A->globalRowStarts[rank];
        Ptmp[cnt].col = T->colMap[T->offd.cols[jj]]; //global id
        Ptmp[cnt].val = Aval*T->offd.vals[jj];
        cnt++;
      }
    }
    //non-local A entries
    start = A->offd.rowStarts[i];
    end   = A->offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->offd.cols[j]-A->NlocalCols;
      const dfloat Aval = -omega*invDi*A->offd.vals[j];

      // entries from recived rows of T
      dlong Tstart = ToffdRowOffsets[col];
      dlong Tend   = ToffdRowOffsets[col+1];
      for (dlong jj=Tstart;jj<Tend;jj++) {
        Ptmp[cnt].row = i + A->globalRowStarts[rank];
        Ptmp[cnt].col = ToffdRows[jj].col; //global id
        Ptmp[cnt].val = Aval*ToffdRows[jj].val;
        cnt++;
      }
    }
  }
  free(ToffdRowOffsets);
  free(ToffdRows);


  //sort entries by the row and col
  qsort(Ptmp, nnzTotal, sizeof(parCOO::nonZero_t), compareNonZeroByRow);

  //count total number of nonzeros;
  dlong nnz =0;
  if (nnzTotal) nnz++;
  for (dlong i=1;i<nnzTotal;i++)
    if ((Ptmp[i].row!=Ptmp[i-1].row)||
        (Ptmp[i].col!=Ptmp[i-1].col)) nnz++;

  parCOO cooP(A->platform, A->comm);

  //copy global partition
  cooP.globalRowStarts = (hlong *) calloc(size+1,sizeof(hlong));
  cooP.globalColStarts = (hlong *) calloc(size+1,sizeof(hlong));
  memcpy(cooP.globalRowStarts, A->globalRowStarts, (size+1)*sizeof(hlong));
  memcpy(cooP.globalColStarts, T->globalColStarts, (size+1)*sizeof(hlong));

  cooP.nnz = nnz;
  cooP.entries = (parCOO::nonZero_t *) calloc(nnz,sizeof(parCOO::nonZero_t));

  //compress nonzeros
  nnz = 0;
  if (nnzTotal) cooP.entries[nnz++] = Ptmp[0];
  for (dlong i=1;i<nnzTotal;i++) {
    if ((Ptmp[i].row!=Ptmp[i-1].row)||
        (Ptmp[i].col!=Ptmp[i-1].col)) {
      cooP.entries[nnz++] = Ptmp[i];
    } else {
      cooP.entries[nnz-1].val += Ptmp[i].val;
    }
  }
  //clean up
  free(Ptmp);

  //build P from coo matrix
  return new parCSR(cooP);
}

} //namespace parAlmond
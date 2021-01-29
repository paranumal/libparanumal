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

  dfloat TOL = 1e-8;
  
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

  int *recvRows = (int *) calloc(A->Ncols-A->NlocalCols, sizeof(int));
  int *sendCounts = (int*) calloc(size, sizeof(int));
  int *recvCounts = (int*) calloc(size, sizeof(int));
  int *sendOffsets = (int*) calloc(size+1, sizeof(int));
  int *recvOffsets = (int*) calloc(size+1, sizeof(int));

  //use the colMap of A to list the needed rows of B
  hlong r=0;
  for (hlong n=A->NlocalCols;n<A->Ncols;n++) {
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

  hlong sendTotal = sendOffsets[size];
  hlong *sendRows = (hlong *) calloc(sendTotal, sizeof(hlong));

  //share the rowIds
  MPI_Alltoallv(recvRows, recvCounts, recvOffsets, MPI_INT,
                sendRows, sendCounts, sendOffsets, MPI_INT,
                B->comm);

  //we now have a list of rows to send, count the nnz to send
  hlong nnzTotal=0;
  for (r=0;r<size;r++) {
    sendCounts[r] =0; //reset
    for (hlong n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      hlong i = (hlong) (sendRows[n]-B->globalRowStarts[rank]); //local row id
      sendCounts[r]+= B->diag.rowStarts[i+1]-B->diag.rowStarts[i]; //count entries in this row
      sendCounts[r]+= B->offd.rowStarts[i+1]-B->offd.rowStarts[i]; //count entries in this row
    }
    nnzTotal += sendCounts[r]; //tally the total
  }

  parCOO::nonZero_t *sendNonZeros = (parCOO::nonZero_t *) calloc(nnzTotal, sizeof(parCOO::nonZero_t));

  nnzTotal=0; //reset
  for (r=0;r<size;r++) {
    for (hlong n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      hlong i = (hlong) (sendRows[n] - B->globalRowStarts[rank]); //local row id
      for (hlong jj=B->diag.rowStarts[i]; jj<B->diag.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = B->diag.cols[jj] + B->globalColStarts[rank];
        sendNonZeros[nnzTotal].val = B->diag.vals[jj];
        nnzTotal++;
      }
      for (hlong jj=B->offd.rowStarts[i]; jj<B->offd.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = B->colMap[B->offd.cols[jj]];
        sendNonZeros[nnzTotal].val = B->offd.vals[jj];
        nnzTotal++;
      }
    }
  }

  MPI_Alltoall(sendCounts, 1, MPI_HLONG,
               recvCounts, 1, MPI_HLONG, A->comm);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }


  hlong Boffdnnz = recvOffsets[size]; //total nonzeros
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
  hlong *BoffdRowOffsets = (hlong *) calloc(A->Ncols-A->NlocalCols+1, sizeof(hlong));

  hlong id=0;
  for (hlong n=0;n<Boffdnnz;n++) {
    hlong row = BoffdRows[n].row;

    while(A->colMap[id+A->NlocalCols]!=row) id++;

    BoffdRowOffsets[id+1]++; //count entry in row
  }

  //cumulative sum
  for (hlong n=0;n<A->Ncols-A->NlocalCols;n++)
    BoffdRowOffsets[n+1] += BoffdRowOffsets[n];


  // The next step to compute C = A*B is to multiply each entry A(i,j) by the
  // row B(j,:), store the all the results, sort them by row+col, and compress
  // the entries

  // Find how big the intermediate form is
  hlong maxrowCount  = 0;
  nnzTotal = 0;
  for (hlong i=0;i<A->Nrows;i++) {
    hlong rowcnt = 0;
    //local entries
    hlong start = A->diag.rowStarts[i];
    hlong end   = A->diag.rowStarts[i+1];
    for (hlong j=start;j<end;j++) {
      const hlong col = A->diag.cols[j];
      const hlong nnzBj =  B->diag.rowStarts[col+1]-B->diag.rowStarts[col]
                        +B->offd.rowStarts[col+1]-B->offd.rowStarts[col];
      nnzTotal += nnzBj;
      rowcnt += nnzBj;
    }
    //non-local entries
    start = A->offd.rowStarts[i];
    end   = A->offd.rowStarts[i+1];
    for (hlong j=start;j<end;j++) {
      const hlong col = A->offd.cols[j]-A->NlocalCols;
      const hlong nnzBj = BoffdRowOffsets[col+1] - BoffdRowOffsets[col];
      nnzTotal += nnzBj;
      rowcnt += nnzBj;
    }
    maxrowCount = mymax(maxrowCount, rowcnt);
  }
  printf("maxrowCount = %lld\n", maxrowCount);
  
  if(nnzTotal<0) printf("OVERFLOW in parAlmondSpMM: A->nnz=%d, B->nnz=%d\n",
			A->diag.nnz, B->diag.nnz);

  parCOO::nonZero_t *Ctmp = (parCOO::nonZero_t *)
                            calloc(nnzTotal, sizeof(parCOO::nonZero_t));

  //  dlong maxrowCount = 100000;
  parCOO::nonZero_t *rowCtmp = (parCOO::nonZero_t *)
    calloc(maxrowCount, sizeof(parCOO::nonZero_t));
  
  printf("sizeof(hlong)=%d\n", sizeof(hlong));
  
  // Fill the intermediate form of C
  hlong cnt = 0;
  for (hlong i=0;i<A->Nrows;i++) {

    dlong rowcnt = 0;

    //local A entries
    hlong start = A->diag.rowStarts[i];
    hlong end   = A->diag.rowStarts[i+1];
    for (hlong j=start;j<end;j++) {
      const hlong col = A->diag.cols[j];
      const dfloat Aval = A->diag.vals[j];
      
      if(col<0) printf("FOUND negative col %lld (1)\n", col);
      
      //local B entries
      hlong Bstart = B->diag.rowStarts[col];
      hlong Bend   = B->diag.rowStarts[col+1];
      for (hlong jj=Bstart;jj<Bend;jj++) {
	dfloat val = Aval*B->diag.vals[jj];
	if(fabs(val)>TOL){
	  rowCtmp[rowcnt].row = i + A->globalRowStarts[rank];
	  rowCtmp[rowcnt].col = B->diag.cols[jj]+B->globalColStarts[rank]; //global id
	  rowCtmp[rowcnt].val = val;
	  rowcnt++;
	}
      }
      //non-local B entries
      Bstart = B->offd.rowStarts[col];
      Bend   = B->offd.rowStarts[col+1];
      for (hlong jj=Bstart;jj<Bend;jj++) {
	dfloat val = Aval*B->offd.vals[jj];
	if(fabs(val)>TOL){
	  rowCtmp[rowcnt].row = i + A->globalRowStarts[rank];
	  rowCtmp[rowcnt].col = B->colMap[B->offd.cols[jj]]; //global id
	  rowCtmp[rowcnt].val = val;
	  rowcnt++;
	}
      }
    }

    if(rowcnt>maxrowCount) printf("WARNING rowcnt=%d, maxrowcnt=%d\n", rowcnt, maxrowCount);
    
    //non-local A entries
    start = A->offd.rowStarts[i];
    end   = A->offd.rowStarts[i+1];
    for (hlong j=start;j<end;j++) {
      const hlong col = A->offd.cols[j]-A->NlocalCols;
      const dfloat Aval = A->offd.vals[j];

      if(col<0) printf("FOUND negative col %lld (2)\n", col);

      // entries from recived rows of B
      hlong Bstart = BoffdRowOffsets[col];
      hlong Bend   = BoffdRowOffsets[col+1];
      for (hlong jj=Bstart;jj<Bend;jj++) {
	dfloat val = Aval*BoffdRows[jj].val;
	if(fabs(val)>TOL){
	  rowCtmp[rowcnt].row = i + A->globalRowStarts[rank];
	  rowCtmp[rowcnt].col = BoffdRows[jj].col; //global id
	  rowCtmp[rowcnt].val = val;
	  rowcnt++;
	}
      }

      if(rowcnt>maxrowCount) printf("WARNING rowcnt=%d, maxrowcnt=%d\n", rowcnt, maxrowCount);
      
    }

    //sort entries by the row and col
    std::sort(rowCtmp, rowCtmp+rowcnt,
	      [](const parCOO::nonZero_t& a, const parCOO::nonZero_t& b) {
		if (a.row < b.row) return true;
		if (a.row > b.row) return false;
		
		return a.col < b.col;
	      });

    //count total number of nonzeros;
    if (rowcnt) Ctmp[cnt++] = rowCtmp[0];
    for (hlong j=1;j<rowcnt;j++) {
      if ((rowCtmp[j].row!=rowCtmp[j-1].row)||
	  (rowCtmp[j].col!=rowCtmp[j-1].col)) {
	Ctmp[cnt++] = rowCtmp[j];
      } else {
	Ctmp[cnt-1].val += rowCtmp[j].val;
      }
    }
  }

  free(BoffdRowOffsets);
  free(BoffdRows);
  free(rowCtmp);
  
  printf("Est nnzTotal=%lld, cnt=%lld\n", nnzTotal, cnt);
  
  nnzTotal = cnt;
  dlong nnz = cnt;
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
  for (hlong i=1;i<nnzTotal;i++) {
    if ((Ctmp[i].row!=Ctmp[i-1].row)||
        (Ctmp[i].col!=Ctmp[i-1].col)) {
      //increment non-zero count if current entry is above tolerance
      if (fabs(cooC.entries[nnz-1].val)>parAlmond::dropTolerance) nnz++;

      cooC.entries[nnz-1] = Ctmp[i];
    } else {
      cooC.entries[nnz-1].val += Ctmp[i].val;
    }
  }
  //clean up
  free(Ctmp);

  //update nnz count
  cooC.nnz = nnz;
  cooC.entries = (parCOO::nonZero_t *)
                    realloc(cooC.entries,
                            cooC.nnz*sizeof(parCOO::nonZero_t));

  //build C from coo matrix
  return new parCSR(cooC);
}

} //namespace parAlmond

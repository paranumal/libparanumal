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

typedef struct {

  hlong row;
  hlong col;
  dfloat val;

} nonzero_t;

static int compareNonZeroByRow(const void *a, const void *b){
  nonzero_t *pa = (nonzero_t *) a;
  nonzero_t *pb = (nonzero_t *) b;

  if (pa->row < pb->row) return -1;
  if (pa->row > pb->row) return +1;

  if (pa->col < pb->col) return -1;
  if (pa->col > pb->col) return +1;

  return 0;
};

parCSR *transpose(parCSR *A){

  // MPI info
  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  hlong *globalRowStarts = A->globalRowStarts;
  hlong *globalColStarts = A->globalColStarts;

  dlong Nrows = (dlong) (globalColStarts[rank+1]-globalColStarts[rank]);
  dlong Ncols = (dlong) (globalRowStarts[rank+1]-globalRowStarts[rank]);

  parCSR *At = new parCSR(Nrows, Ncols, A->platform, A->comm);

  At->globalRowStarts = globalColStarts;
  At->globalColStarts = globalRowStarts;

  At->diag.nnz = A->diag.nnz; //local entries remain local
  At->diag.rowStarts = (dlong *) calloc(At->Nrows+1, sizeof(dlong));

  //start with local entries
  nonzero_t *diagNZs = (nonzero_t *) calloc(A->diag.nnz, sizeof(nonzero_t));

  // Fill nonzero list and count the num of nonzeros per row for transpose
  for(dlong i=0; i<A->Nrows; i++){
    const dlong Jstart = A->diag.rowStarts[i];
    const dlong Jend   = A->diag.rowStarts[i+1];

    for(dlong jj=Jstart; jj<Jend; jj++){
      const dlong row = A->diag.cols[jj];

      diagNZs[jj].row = row;
      diagNZs[jj].col = i;
      diagNZs[jj].val = A->diag.vals[jj];

      At->diag.rowStarts[row+1]++;
    }
  }

  // cumulative sum for rows
  for(dlong i=1; i<=At->Nrows; i++)
    At->diag.rowStarts[i] += At->diag.rowStarts[i-1];

  //sort by row
  qsort(diagNZs, A->diag.nnz, sizeof(nonzero_t), compareNonZeroByRow);

  At->diag.cols = (dlong *)  calloc(At->diag.nnz, sizeof(dlong));
  At->diag.vals = (dfloat *) calloc(At->diag.nnz, sizeof(dfloat));

  for(dlong i=0; i<At->Nrows; i++){
    const dlong Jstart = At->diag.rowStarts[i];
    const dlong Jend   = At->diag.rowStarts[i+1];

    for(dlong jj=Jstart; jj<Jend; jj++){
      At->diag.cols[jj] = diagNZs[jj].col;
      At->diag.vals[jj] = diagNZs[jj].val;
    }
  }
  free(diagNZs);

  // Make the MPI_NONZERO_T data type
  nonzero_t NZ;
  MPI_Datatype MPI_NONZERO_T;
  MPI_Datatype dtype[3] = {MPI_HLONG, MPI_HLONG, MPI_DFLOAT};
  int blength[3] = {1, 1, 1};
  MPI_Aint addr[3], displ[3];
  MPI_Get_address ( &(NZ.row), addr+0);
  MPI_Get_address ( &(NZ.col), addr+1);
  MPI_Get_address ( &(NZ.val), addr+2);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  MPI_Type_create_struct (3, blength, displ, dtype, &MPI_NONZERO_T);
  MPI_Type_commit (&MPI_NONZERO_T);

  nonzero_t *sendNonZeros = (nonzero_t *) calloc(A->offd.nnz, sizeof(nonzero_t));

  // copy data from nonlocal entries into send buffer
  for(dlong i=0;i<A->offd.nzRows;++i){
    const hlong row = A->offd.rows[i] + globalRowStarts[rank]; //global ids
    for (dlong j=A->offd.rowStarts[i];j<A->offd.rowStarts[i+1];j++) {
      const hlong col =  A->colMap[A->offd.cols[j]]; //global ids
      sendNonZeros[j].row = col;
      sendNonZeros[j].col = row;
      sendNonZeros[j].val = A->offd.vals[j];
    }
  }

  //sort by destination row
  qsort(sendNonZeros, A->offd.nnz, sizeof(nonzero_t), compareNonZeroByRow);

  //count number of non-zeros we're sending
  int *sendCounts = (int*) calloc(size, sizeof(int));
  int *recvCounts = (int*) calloc(size, sizeof(int));
  int *sendOffsets = (int*) calloc(size+1, sizeof(int));
  int *recvOffsets = (int*) calloc(size+1, sizeof(int));

  int r=0;
  for (dlong n=0;n<A->offd.nnz;n++) {
    dlong row = sendNonZeros[n].row;
    while(row>=globalColStarts[r+1]) r++;
    sendCounts[r]++;
  }

  MPI_Alltoall(sendCounts, 1, MPI_INT,
               recvCounts, 1, MPI_INT, A->comm);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }
  At->offd.nnz = recvOffsets[size]; //total nonzeros

  nonzero_t *recvNonZeros = (nonzero_t *) calloc(At->offd.nnz, sizeof(nonzero_t));

  MPI_Alltoallv(sendNonZeros, sendCounts, sendOffsets, MPI_NONZERO_T,
                recvNonZeros, recvCounts, recvOffsets, MPI_NONZERO_T,
                A->comm);

  //clean up
  MPI_Barrier(A->comm);
  free(sendNonZeros);
  free(sendCounts);
  free(recvCounts);
  free(sendOffsets);
  free(recvOffsets);

  //sort by row
  qsort(recvNonZeros, At->offd.nnz, sizeof(nonzero_t), compareNonZeroByRow);

  const hlong globalRowOffset = At->globalRowStarts[rank];

  // count how many rows are shared
  At->offd.rowStarts = (dlong *) calloc(At->Nrows+1, sizeof(dlong));
  for (dlong n=0;n<At->offd.nnz;n++) {
    const dlong row = (dlong) (recvNonZeros[n].row - globalRowOffset);
    At->offd.rowStarts[row+1]++;
  }

  At->offd.nzRows=0;
  for(dlong i=0; i<At->Nrows; i++)
    if (At->offd.rowStarts[i+1]>0) At->offd.nzRows++;

  At->offd.rows       = (dlong *) calloc(At->offd.nzRows, sizeof(dlong));
  At->offd.mRowStarts = (dlong *) calloc(At->offd.nzRows+1, sizeof(dlong));

  // cumulative sum
  dlong cnt=0;
  for(dlong i=0; i<At->Nrows; i++) {
    if (At->offd.rowStarts[i+1]>0) {
      At->offd.rows[cnt] = i; //record row id
      At->offd.mRowStarts[cnt+1] = At->offd.mRowStarts[cnt] + At->offd.rowStarts[i+1];
      cnt++;
    }
    At->offd.rowStarts[i+1] += At->offd.rowStarts[i];
  }

  //Halo setup
  hlong *colIds = (hlong *) malloc(At->offd.nnz*sizeof(hlong));
  for (dlong n=0;n<At->offd.nnz;n++) {
    colIds[n] = recvNonZeros[n].col;
  }
  At->haloSetup(colIds); //setup halo, and transform colIds to a local indexing

  //fill the CSR matrix
  At->offd.cols = (dlong *)  calloc(At->offd.nnz, sizeof(dlong));
  At->offd.vals = (dfloat *) calloc(At->offd.nnz, sizeof(dfloat));

  for (dlong n=0;n<At->offd.nnz;n++) {
    At->offd.cols[n] = colIds[n];
    At->offd.vals[n] = recvNonZeros[n].val;
  }

  MPI_Barrier(A->comm);
  free(recvNonZeros);
  free(colIds);

  return At;
}

} //namespace parAlmond
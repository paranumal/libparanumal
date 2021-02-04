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
#include "parAlmond/parAlmondCoarseSolver.hpp"
#include "parAlmond/parAlmondKernels.hpp"

namespace parAlmond {

void oasSolver_t::solve(occa::memory& o_rhs, occa::memory& o_x) {

  A->halo->ExchangeStart(o_rhs, 1, ogs_pfloat);

  //queue local part of gemv
  const pfloat one=1.0;
  const pfloat zero=0.0;
  if (N)
    dGEMVKernel(N,diagTotal,one,o_diagInvAT,o_rhs, zero, o_x);

  A->halo->ExchangeFinish(o_rhs, 1, ogs_pfloat);

  //queue offd part of gemv
  if(offdTotal && N)
    dGEMVKernel(N,offdTotal, one, o_offdInvAT,
                o_rhs+diagTotal*sizeof(pfloat), one, o_x);

  A->halo->Combine(o_x, 1, ogs_pfloat);
}


int oasSolver_t::getTargetSize() {
  MPI_Comm_size(comm, &size);
  return 1000*size;
}

void oasSolver_t::setup(parCSR *_A, bool nullSpace,
                        pfloat *nullVector, pfloat nullSpacePenalty) {

  A = _A;

  comm = A->comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  N = (int) A->Ncols;
  Nrows = A->Nrows;
  Ncols = A->Ncols;

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE")))
  //   {printf("Setting up coarse solver...");fflush(stdout);}

  //make an overlapping patch matrix by collecting the rows
  // corresponding the offd columns

  //need to find where to send local rows
  hlong *recvRows = (hlong *) calloc(A->Ncols-A->Nrows, sizeof(hlong));

  int *sendCounts = (int*) calloc(size, sizeof(int));
  int *recvCounts = (int*) calloc(size, sizeof(int));
  int *sendOffsets = (int*) calloc(size+1, sizeof(int));
  int *recvOffsets = (int*) calloc(size+1, sizeof(int));

  //use the colMap to fill the recv sizes
  int r=0;
  for (int n=A->Nrows;n<A->Ncols;n++) {
    hlong id = A->colMap[n];
    while (id>=A->globalRowStarts[r+1]) r++; //assumes the halo is sorted
    recvCounts[r]++;
    recvRows[n-A->Nrows] = id; //record the row to recv
  }

  //share the counts
  MPI_Alltoall(recvCounts, 1, MPI_INT,
               sendCounts, 1, MPI_INT, comm);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }

  int sendTotal = sendOffsets[size];
  hlong *sendRows = (hlong *) calloc(sendTotal, sizeof(hlong));

  //share the rowIds
  MPI_Alltoallv(recvRows, recvCounts, recvOffsets, MPI_HLONG,
                sendRows, sendCounts, sendOffsets, MPI_HLONG,
                comm);

  //we now have a list of rows to send, count the nnz to send
  dlong nnzTotal=0;
  for (r=0;r<size;r++) {
    sendCounts[r] =0; //reset
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      dlong i = (dlong) (sendRows[n]-A->globalRowStarts[rank]); //local row id
      sendCounts[r]+= A->diag.rowStarts[i+1]-A->diag.rowStarts[i]; //count entries in this row
      sendCounts[r]+= A->offd.rowStarts[i+1]-A->offd.rowStarts[i]; //count entries in this row
    }
    nnzTotal += sendCounts[r]; //tally the total
  }

  parCOO::nonZero_t *sendNonZeros = (parCOO::nonZero_t *) calloc(nnzTotal, sizeof(parCOO::nonZero_t));

  nnzTotal=0; //reset
  for (r=0;r<size;r++) {
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      dlong i = (dlong) (sendRows[n] - A->globalRowStarts[rank]); //local row id
      for (dlong jj=A->diag.rowStarts[i]; jj<A->diag.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = A->diag.cols[jj] + A->globalRowStarts[rank];
        sendNonZeros[nnzTotal].val = A->diag.vals[jj];
        nnzTotal++;
      }
      for (dlong jj=A->offd.rowStarts[i]; jj<A->offd.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = A->colMap[A->offd.cols[jj]];
        sendNonZeros[nnzTotal].val = A->offd.vals[jj];
        nnzTotal++;
      }
    }
  }

  MPI_Alltoall(sendCounts, 1, MPI_INT,
               recvCounts, 1, MPI_INT, comm);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }

  nnzTotal = recvOffsets[size]; //total nonzeros

  parCOO::nonZero_t *recvNonZeros = (parCOO::nonZero_t *) calloc(nnzTotal, sizeof(parCOO::nonZero_t));

  MPI_Alltoallv(sendNonZeros, sendCounts, sendOffsets, MPI_NONZERO_T,
                recvNonZeros, recvCounts, recvOffsets, MPI_NONZERO_T,
                comm);

  //clean up
  MPI_Barrier(comm);
  free(sendNonZeros);
  free(sendCounts);
  free(recvCounts);
  free(sendOffsets);
  free(recvOffsets);

  //we now have all the nonlocal rows (should also be sorted)

  //first re-index the column indices
  dlong id=A->Nrows;
  for (dlong n=0;n<nnzTotal;n++) {
    const hlong row = recvNonZeros[n].row;

    while(A->colMap[id]!=row) id++; //shift along list of recieved columns

    recvNonZeros[n].row = id; //overwrite with new local row id

    //now check the column index
    hlong col = recvNonZeros[n].col;
    if (col >= A->globalRowStarts[rank] && col < A->globalRowStarts[rank+1]) {//local column
      recvNonZeros[n].col = col - A->globalRowStarts[rank];//overwrite with local col id
    } else {
      int flag = 0;
      for (dlong jj=A->Nrows;jj<A->Ncols;jj++) { //look for the right id in the halo
        if (A->colMap[jj]==col) {
          recvNonZeros[n].col = jj;//overwrite with local col id
          flag = 1;
          break;
        }
      }
      if (flag==0) recvNonZeros[n].col = -1; //ignore this entry as its not in our patch the patch
    }
  }

  //assemble the full matrix
  pfloat *coarseA = (pfloat *) calloc(N*N,sizeof(pfloat));
  for (int n=0;n<A->Nrows;n++) {
    const int start = (int) A->diag.rowStarts[n];
    const int end   = (int) A->diag.rowStarts[n+1];
    for (int m=start;m<end;m++) {
      int col = (int) A->diag.cols[m];
      coarseA[n*N+col] = A->diag.vals[m];
    }
  }

  for (int n=0;n<A->offd.nzRows;n++) {
    const int row   = (int) A->offd.rows[n];
    const int start = (int) A->offd.mRowStarts[n];
    const int end   = (int) A->offd.mRowStarts[n+1];
    for (int m=start;m<end;m++) {
      int col = (int) A->offd.cols[m];
      coarseA[row*N+col] = A->offd.vals[m];
    }
  }

  for (int i=0;i<nnzTotal;i++) {
    int n = recvNonZeros[i].row;
    int m = recvNonZeros[i].col;
    if (m>=0)
      coarseA[n*N+m] = recvNonZeros[i].val;
  }

  if (nullSpace) { //A is dense due to nullspace augmentation
    //copy fine nullvector and populate halo
    pfloat *null = (pfloat *) malloc(A->Ncols*sizeof(pfloat));
    memcpy(null, nullVector, A->Nrows*sizeof(pfloat));
    A->halo->Exchange(null, 1, ogs_pfloat);

    for (int n=0;n<N;n++) {
      for (int m=0;m<N;m++) {
        coarseA[n*N+m] += nullSpacePenalty*null[n]*null[m];
      }
    }

    free(null);
  }

  MPI_Barrier(comm);
  free(recvNonZeros);

  matrixInverse(N, coarseA);

  //determine the overlap weighting
  pfloat *weight = (pfloat *) malloc(N*sizeof(pfloat));
  for (int n=0;n<N;n++) weight[n] = 1.0;

  A->halo->Combine(weight, 1, ogs_pfloat);

  for (int n=0;n<N;n++) {
    for (int m=0;m<N;m++) {
      coarseA[n*N+m] *= 1.0/sqrt(weight[n]*weight[m]);
    }
  }
  free(weight);

  //determine size of offd piece
  diagTotal = A->Nrows;
  offdTotal = A->Ncols - A->Nrows;

  //diag piece of invA
  diagInvAT = (pfloat *) calloc(N*diagTotal,sizeof(pfloat));
  for (int n=0;n<N;n++) {
    for (int m=0;m<diagTotal;m++) {
      diagInvAT[n+m*N] = coarseA[n*N+m];
    }
  }

  //offd piece of invA
  offdInvAT = (pfloat *) calloc(N*offdTotal,sizeof(pfloat));
  for (int n=0;n<N;n++) {
    for (int m=0;m<offdTotal;m++) {
      offdInvAT[n+m*N] = coarseA[n*N + m+diagTotal];
    }
  }

  o_diagInvAT = platform.malloc(N*diagTotal*sizeof(pfloat), diagInvAT);
  o_offdInvAT = platform.malloc(N*offdTotal*sizeof(pfloat), offdInvAT);

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE"))) printf("done.\n");
}

void oasSolver_t::syncToDevice() {}

void oasSolver_t::Report(int lev) {

  hlong hNrows = (hlong) N;

  int active = (N>0) ? 1:0;
  int totalActive=0;
  MPI_Allreduce(&active, &totalActive, 1, MPI_INT, MPI_SUM, comm);

  dlong minNrows=0, maxNrows=0;
  hlong totalNrows=0;
  pfloat avgNrows;
  MPI_Allreduce(&N, &maxNrows, 1, MPI_DLONG, MPI_MAX, comm);
  MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, comm);
  avgNrows = (pfloat) totalNrows/totalActive;

  if (N==0) N=maxNrows; //set this so it's ignored for the global min
  MPI_Allreduce(&N, &minNrows, 1, MPI_DLONG, MPI_MIN, comm);

  long long int nnz;
  nnz = A->diag.nnz+A->offd.nnz;

  long long int minNnz=0, maxNnz=0, totalNnz=0;
  MPI_Allreduce(&nnz, &maxNnz,   1, MPI_LONG_LONG_INT, MPI_MAX, comm);
  MPI_Allreduce(&nnz, &totalNnz, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);

  if (nnz==0) nnz = maxNnz; //set this so it's ignored for the global min
  MPI_Allreduce(&nnz, &minNnz, 1, MPI_LONG_LONG_INT, MPI_MIN, comm);

  pfloat nnzPerRow = (Nrows==0) ? 0 : (pfloat) nnz/Nrows;
  pfloat minNnzPerRow=0, maxNnzPerRow=0, avgNnzPerRow=0;
  MPI_Allreduce(&nnzPerRow, &maxNnzPerRow, 1, MPI_PFLOAT, MPI_MAX, comm);
  MPI_Allreduce(&nnzPerRow, &avgNnzPerRow, 1, MPI_PFLOAT, MPI_SUM, comm);
  avgNnzPerRow /= totalActive;

  if (Nrows==0) nnzPerRow = maxNnzPerRow;
  MPI_Allreduce(&nnzPerRow, &minNnzPerRow, 1, MPI_PFLOAT, MPI_MIN, comm);

  std::string name = "OAS             ";

  if (rank==0){
    printf(" %3d  |  parAlmond |  %12lld  |  %12d  | %13d   |   %s|\n", lev, (long long int)totalNrows, minNrows, (int)minNnzPerRow, name.c_str());
    printf("      |            |                |  %12d  | %13d   |                   |\n", maxNrows, (int)maxNnzPerRow);
    printf("      |            |                |  %12d  | %13d   |                   |\n", (int)avgNrows, (int)avgNnzPerRow);
  }
}

oasSolver_t::~oasSolver_t() {
  if (diagInvAT) free(diagInvAT);
  if (offdInvAT) free(offdInvAT);
}

} //namespace parAlmond

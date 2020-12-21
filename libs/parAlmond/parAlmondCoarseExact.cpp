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

//link in the data stream from ogs
namespace ogs {
  extern occa::stream dataStream;
}

namespace parAlmond {

void exactSolver_t::solve(occa::memory& o_rhs, occa::memory& o_x) {

  occa::stream currentStream = platform.device.getStream();

  //queue transfering coarse vector to host for Allgather
  if(N) {
    platform.device.finish();
    platform.device.setStream(ogs::dataStream);
    o_rhs.copyTo(diagRhs, N*sizeof(dfloat), 0, "async: true");
    platform.device.setStream(currentStream);
  }

  //queue local part of gemv
  const dfloat one=1.0;
  const dfloat zero=0.0;
  if (N)
    dGEMVKernel(N,N,one,o_diagInvAT,o_rhs, zero, o_x);

  if(offdTotal) {
    //wait for data to arrive on host
    platform.device.setStream(ogs::dataStream);
    platform.device.finish();


    //gather the offd rhs entries
    MPI_Alltoallv(diagRhs,   sendCounts,   sendOffsets, MPI_DFLOAT,
                  offdRhs, coarseCounts, coarseOffsets, MPI_DFLOAT, comm);

    //queue transfering coarse vector to device
    o_offdRhs.copyFrom(offdRhs, offdTotal*sizeof(dfloat), 0, "async: true");
    platform.device.finish(); //wait for transfer to complete

    platform.device.setStream(currentStream);

    //queue offd part of gemv
    if (N)
      dGEMVKernel(N,offdTotal, one, o_offdInvAT,o_offdRhs, one, o_x);
  }
}


int exactSolver_t::getTargetSize() {
  return 1000;
}

void exactSolver_t::setup(parCSR *_A, bool nullSpace,
                           dfloat *nullVector, dfloat nullSpacePenalty) {

  A = _A;

  comm = A->comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  //copy the global coarse partition as ints
  coarseOffsets = (int* ) calloc(size+1,sizeof(int));
  for (int r=0;r<size+1;r++) coarseOffsets[r] = (int) A->globalRowStarts[r];

  coarseTotal   = coarseOffsets[size];
  coarseOffset  = coarseOffsets[rank];

  N = (int) A->Nrows;
  Nrows = A->Nrows;
  Ncols = A->Ncols;

  coarseCounts = (int*) calloc(size,sizeof(int));

  int sendNNZ = (int) (A->diag.nnz+A->offd.nnz);

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE")))
  //   {printf("Setting up coarse solver...");fflush(stdout);}

  parCOO::nonZero_t *sendNonZeros = (parCOO::nonZero_t *) calloc(sendNNZ, sizeof(parCOO::nonZero_t));

  //populate matrix
  int cnt = 0;
  for (int n=0;n<N;n++) {
    const int start = (int) A->diag.rowStarts[n];
    const int end   = (int) A->diag.rowStarts[n+1];
    for (int m=start;m<end;m++) {
      sendNonZeros[cnt].row = n + coarseOffset;
      sendNonZeros[cnt].col = A->diag.cols[m] + coarseOffset;
      sendNonZeros[cnt].val = A->diag.vals[m];
      cnt++;
    }
  }

  for (int n=0;n<A->offd.nzRows;n++) {
    const int row   = (int) A->offd.rows[n];
    const int start = (int) A->offd.mRowStarts[n];
    const int end   = (int) A->offd.mRowStarts[n+1];
    for (int m=start;m<end;m++) {
      sendNonZeros[cnt].row = row + coarseOffset;
      sendNonZeros[cnt].col = A->colMap[A->offd.cols[m]];
      sendNonZeros[cnt].val = A->offd.vals[m];
      cnt++;
    }
  }

  //get the nonzero counts from all ranks
  int *recvNNZ    = (int*) calloc(size,sizeof(int));
  int *NNZoffsets = (int*) calloc(size+1,sizeof(int));
  MPI_Allgather(&sendNNZ, 1, MPI_INT,
                 recvNNZ, 1, MPI_INT, comm);

  int totalNNZ = 0;
  for (int r=0;r<size;r++) {
    totalNNZ += recvNNZ[r];
    NNZoffsets[r+1] = NNZoffsets[r] + recvNNZ[r];
  }

  parCOO::nonZero_t *recvNonZeros = (parCOO::nonZero_t *) calloc(totalNNZ, sizeof(parCOO::nonZero_t));

  MPI_Allgatherv(sendNonZeros, sendNNZ,             MPI_NONZERO_T,
                 recvNonZeros, recvNNZ, NNZoffsets, MPI_NONZERO_T, comm);

  //gather null vector
  dfloat *nullTotal = (dfloat*) calloc(coarseTotal,sizeof(dfloat));

  for (int r=0;r<size;r++)
    coarseCounts[r] = coarseOffsets[r+1]-coarseOffsets[r];

  MPI_Allgatherv(nullVector,            N,                MPI_DFLOAT,
                  nullTotal, coarseCounts, coarseOffsets, MPI_DFLOAT,
                 comm);

  //clean up
  MPI_Barrier(comm);
  free(sendNonZeros);
  free(NNZoffsets);
  free(recvNNZ);

  //assemble the full matrix
  dfloat *coarseA = (dfloat *) calloc(coarseTotal*coarseTotal,sizeof(dfloat));
  for (int i=0;i<totalNNZ;i++) {
    int n = recvNonZeros[i].row;
    int m = recvNonZeros[i].col;
    coarseA[n*coarseTotal+m] = recvNonZeros[i].val;
  }

  if (nullSpace) { //A is dense due to nullspace augmentation
    for (int n=0;n<coarseTotal;n++) {
      for (int m=0;m<coarseTotal;m++) {
        coarseA[n*coarseTotal+m] += nullSpacePenalty*nullTotal[n]*nullTotal[m];
      }
    }
  }

  free(recvNonZeros);
  free(nullTotal);

  matrixInverse(coarseTotal, coarseA);

  //determine size of offd piece
  offdTotal = coarseTotal - N;

  //shift offsets for MPI_AllGatherv of offd pieces
  for (int r=rank+1;r<=size;r++)
    coarseOffsets[r]-= N;

  //dont copy the local piece in MPI_AllGatherv
  coarseCounts[rank]=0;

  //counts for all-to-all
  sendCounts = (int* ) calloc(size,sizeof(int));
  sendOffsets = (int* ) calloc(size,sizeof(int));
  for (int r=0;r<size;r++) {
    sendCounts[r] = N;
    sendOffsets[r] = 0;
  }
  sendCounts[rank] = 0;

  //diag piece of invA
  diagInvAT = (dfloat *) calloc(N*N,sizeof(dfloat));
  for (int n=0;n<N;n++) {
    for (int m=0;m<N;m++) {
      diagInvAT[n+m*N] = coarseA[(n+coarseOffset)*coarseTotal+(m+coarseOffset)];
    }
  }

  //offd piece of invA
  offdInvAT = (dfloat *) calloc(N*offdTotal,sizeof(dfloat));
  for (int n=0;n<N;n++) {
    for (int m=0;m<coarseOffset;m++) {
      offdInvAT[n+m*N] = coarseA[(n+coarseOffset)*coarseTotal+m];
    }
    for (int m=coarseOffset+N;m<coarseTotal;m++) {
      offdInvAT[n+(m-N)*N] = coarseA[(n+coarseOffset)*coarseTotal+m];
    }
  }

  o_diagInvAT = platform.malloc(N*N*sizeof(dfloat), diagInvAT);
  o_offdInvAT = platform.malloc(N*offdTotal*sizeof(dfloat), offdInvAT);

  diagRhs = (dfloat*) calloc(N,sizeof(dfloat));
  offdRhs = (dfloat*) calloc(offdTotal,sizeof(dfloat));

  o_offdRhs = platform.malloc(offdTotal*sizeof(dfloat));

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE"))) printf("done.\n");
}

void exactSolver_t::syncToDevice() {}

void exactSolver_t::Report(int lev) {

  hlong hNrows = (hlong) N;

  int active = (N>0) ? 1:0;
  int totalActive=0;
  MPI_Allreduce(&active, &totalActive, 1, MPI_INT, MPI_SUM, comm);

  dlong minNrows=0, maxNrows=0;
  hlong totalNrows=0;
  dfloat avgNrows;
  MPI_Allreduce(&N, &maxNrows, 1, MPI_DLONG, MPI_MAX, comm);
  MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, comm);
  avgNrows = (dfloat) totalNrows/totalActive;

  if (N==0) N=maxNrows; //set this so it's ignored for the global min
  MPI_Allreduce(&N, &minNrows, 1, MPI_DLONG, MPI_MIN, comm);

  long long int nnz;
  nnz = A->diag.nnz+A->offd.nnz;

  long long int minNnz=0, maxNnz=0, totalNnz=0;
  MPI_Allreduce(&nnz, &maxNnz,   1, MPI_LONG_LONG_INT, MPI_MAX, A->comm);
  MPI_Allreduce(&nnz, &totalNnz, 1, MPI_LONG_LONG_INT, MPI_SUM, A->comm);

  if (nnz==0) nnz = maxNnz; //set this so it's ignored for the global min
  MPI_Allreduce(&nnz, &minNnz, 1, MPI_LONG_LONG_INT, MPI_MIN, A->comm);

  dfloat nnzPerRow = (Nrows==0) ? 0 : (dfloat) nnz/Nrows;
  dfloat minNnzPerRow=0, maxNnzPerRow=0, avgNnzPerRow=0;
  MPI_Allreduce(&nnzPerRow, &maxNnzPerRow, 1, MPI_DFLOAT, MPI_MAX, A->comm);
  MPI_Allreduce(&nnzPerRow, &avgNnzPerRow, 1, MPI_DFLOAT, MPI_SUM, A->comm);
  avgNnzPerRow /= totalActive;

  if (Nrows==0) nnzPerRow = maxNnzPerRow;
  MPI_Allreduce(&nnzPerRow, &minNnzPerRow, 1, MPI_DFLOAT, MPI_MIN, A->comm);

  std::string name = "Exact Solve     ";

  if (rank==0){
    printf(" %3d  |  parAlmond |  %12lld  |  %12d  | %13d   |   %s|\n", lev, (long long int)totalNrows, minNrows, (int)minNnzPerRow, name.c_str());
    printf("      |            |                |  %12d  | %13d   |                   |\n", maxNrows, (int)maxNnzPerRow);
    printf("      |            |                |  %12d  | %13d   |                   |\n", (int)avgNrows, (int)avgNnzPerRow);
  }
}

exactSolver_t::~exactSolver_t() {
  if (coarseOffsets) free(coarseOffsets);
  if (coarseCounts) free(coarseCounts);
  if (diagInvAT) free(diagInvAT);
  if (offdInvAT) free(offdInvAT);
  if (diagRhs) free(diagRhs);
  if (offdRhs) free(offdRhs);
}

} //namespace parAlmond

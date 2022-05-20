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
#include "parAlmond/parAlmondCoarseSolver.hpp"
#include "parAlmond/parAlmondKernels.hpp"

namespace libp {

namespace parAlmond {

void exactSolver_t::solve(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_x) {

  stream_t currentStream = platform.getStream();

  //queue transfering coarse vector to host for Allgather
  if(N) {
    platform.finish();
    platform.setStream(ogs::ogsBase_t::dataStream);
    o_rhs.copyTo(diagRhs, N, 0, properties_t("async", true));
    platform.setStream(currentStream);
  }

  //queue local part of gemv
  const dfloat one=1.0;
  const dfloat zero=0.0;
  if (N)
    dGEMVKernel(N,N,one,o_diagInvAT,o_rhs, zero, o_x);

  if(offdTotal) {
    //wait for data to arrive on host
    platform.setStream(ogs::ogsBase_t::dataStream);
    platform.finish();


    //gather the offd rhs entries
    comm.Alltoallv(diagRhs,   sendCounts,   sendOffsets,
                   offdRhs, coarseCounts, coarseOffsets);

    //queue transfering coarse vector to device
    o_offdRhs.copyFrom(offdRhs, offdTotal, 0, properties_t("async", true));
    platform.finish(); //wait for transfer to complete

    platform.setStream(currentStream);

    //queue offd part of gemv
    if (N)
      dGEMVKernel(N,offdTotal, one, o_offdInvAT,o_offdRhs, one, o_x);
  }
}


int exactSolver_t::getTargetSize() {
  return 1000;
}

void exactSolver_t::setup(parCSR& _A, bool nullSpace,
                          memory<dfloat> nullVector, dfloat nullSpacePenalty) {

  A = _A;

  comm = A.comm;
  rank = comm.rank();
  size = comm.size();

  //copy the global coarse partition as ints
  coarseOffsets.malloc(size+1);
  for (int r=0;r<size+1;r++) {
    coarseOffsets[r] = static_cast<int>(A.globalRowStarts[r]);
  }

  coarseTotal   = coarseOffsets[size];
  coarseOffset  = coarseOffsets[rank];

  N = static_cast<int>(A.Nrows);
  Nrows = A.Nrows;
  Ncols = A.Ncols;

  coarseCounts.malloc(size,0);

  int sendNNZ = static_cast<int>(A.diag.nnz+A.offd.nnz);

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE")))
  //   {printf("Setting up coarse solver...");fflush(stdout);}

  memory<parCOO::nonZero_t> sendNonZeros(sendNNZ);

  //populate matrix
  int cnt = 0;
  for (int n=0;n<N;n++) {
    const int start = static_cast<int>(A.diag.rowStarts[n]);
    const int end   = static_cast<int>(A.diag.rowStarts[n+1]);
    for (int m=start;m<end;m++) {
      sendNonZeros[cnt].row = n + coarseOffset;
      sendNonZeros[cnt].col = A.diag.cols[m] + coarseOffset;
      sendNonZeros[cnt].val = A.diag.vals[m];
      cnt++;
    }
  }

  for (int n=0;n<A.offd.nzRows;n++) {
    const int row   = static_cast<int>(A.offd.rows[n]);
    const int start = static_cast<int>(A.offd.mRowStarts[n]);
    const int end   = static_cast<int>(A.offd.mRowStarts[n+1]);
    for (int m=start;m<end;m++) {
      sendNonZeros[cnt].row = row + coarseOffset;
      sendNonZeros[cnt].col = A.colMap[A.offd.cols[m]];
      sendNonZeros[cnt].val = A.offd.vals[m];
      cnt++;
    }
  }

  //get the nonzero counts from all ranks
  memory<int> recvNNZ(size);
  memory<int> NNZoffsets(size+1,0);
  comm.Allgather(sendNNZ, recvNNZ);

  int totalNNZ = 0;
  for (int r=0;r<size;r++) {
    totalNNZ += recvNNZ[r];
    NNZoffsets[r+1] = NNZoffsets[r] + recvNNZ[r];
  }

  memory<parCOO::nonZero_t> recvNonZeros(totalNNZ);

  comm.Allgatherv(sendNonZeros, sendNNZ,
                  recvNonZeros, recvNNZ, NNZoffsets);

  //gather null vector
  memory<dfloat> nullTotal(coarseTotal);

  for (int r=0;r<size;r++) {
    coarseCounts[r] = coarseOffsets[r+1]-coarseOffsets[r];
  }

  comm.Allgatherv(nullVector, N,
                  nullTotal, coarseCounts, coarseOffsets);

  //assemble the full matrix
  memory<dfloat> coarseA(coarseTotal*coarseTotal, 0.0);
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

  linAlg_t::matrixInverse(coarseTotal, coarseA);

  //determine size of offd piece
  offdTotal = coarseTotal - N;

  //shift offsets for MPI_AllGatherv of offd pieces
  for (int r=rank+1;r<=size;r++)
    coarseOffsets[r]-= N;

  //dont copy the local piece in MPI_AllGatherv
  coarseCounts[rank]=0;

  //counts for all-to-all
  sendCounts.malloc(size);
  sendOffsets.malloc(size);
  for (int r=0;r<size;r++) {
    sendCounts[r] = N;
    sendOffsets[r] = 0;
  }
  sendCounts[rank] = 0;

  //diag piece of invA
  diagInvAT.malloc(N*N);
  for (int n=0;n<N;n++) {
    for (int m=0;m<N;m++) {
      diagInvAT[n+m*N] = coarseA[(n+coarseOffset)*coarseTotal+(m+coarseOffset)];
    }
  }

  //offd piece of invA
  offdInvAT.malloc(N*offdTotal);
  for (int n=0;n<N;n++) {
    for (int m=0;m<coarseOffset;m++) {
      offdInvAT[n+m*N] = coarseA[(n+coarseOffset)*coarseTotal+m];
    }
    for (int m=coarseOffset+N;m<coarseTotal;m++) {
      offdInvAT[n+(m-N)*N] = coarseA[(n+coarseOffset)*coarseTotal+m];
    }
  }

  o_diagInvAT = platform.malloc<dfloat>(diagInvAT);
  o_offdInvAT = platform.malloc<dfloat>(offdInvAT);

  diagRhs.malloc(N);
  offdRhs.malloc(offdTotal);

  o_offdRhs = platform.malloc<dfloat>(offdTotal);

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE"))) printf("done.\n");
}

void exactSolver_t::syncToDevice() {}

void exactSolver_t::Report(int lev) {

  int totalActive = (N>0) ? 1:0;
  comm.Allreduce(totalActive, Comm::Sum);

  dlong minNrows=N, maxNrows=N;
  hlong totalNrows=N;
  comm.Allreduce(maxNrows, Comm::Max);
  comm.Allreduce(totalNrows, Comm::Sum);
  dfloat avgNrows = (dfloat) totalNrows/totalActive;

  if (N==0) minNrows=maxNrows; //set this so it's ignored for the global min
  comm.Allreduce(minNrows, Comm::Min);

  long long int nnz;
  nnz = A.diag.nnz+A.offd.nnz;

  long long int minNnz=nnz, maxNnz=nnz, totalNnz=nnz;
  comm.Allreduce(maxNnz,   Comm::Max);
  comm.Allreduce(totalNnz, Comm::Sum);

  if (nnz==0) minNnz = maxNnz; //set this so it's ignored for the global min
  comm.Allreduce(minNnz, Comm::Min);

  dfloat nnzPerRow = (Nrows==0) ? 0 : (dfloat) nnz/Nrows;
  dfloat minNnzPerRow=nnzPerRow, maxNnzPerRow=nnzPerRow, avgNnzPerRow=nnzPerRow;
  comm.Allreduce(maxNnzPerRow, Comm::Max);
  comm.Allreduce(avgNnzPerRow, Comm::Sum);
  avgNnzPerRow /= totalActive;

  if (Nrows==0) minNnzPerRow = maxNnzPerRow;
  comm.Allreduce(minNnzPerRow, Comm::Min);

  std::string name = "Exact Solve     ";

  if (rank==0){
    printf(" %3d  |  parAlmond |  %12lld  |  %12d  | %13d   |   %s|\n", lev, (long long int)totalNrows, minNrows, (int)minNnzPerRow, name.c_str());
    printf("      |            |                |  %12d  | %13d   |                   |\n", maxNrows, (int)maxNnzPerRow);
    printf("      |            |                |  %12d  | %13d   |                   |\n", (int)avgNrows, (int)avgNnzPerRow);
  }
}

} //namespace parAlmond

} //namespace libp

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

void oasSolver_t::solve(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_x) {

  A.halo.ExchangeStart(o_rhs, 1);

  //queue local part of gemv
  const dfloat one=1.0;
  const dfloat zero=0.0;
  if (N)
    dGEMVKernel(N,diagTotal,one,o_diagInvAT,o_rhs, zero, o_x);

  A.halo.ExchangeFinish(o_rhs, 1);

  //queue offd part of gemv
  if(offdTotal && N)
    dGEMVKernel(N,offdTotal, one, o_offdInvAT,
                o_rhs+diagTotal, one, o_x);

  A.halo.Combine(o_x, 1);
}


int oasSolver_t::getTargetSize() {
  return 1000*comm.size();
}

void oasSolver_t::setup(parCSR& _A, bool nullSpace,
                        memory<dfloat> nullVector, dfloat nullSpacePenalty) {

  A = _A;

  comm = A.comm;
  rank = comm.rank();
  size = comm.size();

  N = static_cast<int>(A.Ncols);
  Nrows = A.Nrows;
  Ncols = A.Ncols;

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE")))
  //   {printf("Setting up coarse solver...");fflush(stdout);}

  //make an overlapping patch matrix by collecting the rows
  // corresponding the offd columns

  //need to find where to send local rows
  memory<hlong> recvRows(A.Ncols-A.Nrows);

  memory<int> sendCounts(size);
  memory<int> recvCounts(size, 0);
  memory<int> sendOffsets(size+1, 0);
  memory<int> recvOffsets(size+1, 0);

  //use the colMap to fill the recv sizes
  int r=0;
  for (int n=A.Nrows;n<A.Ncols;n++) {
    hlong id = A.colMap[n];
    while (id>=A.globalRowStarts[r+1]) r++; //assumes the halo is sorted
    recvCounts[r]++;
    recvRows[n-A.Nrows] = id; //record the row to recv
  }

  //share the counts
  comm.Alltoall(recvCounts, sendCounts);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }

  int sendTotal = sendOffsets[size];
  memory<hlong> sendRows(sendTotal);

  //share the rowIds
  comm.Alltoallv(recvRows, recvCounts, recvOffsets,
                 sendRows, sendCounts, sendOffsets);

  //we now have a list of rows to send, count the nnz to send
  dlong nnzTotal=0;
  for (r=0;r<size;r++) {
    sendCounts[r] =0; //reset
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      dlong i = static_cast<dlong>(sendRows[n]-A.globalRowStarts[rank]); //local row id
      sendCounts[r]+= A.diag.rowStarts[i+1]-A.diag.rowStarts[i]; //count entries in this row
      sendCounts[r]+= A.offd.rowStarts[i+1]-A.offd.rowStarts[i]; //count entries in this row
    }
    nnzTotal += sendCounts[r]; //tally the total
  }

  memory<parCOO::nonZero_t> sendNonZeros(nnzTotal);

  nnzTotal=0; //reset
  for (r=0;r<size;r++) {
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      dlong i = static_cast<dlong>(sendRows[n] - A.globalRowStarts[rank]); //local row id
      for (dlong jj=A.diag.rowStarts[i]; jj<A.diag.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = A.diag.cols[jj] + A.globalRowStarts[rank];
        sendNonZeros[nnzTotal].val = A.diag.vals[jj];
        nnzTotal++;
      }
      for (dlong jj=A.offd.rowStarts[i]; jj<A.offd.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = A.colMap[A.offd.cols[jj]];
        sendNonZeros[nnzTotal].val = A.offd.vals[jj];
        nnzTotal++;
      }
    }
  }

  comm.Alltoall(sendCounts, recvCounts);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }

  nnzTotal = recvOffsets[size]; //total nonzeros

  memory<parCOO::nonZero_t> recvNonZeros(nnzTotal);

  comm.Alltoallv(sendNonZeros, sendCounts, sendOffsets,
                 recvNonZeros, recvCounts, recvOffsets);

  //we now have all the nonlocal rows (should also be sorted)

  //first re-index the column indices
  dlong id=A.Nrows;
  for (dlong n=0;n<nnzTotal;n++) {
    const hlong row = recvNonZeros[n].row;

    while(A.colMap[id]!=row) id++; //shift along list of recieved columns

    recvNonZeros[n].row = id; //overwrite with new local row id

    //now check the column index
    hlong col = recvNonZeros[n].col;
    if (col >= A.globalRowStarts[rank] && col < A.globalRowStarts[rank+1]) {//local column
      recvNonZeros[n].col = col - A.globalRowStarts[rank];//overwrite with local col id
    } else {
      int flag = 0;
      for (dlong jj=A.Nrows;jj<A.Ncols;jj++) { //look for the right id in the halo
        if (A.colMap[jj]==col) {
          recvNonZeros[n].col = jj;//overwrite with local col id
          flag = 1;
          break;
        }
      }
      if (flag==0) recvNonZeros[n].col = -1; //ignore this entry as its not in our patch the patch
    }
  }

  //assemble the full matrix
  memory<dfloat> coarseA(N*N);
  for (int n=0;n<A.Nrows;n++) {
    const int start = static_cast<int>(A.diag.rowStarts[n]);
    const int end   = static_cast<int>(A.diag.rowStarts[n+1]);
    for (int m=start;m<end;m++) {
      int col = static_cast<int>(A.diag.cols[m]);
      coarseA[n*N+col] = A.diag.vals[m];
    }
  }

  for (int n=0;n<A.offd.nzRows;n++) {
    const int row   = static_cast<int>(A.offd.rows[n]);
    const int start = static_cast<int>(A.offd.mRowStarts[n]);
    const int end   = static_cast<int>(A.offd.mRowStarts[n+1]);
    for (int m=start;m<end;m++) {
      int col = static_cast<int>(A.offd.cols[m]);
      coarseA[row*N+col] = A.offd.vals[m];
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
    memory<dfloat> null(A.Ncols);
    null.copyFrom(nullVector, A.Nrows);
    A.halo.Exchange(null, 1);

    for (int n=0;n<N;n++) {
      for (int m=0;m<N;m++) {
        coarseA[n*N+m] += nullSpacePenalty*null[n]*null[m];
      }
    }
  }

  linAlg_t::matrixInverse(N, coarseA);

  //determine the overlap weighting
  memory<dfloat> weight(N, 1.0);

  A.halo.Combine(weight, 1);

  for (int n=0;n<N;n++) {
    for (int m=0;m<N;m++) {
      coarseA[n*N+m] *= 1.0/sqrt(weight[n]*weight[m]);
    }
  }

  //determine size of offd piece
  diagTotal = A.Nrows;
  offdTotal = A.Ncols - A.Nrows;

  //diag piece of invA
  diagInvAT.malloc(N*diagTotal);
  for (int n=0;n<N;n++) {
    for (int m=0;m<diagTotal;m++) {
      diagInvAT[n+m*N] = coarseA[n*N+m];
    }
  }

  //offd piece of invA
  offdInvAT.malloc(N*offdTotal);
  for (int n=0;n<N;n++) {
    for (int m=0;m<offdTotal;m++) {
      offdInvAT[n+m*N] = coarseA[n*N + m+diagTotal];
    }
  }

  o_diagInvAT = platform.malloc<dfloat>(diagInvAT);
  o_offdInvAT = platform.malloc<dfloat>(offdInvAT);

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE"))) printf("done.\n");
}

void oasSolver_t::syncToDevice() {}

void oasSolver_t::Report(int lev) {

  int totalActive = (N>0) ? 1:0;
  comm.Allreduce(totalActive, Comm::Sum);

  dlong minNrows=N, maxNrows=N;
  hlong totalNrows=N;
  comm.Allreduce(maxNrows, Comm::Max);
  comm.Allreduce(totalNrows, Comm::Sum);
  dfloat avgNrows = static_cast<dfloat>(totalNrows)/totalActive;

  if (N==0) minNrows=maxNrows; //set this so it's ignored for the global min
  comm.Allreduce(minNrows, Comm::Min);

  long long int nnz;
  nnz = A.diag.nnz+A.offd.nnz;

  long long int minNnz=nnz, maxNnz=nnz, totalNnz=nnz;
  comm.Allreduce(maxNnz, Comm::Max);
  comm.Allreduce(totalNnz, Comm::Sum);

  if (nnz==0) minNnz = maxNnz; //set this so it's ignored for the global min
  comm.Allreduce(minNnz, Comm::Min);

  dfloat nnzPerRow = (Nrows==0) ? 0 : static_cast<dfloat>(nnz)/Nrows;
  dfloat minNnzPerRow=nnzPerRow, maxNnzPerRow=nnzPerRow, avgNnzPerRow=nnzPerRow;
  comm.Allreduce(maxNnzPerRow, Comm::Max);
  comm.Allreduce(avgNnzPerRow, Comm::Sum);
  avgNnzPerRow /= totalActive;

  if (Nrows==0) minNnzPerRow = maxNnzPerRow;
  comm.Allreduce(minNnzPerRow, Comm::Min);

  std::string name = "OAS             ";

  if (rank==0){
    printf(" %3d  |  parAlmond |  %12lld  |  %12d  | %13d   |   %s|\n", lev, (long long int)totalNrows, minNrows, (int)minNnzPerRow, name.c_str());
    printf("      |            |                |  %12d  | %13d   |                   |\n", maxNrows, (int)maxNnzPerRow);
    printf("      |            |                |  %12d  | %13d   |                   |\n", (int)avgNrows, (int)avgNnzPerRow);
  }
}

} //namespace parAlmond

} //namespace libp

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

namespace parAlmond {

void coarseSolver_t::solve(occa::memory& o_rhs, occa::memory& o_x) {

  if (gatherLevel) {
    ogs->Gather(o_Gx, o_rhs, ogs_dfloat, ogs_add, ogs_notrans);

    if(N) o_Gx.copyTo(rhsLocal, N*sizeof(dfloat), 0);
  } else {
    if(N) o_rhs.copyTo(rhsLocal, N*sizeof(dfloat), 0);
  }

  //gather the full vector
  MPI_Allgatherv(rhsLocal,             N,                MPI_DFLOAT,
                 rhsCoarse, coarseCounts, coarseOffsets, MPI_DFLOAT, comm);

  //multiply by local part of the exact matrix inverse
  // #pragma omp parallel for
  for (int n=0;n<N;n++) {
    xLocal[n] = 0.;
    for (int m=0;m<coarseTotal;m++) {
      xLocal[n] += invCoarseA[n*coarseTotal+m]*rhsCoarse[m];
    }
  }

  if (gatherLevel) {
    if(N) o_Gx.copyFrom(xLocal, N*sizeof(dfloat), 0);
    ogs->Scatter(o_x, o_Gx, ogs_dfloat, ogs_add, ogs_notrans);
  } else {
    if(N) o_x.copyFrom(xLocal, N*sizeof(dfloat), 0);
  }
}


int coarseSolver_t::getTargetSize() {
  return 1000;
}

typedef struct {

  hlong row;
  hlong col;
  dfloat val;

} nonzero_t;

void coarseSolver_t::setup(parCSR *A, bool nullSpace,
                           dfloat *nullVector, dfloat nullSpacePenalty) {

  comm = A->comm;

  int rank;
  int size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  //copy the global coarse partition as ints
  coarseOffsets = (int* ) calloc(size+1,sizeof(int));
  for (int r=0;r<size+1;r++) coarseOffsets[r] = (int) A->globalRowStarts[r];

  coarseTotal   = coarseOffsets[size];
  coarseOffset  = coarseOffsets[rank];

  N = (int) A->Nrows;

  coarseCounts = (int*) calloc(size,sizeof(int));

  int sendNNZ = (int) (A->diag.nnz+A->offd.nnz);

  if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE")))
    {printf("Setting up coarse solver...");fflush(stdout);}

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

  nonzero_t *sendNonZeros = (nonzero_t *) calloc(sendNNZ, sizeof(nonzero_t));

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

  nonzero_t *recvNonZeros = (nonzero_t *) calloc(totalNNZ, sizeof(nonzero_t));

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
  MPI_Type_free(&MPI_NONZERO_T);
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

  //store only the local rows of the full inverse
  invCoarseA = (dfloat *) calloc(N*coarseTotal,sizeof(dfloat));
  for (int n=0;n<N;n++) {
    for (int m=0;m<coarseTotal;m++) {
      invCoarseA[n*coarseTotal+m] = coarseA[(n+coarseOffset)*coarseTotal+m];
    }
  }

  xLocal   = (dfloat*) calloc(N,sizeof(dfloat));
  rhsLocal = (dfloat*) calloc(N,sizeof(dfloat));

  xCoarse   = (dfloat*) calloc(coarseTotal,sizeof(dfloat));
  rhsCoarse = (dfloat*) calloc(coarseTotal,sizeof(dfloat));

  free(coarseA);

  if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE"))) printf("done.\n");
}

void coarseSolver_t::syncToDevice() {}

coarseSolver_t::~coarseSolver_t() {
  if (coarseOffsets) free(coarseOffsets);
  if (coarseCounts) free(coarseCounts);
  if (invCoarseA) free(invCoarseA);
  if (xLocal) free(xLocal);
  if (rhsLocal) free(rhsLocal);
  if (xCoarse) free(xCoarse);
  if (rhsCoarse) free(rhsCoarse);
}

} //namespace parAlmond

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

namespace parAlmond {

coarseSolver::coarseSolver(settings_t& settings_):
  settings(settings_) {
  gatherLevel = false;
}

int coarseSolver::getTargetSize() {
  return 1000;
}

//set up exact solver using xxt
void coarseSolver::setup(parCSR *A) {

  comm = A->comm;

  int rank, size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  //copy the global coarse partition as ints
  coarseOffsets = (int* ) calloc(size+1,sizeof(int));
  for (int r=0;r<size+1;r++) coarseOffsets[r] = (int) A->globalRowStarts[r];

  coarseTotal   = coarseOffsets[size];
  coarseOffset  = coarseOffsets[rank];

  N = (int) A->Nrows;

  coarseCounts = (int*) calloc(size,sizeof(int));

  // had to move this later
  // if(settings.compareSetting("PARALMOND SMOOTH COARSEST", "TRUE")){
  //   if(rank==0) printf("WARNING !!!!!: not building coarsest level matrix\n");
  //   return; // bail early as this will not get used
  // }



  int sendNNZ = (int) (A->diag->nnz+A->offd->nnz);

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
    int start = (int) A->diag->rowStarts[n];
    int end   = (int) A->diag->rowStarts[n+1];
    for (int m=start;m<end;m++) {
      sendNonZeros[cnt].row = n + coarseOffset;
      sendNonZeros[cnt].col = A->diag->cols[m] + coarseOffset;
      sendNonZeros[cnt].val = A->diag->vals[m];
      cnt++;
    }
    start = (int) A->offd->rowStarts[n];
    end   = (int) A->offd->rowStarts[n+1];
    for (dlong m=start;m<end;m++) {
      sendNonZeros[cnt].row = n + coarseOffset;
      sendNonZeros[cnt].col = A->colMap[A->offd->cols[m]];
      sendNonZeros[cnt].val = A->offd->vals[m];
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

  MPI_Allgatherv(  A->null,          N,                MPI_DFLOAT,
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

  if (A->nullSpace) { //A is dense due to nullspace augmentation
    for (int n=0;n<coarseTotal;n++) {
      for (int m=0;m<coarseTotal;m++) {
        coarseA[n*coarseTotal+m] += A->nullSpacePenalty*nullTotal[n]*nullTotal[m];
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

void coarseSolver::syncToDevice() {}

void coarseSolver::solve(dfloat *rhs, dfloat *x) {

  if (gatherLevel) {
    //weight
    vectorDotStar(ogs->N, 1.0, o_gatherWeight, rhs, 0.0, Sx);
    ogs->Gather(Gx, Sx, ogs_dfloat, ogs_add);

    //gather the full vector
    MPI_Allgatherv(Gx,                  N,                MPI_DFLOAT,
                   rhsCoarse, coarseCounts, coarseOffsets, MPI_DFLOAT, comm);

    //multiply by local part of the exact matrix inverse
    // #pragma omp parallel for

    for (int n=0;n<N;n++) {
#if 1
      xLocal[n] = 0.;
      for (int m=0;m<coarseTotal;m++) {
        xLocal[n] += invCoarseA[n*coarseTotal+m]*rhsCoarse[m];
      }

#else
      xLocal[n] = rhsCoarse[n];
#endif

    }
    ogs->Scatter(x, xLocal, ogs_dfloat, ogs_add);

  } else {
    //gather the full vector
    MPI_Allgatherv(rhs,                  N,                MPI_DFLOAT,
                   rhsCoarse, coarseCounts, coarseOffsets, MPI_DFLOAT, comm);

    printf("HACKING COARSE GRID\n");

    //multiply by local part of the exact matrix inverse
   // #pragma omp parallel for

    for (int n=0;n<N;n++) {
#if 1
      x[n] = 0.;
      for (int m=0;m<coarseTotal;m++) {
        x[n] += invCoarseA[n*coarseTotal+m]*rhsCoarse[m];
      }
#else
      x[n] = rhsCoarse[n];
#endif

    }
  }


}

void coarseSolver::solve(occa::memory o_rhs, occa::memory o_x) {

  if (gatherLevel) {
    //weight
    vectorDotStar(ogs->N, 1.0, o_gatherWeight, o_rhs, 0.0, o_Sx);
    ogs->Gather(o_Gx, o_Sx, ogs_dfloat, ogs_add);

    if(N)
      o_Gx.copyTo(rhsLocal, N*sizeof(dfloat), 0);
  } else {
    if(N)
      o_rhs.copyTo(rhsLocal, N*sizeof(dfloat), 0);
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
    if(N)
      o_Gx.copyFrom(xLocal, N*sizeof(dfloat), 0);
    ogs->Scatter(o_x, o_Gx, ogs_dfloat, ogs_add);
  } else {
    if(N)
      o_x.copyFrom(xLocal, N*sizeof(dfloat), 0);
  }
}

#if 0
//set up exact solver using xxt
void setupExactSolve(parAlmond_t *parAlmond, agmgLevel *level, bool nullSpace, dfloat nullSpacePenalty) {

  int rank, size;
  rank = agmg::rank;
  size = agmg::size;

  int* coarseOffsets = level->globalRowStarts;
  int coarseTotal = coarseOffsets[size];
  int coarseOffset = coarseOffsets[rank];

  int *globalNumbering = (int *) calloc(coarseTotal,sizeof(int));
  for (int n=0;n<coarseTotal;n++)
    globalNumbering[n] = n;

  csr *A = level->A;
  int N = level->Nrows;

  int totalNNZ;
  int *rows;
  int *cols;
  dfloat *vals;

  if(!nullSpace) {
    //if no nullspace, use sparse A
    totalNNZ = A->diagNNZ+A->offdNNZ;
    if (totalNNZ) {
      rows = (int *) calloc(totalNNZ,sizeof(int));
      cols = (int *) calloc(totalNNZ,sizeof(int));
      vals = (dfloat *) calloc(totalNNZ,sizeof(dfloat));
    }

    //populate matrix
    int cnt = 0;
    for (int n=0;n<N;n++) {
      for (int m=A->diagRowStarts[n];m<A->diagRowStarts[n+1];m++) {
        rows[cnt] = n + coarseOffset;
        cols[cnt] = A->diagCols[m] + coarseOffset;
        vals[cnt] = A->diagCoefs[m];
        cnt++;
      }
      for (int m=A->offdRowStarts[n];m<A->offdRowStarts[n+1];m++) {
        rows[cnt] = n + coarseOffset;
        cols[cnt] = A->colMap[A->offdCols[m]];
        vals[cnt] = A->offdCoefs[m];
        cnt++;
      }
    }
  } else {
    totalNNZ = A->Nrows*coarseTotal; //A is dense due to nullspace augmentation
    if (totalNNZ) {
      rows = (int *) calloc(totalNNZ,sizeof(int));
      cols = (int *) calloc(totalNNZ,sizeof(int));
      vals = (dfloat *) calloc(totalNNZ,sizeof(dfloat));
    }

    //gather null vector
    dfloat *nullTotal = (dfloat*) calloc(coarseTotal,sizeof(dfloat));
    int *nullCounts = (int*) calloc(size,sizeof(int));
    for (int r=0;r<size;r++)
      nullCounts[r] = coarseOffsets[r+1]-coarseOffsets[r];

    MPI_Allgatherv(A->null, A->Nrows, MPI_DFLOAT, nullTotal, nullCounts, coarseOffsets, MPI_DFLOAT, agmg::comm);

    //populate matrix
    for (int n=0;n<N;n++) {
      for (int m=0;m<coarseTotal;m++) {
        rows[n*coarseTotal+m] = n + coarseOffset;
        cols[n*coarseTotal+m] = m;
        vals[n*coarseTotal+m] = nullSpacePenalty*nullTotal[n+coarseOffset]*nullTotal[m];
      }
    }

    for (int n=0;n<N;n++) {
      for (int m=A->diagRowStarts[n];m<A->diagRowStarts[n+1];m++) {
        int col = A->diagCols[m] + coarseOffset;
        vals[n*coarseTotal+col] += A->diagCoefs[m];
      }
      for (int m=A->offdRowStarts[n];m<A->offdRowStarts[n+1];m++) {
        int col = A->colMap[A->offdCols[m]];
        vals[n*coarseTotal+col] += A->offdCoefs[m];
      }
    }
  }

  parAlmond->ExactSolve = xxtSetup(A->Nrows,
                                globalNumbering,
                                totalNNZ,
                                rows,
                                cols,
                                vals,
                                0,
                                "int",
                                dfloatString);

  parAlmond->coarseTotal = coarseTotal;
  parAlmond->coarseOffset = coarseOffset;

  parAlmond->xCoarse   = (dfloat*) calloc(coarseTotal,sizeof(dfloat));
  parAlmond->rhsCoarse = (dfloat*) calloc(coarseTotal,sizeof(dfloat));

  free(globalNumbering);
  if (totalNNZ) {
    free(rows);
    free(cols);
    free(vals);
  }

  printf("Done UberCoarse setup\n");
}


void exactCoarseSolve(parAlmond_t *parAlmond, int N, dfloat *rhs, dfloat *x) {

  //use coarse solver
  for (int n=0;n<parAlmond->coarseTotal;n++)
    parAlmond->rhsCoarse[n] =0.;

  for (int n=0;n<N;n++)
    parAlmond->rhsCoarse[n+parAlmond->coarseOffset] = rhs[n];

  xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);

  for (int n=0;n<N;n++)
    x[n] = parAlmond->xCoarse[n+parAlmond->coarseOffset];

}

void device_exactCoarseSolve(parAlmond_t *parAlmond, int N, occa::memory o_rhs, occa::memory o_x) {

  //use coarse solver
  for (int n=0;n<parAlmond->coarseTotal;n++)
    parAlmond->rhsCoarse[n] =0.;

  o_rhs.copyTo(parAlmond->rhsCoarse+parAlmond->coarseOffset);
  xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);
  o_x.copyFrom(parAlmond->xCoarse+parAlmond->coarseOffset,N*sizeof(dfloat));
}
#endif

} //namespace parAlmond

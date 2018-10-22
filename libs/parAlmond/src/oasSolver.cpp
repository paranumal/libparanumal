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

oasSolver::oasSolver(setupAide options_,
                     MPI_Comm comm_) {

  gatherLevel = false;
  options = options_;

  comm = comm_;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
}

int oasSolver::getTargetSize() {
  return 10000*size;
}

void oasSolver::setup(parCSR *A_) {

  A = A_;
  device = A->device;
  N = A->Nrows;

  Nrows = A->Nrows;
  Ncols = A->Ncols;

  //split the MPI comm by rank. All scomms are size=1
  MPI_Comm scomm;
  MPI_Comm_split(comm, rank, 0, &scomm);

  //make an overlapping patch matrix by collecting the offd rows
  parCSR *Ahat = new parCSR(A->Ncols, A->Ncols, scomm, device);

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
                A->comm);

  //we now have a list of rows to send, count the nnz to send
  dlong nnzTotal=0;
  for (r=0;r<size;r++) {
    sendCounts[r] =0; //reset
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      dlong i = (dlong) (sendRows[n]-A->globalRowStarts[rank]); //local row id
      sendCounts[r]+= A->diag->rowStarts[i+1]-A->diag->rowStarts[i]; //count entries in this row
      sendCounts[r]+= A->offd->rowStarts[i+1]-A->offd->rowStarts[i]; //count entries in this row
    }
    nnzTotal += sendCounts[r]; //tally the total
  }

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

  nonzero_t *sendNonZeros = (nonzero_t *) calloc(nnzTotal, sizeof(nonzero_t));

  nnzTotal=0; //reset
  for (r=0;r<size;r++) {
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      dlong i = (dlong) (sendRows[n] - A->globalRowStarts[rank]); //local row id
      for (dlong jj=A->diag->rowStarts[i]; jj<A->diag->rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = A->diag->cols[jj] + A->globalRowStarts[rank];
        sendNonZeros[nnzTotal].val = A->diag->vals[jj];
        nnzTotal++;
      }
      for (dlong jj=A->offd->rowStarts[i]; jj<A->offd->rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = A->colMap[A->offd->cols[jj]];
        sendNonZeros[nnzTotal].val = A->offd->vals[jj];
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

  nnzTotal = recvOffsets[size]; //total nonzeros

  nonzero_t *recvNonZeros = (nonzero_t *) calloc(nnzTotal, sizeof(nonzero_t));

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

  //we now have all the nonlocal rows (should also be sorted)

  //fill in the additional entries for the overlapping patch matrix
  Ahat->diag->rowStarts = (dlong *) calloc(Ahat->Nrows+1, sizeof(dlong));
  Ahat->offd->rowStarts = (dlong *) calloc(Ahat->Nrows+1, sizeof(dlong));

  for (dlong i=0;i<A->Nrows;i++) {
    Ahat->diag->rowStarts[i+1] = A->diag->rowStarts[i+1]-A->diag->rowStarts[i]
                                +A->offd->rowStarts[i+1]-A->offd->rowStarts[i];
  }

  dlong id=A->Nrows;
  for (dlong n=0;n<nnzTotal;n++) {
    hlong row = recvNonZeros[n].row;

    while(A->colMap[id]!=row) id++;

    recvNonZeros[n].row = id; //overwrite with local row id

    hlong col = recvNonZeros[n].col;
    if (col >= A->globalRowStarts[rank] && col < A->globalRowStarts[rank+1]) {
      Ahat->diag->rowStarts[id+1]++;
      recvNonZeros[n].col = col - A->globalRowStarts[rank];//overwrite with local col id
    } else {
      int flag = 0;
      for (dlong jj=A->Nrows;jj<A->Ncols;jj++) { //look for the right id in the halo
        if (A->colMap[jj]==col) {
          Ahat->diag->rowStarts[id+1]++;
          recvNonZeros[n].col = jj;//overwrite with local col id
          flag = 1;
          break;
        }
      }
      if (flag==0) recvNonZeros[n].col = -1; //ignore this entry as its not local to the patch
    }
  }

  // cumulative sum
  for(dlong i=0; i<Ahat->Nrows; i++)
    Ahat->diag->rowStarts[i+1] += Ahat->diag->rowStarts[i];

  Ahat->diag->nnz = Ahat->diag->rowStarts[Ahat->Nrows];
  Ahat->offd->nnz = Ahat->offd->rowStarts[Ahat->Nrows];

  Ahat->globalRowStarts = (hlong *) calloc(2,sizeof(hlong));
  Ahat->globalRowStarts[1] = Ahat->Nrows;
  Ahat->globalColStarts = Ahat->globalRowStarts;
  Ahat->haloSetup(NULL);

  //fill the CSR matrices
  Ahat->diagA   = (dfloat *) calloc(Ahat->Ncols, sizeof(dfloat));
  Ahat->diagInv = (dfloat *) calloc(Ahat->Ncols, sizeof(dfloat));
  Ahat->diag->cols = (dlong *)  calloc(Ahat->diag->nnz, sizeof(dlong));
  Ahat->diag->vals = (dfloat *) calloc(Ahat->diag->nnz, sizeof(dfloat));

  dlong cnt=0;
  for (dlong i=0;i<A->Nrows;i++) {
    for (dlong jj=A->diag->rowStarts[i]; jj<A->diag->rowStarts[i+1];jj++){
      Ahat->diag->cols[cnt] = A->diag->cols[jj];
      Ahat->diag->vals[cnt] = A->diag->vals[jj];
      cnt++;
    }
    for (dlong jj=A->offd->rowStarts[i]; jj<A->offd->rowStarts[i+1];jj++){
      Ahat->diag->cols[cnt] = A->offd->cols[jj];
      Ahat->diag->vals[cnt] = A->offd->vals[jj];
      cnt++;
    }
  }

  for (dlong n=0;n<nnzTotal;n++) {
    dlong col = (dlong) recvNonZeros[n].col;
    dfloat val = recvNonZeros[n].val;

    if (col<0) continue;

    Ahat->diag->cols[cnt] = col;
    Ahat->diag->vals[cnt] = val;
    cnt++;
  }

  //copy the diagonal info directly from the original parCSR
  memcpy(Ahat->diagA, A->diagA, A->Ncols*sizeof(dfloat));
  memcpy(Ahat->diagInv, A->diagInv, A->Ncols*sizeof(dfloat));


  Ahat->nullSpace = A->nullSpace;
  Ahat->nullSpacePenalty = A->nullSpacePenalty;

  Ahat->null = (dfloat*) calloc(Ahat->Nrows,sizeof(dfloat));
  memcpy(Ahat->null, A->null, A->Nrows*sizeof(dfloat));
  ogsGatherScatter(Ahat->null, ogsDfloat, ogsAdd, A->ogs);


  MPI_Barrier(A->comm);
  MPI_Type_free(&MPI_NONZERO_T);
  free(recvNonZeros);

  //set up an AMG solver for this patch
  CoarseType ctype = EXACTSOLVER;
  patchSolver = new solver_t(device, scomm, options, ctype);
  patchSolver->AMGSetup(Ahat);

  /* ----- now make the uber coarse problem -------*/
  hlong *FineToCoarse = (hlong *) malloc(A->Nrows*sizeof(hlong));
  hlong *globalAggStarts = (hlong *) calloc(size+1,sizeof(hlong));

  //all dofs on this rank coarsen to a single agglomerate
  for (int n=0;n<A->Nrows;n++) FineToCoarse[n] = rank;

  // this assumes there is at least one dof per rank before coarsening
  for (int r=0;r<size+1;r++) globalAggStarts[r] = r;

  dfloat *nullCoarseA;
  coarseP = constructProlongation(A, FineToCoarse, globalAggStarts, &nullCoarseA);
  coarseR = transpose(coarseP);
  coarseA = galerkinProd(A, coarseP);
  coarseA->null = nullCoarseA;

  //set up an AMG solver for the ubercoarse problem
  uberCoarseSolver = new exactSolver(options, comm);
  uberCoarseSolver->setup(coarseA);

  xCoarse   = (dfloat*) calloc(1,sizeof(dfloat));
  rhsCoarse = (dfloat*) calloc(1,sizeof(dfloat));

  if (A->Nshared)
    A->o_haloIds = device.malloc(A->Nshared*sizeof(dlong), A->haloIds);
}

void oasSolver::Report(int base) {

  int numLevels = patchSolver->numLevels;
  int maxNumLevels=0;

  MPI_Allreduce(&numLevels, &maxNumLevels, 1, MPI_INT, MPI_MAX, comm);

  for(int lev=0; lev<maxNumLevels-1; lev++) {

    dlong Nrows = (lev<numLevels-1) ? patchSolver->levels[lev]->Nrows : 0;
    hlong hNrows = (hlong) Nrows;

    int active = (Nrows>0) ? 1:0;
    int totalActive=0;
    MPI_Allreduce(&active, &totalActive, 1, MPI_INT, MPI_SUM, comm);

    dlong minNrows=0, maxNrows=0;
    hlong totalNrows=0;
    dfloat avgNrows;
    MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_DLONG, MPI_MAX, comm);
    MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, comm);
    avgNrows = (dfloat) totalNrows/totalActive;

    if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
    MPI_Allreduce(&Nrows, &minNrows, 1, MPI_DLONG, MPI_MIN, comm);


    long long int nnz;
    agmgLevel *level = (lev<numLevels-1) ? (agmgLevel*) patchSolver->levels[lev] : NULL;
    nnz = (lev<numLevels-1) ? level->A->diag->nnz+A->offd->nnz : 0;

    long long int minNnz=0, maxNnz=0, totalNnz=0;
    dfloat avgNnz;
    MPI_Allreduce(&nnz, &maxNnz,   1, MPI_LONG_LONG_INT, MPI_MAX, comm);
    MPI_Allreduce(&nnz, &totalNnz, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);
    avgNnz = (dfloat) totalNnz/totalActive;

    if (nnz==0) nnz = maxNnz; //set this so it's ignored for the global min
    MPI_Allreduce(&nnz, &minNnz, 1, MPI_LONG_LONG_INT, MPI_MIN, comm);

    dfloat nnzPerRow = (Nrows==0) ? 0 : (dfloat) nnz/Nrows;
    dfloat minNnzPerRow=0, maxNnzPerRow=0, avgNnzPerRow=0;
    MPI_Allreduce(&nnzPerRow, &maxNnzPerRow, 1, MPI_DFLOAT, MPI_MAX, comm);
    MPI_Allreduce(&nnzPerRow, &avgNnzPerRow, 1, MPI_DFLOAT, MPI_SUM, comm);
    avgNnzPerRow /= totalActive;

    if (Nrows==0) nnzPerRow = maxNnzPerRow;
    MPI_Allreduce(&nnzPerRow, &minNnzPerRow, 1, MPI_DFLOAT, MPI_MIN, comm);

    char smootherString[BUFSIZ];
    if (patchSolver->stype==DAMPED_JACOBI)
      strcpy(smootherString, "Damped Jacobi   ");
    else if (patchSolver->stype==CHEBYSHEV)
      strcpy(smootherString, "Chebyshev       ");

    if (rank==0){
      printf(" %3d  |  parAlmond |  %12d  | %13d   |   %s|\n", lev+base,minNrows, (int)minNnzPerRow, smootherString);
      printf("      |  OAS       |  %12d  | %13d   |                   |\n", maxNrows, (int)maxNnzPerRow);
      printf("      |            |  %12d  | %13d   |                   |\n", (int)avgNrows, (int)avgNnzPerRow);
    }
  }

  dlong Nrows = patchSolver->levels[numLevels-1]->Nrows;
  hlong hNrows = (hlong) Nrows;

  int active = (Nrows>0) ? 1:0;
  int totalActive=0;
  MPI_Allreduce(&active, &totalActive, 1, MPI_INT, MPI_SUM, comm);

  dlong minNrows=0, maxNrows=0;
  hlong totalNrows=0;
  dfloat avgNrows;
  MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_DLONG, MPI_MAX, comm);
  MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, comm);
  avgNrows = (dfloat) totalNrows/totalActive;

  if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
  MPI_Allreduce(&Nrows, &minNrows, 1, MPI_DLONG, MPI_MIN, comm);

  if (rank==0){
    printf(" %3d  |   Exact    |  %12d  |-------------------------------------|\n", base+maxNumLevels-1, minNrows);
    printf("      |   Solves   |  %12d  |                                     |\n", maxNrows);
    printf("      |            |  %12d  |-------------------------------------|\n", (int)avgNrows);

    printf(" %3d  |   Coarse   |  %12d  |-------------------------------------|\n", base+maxNumLevels, 1);
    printf("      |   Solve    |  %12d  | Total Size:   %5d  x%5d         |\n", 1, size, size);
    printf("      |            |  %12d  |-------------------------------------|\n", (int)1);
  }
}

void oasSolver::solve(dfloat *rhs, dfloat *x) {
  printf("ERROR: Calling OAS coarse solver with host memory not supported.\n");
  exit(-1);
}

void oasSolver::solve(occa::memory o_rhs, occa::memory o_x) {

  if (gatherLevel) {
    //gather to global indexing
    ogsGather(o_Grhs, o_rhs, ogsDfloat, ogsAdd, ogs);

    //agglomerate to one dof (just a scaled sum)
    A->haloExchangeStart(o_Grhs);
    dlong numBlocks = (N < NBLOCKS) ? N : NBLOCKS;
    vectorSumKernel(numBlocks, N, o_Grhs, o_reductionScratch);
    A->haloExchangeFinish(o_Grhs);
    o_reductionScratch.copyTo(rhsCoarse, 1*sizeof(dfloat));
    rhsCoarse[0] *= coarseR->diag->vals[0];

    patchSolver->levels[0]->o_x   = o_Gx;
    patchSolver->levels[0]->o_rhs = o_Grhs;

    //queue up the vcycle on the device patch
    patchSolver->device_vcycle(0);

    //while the device is busy, collect and solve the ubercoarse problem
    uberCoarseSolver->solve(rhsCoarse, xCoarse);

    //add the contributions from the patches together and scale by overlap degree
    ogsGatherScatter(o_Gx, ogsDfloat, ogsAdd, A->ogs); //this syncs host+device
    vectorDotStar(N, A->ogs->o_invDegree, o_Gx);

    //augment the patch solution with the ubercoarse solve
    vectorAddScalar(N, coarseP->diag->vals[0]*xCoarse[0], o_Gx);

    //scatter back to local indexing
    ogsScatter(o_x, o_Gx, ogsDfloat, ogsAdd, ogs);
  } else {
    //agglomerate to one dof (just a scaled sum)
    A->haloExchangeStart(o_rhs);
    dlong numBlocks = (N < NBLOCKS) ? N : NBLOCKS;
    vectorSumKernel(numBlocks, N, o_rhs, o_reductionScratch);
    A->haloExchangeFinish(o_rhs);
    o_reductionScratch.copyTo(rhsCoarse, 1*sizeof(dfloat));
    rhsCoarse[0] *= coarseR->diag->vals[0];

    patchSolver->levels[0]->o_x   = o_x;
    patchSolver->levels[0]->o_rhs = o_rhs;

    //queue up the vcycle on the device patch
    patchSolver->device_vcycle(0);

    //while the device is busy, collect and solve the ubercoarse problem
    uberCoarseSolver->solve(rhsCoarse, xCoarse);

    //add the contributions from the patches together and scale by overlap degree
    ogsGatherScatter(o_x, ogsDfloat, ogsAdd, A->ogs); //this syncs host+device
    vectorDotStar(N, A->ogs->o_invDegree, o_x);

    //augment the patch solution with the ubercoarse solve
    vectorAddScalar(N, coarseP->diag->vals[0]*xCoarse[0], o_x);
  }
}

} //namespace parAlmond
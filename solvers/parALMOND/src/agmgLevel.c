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

#include "agmg.h"

// parAlmond's function call-backs
void agmgAx(void **args, dfloat *x, dfloat *Ax){
  parAlmond_t *parAlmond = (parAlmond_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  axpy(level->A, 1.0, x, 0.0, Ax,parAlmond->nullSpace,parAlmond->nullSpacePenalty);
}

void agmgCoarsen(void **args, dfloat *r, dfloat *Rr){
  // parAlmond_t *parAlmond = (parAlmond_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  axpy(level->R, 1.0, r, 0.0, Rr,false,0.);
}

void agmgProlongate(void **args, dfloat *x, dfloat *Px){
  // parAlmond_t *parAlmond = (parAlmond_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  axpy(level->P, 1.0, x, 1.0, Px,false,0.);
}

void agmgSmooth(void **args, dfloat *rhs, dfloat *x, bool x_is_zero){
  parAlmond_t *parAlmond = (parAlmond_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  if(level->stype == JACOBI){
    smoothJacobi(parAlmond, level, level->A, rhs, x, x_is_zero);
  } else if(level->stype == DAMPED_JACOBI){
    smoothDampedJacobi(parAlmond, level, level->A, rhs, x, x_is_zero);
  } else if(level->stype == CHEBYSHEV){
    smoothChebyshev(parAlmond, level, level->A, rhs, x, x_is_zero);
  }
}

void device_agmgAx(void **args, occa::memory &o_x, occa::memory &o_Ax){
  parAlmond_t *parAlmond = (parAlmond_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  axpy(parAlmond,level->deviceA, 1.0, o_x, 0.0, o_Ax,parAlmond->nullSpace,parAlmond->nullSpacePenalty);
}

void device_agmgCoarsen(void **args, occa::memory &o_r, occa::memory &o_Rr){
  parAlmond_t *parAlmond = (parAlmond_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  axpy(parAlmond, level->deviceR, 1.0, o_r, 0.0, o_Rr,false,0.);
}

void device_agmgProlongate(void **args, occa::memory &o_x, occa::memory &o_Px){
  parAlmond_t *parAlmond = (parAlmond_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  axpy(parAlmond, level->dcsrP, 1.0, o_x, 1.0, o_Px);
}

void device_agmgSmooth(void **args, occa::memory &o_rhs, occa::memory &o_x, bool x_is_zero){
  parAlmond_t *parAlmond = (parAlmond_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  if(level->stype == JACOBI){
    smoothJacobi(parAlmond, level, level->deviceA, o_rhs, o_x, x_is_zero);
  } else if(level->stype == DAMPED_JACOBI){
    smoothDampedJacobi(parAlmond, level, level->deviceA, o_rhs, o_x, x_is_zero);
  } else if(level->stype == CHEBYSHEV){
    smoothChebyshev(parAlmond, level, level->deviceA, o_rhs, o_x, x_is_zero);
  }
}

dfloat rhoDinvA(parAlmond_t *parAlmond, csr *A, dfloat *invD);

void setupSmoother(parAlmond_t *parAlmond, agmgLevel *level, SmoothType s){

  level->stype = s;

  if((s == DAMPED_JACOBI)||(s == CHEBYSHEV)){
    // estimate rho(invD * A)
    dfloat rho=0;

    if(level->A->Nrows)
      level->A->diagInv = (dfloat *) calloc(level->A->Nrows, sizeof(dfloat));

    for (dlong i=0;i<level->A->Nrows;i++) {
      dfloat diag = level->A->diagCoefs[level->A->diagRowStarts[i]];
      if (parAlmond->nullSpace) {
        diag += parAlmond->nullSpacePenalty*level->A->null[i]*level->A->null[i];
      }
      level->A->diagInv[i] = 1.0/diag;
    }

    rho = rhoDinvA(parAlmond, level->A, level->A->diagInv);

    if (s == DAMPED_JACOBI) {

      level->smoother_params = (dfloat *) calloc(1,sizeof(dfloat));

      level->smoother_params[0] = (4./3.)/rho;

      //temp storage for smoothing
      if (level->Ncols) level->smootherResidual = (dfloat *) calloc(level->Ncols,sizeof(dfloat));
      if (level->Ncols) level->o_smootherResidual = parAlmond->device.malloc(level->Ncols*sizeof(dfloat),level->smootherResidual);

    } else if (s == CHEBYSHEV) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

      //temp storage for smoothing
      if (level->Ncols) level->smootherResidual = (dfloat *) calloc(level->Ncols,sizeof(dfloat));
      if (level->Ncols) level->smootherResidual2 = (dfloat *) calloc(level->Ncols,sizeof(dfloat));
      if (level->Ncols) level->smootherUpdate   = (dfloat *) calloc(level->Ncols,sizeof(dfloat));
      if (level->Ncols) level->o_smootherResidual  = parAlmond->device.malloc(level->Ncols*sizeof(dfloat),level->smootherResidual);
      if (level->Ncols) level->o_smootherResidual2 = parAlmond->device.malloc(level->Ncols*sizeof(dfloat),level->smootherResidual);
      if (level->Ncols) level->o_smootherUpdate    = parAlmond->device.malloc(level->Ncols*sizeof(dfloat),level->smootherUpdate);
    }
  }
}

extern "C"{
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
  double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
}


static void eig(const int Nrows, double *A, double *WR,
    double *WI){

  if(Nrows){
  int NB  = 256;
  char JOBVL  = 'V';
  char JOBVR  = 'V';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK  = (NB+2)*N;

  double *WORK  = new double[LWORK];
  double *VL  = new double[Nrows*Nrows];
  double *VR  = new double[Nrows*Nrows];

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
    VL, &LDA, VR, &LDA, WORK, &LWORK, &INFO);


  assert(INFO == 0);

  delete [] VL;
  delete [] VR;
  delete [] WORK;
  }
}

dfloat rhoDinvA(parAlmond_t* parAlmond,csr *A, dfloat *invD){

  const dlong N = A->Nrows;
  const dlong M = A->Ncols;

  int k = 10;

  int rank, size;
  rank = agmg::rank;
  size = agmg::size;

  hlong Nlocal = (hlong) N;
  hlong Ntotal = 0;
  MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_HLONG, MPI_SUM, agmg::comm);
  if(k > Ntotal)
    k = (int) Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double *H = (double *) calloc(k*k,sizeof(double));

  // allocate memory for basis
  dfloat **V = (dfloat **) calloc(k+1, sizeof(dfloat *));
  dfloat *Vx = (dfloat *) calloc(M, sizeof(dfloat));

  for(int i=0; i<=k; i++)
    V[i] = (dfloat *) calloc(N, sizeof(dfloat));

  // generate a random vector for initial basis vector
  for (dlong i=0;i<N;i++)
    Vx[i] = (dfloat) drand48();

  dfloat norm_vo = 0.;
  for (dlong i=0;i<N;i++)
    norm_vo += Vx[i]*Vx[i];

  dfloat gNorm_vo = 0;
  MPI_Allreduce(&norm_vo, &gNorm_vo, 1, MPI_DFLOAT, MPI_SUM, agmg::comm);
  gNorm_vo = sqrt(gNorm_vo);

  for (dlong i=0;i<N;i++)
    Vx[i] /= gNorm_vo;

  for (dlong i=0;i<N;i++)
    V[0][i] = Vx[i];

  for(int j=0; j<k; j++){

    for (dlong i=0;i<N;i++)
      Vx[i] = V[j][i];

    // v[j+1] = invD*(A*v[j])
    axpy(A, 1.0, Vx, 0., V[j+1],parAlmond->nullSpace,parAlmond->nullSpacePenalty);

    dotStar(N, invD, V[j+1]);

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij = innerProd(N, V[i], V[j+1]);
      dfloat ghij = 0;
      MPI_Allreduce(&hij, &ghij, 1, MPI_DFLOAT, MPI_SUM, agmg::comm);

      // v[j+1] = v[j+1] - hij*v[i]
      vectorAdd(N,-ghij, V[i], 1.0, V[j+1]);

      H[i + j*k] = (double) ghij;
    }

    if(j+1 < k){

      dfloat norm_vj = 0.;
      for (dlong i=0;i<N;i++)
        norm_vj += V[j+1][i]*V[j+1][i];

      dfloat gNorm_vj;
      MPI_Allreduce(&norm_vj, &gNorm_vj, 1, MPI_DFLOAT, MPI_SUM, agmg::comm);
      gNorm_vj = sqrt(gNorm_vj);

      H[j+1+ j*k] = (double) gNorm_vj;

      scaleVector(N,V[j+1], 1./H[j+1 + j*k]);
    }
  }

  double *WR = (double *) calloc(k,sizeof(double));
  double *WI = (double *) calloc(k,sizeof(double));

  eig(k, H, WR, WI);

  double rho = 0.;

  for(int i=0; i<k; i++){
    double rho_i  = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);

    if(rho < rho_i) {
      rho = rho_i;
    }
  }

  free(H);
  free(WR);
  free(WI);

  // free memory
  for(int i=0; i<=k; i++){
    free(V[i]);
  }

  if ((rank==0)&& (parAlmond->options.compareArgs("VERBOSE","TRUE"))) printf("weight = %g \n", rho);

  return rho;
}

void matrixInverse(int N, dfloat *A);

//set up exact solver using xxt
void setupExactSolve(parAlmond_t *parAlmond, agmgLevel *level, bool nullSpace, dfloat nullSpacePenalty) {

  int rank, size;
  rank = agmg::rank;
  size = agmg::size;

  //copy the global coarse partition as ints
  int *coarseOffsets = (int* ) calloc(size+1,sizeof(int));
  for (int r=0;r<size+1;r++) coarseOffsets[r] = (int) level->globalRowStarts[r];
  
  int  coarseTotal   = coarseOffsets[size];
  int  coarseOffset  = coarseOffsets[rank];

  csr *A = level->A;
  int N = (int) level->Nrows;

  int localNNZ;
  int *rows;
  int *cols;
  dfloat *vals;

  if((rank==0)&&(parAlmond->options.compareArgs("VERBOSE","TRUE"))) printf("Setting up coarse solver...");fflush(stdout);

  if(!nullSpace) {
    //if no nullspace, use sparse A
    localNNZ = (int) (A->diagNNZ+A->offdNNZ);
    
    if (localNNZ) {
      rows = (int *) calloc(localNNZ,sizeof(int));
      cols = (int *) calloc(localNNZ,sizeof(int));
      vals = (dfloat *) calloc(localNNZ,sizeof(dfloat));
    }

    //populate matrix
    int cnt = 0;
    for (int n=0;n<N;n++) {
      int start = (int) A->diagRowStarts[n];
      int end   = (int) A->diagRowStarts[n+1];
      for (int m=start;m<end;m++) {
        rows[cnt] = n + coarseOffset;
        cols[cnt] = (int) (A->diagCols[m] + coarseOffset);
        vals[cnt] = A->diagCoefs[m];
        cnt++;
      }
      start = (int) A->offdRowStarts[n];
      end   = (int) A->offdRowStarts[n+1];
      for (dlong m=A->offdRowStarts[n];m<A->offdRowStarts[n+1];m++) {
        rows[cnt] = n + coarseOffset;
        cols[cnt] = (int) A->colMap[A->offdCols[m]];
        vals[cnt] = A->offdCoefs[m];
        cnt++;
      }
    }
  } else {
    localNNZ = (int) (A->Nrows*coarseTotal); //A is dense due to nullspace augmentation

    if (localNNZ) {
      rows = (int *) calloc(localNNZ,sizeof(int));
      cols = (int *) calloc(localNNZ,sizeof(int));
      vals = (dfloat *) calloc(localNNZ,sizeof(dfloat));
    }

    //gather null vector
    dfloat *nullTotal = (dfloat*) calloc(coarseTotal,sizeof(dfloat));
    int *nullCounts = (int*) calloc(size,sizeof(int));
    for (int r=0;r<size;r++) 
      nullCounts[r] = coarseOffsets[r+1]-coarseOffsets[r];
    
    MPI_Allgatherv(A->null, N, MPI_DFLOAT, nullTotal, nullCounts, coarseOffsets, MPI_DFLOAT, agmg::comm);

    //populate matrix
    for (int n=0;n<N;n++) {
      for (int m=0;m<coarseTotal;m++) {    
        rows[n*coarseTotal+m] = n + coarseOffset;
        cols[n*coarseTotal+m] = m;
        vals[n*coarseTotal+m] = nullSpacePenalty*nullTotal[n+coarseOffset]*nullTotal[m];
      }
    }

    for (int n=0;n<N;n++) {
      int start = (int) A->diagRowStarts[n];
      int end   = (int) A->diagRowStarts[n+1];
      for (int m=start;m<end;m++) {
        int col = (int) (A->diagCols[m] + coarseOffset);
        vals[n*coarseTotal+col] += A->diagCoefs[m];
      }
      start = (int) A->offdRowStarts[n];
      end   = (int) A->offdRowStarts[n+1];
      for (int m=start;m<end;m++) {
        int col = (int) A->colMap[A->offdCols[m]];
        vals[n*coarseTotal+col] += A->offdCoefs[m];
      }
    }
  }

  //ge the nonzero counts from all ranks
  int *NNZ = (int*) calloc(size,sizeof(int));  
  int *NNZoffsets = (int*) calloc(size+1,sizeof(int));  
  MPI_Allgather(&localNNZ, 1, MPI_INT, NNZ, 1, MPI_INT, agmg::comm);

  int totalNNZ = 0;
  for (int r=0;r<size;r++) {
    totalNNZ += NNZ[r];
    NNZoffsets[r+1] = NNZoffsets[r] + NNZ[r];
  }

  int *Arows = (int *) calloc(totalNNZ,sizeof(int));
  int *Acols = (int *) calloc(totalNNZ,sizeof(int));
  dfloat *Avals = (dfloat *) calloc(totalNNZ,sizeof(dfloat));

  MPI_Allgatherv(rows, localNNZ, MPI_INT, Arows, NNZ, NNZoffsets, MPI_INT, agmg::comm);
  MPI_Allgatherv(cols, localNNZ, MPI_INT, Acols, NNZ, NNZoffsets, MPI_INT, agmg::comm);
  MPI_Allgatherv(vals, localNNZ, MPI_DFLOAT, Avals, NNZ, NNZoffsets, MPI_DFLOAT, agmg::comm);

  //assemble the full matrix
  dfloat *coarseA = (dfloat *) calloc(coarseTotal*coarseTotal,sizeof(dfloat));
  for (int i=0;i<totalNNZ;i++) {
    int n = Arows[i];
    int m = Acols[i];
    coarseA[n*coarseTotal+m] = Avals[i];
  }

  matrixInverse(coarseTotal, coarseA);

  //store only the local rows of the full inverse
  parAlmond->invCoarseA = (dfloat *) calloc(A->Nrows*coarseTotal,sizeof(dfloat));
  for (int n=0;n<N;n++) {
    for (int m=0;m<coarseTotal;m++) {
      parAlmond->invCoarseA[n*coarseTotal+m] = coarseA[(n+coarseOffset)*coarseTotal+m];
    } 
  }

  parAlmond->coarseTotal = coarseTotal;
  parAlmond->coarseOffset = coarseOffset;
  parAlmond->coarseOffsets = coarseOffsets;
  parAlmond->coarseCounts = (int*) calloc(size,sizeof(int));
    for (int r=0;r<size;r++) 
      parAlmond->coarseCounts[r] = coarseOffsets[r+1]-coarseOffsets[r];

  parAlmond->xCoarse   = (dfloat*) calloc(coarseTotal,sizeof(dfloat));
  parAlmond->rhsCoarse = (dfloat*) calloc(coarseTotal,sizeof(dfloat));

  if (localNNZ) {
    free(rows);
    free(cols);
    free(vals);
  }

  if (totalNNZ) {
    free(Arows);
    free(Acols);
    free(Avals);
  }

  if(coarseTotal) {
    free(coarseA);
  }

  if((rank==0)&&(parAlmond->options.compareArgs("VERBOSE","TRUE"))) printf("done.\n");
}


void exactCoarseSolve(parAlmond_t *parAlmond, int N, dfloat *rhs, dfloat *x) {

  //gather the full vector
  MPI_Allgatherv(rhs, N, MPI_DFLOAT, parAlmond->rhsCoarse, parAlmond->coarseCounts, parAlmond->coarseOffsets, MPI_DFLOAT, agmg::comm);

  //multiply by local part of the exact matrix inverse
  #pragma omp parallel for
  for (int n=0;n<N;n++) {
    x[n] = 0.;
    for (int m=0;m<parAlmond->coarseTotal;m++) {
      x[n] += parAlmond->invCoarseA[n*parAlmond->coarseTotal+m]*parAlmond->rhsCoarse[m];
    }
  }
}

void device_exactCoarseSolve(parAlmond_t *parAlmond, int N, occa::memory o_rhs, occa::memory o_x) {

  dfloat *rhs = parAlmond->levels[parAlmond->numLevels-1]->rhs;
  dfloat *x = parAlmond->levels[parAlmond->numLevels-1]->x;

  //use coarse solver
  o_rhs.copyTo(rhs);
  //gather the full vector
  MPI_Allgatherv(rhs, N, MPI_DFLOAT, parAlmond->rhsCoarse, parAlmond->coarseCounts, parAlmond->coarseOffsets, MPI_DFLOAT, agmg::comm);

  //multiply by local part of the exact matrix inverse
  #pragma omp parallel for
  for (int n=0;n<N;n++) {
    x[n] = 0.;
    for (int m=0;m<parAlmond->coarseTotal;m++) {
      x[n] += parAlmond->invCoarseA[n*parAlmond->coarseTotal+m]*parAlmond->rhsCoarse[m];
    }
  }

  o_x.copyFrom(x);
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

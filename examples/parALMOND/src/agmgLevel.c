#include "agmg.h"

// parAlmond's function call-backs
void agmgAx(void **args, dfloat *x, dfloat *Ax){
  parAlmond_t *parAlmond = (parAlmond_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  axpy(level->A, 1.0, x, 0.0, Ax,parAlmond->nullSpace,parAlmond->nullSpacePenalty);
}

void agmgCoarsen(void **args, dfloat *r, dfloat *Rr){
  parAlmond_t *parAlmond = (parAlmond_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  axpy(level->R, 1.0, r, 0.0, Rr,false,0.);
}

void agmgProlongate(void **args, dfloat *x, dfloat *Px){
  parAlmond_t *parAlmond = (parAlmond_t *) args[0];
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

    for (iint i=0;i<level->A->Nrows;i++) {
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

      printf("weight = %g \n", level->smoother_params[0]);

      //temp storage for smoothing
      if (level->Ncols) level->smootherResidual = (dfloat *) calloc(level->Ncols,sizeof(dfloat));
      if (level->Ncols) level->o_smootherResidual = parAlmond->device.malloc(level->Ncols*sizeof(dfloat),level->smootherResidual);

    } else if (s == CHEBYSHEV) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

      printf("weight = %g \n", level->smoother_params[0]);

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

dfloat rhoDinvA(parAlmond_t* parAlmond,csr *A, dfloat *invD){

  const iint N = A->Nrows;
  const iint M = A->Ncols;

  int k = 10;

  iint rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint Ntotal=0;
  MPI_Allreduce(&N, &Ntotal, 1, MPI_IINT, MPI_SUM, MPI_COMM_WORLD);
  if(k > Ntotal)
    k = Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double *H = (double *) calloc(k*k,sizeof(double));

  // allocate memory for basis
  dfloat **V = (dfloat **) calloc(k+1, sizeof(dfloat *));
  dfloat *Vx = (dfloat *) calloc(M, sizeof(dfloat));

  for(int i=0; i<=k; i++)
    V[i] = (dfloat *) calloc(N, sizeof(dfloat));

  // generate a random vector for initial basis vector
  for (iint i=0;i<N;i++)
    Vx[i] = (dfloat) drand48();

  dfloat norm_vo = 0.;
  for (iint i=0;i<N;i++)
    norm_vo += Vx[i]*Vx[i];

  dfloat gNorm_vo = 0;
  MPI_Allreduce(&norm_vo, &gNorm_vo, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
  gNorm_vo = sqrt(gNorm_vo);

  for (iint i=0;i<N;i++)
    Vx[i] /= gNorm_vo;

  for (iint i=0;i<N;i++)
    V[0][i] = Vx[i];

  for(int j=0; j<k; j++){

    for (iint i=0;i<N;i++)
      Vx[i] = V[j][i];

    // v[j+1] = invD*(A*v[j])
    axpy(A, 1.0, Vx, 0., V[j+1],parAlmond->nullSpace,parAlmond->nullSpacePenalty);

    dotStar(N, invD, V[j+1]);

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij = innerProd(N, V[i], V[j+1]);
      dfloat ghij = 0;
      MPI_Allreduce(&hij, &ghij, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

      // v[j+1] = v[j+1] - hij*v[i]
      vectorAdd(N,-ghij, V[i], 1.0, V[j+1]);

      H[i + j*k] = (double) ghij;
    }

    if(j+1 < k){

      dfloat norm_vj = 0.;
      for (iint i=0;i<N;i++)
        norm_vj += V[j+1][i]*V[j+1][i];

      dfloat gNorm_vj;
      MPI_Allreduce(&norm_vj, &gNorm_vj, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
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

  return rho;
}

void matrixInverse(int N, dfloat *A);

//set up exact solver using xxt
void setupExactSolve(parAlmond_t *parAlmond, agmgLevel *level, bool nullSpace, dfloat nullSpacePenalty) {

  iint rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint* coarseOffsets = level->globalRowStarts;
  iint coarseTotal = coarseOffsets[size];
  iint coarseOffset = coarseOffsets[rank];

  csr *A = level->A;
  iint N = level->Nrows;

  iint localNNZ;
  iint *rows;
  iint *cols;
  dfloat *vals;

  if(!nullSpace) {
    //if no nullspace, use sparse A
    localNNZ = A->diagNNZ+A->offdNNZ;
    if (localNNZ) {
      rows = (iint *) calloc(localNNZ,sizeof(iint));
      cols = (iint *) calloc(localNNZ,sizeof(iint));
      vals = (dfloat *) calloc(localNNZ,sizeof(dfloat));
    }

    //populate matrix
    int cnt = 0;
    for (iint n=0;n<N;n++) {
      for (iint m=A->diagRowStarts[n];m<A->diagRowStarts[n+1];m++) {
        rows[cnt] = n + coarseOffset;
        cols[cnt] = A->diagCols[m] + coarseOffset;
        vals[cnt] = A->diagCoefs[m];
        cnt++;
      }
      for (iint m=A->offdRowStarts[n];m<A->offdRowStarts[n+1];m++) {
        rows[cnt] = n + coarseOffset;
        cols[cnt] = A->colMap[A->offdCols[m]];
        vals[cnt] = A->offdCoefs[m];
        cnt++;
      }
    }
  } else {
    localNNZ = A->Nrows*coarseTotal; //A is dense due to nullspace augmentation
    if (localNNZ) {
      rows = (iint *) calloc(localNNZ,sizeof(iint));
      cols = (iint *) calloc(localNNZ,sizeof(iint));
      vals = (dfloat *) calloc(localNNZ,sizeof(dfloat));
    }

    //gather null vector
    dfloat *nullTotal = (dfloat*) calloc(coarseTotal,sizeof(dfloat));
    int *nullCounts = (int*) calloc(size,sizeof(int));
    for (int r=0;r<size;r++) 
      nullCounts[r] = coarseOffsets[r+1]-coarseOffsets[r];
    
    MPI_Allgatherv(A->null, A->Nrows, MPI_DFLOAT, nullTotal, nullCounts, coarseOffsets, MPI_DFLOAT, MPI_COMM_WORLD);

    //populate matrix
    for (iint n=0;n<N;n++) {
      for (iint m=0;m<coarseTotal;m++) {    
        rows[n*coarseTotal+m] = n + coarseOffset;
        cols[n*coarseTotal+m] = m;
        vals[n*coarseTotal+m] = nullSpacePenalty*nullTotal[n+coarseOffset]*nullTotal[m];
      }
    }

    for (iint n=0;n<N;n++) {
      for (iint m=A->diagRowStarts[n];m<A->diagRowStarts[n+1];m++) {
        int col = A->diagCols[m] + coarseOffset;
        vals[n*coarseTotal+col] += A->diagCoefs[m];
      }
      for (iint m=A->offdRowStarts[n];m<A->offdRowStarts[n+1];m++) {
        int col = A->colMap[A->offdCols[m]];
        vals[n*coarseTotal+col] += A->offdCoefs[m];
      }
    }
  }

  //ge the nonzero counts from all ranks
  int *NNZ = (int*) calloc(size,sizeof(int));  
  int *NNZoffsets = (int*) calloc(size+1,sizeof(int));  
  MPI_Allgather(&localNNZ, 1, MPI_INT, NNZ, 1, MPI_INT, MPI_COMM_WORLD);

  int totalNNZ = 0;
  for (int r=0;r<size;r++) {
    totalNNZ += NNZ[r];
    NNZoffsets[r+1] = NNZoffsets[r] + NNZ[r];
  }

  int *Arows = (iint *) calloc(totalNNZ,sizeof(iint));
  int *Acols = (iint *) calloc(totalNNZ,sizeof(iint));
  dfloat *Avals = (dfloat *) calloc(totalNNZ,sizeof(dfloat));

  MPI_Allgatherv(rows, localNNZ, MPI_INT, Arows, NNZ, NNZoffsets, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgatherv(cols, localNNZ, MPI_INT, Acols, NNZ, NNZoffsets, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgatherv(vals, localNNZ, MPI_DFLOAT, Avals, NNZ, NNZoffsets, MPI_DFLOAT, MPI_COMM_WORLD);

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
  for (int n=0;n<A->Nrows;n++) {
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

  printf("Done UberCoarse setup\n");
}


void exactCoarseSolve(parAlmond_t *parAlmond, int N, dfloat *rhs, dfloat *x) {

  //gather the full vector
  MPI_Allgatherv(rhs, N, MPI_DFLOAT, parAlmond->rhsCoarse, parAlmond->coarseCounts, parAlmond->coarseOffsets, MPI_DFLOAT, MPI_COMM_WORLD);

  //multiply by local part of the exact matrix inverse
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
  MPI_Allgatherv(rhs, N, MPI_DFLOAT, parAlmond->rhsCoarse, parAlmond->coarseCounts, parAlmond->coarseOffsets, MPI_DFLOAT, MPI_COMM_WORLD);

  //multiply by local part of the exact matrix inverse
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

  iint rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint* coarseOffsets = level->globalRowStarts;
  iint coarseTotal = coarseOffsets[size];
  iint coarseOffset = coarseOffsets[rank];

  iint *globalNumbering = (iint *) calloc(coarseTotal,sizeof(iint));
  for (iint n=0;n<coarseTotal;n++)
    globalNumbering[n] = n;

  csr *A = level->A;
  iint N = level->Nrows;

  iint totalNNZ;
  iint *rows;
  iint *cols;
  dfloat *vals;

  if(!nullSpace) {
    //if no nullspace, use sparse A
    totalNNZ = A->diagNNZ+A->offdNNZ;
    if (totalNNZ) {
      rows = (iint *) calloc(totalNNZ,sizeof(iint));
      cols = (iint *) calloc(totalNNZ,sizeof(iint));
      vals = (dfloat *) calloc(totalNNZ,sizeof(dfloat));
    }

    //populate matrix
    int cnt = 0;
    for (iint n=0;n<N;n++) {
      for (iint m=A->diagRowStarts[n];m<A->diagRowStarts[n+1];m++) {
        rows[cnt] = n + coarseOffset;
        cols[cnt] = A->diagCols[m] + coarseOffset;
        vals[cnt] = A->diagCoefs[m];
        cnt++;
      }
      for (iint m=A->offdRowStarts[n];m<A->offdRowStarts[n+1];m++) {
        rows[cnt] = n + coarseOffset;
        cols[cnt] = A->colMap[A->offdCols[m]];
        vals[cnt] = A->offdCoefs[m];
        cnt++;
      }
    }
  } else {
    totalNNZ = A->Nrows*coarseTotal; //A is dense due to nullspace augmentation
    if (totalNNZ) {
      rows = (iint *) calloc(totalNNZ,sizeof(iint));
      cols = (iint *) calloc(totalNNZ,sizeof(iint));
      vals = (dfloat *) calloc(totalNNZ,sizeof(dfloat));
    }

    //gather null vector
    dfloat *nullTotal = (dfloat*) calloc(coarseTotal,sizeof(dfloat));
    int *nullCounts = (int*) calloc(size,sizeof(int));
    for (int r=0;r<size;r++) 
      nullCounts[r] = coarseOffsets[r+1]-coarseOffsets[r];
    
    MPI_Allgatherv(A->null, A->Nrows, MPI_DFLOAT, nullTotal, nullCounts, coarseOffsets, MPI_DFLOAT, MPI_COMM_WORLD);

    //populate matrix
    for (iint n=0;n<N;n++) {
      for (iint m=0;m<coarseTotal;m++) {    
        rows[n*coarseTotal+m] = n + coarseOffset;
        cols[n*coarseTotal+m] = m;
        vals[n*coarseTotal+m] = nullSpacePenalty*nullTotal[n+coarseOffset]*nullTotal[m];
      }
    }

    for (iint n=0;n<N;n++) {
      for (iint m=A->diagRowStarts[n];m<A->diagRowStarts[n+1];m++) {
        int col = A->diagCols[m] + coarseOffset;
        vals[n*coarseTotal+col] += A->diagCoefs[m];
      }
      for (iint m=A->offdRowStarts[n];m<A->offdRowStarts[n+1];m++) {
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
                                iintString,
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
  for (iint n=0;n<parAlmond->coarseTotal;n++)
    parAlmond->rhsCoarse[n] =0.;

  for (iint n=0;n<N;n++)
    parAlmond->rhsCoarse[n+parAlmond->coarseOffset] = rhs[n];

  xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);

  for (iint n=0;n<N;n++)
    x[n] = parAlmond->xCoarse[n+parAlmond->coarseOffset];

}

void device_exactCoarseSolve(parAlmond_t *parAlmond, int N, occa::memory o_rhs, occa::memory o_x) {

  //use coarse solver
  for (iint n=0;n<parAlmond->coarseTotal;n++)
    parAlmond->rhsCoarse[n] =0.;

  o_rhs.copyTo(parAlmond->rhsCoarse+parAlmond->coarseOffset);
  xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);
  o_x.copyFrom(parAlmond->xCoarse+parAlmond->coarseOffset,N*sizeof(dfloat));
}
#endif
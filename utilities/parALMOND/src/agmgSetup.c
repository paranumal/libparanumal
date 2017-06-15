#include "parAlmond.h"

parAlmond_t * agmgSetup(csr *A, dfloat *nullA, iint *globalRowStarts, const char* options){
  iint rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  parAlmond_t *parAlmond = (parAlmond_t *) calloc(1,sizeof(parAlmond_t));

  const iint coarseSize = 10;

  double seed = 1.0;//MPI_Wtime();
  double gSeed;
  MPI_Allreduce(&seed, &gSeed, 1, MPI_LONG, MPI_BXOR, MPI_COMM_WORLD);
  srand48(gSeed);

  agmgLevel **levels = (agmgLevel **) calloc(MAX_LEVELS,sizeof(agmgLevel *));

  levels[0] = (agmgLevel *) calloc(1,sizeof(agmgLevel));

  //copy A matrix
  levels[0]->A = (csr *) calloc(1,sizeof(csr));

  levels[0]->A->Nrows = A->Nrows;
  levels[0]->A->Ncols = A->Ncols;
  levels[0]->A->nnz = A->nnz;

  levels[0]->A->rowStarts = (iint *) calloc(A->Nrows+1,sizeof(iint));
  levels[0]->A->cols      = (iint *) calloc(A->nnz, sizeof(iint));
  levels[0]->A->coefs     = (dfloat *) calloc(A->nnz, sizeof(dfloat));
  for (iint i=0; i<A->Nrows+1; i++)
    levels[0]->A->rowStarts[i] = A->rowStarts[i];

  for (iint i=0; i<A->nnz; i++) {
    levels[0]->A->cols[i] = A->cols[i];   
    levels[0]->A->coefs[i] = A->coefs[i];   
  }
  levels[0]->A->NsendTotal = A->NsendTotal;
  levels[0]->A->NrecvTotal = A->NrecvTotal;
  levels[0]->A->NHalo      = A->NHalo;

  //set up level size
  levels[0]->nullA = nullA;
  levels[0]->Nrows = A->Nrows;
  levels[0]->Ncols = A->Ncols;

  if (globalRowStarts) {
    levels[0]->globalRowStarts = (iint *) calloc(size+1,sizeof(iint));
    for (iint r=0;r<size+1;r++) 
      levels[0]->globalRowStarts[r] = globalRowStarts[r];
  }

  int numLevels = 1;
  int lev =0;

  bool done = false;
  while(!done){
    const iint dim = levels[lev]->A->Nrows;
    csr *coarseA = (csr *) calloc(1,sizeof(csr));
    dfloat *nullCoarseA;

    coarsen(levels[lev], &coarseA, &nullCoarseA); 

    const iint coarseDim = coarseA->Nrows;

    // allocate vectors required
    //allocate(levels[lev]);

    SmoothType s = DAMPED_JACOBI;
    //SmoothType s = JACOBI;

    setup_smoother(levels[lev], s);

    numLevels++;

    levels[lev+1] = (agmgLevel *) calloc(1,sizeof(agmgLevel));
    levels[lev+1]->A = coarseA;
    levels[lev+1]->nullA = nullCoarseA;
    levels[lev+1]->Nrows = coarseA->Nrows;
    levels[lev+1]->Ncols = coarseA->Ncols;

    if (globalRowStarts) {
      levels[lev+1]->globalRowStarts = (iint *) calloc(size+1,sizeof(iint));

      //figure out global partitioning for this level
      iint chunk = coarseA->Nrows/size;
      iint remainder = coarseA->Nrows - chunk*size;

      for (iint r=0;r<size+1;r++)
        if (globalRowStarts)
          levels[lev+1]->globalRowStarts[r] = r*chunk + (r<remainder ? r : remainder);
    }

    if(coarseA->Nrows <= coarseSize || dim < 2*coarseDim){
      //allocate(levels[lev+1]);
      setup_smoother(levels[lev+1],JACOBI);
      break;
    }
    lev++;
  }

  parAlmond->ktype = PCG;


  //Now that AGMG is setup, distribute the operators between the processors and set up the halo
  if (globalRowStarts) {
    for (int n=0;n<numLevels-1;n++) {

      levels[n]->A = distribute(levels[n]->A,
                                    levels[n]->globalRowStarts,
                                    levels[n]->globalRowStarts);
      levels[n]->P = distribute(levels[n]->P,
                                    levels[n]->globalRowStarts,
                                    levels[n+1]->globalRowStarts);
      levels[n]->R = distribute(levels[n]->R,
                                    levels[n+1]->globalRowStarts,
                                    levels[n]->globalRowStarts);
      
      iint M    = levels[n]->A->Nrows;
      iint Nmax = levels[n]->A->Ncols;

      Nmax = levels[n]->R->Ncols > Nmax ? levels[n]->R->Ncols : Nmax;
      if (n>0) Nmax = levels[n-1]->P->Ncols > Nmax ? levels[n-1]->P->Ncols : Nmax;

      levels[n]->Nrows = M;
      levels[n]->Ncols = Nmax;
    }
    levels[numLevels-1]->A = distribute(levels[numLevels-1]->A,
                                  levels[numLevels-1]->globalRowStarts,
                                  levels[numLevels-1]->globalRowStarts);

    iint M    = levels[numLevels-1]->A->Nrows;
    iint Nmax = levels[numLevels-1]->A->Ncols;

    if (numLevels>1) Nmax = levels[numLevels-2]->P->Ncols > Nmax ? levels[numLevels-2]->P->Ncols : Nmax;

    levels[numLevels-1]->Nrows = M;
    levels[numLevels-1]->Ncols = Nmax;
  }

  //allocate vectors required
  for (int n=0;n<numLevels;n++) {
    iint M = levels[n]->Nrows;
    iint N = levels[n]->Ncols;

    if ((n>0)&&(n<numLevels-1)) { //kcycle vectors
      levels[n]->ckp1 = (dfloat *) calloc(N,sizeof(dfloat)); 
      levels[n]->vkp1 = (dfloat *) calloc(M,sizeof(dfloat)); 
      levels[n]->wkp1 = (dfloat *) calloc(M,sizeof(dfloat));
    }
    levels[n]->x    = (dfloat *) calloc(N,sizeof(dfloat));
    levels[n]->rhs  = (dfloat *) calloc(M,sizeof(dfloat));
    levels[n]->res  = (dfloat *) calloc(N,sizeof(dfloat));
  }

  //set up base solver using xxt
  if (strstr(options,"UBERGRID")) {
    iint N = levels[numLevels-1]->Nrows;

    iint* coarseN       = (iint *) calloc(size,sizeof(iint));
    iint* coarseOffsets = (iint *) calloc(size+1,sizeof(iint));

    MPI_Allgather(&N, 1, MPI_IINT, coarseN, 1, MPI_IINT, MPI_COMM_WORLD);

    coarseOffsets[0] = 0;
    for (iint r=0;r<size;r++)
      coarseOffsets[r+1] = coarseOffsets[r] + coarseN[r];

    iint coarseTotal = coarseOffsets[size];
    iint coarseOffset = coarseOffsets[rank];

    iint *globalNumbering = (iint *) calloc(coarseTotal,sizeof(iint));
    for (iint n=0;n<coarseTotal;n++)
      globalNumbering[n] = n;

    iint nnz = levels[numLevels-1]->A->nnz;
    iint *rows;
    iint *cols;
    dfloat *vals;
    if (nnz) {
      rows = (iint *) calloc(nnz,sizeof(iint));
      cols = (iint *) calloc(nnz,sizeof(iint));
      vals = (dfloat *) calloc(nnz,sizeof(dfloat));
    }

    //populate A matrix
    for (iint n=0;n<N;n++) {
      for (iint m=levels[numLevels-1]->A->rowStarts[n];
               m<levels[numLevels-1]->A->rowStarts[n+1];m++) {
        rows[m]  = n + parAlmond->coarseOffset;
        iint col = levels[numLevels-1]->A->cols[m];
        cols[m]  = levels[numLevels-1]->A->colMap[col];
        vals[m]  = levels[numLevels-1]->A->coefs[m];
      }
    }

    // need to create numbering for really coarse grid on each process for xxt
    parAlmond->Acoarse = xxtSetup(coarseTotal,
                                  globalNumbering,
                                  nnz,
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

    free(coarseN);
    free(coarseOffsets);
    free(globalNumbering);
    if (nnz) {
      free(rows);
      free(cols);
      free(vals);
    }

    printf("Done UberCoarse setup\n"); 
  }

  parAlmond->levels = levels;
  parAlmond->numLevels = numLevels;

  return parAlmond;
}


void sync_setup_on_device(parAlmond_t *parAlmond, occa::device dev){
  //set occa device pointer
  parAlmond->device = dev;
  buildAlmondKernels(parAlmond);

  for(int i=0; i<parAlmond->numLevels; i++){
    iint N = parAlmond->levels[i]->Ncols;
    iint M = parAlmond->levels[i]->Nrows;
    
    parAlmond->levels[i]->deviceA = newHYB(parAlmond, parAlmond->levels[i]->A);
    if (i < parAlmond->numLevels-1) {
      parAlmond->levels[i]->dcsrP   = newDCSR(parAlmond, parAlmond->levels[i]->P);
      parAlmond->levels[i]->deviceR = newHYB(parAlmond, parAlmond->levels[i]->R);
    }

    if (N) parAlmond->levels[i]->o_x   = parAlmond->device.malloc(N*sizeof(dfloat), parAlmond->levels[i]->x);
    if (N) parAlmond->levels[i]->o_res = parAlmond->device.malloc(N*sizeof(dfloat), parAlmond->levels[i]->res);
    if (M) parAlmond->levels[i]->o_rhs = parAlmond->device.malloc(M*sizeof(dfloat), parAlmond->levels[i]->rhs);

    if(i > 0){
      if (N) parAlmond->levels[i]->o_ckp1 = parAlmond->device.malloc(N*sizeof(dfloat), parAlmond->levels[i]->x);
      if (M) parAlmond->levels[i]->o_wkp1 = parAlmond->device.malloc(M*sizeof(dfloat), parAlmond->levels[i]->x);
      if (M) parAlmond->levels[i]->o_vkp1 = parAlmond->device.malloc(M*sizeof(dfloat), parAlmond->levels[i]->x);
    }
  }

  //buffer for innerproducts in kcycle
  dfloat dummy[3];
  parAlmond->o_rho  = parAlmond->device.malloc(3*sizeof(dfloat), dummy);

  //if using matrix-free A action, free unnecessary buffers
  if (strstr(parAlmond->options,"MATRIXFREE")) {
    parAlmond->levels[0]->deviceA->E->o_cols.free();
    parAlmond->levels[0]->deviceA->E->o_coefs.free();
    if (parAlmond->levels[0]->deviceA->C->nnz) {
      parAlmond->levels[0]->deviceA->C->o_offsets.free();
      parAlmond->levels[0]->deviceA->C->o_cols.free();
      parAlmond->levels[0]->deviceA->C->o_coefs.free();
    }
    if(parAlmond->levels[0]->deviceA->NsendTotal) {
      parAlmond->levels[0]->deviceA->o_haloElementList.free();
      parAlmond->levels[0]->deviceA->o_haloBuffer.free();
    }
  }
}
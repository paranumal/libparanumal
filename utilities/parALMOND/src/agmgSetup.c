#include "parAlmond.h"

parAlmond_t * agmgSetup(csr *A, dfloat *nullA, iint *globalRowStarts, const char* options){
  iint rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  parAlmond_t *parAlmond = (parAlmond_t *) calloc(1,sizeof(parAlmond_t));

  // approximate Nrows at coarsest level
  const iint coarseSize = 10;

  double seed = (double) rank;//MPI_Wtime();
  srand48(seed);

  agmgLevel **levels = (agmgLevel **) calloc(MAX_LEVELS,sizeof(agmgLevel *));

  levels[0] = (agmgLevel *) calloc(1,sizeof(agmgLevel));

  //copy A matrix and null vector
  levels[0]->A = A;
  levels[0]->nullA = nullA;

  //set up level size
  levels[0]->Nrows = A->Nrows;
  levels[0]->Ncols = A->Ncols;

  //copy global partiton
  levels[0]->globalRowStarts = (iint *) calloc(size+1,sizeof(iint));
  for (iint r=0;r<size+1;r++)
      levels[0]->globalRowStarts[r] = globalRowStarts[r];

  int numLevels = 1;
  int lev =0;

  bool done = false;
  while(!done){
    const iint dim = levels[lev]->A->Nrows;
    csr *coarseA = (csr *) calloc(1,sizeof(csr));
    dfloat *nullCoarseA;

    coarsen(levels[lev], &coarseA, &nullCoarseA);

    const iint coarseDim = coarseA->Nrows;

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

//create coarsened problem
void coarsen(agmgLevel *level, csr **coarseA, dfloat **nullCoarseA){

  // establish the graph of strong connections
  level->threshold = 0.5;

  csr *C = strong_graph(level->A, level->threshold);

  iint *FineToCoarse = form_aggregates(level, C);

  //find_aggregate_owners(level,FineToCoarse);

  construct_interpolator(level, FineToCoarse, nullCoarseA);

  *coarseA = galerkinProd(level);
}

csr * strong_graph(csr *A, dfloat threshold){

  const iint N = A->Nrows;
  const iint M = A->Ncols;

  csr *C = (csr *) calloc(1, sizeof(csr));

  C->Nrows = N;
  C->Ncols = M;

  C->diagRowStarts = (iint *) calloc(N+1,sizeof(iint));
  if(A->offdNNZ) C->offdRowStarts = (iint *) calloc(N+1,sizeof(iint));

  dfloat *maxOD;
  if (N) maxOD = (dfloat *) calloc(N,sizeof(dfloat));

  //store the diagonal of A for all needed columns
  dfloat *diagA = (dfloat *) calloc(M,sizeof(dfloat));
  for (iint i=0;i<N;i++)
    diagA[i] = A->diagCoefs[A->diagRowStarts[i]];
  crsHaloExchange(A, sizeof(dfloat), diagA, A->sendBuffer, diagA);

  for(iint i=0; i<N; i++){
    dfloat sign = (diagA[i] >= 0) ? 1:-1;
    dfloat Aii = fabs(diagA[i]]);

    //find maxOD
    //local entries
    iint Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];
    for(iint jj= Jstart+1; jj<Jend; jj++){
      iint col = A->diagCols[jj];
      dfloat Ajj = fabs(diagA[col]]);
      dfloat OD = -sign*A->diagCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > maxOD[i]) maxOD[i] = OD;
    }
    //non-local entries
    Jstart = A->offdRowStarts[i], Jend = A->offdRowStarts[i+1];
    for(iint jj= Jstart; jj<Jend; jj++){
      iint col = A->offdCols[jj];
      dfloat Ajj = fabs(diagA[col]]);
      dfloat OD = -sign*A->offdCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > maxOD[i]) maxOD[i] = OD;
    }

    iint diag_strong_per_row = 1; // diagonal entry
    //local entries
    Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];
    for(iint jj = Jstart+1; jj<Jend; jj++){
      iint col = A->diagCols[jj];
      dfloat Ajj = fabs(diagA[col]]);
      dfloat OD = -sign*A->diagCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD) diag_strong_per_row++;
    }
    iint offd_strong_per_row = 0;
    //non-local entries
    Jstart = A->offdRowStarts[i], Jend = A->offdRowStarts[i+1];
    for(iint jj= Jstart; jj<Jend; jj++){
      iint col = A->offdCols[jj];
      dfloat Ajj = fabs(diagA[col]]);
      dfloat OD = -sign*A->offdCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD) offd_strong_per_row++;
    }

    C->diagRowStarts[i+1] = diag_strong_per_row;
    C->offdRowStarts[i+1] = offd_strong_per_row;
  }

  // cumulative sum
  for(iint i=1; i<N+1 ; i++) {
    C->diagRowStarts[i] += C->diagRowStarts[i-1];
    C->offdRowStarts[i] += C->offdRowStarts[i-1];
  }

  C->diagNNZ = C->diagRowStarts[m];
  C->offdNNZ = C->offdRowStarts[m];

  if (C->diagNNZ) C->diagCols = (iint *) calloc(C->diagNNZ, sizeof(iint));
  if (C->offdNNZ) C->offdCols = (iint *) calloc(C->offdNNZ, sizeof(iint));

  // fill in the columns for strong connections
  for(iint i=0; i<N; i++){
    dfloat sign = (diagA[i] >= 0) ? 1:-1;
    dfloat Aii = fabs(diagA[i]);

    iint diagCounter = C->diagRowStarts[i];
    iint offdCounter = C->offdRowStarts[i];

    //local entries
    C->diagCols[diagCounter++] = i;// diag entry
    iint Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];
    for(jj = Jstart+1; jj<Jend; jj++){
      iint col = A->diagCols[jj];
      dfloat Ajj = fabs(diagA[col]]);
      dfloat OD = -sign*A->diagCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD[i])
        C->diagCols[diagCounter++] = A->cols[jj];
    }
    Jstart = A->offdRowStarts[i], Jend = A->offdRowStarts[i+1];
    for(jj = Jstart; jj<Jend; jj++){
      iint col = A->offdCols[jj];
      dfloat Ajj = fabs(diagA[col]]);
      dfloat OD = -sign*A->offdCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD[i])
        C->offdCols[offdCounter++] = A->cols[jj];
    }
  }
  if(N) free(maxOD);

  return C;
}

bool customLess(iint smax, dfloat rmax, iint imax, iint s, dfloat r, iint i){

  if(s > smax) return true;
  if(smax > s) return false;

  if(r > rmax) return true;
  if(rmax > r) return false;

  if(i > imax) return true;
  if(i < imax) return false;

  return false;
}

iint * form_aggregates(agmgLevel *level, csr *C){

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const iint N   = C->Nrows;
  const iint M   = C->Ncols;
  const iint diagNNZ = C->diagNNZ;
  const iint offdNNZ = C->offdNNZ;

  iint *FineToCoarse = (iint *) calloc(M, sizeof(iint));
  for (iint i =0;i<M;i++) FineToCoarse[i] = -1;

  dfloat *rands  = (dfloat *) calloc(M, sizeof(dfloat));
  iint   *states = (iint *)   calloc(M, sizeof(iint));

  dfloat *Tr = (dfloat *) calloc(M, sizeof(dfloat));
  iint   *Ts = (iint *)   calloc(M, sizeof(iint));
  iint   *Ti = (iint *)   calloc(M, sizeof(iint));
  iint   *Tc = (iint *)   calloc(M, sizeof(iint));

  for(iint i=0; i<M; i++){
    rands[i] = (dfloat) drand48();
    states[i] = 0;
  }

  // add the number of strong connections
  for(iint i=0; i<diagNNZ; i++)
    rands[C->diagCols[i]] += 1.;
  for(iint i=0; i<offdNNZ; i++)
    rands[C->offdCols[i]] += 1.;

  csr *A = level->A;
  iint *globalRowStarts = level->globalRowStarts;

  iint *iintSendBuffer;
  dfloat *dfloatSendBuffer;
  if (level->A->NsendTotal) {
    iintSendBuffer = (iint *) calloc(A->NsendTotal,sizeof(iint));
    dfloatSendBuffer = (dfloat *) calloc(A->NsendTotal,sizeof(dfloat));
  }
  //share randomizer values
  crsHaloExchange(A, sizeof(dfloat), rands, dfloatSendBuffer, rands);

  bool done = false;
  while(!done){
    // first neighbours
    for(iint i=0; i<N; i++){

      iint smax = states[i];
      dfloat rmax = rands[i];
      iint imax = i + globalRowStarts[rank];

      if(smax != 1){
        //local entries
        for(iint jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
          const iint col = C->diagCols[jj];
          if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
            smax = states[col];
            rmax = rands[col];
            imax = col + globalRowStarts[rank];
          }
        }
        //nonlocal entries
        for(iint jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
          const iint col = C->offdCols[jj];
          if(customLess(smax, rmax, imax, states[col], rands[col], A->offdColMap[col])) {
            smax = states[col];
            rmax = rands[col];
            imax = A->offdColMap[col];
          }
        }
      }
      Ts[i] = smax;
      Tr[i] = rmax;
      Ti[i] = imax;
    }

    //share results
    crsHaloExchange(A, sizeof(dfloat), Tr, dfloatSendBuffer, Tr);
    crsHaloExchange(A, sizeof(iint), Ts, iintSendBuffer, Ts);
    crsHaloExchange(A, sizeof(iint), Ti, iintSendBuffer, Ti);

    // second neighbours
    for(iint i=0; i<N; i++){
      iint   smax = Ts[i];
      dfloat rmax = Tr[i];
      iint   imax = Ti[i];

      //local entries
      for(iint jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
        const iint col = C->diagCols[jj];
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }
      //nonlocal entries
      for(iint jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
        const iint col = C->offdCols[jj];
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }

      // if I am the strongest among all the 1 and 2 ring neighbours
      // I am an MIS node
      if((states[i] == 0) && (imax == (i + globalRowStarts[rank])))
        states[i] = 1;

      // if there is an MIS node within distance 2, I am removed
      if((states[i] == 0) && (smax == 1))
        states[i] = -1;
    }

    crsHaloExchange(A, sizeof(iint), Ts, iintSendBuffer, Ts);

    // if number of undecided nodes = 0, algorithm terminates
    iint cnt = std::count(states, states+N, 0);
    MPI_Allreduce(&cnt,&done,1,MPI_IINT, MPI_SUM,MPI_COMM_WORLD);
    done = (done == 0);
  }

  level->numAggs = 0;
  level->globalAggOffset = (iint *) calloc(size+1,sizeof(iint));
  // count the coarse nodes/aggregates
  for(iint i=0; i<N; i++)
    if(states[i] == 1) level->numAggs++;

  MPI_Allgather(&(level->numAggs),1,MPI_IINT,level->globalAggOffset,1,MPI_IINT,MPI_COMM_WORLD);

  for (iint r=0;r<size;r++)
    level->globalAggOffset[r+1] += level->globalAggOffset[r];
  level->gNumAggs = level->globalAggOffset[size];

  level->numAggs = 0;
  // enumerate the coarse nodes/aggregates
  for(iint i=0; i<N; i++)
    if(states[i] == 1)
      FineToCoarse[i] = aggOffset[rank] + level->numAggs++;

  //share the initial aggregate flags
  crsHaloExchange(A, sizeof(iint), FineToCoarse, iintSendBuffer, FineToCoarse);

  // form the aggregates
  for(iint i=0; i<N; i++){
    iint   smax = states[i];
    dfloat rmax = rands[i];
    iint   imax = i + globalRowStarts[rank];
    iint   cmax = FineToCoarse[i];

    if(smax != 1){
      //local entries
      for(iint jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
        const iint col = C->diagCols[jj];
        if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
          smax = states[col];
          rmax = rands[col];
          imax = col + globalRowStarts[rank];
          cmax = FineToCoarse[col];
        }
      }
      //nonlocal entries
      for(iint jj=C->offdRowStarts[i]+1;jj<C->offdRowStarts[i+1];jj++){
        const iint col = C->offdCols[jj];
        if(customLess(smax, rmax, imax, states[col], rands[col], A->offdColMap[col])){
          smax = states[col];
          rmax = rands[col];
          imax = A->offdColMap[col];
          cmax = FineToCoarse[col];
        }
      }
    }
    Ts[i] = smax;
    Tr[i] = rmax;
    Ti[i] = imax;
    Tc[i] = cmax;

    if((states[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  crsHaloExchange(A, sizeof(iint), FineToCoarse, iintSendBuffer, FineToCoarse);
  crsHaloExchange(A, sizeof(dfloat), Tr, dfloatSendBuffer, Tr);
  crsHaloExchange(A, sizeof(iint), Ts, iintSendBuffer, Ts);
  crsHaloExchange(A, sizeof(iint), Ti, iintSendBuffer, Ti);
  crsHaloExchange(A, sizeof(iint), Tc, iintSendBuffer, Tc);

  // second neighbours
  for(iint i=0; i<N; i++){
    iint   smax = Ts[i];
    dfloat rmax = Tr[i];
    iint   imax = Ti[i];
    iint   cmax = Tc[i];

    //local entries
    for(iint jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
      const iint col = C->diagCols[jj];
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }
    //nonlocal entries
    for(iint jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
      const iint col = C->offdCols[jj];
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }

    if((states[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  free(rands);
  free(states);
  free(Tr);
  free(Ts);
  free(Ti);
  free(Tc);
  if (level->A->NsendTotal) {
    free(iintSendBuffer);
    free(dfloatSendBuffer);
  }

  //TODO maybe free C here?

  return FineToCoarse;
}

struct key_value_pair1{
  long key;
  long value;
};

int compare_key1(const void *a, const void *b){
  struct key_value_pair1 *pa = (struct key_value_pair1 *) a;
  struct key_value_pair1 *pb = (struct key_value_pair1 *) b;

  if (pa->key < pb->key) return -1;
  if (pa->key > pb->key) return +1;

  return 0;
};

typedef struct {

  iint fineId;
  iint coarseId;
  iint newCoarseId;

  iint orginRank;
  iint ownerRank;

} parallelAggregate_t;

int compareOwner(const void *a, const void *b){
  parallelAggregate_t *pa = (parallelAggregate_t *) a;
  parallelAggregate_t *pb = (parallelAggregate_t *) b;

  if (pa->ownerRank < pb->ownerRank) return -1;
  if (pa->ownerRank > pb->ownerRank) return +1;

  return 0;
};

int compareAgg(const void *a, const void *b){
  parallelAggregate_t *pa = (parallelAggregate_t *) a;
  parallelAggregate_t *pb = (parallelAggregate_t *) b;

  if (pa->coarseId < pb->coarseId) return -1;
  if (pa->coarseId > pb->coarseId) return +1;

  if (pa->orginRank < pb->orginRank) return -1;
  if (pa->orginRank > pb->orginRank) return +1;

  return 0;
};

void find_aggregate_owners(agmgLevel *level, iint* FineToCoarse) {
  // MPI info
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //Need to establish 'ownership' of aggregates
  //populate aggregate array
  iint gNumAggs = level->gNumAggs; //total number of aggregates
  parallelAggregate_t *sendAggs = (parallelAggregate_t *) calloc(m,sizeof(parallelAggregate_t));
  for (iint i=0;i<m;i++) {
    sendAggs[i].localId = i;
    sendAggs[i].orginRank = rank;

    sendAggs[i].coarseId = FineToCoarse[i];
    //set a temporary owner. Evenly distibute aggregates amoungst ranks
    sendAggs.ownerRank[i] = (FineToCoarse[i]*size)/gNumAggs;
  }

  //sort by owning rank for all_reduce
  qsort(sendAggs, m, sizeof(parallelAggregate_t), compareOwner);

  iint *sendCounts = (iint *) calloc(size,sizeof(iint));
  iint *recvCounts = (iint *) calloc(size,sizeof(iint));
  iint *sendOffsets = (iint *) calloc(size+1,sizeof(iint));
  iint *recvOffsets = (iint *) calloc(size+1,sizeof(iint));

  for(iint i=0;i<m;++n)
    sendCounts[sendAggs[i].ownerRank] += sizeof(parallelAggregate_t);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_IINT, recvCounts, 1, MPI_IINT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  iint recvNtotal = 0;
  for(iint r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r]/sizeof(parallelAggregate_t);
  }
  parallelAggregate_t *recvAggs = (parallelAggregate_t *) calloc(recvNtotal,sizeof(parallelAggregate_t));

  MPI_Alltoallv(sendAggs, sendCounts, sendOffsets, MPI_CHAR,
                recvAggs, recvCounts, recvOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  //sort by coarse aggregate number, and then by original rank
  qsort(recvAggs, recvNtotal, sizeof(parallelAggregate_t), compareAgg);

  //count the number of unique aggregates here
  iint numUniqueAggs =0;
  if (recvNtotal) numAggregates++;
  for (iint i=1;i<recvNtotal;i++)
    if(recvAggs[i].coarseId!=recvAggs[i-1].coarseId) numUniqueAggs++;

  //get their locations in the array
  iint *aggStarts;
  if (NumUniqueAggs)
    aggStarts = (iint *) calloc(NumUniqueAggs+1,sizeof(iint));
  iint cnt = 1;
  for (iint i=1;i<recvNtotal;i++)
    if(recvAggs[i].coarseId!=recvAggs[i-1].coarseId) aggStarts[cnt++]=i;
  aggStarts[NumUniqueAggs] = recvNtotal;

  //use a random dfloat for each rank to break ties.
  dfloat rand = (dfloat) drand48();
  dfloat *gRands = (dfloat *) calloc(size,sizeof(dfloat));
  MPI_Allgather(&rand, 1, MPI_DFLOAT, gRands, 1, MPI_DFLOAT, MPI_COMM_WORLD);

  //determine the aggregates majority owner
  dfloat *rankCounts = (dfloat *) calloc(size,sizeof(dfloat));
  for (iint n=0;n<NumUniqueAggs) {
    //populate randomizer
    for (iint r=0;r>size;r++)
      rankCounts[r] = gRands[r];

    //count the number of contributions to the aggregate from the separate ranks
    for (iint i=aggStarts[n];i<aggStarts[n+1];i++)
      rankCounts[recvAggs[i].orginRank]++;

    //find which rank is contributing the most to this aggregate
    iint owningRank = 0;
    dfloat maxEntries = rankCounts[0];
    for (iint r=1;r<size;r++) {
      if (rankCounts[r]>maxEntries) {
        owningRank = r;
        maxEntries = rankCounts[r];
      }
    }

    //set this aggregate's owner
    for (iint i=aggStarts[n];i<aggStarts[n+1];i++)
      recvAggs[i].owningRank = owningRank;
  }
  free(gRands); free(rankCounts);
  free(aggStarts);

  //sort by owning rank
  qsort(recvAggs, recvNtotal, sizeof(parallelAggregate_t), compareOwner);

  iint *newSendCounts = (iint *) calloc(size,sizeof(iint));
  iint *newRecvCounts = (iint *) calloc(size,sizeof(iint));
  iint *newSendOffsets = (iint *) calloc(size+1,sizeof(iint));
  iint *newRecvOffsets = (iint *) calloc(size+1,sizeof(iint));

  for(iint i=0;i<m;++n)
    newSendCounts[sendAggs[i].ownerRank] += sizeof(parallelAggregate_t);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(newSendCounts, 1, MPI_IINT, newRecvCounts, 1, MPI_IINT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  iint newRecvNtotal = 0;
  for(iint r=0;r<size;++r){
    newSendOffsets[r+1] = newSendOffsets[r] + newSendCounts[r];
    newRecvOffsets[r+1] = newRecvOffsets[r] + newRecvCounts[r];
    newRecvNtotal += newRecvCounts[r]/sizeof(parallelAggregate_t);
  }
  parallelAggregate_t *newRecvAggs = (parallelAggregate_t *) calloc(newRecvNtotal,sizeof(parallelAggregate_t));

  MPI_Alltoallv(   recvAggs, newSendCounts, newSendOffsets, MPI_CHAR,
                newRecvAggs, newRecvCounts, newRecvOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  //sort by coarse aggregate number, and then by original rank
  qsort(newRecvAggs, newRecvNtotal, sizeof(parallelAggregate_t), compareAgg);

  //count the number of unique aggregates this rank owns
  level->numAggs = 0;
  if (newRecvNtotal) level->numAggs++;
  for (iint i=1;i<newRecvNtotal;i++)
    if(newRecvAggs[i].coarseId!=newRecvAggs[i-1].coarseId) level->numAggs++;

  //determine a global numbering of the aggregates
  MPI_Allgather(&(level->numAggs), 1, MPI_IINT, level->globalAggOffset, 1, MPI_IINT, MPI_COMM_WORLD);

  for (iint r=0;r<size;r++)
    level->globalAggOffset[r+1] += level->globalAggOffset[r];
  level->gNumAggs = level->globalAggOffset[size];

  //set the new global coarse index
  iint cnt = level->globalAggOffset[rank];
  if (newRecvNtotal) newRecvAggs[0].newCoarseId = cnt;
  for (iint i=1;i<newRecvNtotal;i++) {
    if(newRecvAggs[i].coarseId!=newRecvAggs[i-1].coarseId) cnt++;

    newRecvAggs[i].newCoarseId = cnt;
  }

  //send the aggregate data back
  MPI_Alltoallv(newRecvAggs, newRecvCounts, newRecvOffsets, MPI_CHAR,
                   recvAggs, newSendCounts, newSendOffsets, MPI_CHAR,
                MPI_COMM_WORLD);
  MPI_Alltoallv(recvAggs, recvCounts, recvOffsets, MPI_CHAR,
                sendAggs, sendCounts, sendOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  free(recvAggs);
  free(sendCounts);  free(recvCounts);
  free(sendOffsets); free(recvOffsets);
  free(newRecvAggs);
  free(newSendCounts);  free(newRecvCounts);
  free(newSendOffsets); free(newRecvOffsets);

  //record the new FineToCoarse map
  for (iint i=0;i<m;i++)
    FineToCoarse[sendAggs[i].localId] = sendAggs[i].newCoarseId;
}


void construct_interpolator(agmgLevel *level, iint *FineToCoarse, dfloat **nullCoarseA){
  // MPI info
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const iint N = level->A->Nrows;
  const iint M = level->A->Ncols;
  const iint NCoarse = level->numAggs; //local num agg

  iint *globalAggOffset = level->globalAggOffset;
  const iint globalAggStart = level->globalAggOffset[rank];

  csr* P = (csr *) calloc(1, sizeof(csr));

  P->Nrows = N;
  P->Ncols = NCoarse;
  P->Nhalo = 0;

  P->diagRowStarts = (iint *) calloc(N+1, sizeof(iint));
  P->offdRowStarts = (iint *) calloc(N+1, sizeof(iint));

  // each row has exactly one nonzero per row
  P->diagNNZ =0;
  P->offdNNZ =0;
  for(iint i=0; i<N; i++) {
    iint col = FineToCoarse[i];
    if ((col>globalAggOffset[rank]-1)||(col<globalAggOffset[rank+1])) {
      P->diagNNZ++;
      P->diagRowStarts[i]++;
    } else {
      P->offdNNZ++;
      P->offdRowStarts[i]++;
    }
  }
  for(iint i=0; i<N; i++) {
    P->diagRowStarts[i+1] += P->diagRowStarts[i];
    P->offdRowStarts[i+1] += P->offdRowStarts[i];
  }

  if (P->diagNNZ) {
    P->diagCols  = (iint *)   calloc(P->diagNNZ, sizeof(iint));
    P->diagCoefs = (dfloat *) calloc(P->diagNNZ, sizeof(dfloat));
  }
  if (P->offdNNZ) {
    P->offdCols  = (iint *)   calloc(P->offdNNZ, sizeof(iint));
    P->offdCoefs = (dfloat *) calloc(P->offdNNZ, sizeof(dfloat));
  }

  iint diagCnt = 0;
  iint offdCnt = 0;
  for(iint i=0; i<N; i++) {
    iint col = FineToCoarse[i];
    if ((col>globalAggOffset[rank]-1)||(col<globalAggOffset[rank+1])) {
      P->diagCols[diagCnt] = col - globalAggStart; //local index
      p->diagCoefs[diagCnt++] = level->nullA[i];
    } else {
      P->offdCols[offdCnt] = col;
      P->offdCoefs[offdCnt++] = level->nullA[i];
    }
  }

  if (P->offdNNZ) {
    //we now need to reorder the x vector for the halo, and shift the column indices
    iint *col = (iint *) calloc(P->offdNNZ,sizeof(iint));
    for (iint i=0;i<P->offdNNZ;i++)
      col[i] = offdCols[i]; //copy non-local column global ids

    //sort by global index
    std::sort(col,col+P->offdNNZ);

    //count unique non-local column ids
    P->Nhalo = 0;
    for (iint i=1;i<P->offdNNZ;i++)
      if (col[i]!=col[i-1])  col[++P->Nhalo] = col[i];
    P->Nhalo++; //number of unique columns

    P->Ncols += P->NHalo;

    //save global column ids in colMap
    P->offdColMap    = (iint *)   calloc(P->Nhalo, sizeof(iint));
    for (iint i=0; i<P->Nhalo; i++)
      P->offdColMap[i] = col[i];
    free(col);

    //shift the column indices to local indexing
    for (iint i=0;i<P->offdNNZ;i++) {
      iint gcol = P->offdCols[i];
      for (iint m=0;m<P->Nhalo;m++) {
        if (gcol == P->colMap[m])
          P->offdCols[i] = m;
      }
    }
  }

  csrHaloSetup(P,globalAggOffset);

  // normalize the columns of P
  *nullCoarseA = (dfloat *) calloc(P->Ncols,sizeof(dfloat));

  for(iint i=0; i<diagNNZ; i++)
    (*nullCoarseA)[P->diagCols[i]] += P->diagCoefs[i] * P->diagCoefs[i];
  for(iint i=0; i<offdNNZ; i++)
    (*nullCoarseA)[P->offdCols[i]] += P->offdCoefs[i] * P->offdCoefs[i];

  for(iint i=0; i<NCoarse; i++)
    (*nullCoarseA)[i] = sqrt((*nullCoarseA)[i]);

  crsHaloExchange(P, sizeof(dfloat), *nullCoarseA, P->sendBuffer, *nullCoarseA);

  for(iint i=0; i<diagNNZ; i++)
    P->diagCoefs[i] /= (*nullCoarseA)[P->diagCols[i]];
  for(iint i=0; i<offdNNZ; i++)
    P->offdCoefs[i] /= (*nullCoarseA)[P->offdCols[i]];

  level->R = transpose(level, P, level->globalRowStarts, globalAggOffset);
  level->P = P;
}

typedef struct {

  iint row;
  iint col;
  dfloat val;
  iint owner;

} nonzero_t;

int compareNonZero(const void *a, const void *b){
  nonzero_t *pa = (nonzero_t *) a;
  nonzero_t *pb = (nonzero_t *) b;

  if (pa->owner < pb->owner) return -1;
  if (pa->owner > pb->owner) return +1;

  if (pa->row < pb->row) return -1;
  if (pa->row > pb->row) return +1;

  if (pa->col < pb->col) return -1;
  if (pa->col > pb->col) return +1;

  return 0;
};

csr * transpose(agmgLevel* level, csr *A
                iint *globalRowStarts, iint *globalColStarts){

  // MPI info
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  csr *At = (csr *) calloc(1,sizeof(csr));

  At->Nrows = A->Ncols-A->Nhalo;
  At->Ncols = A->Nrows;
  At->diagNNZ   = A->diagNNZ; //local entries remain local

  At->diagRowStarts = (iint *)   calloc(At->Nrows+1, sizeof(iint));
  At->offdRowStarts = (iint *)   calloc(At->Nrows+1, sizeof(iint));

  //start with local entries
  if (A->diagNNZ) {
    At->diagCols      = (iint *)   calloc(At->diagNNZ, sizeof(iint));
    At->diagCoefs     = (dfloat *) calloc(At->diagNNZ, sizeof(dfloat));
  }

  // count the no of nonzeros per row for transpose
  for(iint i=0; i<A->diagNNZ; i++){
    iint row = A->diagCols[i];
    At->diagRowStarts[row+1]++;
  }

  // cumulative sum for rows
  for(iint i=1; i<=At->Nrows; i++)
    At->diagRowStarts[i] += At->diagRowStarts[i-1];

  iint *counter = (iint *) calloc(At->Nrows+1,sizeof(iint));
  for (iint i=0; i<At->Nrows+1; i++)
    counter[i] = At->rowStarts[i];

  for(iint i=0; i<A->Nrows; i++){
    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    for(iint jj=Jstart; jj<Jend; jj++){
      iint row = A->diagCols[jj];
      At->diagCols[counter[row]]  = i;
      At->diagCoefs[counter[row]] = A->coefs[jj];

      counter[row]++;
    }
  }
  free(counter);

  //now the nonlocal entries. Need to reverse the halo exchange to send the nonzeros
  // count outgoing and incoming meshes
  iint tag = 999;

  nonzero_t *sendNonZeros;
  if (A->offdNNZ)
    sendNonZeros = (nonzero_t *) calloc(A->offdNNZ,sizeof(nonzero_t));

  iint *Nsend = (iint*) calloc(size, sizeof(iint));
  iint *Nrecv = (iint*) calloc(size, sizeof(iint));

  // copy data from nonlocal entries into send buffer
  for(iint i=0;i<A->Nrows+1;++i){
    for (iint j=A->offdRowStarts[i];j<A->offdRowStarts[i+1];j++) {
      iint col =  A->offdColMap[A->offdCols[j]]; //global ids
      for (iint r=0;r<size;r++) { //find owner's rank
        if (globalColStarts[r]-1<id && id < globalColStarts[r+1]) {
          Nsend[r]++;
          sendNonZeros[j].owner = r;
        }
      }
      sendNonZeros[j].row = col;
      sendNonZeros[j].col = i + globalRowStarts[rank];     //global ids
      sendNonZeros[j].val = A->offdCoefs[j];
    }
  }

  //sort outgoing nonzeros by owner, then row and col
  qsort(sendNonZeros, A->offdNNZ, sizeof(nonzero_t), compareNonZero);

  MPI_Alltoall(Nsend, 1, MPI_IINT, Nrecv, 1, MPI_IINT, MPI_COMM_WORLD);

  //count incoming nonzeros
  At->offdNNZ = 0;
  for (iint r=0;r<size;r++)
    At->offdNNZ += Nrecv[r];

  nonzero_t *recvNonZeros;
  if (At->offdNNZ)
    recvNonZeros = (nonzero_t *) calloc(A->offdNNZ,sizeof(nonzero_t));

  // initiate immediate send and receives to each other process as needed
  iint recvOffset = 0;
  iint sendOffset = 0;
  iint sendMessage = 0, recvMessage = 0;
  for(iint r=0;r<size;++r){
    if (At->offdNNZ) {
      if(Nrecv[r]) {
        MPI_Irecv(((char*)recvNonZeros)+recvOffset, Nrecv[r]*sizeof(nonzero_t),
                      MPI_CHAR, r, tag, MPI_COMM_WORLD,
                      (MPI_Request*)A->haloSendRequests+recvMessage);
        recvOffset += Nrecv[r]*sizeof(nonzero_t);
        ++recvMessage;
      }
    }
    if (A->offdNNZ) {
      if(Nsend[r]){
        MPI_Isend(((char*)sendNonZeros)+sendOffset, Nsend[r]*sizeof(nonzero_t),
                      MPI_CHAR, r, tag, MPI_COMM_WORLD,
                      (MPI_Request*)A->haloRecvRequests+sendMessage);
        sendOffset += Nsend[r]*sizeof(nonzero_t);
        ++sendMessage;
      }
    }
  }

  // Wait for all sent messages to have left and received messages to have arrived
  if (A->offdNNZ) {
    MPI_Status *sendStatus = (MPI_Status*) calloc(sendMessage-1, sizeof(MPI_Status));
    MPI_Waitall(sendMessage-1, (MPI_Request*)A->haloRecvRequests, sendStatus);
    free(sendStatus);
  }
  if (At->offdNNZ) {
    MPI_Status *recvStatus = (MPI_Status*) calloc(recvMessage-1, sizeof(MPI_Status));
    MPI_Waitall(recvMessage-1, (MPI_Request*)A->haloSendRequests, recvStatus);
    free(recvStatus);
  }
  if (A->offdNNZ) free(sendNonZeros);
  free(Nsend); free(Nrecv);

  //sort recieved nonzeros by row and col
  qsort(recvNonZeros, At->offdNNZ, sizeof(nonzero_t), compareNonZero);

  if (At->offdNNZ) {
    At->offdCols  = (iint *)   calloc(At->offdNNZ,sizeof(iint));
    At->offdCoefs = (dfloat *) calloc(At->offdNNZ, sizeof(dfloat));

    //find row starts
    iint cnt =0; //row start counter
    for(iint n=0;n<At->offdNNZ;++n)
      if(cnt==0 || (recvNonZeros[n].row!=recvNonZeros[n-1].row))
        At->offdRowStarts[cnt++] = n;
    At->offdRowStarts[cnt] = At->offdNNZ;

    //fill cols and coefs
    for (iint i=0; i<At->Nrows; i++) {
      for (iint j=At->offdRowStarts[i]; j<At->offdRowStarts[i+1]; j++) {
        A->offdCols[j]  = recvNonZeros[j].col;
        A->offdCoefs[j] = recvNonZeros[j].val;
      }
    }
    free(recvNonZeros);

    //we now need to reorder the x vector for the halo, and shift the column indices
    iint *col = (iint *) calloc(At->offdNNZ,sizeof(iint));
    for (iint n=0;n<At->offdNNZ;n++)
      col[n] = At->offdCols[n]; //copy non-local column global ids

    //sort by global index
    std::sort(col,col+At->offdNNZ);

    //count unique non-local column ids
    At->Nhalo = 0;
    for (iint n=1;n<At->offdNNZ;n++)
      if (col[n]!=col[n-1])  col[++At->Nhalo] = col[n];
    At->Nhalo++; //number of unique columns

    At->Ncols += At->NHalo;

    //save global column ids in colMap
    At->offdColMap = (iint *)   calloc(At->Nhalo, sizeof(iint));
    for (iint n=0; n<At->Nhalo; n++)
      At->offdColMap[n] = col[n];
    free(col);

    //shift the column indices to local indexing
    for (iint n=0;n<At->offdNNZ;n++) {
      iint gcol = At->offdCols[n];
      for (iint m=0;m<At->Nhalo;m++) {
        if (gcol == At->offdColMap[m])
          At->offdCols[n] = m;
      }
    }
  }

  csrHaloSetup(At,globalRowStarts);

  return At;
}

csr *galerkinProd(agmgLevel *level){

  iint numAggs = level->numAggss; //local number of aggregates

  csr *A = level->A;
  csr *P = level->P;
  csr *R = level->R;

  csr *RAP = (csr*) calloc(1,sizeof(csr));

  RAP->Nrows = numAggs;
  RAP->Ncols = numAggs;

  RAP->diagRowStarts = (iint *) calloc(numAggs+1, sizeof(iint));

  dfloat *diagDummyCoefs = (dfloat*) calloc(A->diagNNZ,sizeof(dfloat));
  for (iint i=0; i<A->diagNNZ;i++) //copy A coefs
    diagDummyCoefs[i] = A->diagCoefs[i];

  long base = numAggs+1;
  struct key_value_pair1 *pair = (struct key_value_pair1 *)
                  calloc(A->nnz,sizeof(struct key_value_pair1));

  for(iint i=0; i<A->Nrows; i++){
    for(iint jj=A->rowStarts[i]; jj<A->rowStarts[i+1]; jj++){
      iint j = A->cols[jj];

      diagDummyCoefs[jj] *= (P->coefs[i] * P->coefs[j]);

      iint I = P->cols[i];
      iint J = P->cols[j];

      pair[jj].value = jj;

      if (I == J)
        pair[jj].key = ((long)I) * base;
      else
        pair[jj].key = ((long)I) * base + ((long) (J+1));

    }
  }

  qsort(pair, A->nnz, sizeof(struct key_value_pair1), compare_key1);

  iint nnz = 1;
  for(iint i=1; i<A->nnz; i++)
    if(pair[i].key > pair[i-1].key) nnz++;

  RAP->nnz = nnz;

  RAP->cols  = (iint *) calloc(nnz, sizeof(iint));
  RAP->coefs = (dfloat *) calloc(nnz, sizeof(dfloat));

  iint count = 0;
  RAP->coefs[count] = dummyCoefs[pair[0].value];
  RAP->cols[count] = 0; // This assumes that first entry is diagonal
  RAP->rowStarts[1]++;
  for(iint i=1; i<A->nnz; i++){
    iint id = pair[i].value;

    if(pair[i].key == pair[i-1].key)
      RAP->coefs[count] += dummyCoefs[id];
    else {
      RAP->coefs[count+1] = dummyCoefs[id];

      iint J = pair[i].key % base;
      iint I = pair[i].key / base;
      J = (J==0) ? I : J-1;

      RAP->cols[count+1] = J;
      RAP->rowStarts[I+1]++;

      count++;
    }
  }
  free(dummyCoefs);

  // cumulative sum
  for(iint i=1; i<=numAggs; i++)
    RAP->rowStarts[i] += RAP->rowStarts[i-1];

  //MPI comms
  RAP->NrecvTotal =0;
  RAP->NsendTotal =0;

  return RAP;
}


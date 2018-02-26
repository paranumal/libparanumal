#include "agmg.h"

csr *strong_graph(csr *A, dfloat threshold);
bool customLess(int smax, dfloat rmax, int imax, int s, dfloat r, int i);
int *form_aggregates(agmgLevel *level, csr *C);
void find_aggregate_owners(agmgLevel *level, int* FineToCoarse, const char *options);
csr *construct_interpolator(agmgLevel *level, int *FineToCoarse, dfloat **nullCoarseA);
csr *transpose(agmgLevel* level, csr *A, int *globalRowStarts, int *globalColStarts);
csr *galerkinProd(agmgLevel *level, csr *R, csr *A, csr *P);
void coarsenAgmgLevel(agmgLevel *level, csr **coarseA, csr **P, csr **R, dfloat **nullCoarseA, const char *options);


void agmgSetup(parAlmond_t *parAlmond, csr *A, dfloat *nullA, int *globalRowStarts, const char* options){

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // approximate Nrows at coarsest level
  int gCoarseSize = 1000;

  double seed = (double) rank;
  srand48(seed);

  agmgLevel **levels = parAlmond->levels;

  int lev = parAlmond->numLevels; //add this level to the end of the chain

  levels[lev] = (agmgLevel *) calloc(1,sizeof(agmgLevel));
  levels[lev]->gatherLevel = false;
  levels[lev]->weightedInnerProds = false;
  parAlmond->numLevels++;

  //copy A matrix and null vector
  levels[lev]->A = A;
  levels[lev]->A->null = nullA;

  levels[lev]->Nrows = A->Nrows;
  levels[lev]->Ncols = A->Ncols;

  
  SmoothType smoothType;
  if (strstr(options,"CHEBYSHEV")) {
    smoothType = CHEBYSHEV;
  } else { //default to DAMPED_JACOBI
    smoothType = DAMPED_JACOBI;
  }
  setupSmoother(parAlmond, levels[lev], smoothType);

  levels[lev]->deviceA = newHYB(parAlmond, levels[lev]->A);

  //set operator callback
  void **args = (void **) calloc(2,sizeof(void*));
  args[0] = (void *) parAlmond;
  args[1] = (void *) levels[lev];

  levels[lev]->AxArgs = args;
  levels[lev]->smoothArgs = args;
  levels[lev]->Ax = agmgAx;
  levels[lev]->smooth = agmgSmooth;
  levels[lev]->device_Ax = device_agmgAx;
  levels[lev]->device_smooth = device_agmgSmooth;

  //copy global partiton
  levels[lev]->globalRowStarts = (int *) calloc(size+1,sizeof(int));
  for (int r=0;r<size+1;r++)
      levels[lev]->globalRowStarts[r] = globalRowStarts[r];

  int localSize = levels[lev]->A->Nrows;
  int globalSize = 0;
  MPI_Allreduce(&localSize, &globalSize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  //if the system if already small, dont create MG levels
  bool done = false;
  if(globalSize <= gCoarseSize){
    setupExactSolve(parAlmond, levels[lev],parAlmond->nullSpace,parAlmond->nullSpacePenalty);
    //setupSmoother(parAlmond, levels[lev], smoothType);
    done = true;
  }
  while(!done){
    // create coarse MG level
    levels[lev+1] = (agmgLevel *) calloc(1,sizeof(agmgLevel));
    dfloat *nullCoarseA;

    //printf("Setting up coarse level %d\n", lev+1);

    coarsenAgmgLevel(levels[lev], &(levels[lev+1]->A), &(levels[lev+1]->P),
                                  &(levels[lev+1]->R), &nullCoarseA, parAlmond->options);

    //set dimensions of the fine level (max among the A,R ops)
    levels[lev]->Ncols = mymax(levels[lev]->Ncols, levels[lev+1]->R->Ncols);

    parAlmond->numLevels++;

    levels[lev+1]->A->null = nullCoarseA;
    levels[lev+1]->Nrows = levels[lev+1]->A->Nrows;
    levels[lev+1]->Ncols = mymax(levels[lev+1]->A->Ncols, levels[lev+1]->P->Ncols);
    levels[lev+1]->globalRowStarts = levels[lev]->globalAggStarts;
    
    setupSmoother(parAlmond, levels[lev+1], smoothType);

    levels[lev+1]->deviceA = newHYB (parAlmond, levels[lev+1]->A);
    levels[lev+1]->deviceR = newHYB (parAlmond, levels[lev+1]->R);
    levels[lev+1]->dcsrP   = newDCOO(parAlmond, levels[lev+1]->P);

    //set operator callback
    void **args = (void **) calloc(2,sizeof(void*));
    args[0] = (void *) parAlmond;
    args[1] = (void *) levels[lev+1];

    levels[lev+1]->AxArgs = args;
    levels[lev+1]->coarsenArgs = args;
    levels[lev+1]->prolongateArgs = args;
    levels[lev+1]->smoothArgs = args;

    levels[lev+1]->Ax = agmgAx;
    levels[lev+1]->coarsen = agmgCoarsen;
    levels[lev+1]->prolongate = agmgProlongate;
    levels[lev+1]->smooth = agmgSmooth;

    levels[lev+1]->device_Ax = device_agmgAx;
    levels[lev+1]->device_coarsen = device_agmgCoarsen;
    levels[lev+1]->device_prolongate = device_agmgProlongate;
    levels[lev+1]->device_smooth = device_agmgSmooth;

    const int localCoarseDim = levels[lev+1]->A->Nrows;
    int globalCoarseSize;
    MPI_Allreduce(&localCoarseDim, &globalCoarseSize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(globalCoarseSize <= gCoarseSize || globalSize < 2*globalCoarseSize){
      setupExactSolve(parAlmond, levels[lev+1],parAlmond->nullSpace,parAlmond->nullSpacePenalty);
      //setupSmoother(parAlmond, levels[lev+1], smoothType);
      break;
    }

    globalSize = globalCoarseSize;
    lev++;
  } 
  
  //allocate vectors required
  occa::device device = parAlmond->device;
  for (int n=0;n<parAlmond->numLevels;n++) {
    int N = levels[n]->Nrows;
    int M = levels[n]->Ncols;

    if ((n>0)&&(n<parAlmond->numLevels)) { //kcycle vectors
      if (M) levels[n]->ckp1 = (dfloat *) calloc(M,sizeof(dfloat));
      if (N) levels[n]->vkp1 = (dfloat *) calloc(N,sizeof(dfloat));
      if (N) levels[n]->wkp1 = (dfloat *) calloc(N,sizeof(dfloat));

      if (M) levels[n]->o_ckp1 = device.malloc(M*sizeof(dfloat),levels[n]->ckp1);
      if (N) levels[n]->o_vkp1 = device.malloc(N*sizeof(dfloat),levels[n]->vkp1);
      if (N) levels[n]->o_wkp1 = device.malloc(N*sizeof(dfloat),levels[n]->wkp1);
    }
    if (M) levels[n]->x    = (dfloat *) calloc(M,sizeof(dfloat));
    if (M) levels[n]->res  = (dfloat *) calloc(M,sizeof(dfloat));
    if (N) levels[n]->rhs  = (dfloat *) calloc(N,sizeof(dfloat));

    if (M) levels[n]->o_x   = device.malloc(M*sizeof(dfloat),levels[n]->x);
    if (M) levels[n]->o_res = device.malloc(M*sizeof(dfloat),levels[n]->res);
    if (N) levels[n]->o_rhs = device.malloc(N*sizeof(dfloat),levels[n]->rhs);
  }
}

void parAlmondReport(parAlmond_t *parAlmond) {

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0) {
    printf("------------------ParAlmond Report-----------------------------------\n");
    printf("---------------------------------------------------------------------\n");
    printf("level| active ranks |   dimension   |  nnzs         |  nnz/row      |\n");
    printf("     |              | (min,max,avg) | (min,max,avg) | (min,max,avg) |\n");
    printf("---------------------------------------------------------------------\n");
  }

  for(int lev=0; lev<parAlmond->numLevels; lev++){

    int Nrows = parAlmond->levels[lev]->Nrows;

    int active = (Nrows>0) ? 1:0;
    int totalActive=0;
    MPI_Allreduce(&active, &totalActive, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int minNrows=0, maxNrows=0, totalNrows=0;
    dfloat avgNrows;
    MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&Nrows, &totalNrows, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    avgNrows = (dfloat) totalNrows/totalActive;

    if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
    MPI_Allreduce(&Nrows, &minNrows, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);


    int nnz;
    if (parAlmond->levels[lev]->A)
      nnz = parAlmond->levels[lev]->A->diagNNZ+parAlmond->levels[lev]->A->offdNNZ;
    else
      nnz =0;
    int minNnz=0, maxNnz=0, totalNnz=0;
    dfloat avgNnz;
    MPI_Allreduce(&nnz, &maxNnz, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&nnz, &totalNnz, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    avgNnz = (dfloat) totalNnz/totalActive;

    if (nnz==0) nnz = maxNnz; //set this so it's ignored for the global min
    MPI_Allreduce(&nnz, &minNnz, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    Nrows = parAlmond->levels[lev]->Nrows;
    dfloat nnzPerRow = (Nrows==0) ? 0 : (dfloat) nnz/Nrows;
    dfloat minNnzPerRow=0, maxNnzPerRow=0, avgNnzPerRow=0;
    MPI_Allreduce(&nnzPerRow, &maxNnzPerRow, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&nnzPerRow, &avgNnzPerRow, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
    avgNnzPerRow /= totalActive;

    if (Nrows==0) nnzPerRow = maxNnzPerRow;
    MPI_Allreduce(&nnzPerRow, &minNnzPerRow, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

    if (rank==0){
      printf(" %3d |        %4d  |   %10.2f  |   %10.2f  |   %10.2f  |\n",
        lev, totalActive, (dfloat)minNrows, (dfloat)minNnz, minNnzPerRow);
      printf("     |              |   %10.2f  |   %10.2f  |   %10.2f  |\n",
        (dfloat)maxNrows, (dfloat)maxNnz, maxNnzPerRow);
      printf("     |              |   %10.2f  |   %10.2f  |   %10.2f  |\n",
        avgNrows, avgNnz, avgNnzPerRow);
    }
  }
  if(rank==0)
    printf("---------------------------------------------------------------------\n");
}


//create coarsened problem
void coarsenAgmgLevel(agmgLevel *level, csr **coarseA, csr **P, csr **R, dfloat **nullCoarseA, const char* options){

  // establish the graph of strong connections
  level->threshold = 0.5;

  csr *C = strong_graph(level->A, level->threshold);

  int *FineToCoarse = form_aggregates(level, C);

  find_aggregate_owners(level,FineToCoarse,options);

  *P = construct_interpolator(level, FineToCoarse, nullCoarseA);
  *R = transpose(level, *P, level->globalRowStarts, level->globalAggStarts);
  *coarseA = galerkinProd(level, *R, level->A, *P);
}

csr * strong_graph(csr *A, dfloat threshold){

  const int N = A->Nrows;
  const int M = A->Ncols;

  csr *C = (csr *) calloc(1, sizeof(csr));

  C->Nrows = N;
  C->Ncols = M;

  C->diagRowStarts = (int *) calloc(N+1,sizeof(int));
  C->offdRowStarts = (int *) calloc(N+1,sizeof(int));

  dfloat *maxOD;
  if (N) maxOD = (dfloat *) calloc(N,sizeof(dfloat));

  //store the diagonal of A for all needed columns
  dfloat *diagA = (dfloat *) calloc(M,sizeof(dfloat));
  for (int i=0;i<N;i++)
    diagA[i] = A->diagCoefs[A->diagRowStarts[i]];
  csrHaloExchange(A, sizeof(dfloat), diagA, A->sendBuffer, diagA+A->NlocalCols);

  #pragma omp parallel for
  for(int i=0; i<N; i++){
    dfloat sign = (diagA[i] >= 0) ? 1:-1;
    dfloat Aii = fabs(diagA[i]);

    //find maxOD
    //local entries
    int Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];
    for(int jj= Jstart+1; jj<Jend; jj++){
      int col = A->diagCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->diagCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > maxOD[i]) maxOD[i] = OD;
    }
    //non-local entries
    Jstart = A->offdRowStarts[i], Jend = A->offdRowStarts[i+1];
    for(int jj= Jstart; jj<Jend; jj++){
      int col = A->offdCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->offdCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > maxOD[i]) maxOD[i] = OD;
    }

    int diag_strong_per_row = 1; // diagonal entry
    //local entries
    Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];
    for(int jj = Jstart+1; jj<Jend; jj++){
      int col = A->diagCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->diagCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD[i]) diag_strong_per_row++;
    }
    int offd_strong_per_row = 0;
    //non-local entries
    Jstart = A->offdRowStarts[i], Jend = A->offdRowStarts[i+1];
    for(int jj= Jstart; jj<Jend; jj++){
      int col = A->offdCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->offdCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD[i]) offd_strong_per_row++;
    }

    C->diagRowStarts[i+1] = diag_strong_per_row;
    C->offdRowStarts[i+1] = offd_strong_per_row;
  }

  // cumulative sum
  for(int i=1; i<N+1 ; i++) {
    C->diagRowStarts[i] += C->diagRowStarts[i-1];
    C->offdRowStarts[i] += C->offdRowStarts[i-1];
  }

  C->diagNNZ = C->diagRowStarts[N];
  C->offdNNZ = C->offdRowStarts[N];

  if (C->diagNNZ) C->diagCols = (int *) calloc(C->diagNNZ, sizeof(int));
  if (C->offdNNZ) C->offdCols = (int *) calloc(C->offdNNZ, sizeof(int));

  // fill in the columns for strong connections
  #pragma omp parallel for
  for(int i=0; i<N; i++){
    dfloat sign = (diagA[i] >= 0) ? 1:-1;
    dfloat Aii = fabs(diagA[i]);

    int diagCounter = C->diagRowStarts[i];
    int offdCounter = C->offdRowStarts[i];

    //local entries
    C->diagCols[diagCounter++] = i;// diag entry
    int Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];
    for(int jj = Jstart+1; jj<Jend; jj++){
      int col = A->diagCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->diagCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD[i])
        C->diagCols[diagCounter++] = A->diagCols[jj];
    }
    Jstart = A->offdRowStarts[i], Jend = A->offdRowStarts[i+1];
    for(int jj = Jstart; jj<Jend; jj++){
      int col = A->offdCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->offdCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD[i])
        C->offdCols[offdCounter++] = A->offdCols[jj];
    }
  }
  if(N) free(maxOD);

  return C;
}

bool customLess(int smax, dfloat rmax, int imax, int s, dfloat r, int i){

  if(s > smax) return true;
  if(smax > s) return false;

  if(r > rmax) return true;
  if(rmax > r) return false;

  if(i > imax) return true;
  if(i < imax) return false;

  return false;
}

int * form_aggregates(agmgLevel *level, csr *C){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int N   = C->Nrows;
  const int M   = C->Ncols;
  const int diagNNZ = C->diagNNZ;
  const int offdNNZ = C->offdNNZ;

  int *FineToCoarse = (int *) calloc(M, sizeof(int));
  for (int i =0;i<M;i++) FineToCoarse[i] = -1;

  dfloat *rands  = (dfloat *) calloc(M, sizeof(dfloat));
  int   *states = (int *)   calloc(M, sizeof(int));

  dfloat *Tr = (dfloat *) calloc(M, sizeof(dfloat));
  int   *Ts = (int *)   calloc(M, sizeof(int));
  int   *Ti = (int *)   calloc(M, sizeof(int));
  int   *Tc = (int *)   calloc(M, sizeof(int));

  csr *A = level->A;
  int *globalRowStarts = level->globalRowStarts;

  int *intSendBuffer;
  dfloat *dfloatSendBuffer;
  if (level->A->NsendTotal) {
    intSendBuffer = (int *) calloc(A->NsendTotal,sizeof(int));
    dfloatSendBuffer = (dfloat *) calloc(A->NsendTotal,sizeof(dfloat));
  }

  for(int i=0; i<N; i++)
    rands[i] = (dfloat) drand48();

  for(int i=0; i<N; i++)
    states[i] = 0;

  // add the number of non-zeros in each column
  //local non-zeros
  for(int i=0; i<diagNNZ; i++)
    rands[C->diagCols[i]] += 1.;

  int *nnzCnt, *recvNnzCnt;
  if (A->NHalo) nnzCnt = (int *) calloc(A->NHalo,sizeof(int));
  if (A->NsendTotal) recvNnzCnt = (int *) calloc(A->NsendTotal,sizeof(int));

  //count the non-local non-zeros
  for (int i=0;i<offdNNZ;i++)
    nnzCnt[C->offdCols[i]-A->NlocalCols]++;

  //do a reverse halo exchange
  int tag = 999;

  // initiate immediate send  and receives to each other process as needed
  int recvOffset = 0;
  int sendOffset = 0;
  int sendMessage = 0, recvMessage = 0;
  for(int r=0;r<size;++r){
    if (A->NsendTotal) {
      if(A->NsendPairs[r]) {
        MPI_Irecv(recvNnzCnt+sendOffset, A->NsendPairs[r], MPI_INT, r, tag,
            MPI_COMM_WORLD, (MPI_Request*)A->haloSendRequests+sendMessage);
        sendOffset += A->NsendPairs[r];
        ++sendMessage;
      }
    }
    if (A->NrecvTotal) {
      if(A->NrecvPairs[r]){
        MPI_Isend(nnzCnt+recvOffset, A->NrecvPairs[r], MPI_INT, r, tag,
            MPI_COMM_WORLD, (MPI_Request*)A->haloRecvRequests+recvMessage);
        recvOffset += A->NrecvPairs[r];
        ++recvMessage;
      }
    }
  }

  // Wait for all sent messages to have left and received messages to have arrived
  if (A->NrecvTotal) {
    MPI_Status *sendStatus = (MPI_Status*) calloc(A->NsendMessages, sizeof(MPI_Status));
    MPI_Waitall(A->NsendMessages, (MPI_Request*)A->haloSendRequests, sendStatus);
    free(sendStatus);
  }
  if (A->NsendTotal) {
    MPI_Status *recvStatus = (MPI_Status*) calloc(A->NrecvMessages, sizeof(MPI_Status));
    MPI_Waitall(A->NrecvMessages, (MPI_Request*)A->haloRecvRequests, recvStatus);
    free(recvStatus);
  }

  for(int i=0;i<A->NsendTotal;++i){
    // local index of outgoing element in halo exchange
    int id = A->haloElementList[i];

    rands[id] += recvNnzCnt[i];
  }

  if (A->NHalo) free(nnzCnt);
  if (A->NsendTotal) free(recvNnzCnt);

  //share randomizer values
  csrHaloExchange(A, sizeof(dfloat), rands, dfloatSendBuffer, rands+A->NlocalCols);



  int done = 0;
  while(!done){
    // first neighbours
    #pragma omp parallel for
    for(int i=0; i<N; i++){

      int smax = states[i];
      dfloat rmax = rands[i];
      int imax = i + globalRowStarts[rank];

      if(smax != 1){
        //local entries
        for(int jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
          const int col = C->diagCols[jj];
          if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
            smax = states[col];
            rmax = rands[col];
            imax = col + globalRowStarts[rank];
          }
        }
        //nonlocal entries
        for(int jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
          const int col = C->offdCols[jj];
          if(customLess(smax, rmax, imax, states[col], rands[col], A->colMap[col])) {
            smax = states[col];
            rmax = rands[col];
            imax = A->colMap[col];
          }
        }
      }
      Ts[i] = smax;
      Tr[i] = rmax;
      Ti[i] = imax;
    }

    //share results
    csrHaloExchange(A, sizeof(dfloat), Tr, dfloatSendBuffer, Tr+A->NlocalCols);
    csrHaloExchange(A, sizeof(int), Ts, intSendBuffer, Ts+A->NlocalCols);
    csrHaloExchange(A, sizeof(int), Ti, intSendBuffer, Ti+A->NlocalCols);

    // second neighbours
    #pragma omp parallel for
    for(int i=0; i<N; i++){
      int   smax = Ts[i];
      dfloat rmax = Tr[i];
      int   imax = Ti[i];

      //local entries
      for(int jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
        const int col = C->diagCols[jj];
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }
      //nonlocal entries
      for(int jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
        const int col = C->offdCols[jj];
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

    csrHaloExchange(A, sizeof(int), states, intSendBuffer, states+A->NlocalCols);

    // if number of undecided nodes = 0, algorithm terminates
    int cnt = std::count(states, states+N, 0);
    MPI_Allreduce(&cnt,&done,1,MPI_INT, MPI_SUM,MPI_COMM_WORLD);
    done = (done == 0) ? 1 : 0;
  }

  int numAggs = 0;
  level->globalAggStarts = (int *) calloc(size+1,sizeof(int));
  // count the coarse nodes/aggregates
  for(int i=0; i<N; i++)
    if(states[i] == 1) numAggs++;

  MPI_Allgather(&numAggs,1,MPI_INT,level->globalAggStarts+1,1,MPI_INT,MPI_COMM_WORLD);

  for (int r=0;r<size;r++)
    level->globalAggStarts[r+1] += level->globalAggStarts[r];

  numAggs = 0;
  // enumerate the coarse nodes/aggregates
  for(int i=0; i<N; i++)
    if(states[i] == 1)
      FineToCoarse[i] = level->globalAggStarts[rank] + numAggs++;

  //share the initial aggregate flags
  csrHaloExchange(A, sizeof(int), FineToCoarse, intSendBuffer, FineToCoarse+A->NlocalCols);

  // form the aggregates
  #pragma omp parallel for
  for(int i=0; i<N; i++){
    int   smax = states[i];
    dfloat rmax = rands[i];
    int   imax = i + globalRowStarts[rank];
    int   cmax = FineToCoarse[i];

    if(smax != 1){
      //local entries
      for(int jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
        const int col = C->diagCols[jj];
        if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
          smax = states[col];
          rmax = rands[col];
          imax = col + globalRowStarts[rank];
          cmax = FineToCoarse[col];
        }
      }
      //nonlocal entries
      for(int jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
        const int col = C->offdCols[jj];
        if(customLess(smax, rmax, imax, states[col], rands[col], A->colMap[col])){
          smax = states[col];
          rmax = rands[col];
          imax = A->colMap[col];
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

  csrHaloExchange(A, sizeof(int), FineToCoarse, intSendBuffer, FineToCoarse+A->NlocalCols);
  csrHaloExchange(A, sizeof(dfloat), Tr, dfloatSendBuffer, Tr+A->NlocalCols);
  csrHaloExchange(A, sizeof(int), Ts, intSendBuffer, Ts+A->NlocalCols);
  csrHaloExchange(A, sizeof(int), Ti, intSendBuffer, Ti+A->NlocalCols);
  csrHaloExchange(A, sizeof(int), Tc, intSendBuffer, Tc+A->NlocalCols);

  // second neighbours
  #pragma omp parallel for
  for(int i=0; i<N; i++){
    int   smax = Ts[i];
    dfloat rmax = Tr[i];
    int   imax = Ti[i];
    int   cmax = Tc[i];

    //local entries
    for(int jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
      const int col = C->diagCols[jj];
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }
    //nonlocal entries
    for(int jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
      const int col = C->offdCols[jj];
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

  csrHaloExchange(A, sizeof(int), FineToCoarse, intSendBuffer, FineToCoarse+A->NlocalCols);

  free(rands);
  free(states);
  free(Tr);
  free(Ts);
  free(Ti);
  free(Tc);
  if (level->A->NsendTotal) {
    free(intSendBuffer);
    free(dfloatSendBuffer);
  }

  //TODO maybe free C here?

  return FineToCoarse;
}

typedef struct {

  int fineId;
  int coarseId;
  int newCoarseId;

  int originRank;
  int ownerRank;

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

  if (pa->originRank < pb->originRank) return -1;
  if (pa->originRank > pb->originRank) return +1;

  return 0;
};

int compareOrigin(const void *a, const void *b){
  parallelAggregate_t *pa = (parallelAggregate_t *) a;
  parallelAggregate_t *pb = (parallelAggregate_t *) b;

  if (pa->originRank < pb->originRank) return -1;
  if (pa->originRank > pb->originRank) return +1;

  return 0;
};

void find_aggregate_owners(agmgLevel *level, int* FineToCoarse, const char* options) {
  // MPI info
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int N = level->A->Nrows;

  //Need to establish 'ownership' of aggregates
  
  //Keep the current partitioning for STRONGNODES. 
  // The rank that had the strong node for each aggregate owns the aggregate
  if (strstr(options,"STRONGNODES")) return;

  /* Setup description of the MPI_PARALLEL_AGGREGATE struct */
  MPI_Datatype MPI_PARALLEL_AGGREGATE;
  MPI_Datatype oldtypes[1] = {MPI_INT};
  int blockcounts[1] = {5};
  MPI_Aint  entryoffsets[1] = {0};

  /* Now define structured type and commit it */
  MPI_Type_struct(1, blockcounts, entryoffsets, oldtypes, &MPI_PARALLEL_AGGREGATE);
  MPI_Type_commit(&MPI_PARALLEL_AGGREGATE);

  //populate aggregate array
  int gNumAggs = level->globalAggStarts[size]; //total number of aggregates
  
  parallelAggregate_t *sendAggs;
  if (N) sendAggs = (parallelAggregate_t *) calloc(N,sizeof(parallelAggregate_t));
  for (int i=0;i<N;i++) {
    sendAggs[i].fineId = i;
    sendAggs[i].originRank = rank;

    sendAggs[i].coarseId = FineToCoarse[i];

    //set a temporary owner. Evenly distibute aggregates amoungst ranks
    sendAggs[i].ownerRank = (FineToCoarse[i]*size)/gNumAggs;
  }

  //sort by owning rank for all_reduce
  qsort(sendAggs, N, sizeof(parallelAggregate_t), compareOwner);

  int *sendCounts = (int *) calloc(size,sizeof(int));
  int *recvCounts = (int *) calloc(size,sizeof(int));
  int *sendOffsets = (int *) calloc(size+1,sizeof(int));
  int *recvOffsets = (int *) calloc(size+1,sizeof(int));

  for(int i=0;i<N;++i)
    sendCounts[sendAggs[i].ownerRank]++;

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  int recvNtotal = 0;
  for(int r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r];
  }
  parallelAggregate_t *recvAggs = (parallelAggregate_t *) calloc(recvNtotal,sizeof(parallelAggregate_t));

  MPI_Alltoallv(sendAggs, sendCounts, sendOffsets, MPI_PARALLEL_AGGREGATE,
                recvAggs, recvCounts, recvOffsets, MPI_PARALLEL_AGGREGATE,
                MPI_COMM_WORLD);

  //sort by coarse aggregate number, and then by original rank
  qsort(recvAggs, recvNtotal, sizeof(parallelAggregate_t), compareAgg);

  //count the number of unique aggregates here
  int NumUniqueAggs =0;
  if (recvNtotal) NumUniqueAggs++;
  for (int i=1;i<recvNtotal;i++)
    if(recvAggs[i].coarseId!=recvAggs[i-1].coarseId) NumUniqueAggs++;

  //get their locations in the array
  int *aggStarts;
  if (NumUniqueAggs)
    aggStarts = (int *) calloc(NumUniqueAggs+1,sizeof(int));
  int cnt = 1;
  for (int i=1;i<recvNtotal;i++)
    if(recvAggs[i].coarseId!=recvAggs[i-1].coarseId) aggStarts[cnt++] = i;
  aggStarts[NumUniqueAggs] = recvNtotal;


  if (strstr(options,"DISTRIBUTED")) { //rank that contributes most to the aggregate ownes it
    //use a random dfloat for each rank to break ties.
    dfloat rand = (dfloat) drand48();
    dfloat *gRands = (dfloat *) calloc(size,sizeof(dfloat));
    MPI_Allgather(&rand, 1, MPI_DFLOAT, gRands, 1, MPI_DFLOAT, MPI_COMM_WORLD);

    //determine the aggregates majority owner
    dfloat *rankCounts = (dfloat *) calloc(size,sizeof(dfloat));
    for (int n=0;n<NumUniqueAggs;n++) {
      //populate randomizer
      for (int r=0;r<size;r++)
        rankCounts[r] = gRands[r];

      //count the number of contributions to the aggregate from the separate ranks
      for (int i=aggStarts[n];i<aggStarts[n+1];i++)
        rankCounts[recvAggs[i].originRank]++;

      //find which rank is contributing the most to this aggregate
      int ownerRank = 0;
      dfloat maxEntries = rankCounts[0];
      for (int r=1;r<size;r++) {
        if (rankCounts[r]>maxEntries) {
          ownerRank = r;
          maxEntries = rankCounts[r];
        }
      }

      //set this aggregate's owner
      for (int i=aggStarts[n];i<aggStarts[n+1];i++)
        recvAggs[i].ownerRank = ownerRank;
    }
    free(gRands); free(rankCounts);
  } else { //default SATURATE: always choose the lowest rank to own the aggregate
    for (int n=0;n<NumUniqueAggs;n++) {
      
      int minrank = size;

      //count the number of contributions to the aggregate from the separate ranks
      for (int i=aggStarts[n];i<aggStarts[n+1];i++){

        minrank = (recvAggs[i].originRank<minrank) ? recvAggs[i].originRank : minrank;
      }

      //set this aggregate's owner
      for (int i=aggStarts[n];i<aggStarts[n+1];i++)
        recvAggs[i].ownerRank = minrank;
    }
  }
  free(aggStarts);

  //sort by owning rank
  qsort(recvAggs, recvNtotal, sizeof(parallelAggregate_t), compareOwner);

  int *newSendCounts = (int *) calloc(size,sizeof(int));
  int *newRecvCounts = (int *) calloc(size,sizeof(int));
  int *newSendOffsets = (int *) calloc(size+1,sizeof(int));
  int *newRecvOffsets = (int *) calloc(size+1,sizeof(int));

  for(int i=0;i<recvNtotal;++i)
    newSendCounts[recvAggs[i].ownerRank]++;

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(newSendCounts, 1, MPI_INT, newRecvCounts, 1, MPI_INT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  int newRecvNtotal = 0;
  for(int r=0;r<size;++r){
    newSendOffsets[r+1] = newSendOffsets[r] + newSendCounts[r];
    newRecvOffsets[r+1] = newRecvOffsets[r] + newRecvCounts[r];
    newRecvNtotal += newRecvCounts[r];
  }
  parallelAggregate_t *newRecvAggs = (parallelAggregate_t *) calloc(newRecvNtotal,sizeof(parallelAggregate_t));

  MPI_Alltoallv(   recvAggs, newSendCounts, newSendOffsets, MPI_PARALLEL_AGGREGATE,
                newRecvAggs, newRecvCounts, newRecvOffsets, MPI_PARALLEL_AGGREGATE,
                MPI_COMM_WORLD);

  //sort by coarse aggregate number, and then by original rank
  qsort(newRecvAggs, newRecvNtotal, sizeof(parallelAggregate_t), compareAgg);

  //count the number of unique aggregates this rank owns
  int numAggs = 0;
  if (newRecvNtotal) numAggs++;
  for (int i=1;i<newRecvNtotal;i++)
    if(newRecvAggs[i].coarseId!=newRecvAggs[i-1].coarseId) numAggs++;

  //determine a global numbering of the aggregates
  MPI_Allgather(&numAggs, 1, MPI_INT, level->globalAggStarts+1, 1, MPI_INT, MPI_COMM_WORLD);

  for (int r=0;r<size;r++)
    level->globalAggStarts[r+1] += level->globalAggStarts[r];

  //set the new global coarse index
  cnt = level->globalAggStarts[rank];
  if (newRecvNtotal) newRecvAggs[0].newCoarseId = cnt;
  for (int i=1;i<newRecvNtotal;i++) {
    if(newRecvAggs[i].coarseId!=newRecvAggs[i-1].coarseId) cnt++;

    newRecvAggs[i].newCoarseId = cnt;
  }

  //sort by owning rank
  qsort(newRecvAggs, newRecvNtotal, sizeof(parallelAggregate_t), compareOrigin);

  for(int r=0;r<size;r++) sendCounts[r] = 0;  
  for(int r=0;r<=size;r++) {
    sendOffsets[r] = 0;
    recvOffsets[r] = 0;
  }

  for(int i=0;i<newRecvNtotal;++i)
    sendCounts[newRecvAggs[i].originRank]++;

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  recvNtotal = 0;
  for(int r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r];
  }

  //send the aggregate data back
  MPI_Alltoallv(newRecvAggs, sendCounts, sendOffsets, MPI_PARALLEL_AGGREGATE,
                   sendAggs, recvCounts, recvOffsets, MPI_PARALLEL_AGGREGATE,
                MPI_COMM_WORLD);

  //clean up
  MPI_Type_free(&MPI_PARALLEL_AGGREGATE);

  free(recvAggs);
  free(sendCounts);  free(recvCounts);
  free(sendOffsets); free(recvOffsets);
  free(newRecvAggs);
  free(newSendCounts);  free(newRecvCounts);
  free(newSendOffsets); free(newRecvOffsets);

  //record the new FineToCoarse map
  for (int i=0;i<N;i++)
    FineToCoarse[sendAggs[i].fineId] = sendAggs[i].newCoarseId;

  if (N) free(sendAggs);
}


csr *construct_interpolator(agmgLevel *level, int *FineToCoarse, dfloat **nullCoarseA){
  // MPI info
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int N = level->A->Nrows;
  const int M = level->A->Ncols;

  int *globalAggStarts = level->globalAggStarts;

  const int globalAggOffset = level->globalAggStarts[rank];
  const int NCoarse = globalAggStarts[rank+1]-globalAggStarts[rank]; //local num agg

  csr* P = (csr *) calloc(1, sizeof(csr));

  P->Nrows = N;
  P->Ncols = NCoarse;

  P->NlocalCols = NCoarse;
  P->NHalo = 0;

  P->diagRowStarts = (int *) calloc(N+1, sizeof(int));
  P->offdRowStarts = (int *) calloc(N+1, sizeof(int));

  // each row has exactly one nonzero per row
  P->diagNNZ =0;
  P->offdNNZ =0;
  for(int i=0; i<N; i++) {
    int col = FineToCoarse[i];
    if ((col>globalAggOffset-1)&&(col<globalAggOffset+NCoarse)) {
      P->diagNNZ++;
      P->diagRowStarts[i+1]++;
    } else {
      P->offdNNZ++;
      P->offdRowStarts[i+1]++;
    }
  }
  for(int i=0; i<N; i++) {
    P->diagRowStarts[i+1] += P->diagRowStarts[i];
    P->offdRowStarts[i+1] += P->offdRowStarts[i];
  }

  if (P->diagNNZ) {
    P->diagCols  = (int *)   calloc(P->diagNNZ, sizeof(int));
    P->diagCoefs = (dfloat *) calloc(P->diagNNZ, sizeof(dfloat));
  }
  if (P->offdNNZ) {
    P->offdCols  = (int *)   calloc(P->offdNNZ, sizeof(int));
    P->offdCoefs = (dfloat *) calloc(P->offdNNZ, sizeof(dfloat));
  }

  int diagCnt = 0;
  int offdCnt = 0;
  for(int i=0; i<N; i++) {
    int col = FineToCoarse[i];
    if ((col>globalAggStarts[rank]-1)&&(col<globalAggStarts[rank+1])) {
      P->diagCols[diagCnt] = col - globalAggOffset; //local index
      P->diagCoefs[diagCnt++] = level->A->null[i];
    } else {
      P->offdCols[offdCnt] = col;
      P->offdCoefs[offdCnt++] = level->A->null[i];
    }
  }

  //record global indexing of columns
  P->colMap = (int *)   calloc(P->Ncols, sizeof(int));
  for (int i=0;i<P->Ncols;i++)
    P->colMap[i] = i + globalAggOffset;

  if (P->offdNNZ) {
    //we now need to reorder the x vector for the halo, and shift the column indices
    int *col = (int *) calloc(P->offdNNZ,sizeof(int));
    for (int i=0;i<P->offdNNZ;i++)
      col[i] = P->offdCols[i]; //copy non-local column global ids

    //sort by global index
    std::sort(col,col+P->offdNNZ);

    //count unique non-local column ids
    P->NHalo = 0;
    for (int i=1;i<P->offdNNZ;i++)
      if (col[i]!=col[i-1])  col[++P->NHalo] = col[i];
    P->NHalo++; //number of unique columns

    P->Ncols += P->NHalo;

    //save global column ids in colMap
    P->colMap = (int *) realloc(P->colMap, P->Ncols*sizeof(int));
    for (int i=0; i<P->NHalo; i++)
      P->colMap[i+P->NlocalCols] = col[i];
    free(col);

    //shift the column indices to local indexing
    for (int i=0;i<P->offdNNZ;i++) {
      int gcol = P->offdCols[i];
      for (int m=P->NlocalCols;m<P->Ncols;m++) {
        if (gcol == P->colMap[m])
          P->offdCols[i] = m;
      }
    }
  }

  csrHaloSetup(P,globalAggStarts);

  // normalize the columns of P
  *nullCoarseA = (dfloat *) calloc(P->Ncols,sizeof(dfloat));

  //add local nonzeros
  for(int i=0; i<P->diagNNZ; i++)
    (*nullCoarseA)[P->diagCols[i]] += P->diagCoefs[i] * P->diagCoefs[i];

  dfloat *nnzSum, *recvNnzSum;
  if (P->NHalo) nnzSum = (dfloat *) calloc(P->NHalo,sizeof(dfloat));
  if (P->NsendTotal) recvNnzSum = (dfloat *) calloc(P->NsendTotal,sizeof(dfloat));

  //add the non-local non-zeros
  for (int i=0;i<P->offdNNZ;i++)
    nnzSum[P->offdCols[i]-P->NlocalCols] += P->offdCoefs[i] * P->offdCoefs[i];

  //do a reverse halo exchange
  int tag = 999;

  // initiate immediate send  and receives to each other process as needed
  int recvOffset = 0;
  int sendOffset = 0;
  int sendMessage = 0, recvMessage = 0;
  for(int r=0;r<size;++r){
    if (P->NsendTotal) {
      if(P->NsendPairs[r]) {
        MPI_Irecv(recvNnzSum+sendOffset, P->NsendPairs[r], MPI_DFLOAT, r, tag,
            MPI_COMM_WORLD, (MPI_Request*)P->haloSendRequests+sendMessage);
        sendOffset += P->NsendPairs[r];
        ++sendMessage;
      }
    }
    if (P->NrecvTotal) {
      if(P->NrecvPairs[r]){
        MPI_Isend(nnzSum+recvOffset, P->NrecvPairs[r], MPI_DFLOAT, r, tag,
            MPI_COMM_WORLD, (MPI_Request*)P->haloRecvRequests+recvMessage);
        recvOffset += P->NrecvPairs[r];
        ++recvMessage;
      }
    }
  }

  // Wait for all sent messages to have left and received messages to have arrived
  if (P->NrecvTotal) {
    MPI_Status *sendStatus = (MPI_Status*) calloc(P->NsendMessages, sizeof(MPI_Status));
    MPI_Waitall(P->NsendMessages, (MPI_Request*)P->haloSendRequests, sendStatus);
    free(sendStatus);
  }
  if (P->NsendTotal) {
    MPI_Status *recvStatus = (MPI_Status*) calloc(P->NrecvMessages, sizeof(MPI_Status));
    MPI_Waitall(P->NrecvMessages, (MPI_Request*)P->haloRecvRequests, recvStatus);
    free(recvStatus);
  }

  for(int i=0;i<P->NsendTotal;++i){
    // local index of outgoing element in halo exchange
    int id = P->haloElementList[i];

    (*nullCoarseA)[id] += recvNnzSum[i];
  }

  if (P->NHalo) free(nnzSum);

  for(int i=0; i<NCoarse; i++)
    (*nullCoarseA)[i] = sqrt((*nullCoarseA)[i]);

  csrHaloExchange(P, sizeof(dfloat), *nullCoarseA, P->sendBuffer, *nullCoarseA+P->NlocalCols);

  for(int i=0; i<P->diagNNZ; i++)
    P->diagCoefs[i] /= (*nullCoarseA)[P->diagCols[i]];
  for(int i=0; i<P->offdNNZ; i++)
    P->offdCoefs[i] /= (*nullCoarseA)[P->offdCols[i]];

  MPI_Barrier(MPI_COMM_WORLD);
  if (P->NsendTotal) free(recvNnzSum);

  return P;
}

typedef struct {

  int row;
  int col;
  dfloat val;
  int owner;

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

csr * transpose(agmgLevel* level, csr *A,
                int *globalRowStarts, int *globalColStarts){

  // MPI info
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  csr *At = (csr *) calloc(1,sizeof(csr));

  At->Nrows = A->Ncols-A->NHalo;
  At->Ncols = A->Nrows;
  At->diagNNZ   = A->diagNNZ; //local entries remain local

  At->NlocalCols = At->Ncols;

  At->diagRowStarts = (int *)   calloc(At->Nrows+1, sizeof(int));
  At->offdRowStarts = (int *)   calloc(At->Nrows+1, sizeof(int));

  //start with local entries
  if (A->diagNNZ) {
    At->diagCols      = (int *)   calloc(At->diagNNZ, sizeof(int));
    At->diagCoefs     = (dfloat *) calloc(At->diagNNZ, sizeof(dfloat));
  }

  // count the num of nonzeros per row for transpose
  for(int i=0; i<A->diagNNZ; i++){
    int row = A->diagCols[i];
    At->diagRowStarts[row+1]++;
  }

  // cumulative sum for rows
  for(int i=1; i<=At->Nrows; i++)
    At->diagRowStarts[i] += At->diagRowStarts[i-1];

  int *counter = (int *) calloc(At->Nrows+1,sizeof(int));
  for (int i=0; i<At->Nrows+1; i++)
    counter[i] = At->diagRowStarts[i];

  for(int i=0; i<A->Nrows; i++){
    const int Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];

    for(int jj=Jstart; jj<Jend; jj++){
      int row = A->diagCols[jj];
      At->diagCols[counter[row]]  = i;
      At->diagCoefs[counter[row]] = A->diagCoefs[jj];

      counter[row]++;
    }
  }
  free(counter);

  //record global indexing of columns
  At->colMap = (int *)   calloc(At->Ncols, sizeof(int));
  for (int i=0;i<At->Ncols;i++)
    At->colMap[i] = i + globalRowStarts[rank];

  //now the nonlocal entries. Need to reverse the halo exchange to send the nonzeros
  int tag = 999;

  nonzero_t *sendNonZeros;
  if (A->offdNNZ)
    sendNonZeros = (nonzero_t *) calloc(A->offdNNZ,sizeof(nonzero_t));

  int *Nsend = (int*) calloc(size, sizeof(int));
  int *Nrecv = (int*) calloc(size, sizeof(int));

  for(int r=0;r<size;r++) {
    Nsend[r] =0;
    Nrecv[r] =0;
  }

  // copy data from nonlocal entries into send buffer
  for(int i=0;i<A->Nrows;++i){
    for (int j=A->offdRowStarts[i];j<A->offdRowStarts[i+1];j++) {
      int col =  A->colMap[A->offdCols[j]]; //global ids
      for (int r=0;r<size;r++) { //find owner's rank
        if ((globalColStarts[r]-1<col) && (col < globalColStarts[r+1])) {
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
  if (A->offdNNZ)
    qsort(sendNonZeros, A->offdNNZ, sizeof(nonzero_t), compareNonZero);

  MPI_Alltoall(Nsend, 1, MPI_INT, Nrecv, 1, MPI_INT, MPI_COMM_WORLD);

  //count incoming nonzeros
  At->offdNNZ = 0;
  for (int r=0;r<size;r++)
    At->offdNNZ += Nrecv[r];

  nonzero_t *recvNonZeros;
  if (At->offdNNZ)
    recvNonZeros = (nonzero_t *) calloc(At->offdNNZ,sizeof(nonzero_t));

  // initiate immediate send and receives to each other process as needed
  int recvOffset = 0;
  int sendOffset = 0;
  int sendMessage = 0, recvMessage = 0;
  for(int r=0;r<size;++r){
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
    MPI_Status *sendStatus = (MPI_Status*) calloc(sendMessage, sizeof(MPI_Status));
    MPI_Waitall(sendMessage, (MPI_Request*)A->haloRecvRequests, sendStatus);
    free(sendStatus);
  }
  if (At->offdNNZ) {
    MPI_Status *recvStatus = (MPI_Status*) calloc(recvMessage, sizeof(MPI_Status));
    MPI_Waitall(recvMessage, (MPI_Request*)A->haloSendRequests, recvStatus);
    free(recvStatus);
  }
  if (A->offdNNZ) free(sendNonZeros);

  //free(Nsend); free(Nrecv);

  if (At->offdNNZ) {
    //sort recieved nonzeros by row and col
    qsort(recvNonZeros, At->offdNNZ, sizeof(nonzero_t), compareNonZero);

    At->offdCols  = (int *)   calloc(At->offdNNZ,sizeof(int));
    At->offdCoefs = (dfloat *) calloc(At->offdNNZ, sizeof(dfloat));

    //find row starts
    for(int n=0;n<At->offdNNZ;++n) {
      int row = recvNonZeros[n].row - globalColStarts[rank];
      At->offdRowStarts[row+1]++;
    }
    //cumulative sum
    for (int i=0;i<At->Nrows;i++)
      At->offdRowStarts[i+1] += At->offdRowStarts[i];

    //fill cols and coefs
    for (int i=0; i<At->Nrows; i++) {
      for (int j=At->offdRowStarts[i]; j<At->offdRowStarts[i+1]; j++) {
        At->offdCols[j]  = recvNonZeros[j].col;
        At->offdCoefs[j] = recvNonZeros[j].val;
      }
    }
    free(recvNonZeros);

    //we now need to reorder the x vector for the halo, and shift the column indices
    int *col = (int *) calloc(At->offdNNZ,sizeof(int));
    for (int n=0;n<At->offdNNZ;n++)
      col[n] = At->offdCols[n]; //copy non-local column global ids

    //sort by global index
    std::sort(col,col+At->offdNNZ);

    //count unique non-local column ids
    At->NHalo = 0;
    for (int n=1;n<At->offdNNZ;n++)
      if (col[n]!=col[n-1])  col[++At->NHalo] = col[n];
    At->NHalo++; //number of unique columns

    At->Ncols += At->NHalo;

    //save global column ids in colMap
    At->colMap = (int *) realloc(At->colMap,At->Ncols*sizeof(int));
    for (int n=0; n<At->NHalo; n++)
      At->colMap[n+At->NlocalCols] = col[n];
    free(col);

    //shift the column indices to local indexing
    for (int n=0;n<At->offdNNZ;n++) {
      int gcol = At->offdCols[n];
      for (int m=At->NlocalCols;m<At->Ncols;m++) {
        if (gcol == At->colMap[m])
          At->offdCols[n] = m;
      }
    }
  }

  csrHaloSetup(At,globalRowStarts);

  return At;
}

typedef struct {

  int coarseId;
  dfloat coef;

} pEntry_t;

typedef struct {

  int I;
  int J;
  dfloat coef;

} rapEntry_t;

int compareRAPEntries(const void *a, const void *b){
  rapEntry_t *pa = (rapEntry_t *) a;
  rapEntry_t *pb = (rapEntry_t *) b;

  if (pa->I < pb->I) return -1;
  if (pa->I > pb->I) return +1;

  if (pa->J < pb->J) return -1;
  if (pa->J > pb->J) return +1;

  return 0;
};

csr *galerkinProd(agmgLevel *level, csr *R, csr *A, csr *P){

  // MPI info
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int *globalAggStarts = level->globalAggStarts;
  int *globalRowStarts = level->globalRowStarts;

  int globalAggOffset = globalAggStarts[rank];

  /* Setup description of the MPI_RAPENTRY_T struct */
  MPI_Datatype MPI_RAPENTRY_T;
  MPI_Datatype oldtypes[2] = {MPI_INT, MPI_DFLOAT};
  int blockcounts[2] = {2, 1};

  MPI_Aint intext;
  MPI_Type_extent(MPI_INT, &intext);
  MPI_Aint  rapEntryoffsets[2] = {0, 2*intext};

  /* Now define structured type and commit it */
  MPI_Type_struct(2, blockcounts, rapEntryoffsets, oldtypes, &MPI_RAPENTRY_T);
  MPI_Type_commit(&MPI_RAPENTRY_T);

  //The galerkin product can be computed as
  // (RAP)_IJ = sum_{i in Agg_I} sum_{j in Agg_j} P_iI A_ij P_jJ
  // Since each row of P has only one entry, we can share the ncessary
  // P entries, form the products, and send them to their destination rank

  int N = A->Nrows;
  int M = A->Ncols;

  //printf("Level has %d rows, and is making %d aggregates\n", N, globalAggStarts[rank+1]-globalAggStarts[rank]);

  pEntry_t *PEntries;
  if (M) PEntries = (pEntry_t *) calloc(M,sizeof(pEntry_t));

  //record the entries of P that this rank has
  int cnt =0;
  for (int i=0;i<N;i++) {
    for (int j=P->diagRowStarts[i];j<P->diagRowStarts[i+1];j++) {
      PEntries[cnt].coarseId = P->diagCols[j] + globalAggOffset; //global ID
      PEntries[cnt].coef = P->diagCoefs[j];
      cnt++;
    }
    for (int j=P->offdRowStarts[i];j<P->offdRowStarts[i+1];j++) {
      PEntries[cnt].coarseId = P->colMap[P->offdCols[j]]; //global ID
      PEntries[cnt].coef = P->offdCoefs[j];
      cnt++;
    }
  }

  pEntry_t *entrySendBuffer;
  if (A->NsendTotal)
    entrySendBuffer = (pEntry_t *) calloc(A->NsendTotal,sizeof(pEntry_t));

  //fill in the entires of P needed in the halo
  csrHaloExchange(A, sizeof(pEntry_t), PEntries, entrySendBuffer, PEntries+A->NlocalCols);
  if (A->NsendTotal) free(entrySendBuffer);

  rapEntry_t *RAPEntries;
  int totalNNZ = A->diagNNZ+A->offdNNZ;
  if (totalNNZ) {
    RAPEntries = (rapEntry_t *) calloc(totalNNZ,sizeof(rapEntry_t));
  } else {
    RAPEntries = (rapEntry_t *) calloc(1,sizeof(rapEntry_t)); //MPI_AlltoAll doesnt like null pointers
  }

  //for the RAP products
  cnt =0;
  for (int i=0;i<N;i++) {
    for (int j=A->diagRowStarts[i];j<A->diagRowStarts[i+1];j++) {
      int col  = A->diagCols[j];
      dfloat coef = A->diagCoefs[j];

      RAPEntries[cnt].I = PEntries[i].coarseId;
      RAPEntries[cnt].J = PEntries[col].coarseId;
      RAPEntries[cnt].coef = coef*PEntries[i].coef*PEntries[col].coef;
      cnt++;
    }
  }
  for (int i=0;i<N;i++) {
    for (int j=A->offdRowStarts[i];j<A->offdRowStarts[i+1];j++) {
      int col  = A->offdCols[j];
      dfloat coef = A->offdCoefs[j];

      RAPEntries[cnt].I = PEntries[i].coarseId;
      RAPEntries[cnt].J = PEntries[col].coarseId;
      RAPEntries[cnt].coef = PEntries[i].coef*coef*PEntries[col].coef;
      cnt++;
    }
  }

  //sort entries by the coarse row and col
  if (totalNNZ) qsort(RAPEntries, totalNNZ, sizeof(rapEntry_t), compareRAPEntries);

  int *sendCounts = (int *) calloc(size,sizeof(int));
  int *recvCounts = (int *) calloc(size,sizeof(int));
  int *sendOffsets = (int *) calloc(size+1,sizeof(int));
  int *recvOffsets = (int *) calloc(size+1,sizeof(int));

  for(int i=0;i<totalNNZ;++i) {
    int id = RAPEntries[i].I;
    for (int r=0;r<size;r++) {
      if (globalAggStarts[r]-1<id && id < globalAggStarts[r+1])
        sendCounts[r]++;
    }
  }

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  int recvNtotal = 0;
  for(int r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r];
  }
  rapEntry_t *recvRAPEntries;
  if (recvNtotal) {
    recvRAPEntries = (rapEntry_t *) calloc(recvNtotal,sizeof(rapEntry_t));
  } else {
    recvRAPEntries = (rapEntry_t *) calloc(1,sizeof(rapEntry_t));//MPI_AlltoAll doesnt like null pointers
  }

  MPI_Alltoallv(RAPEntries, sendCounts, sendOffsets, MPI_RAPENTRY_T,
                recvRAPEntries, recvCounts, recvOffsets, MPI_RAPENTRY_T,
                MPI_COMM_WORLD);

  //sort entries by the coarse row and col
  if (recvNtotal) qsort(recvRAPEntries, recvNtotal, sizeof(rapEntry_t), compareRAPEntries);

  //count total number of nonzeros;
  int nnz =0;
  if (recvNtotal) nnz++;
  for (int i=1;i<recvNtotal;i++)
    if ((recvRAPEntries[i].I!=recvRAPEntries[i-1].I)||
          (recvRAPEntries[i].J!=recvRAPEntries[i-1].J)) nnz++;

  rapEntry_t *newRAPEntries;
  if (nnz) {
    newRAPEntries = (rapEntry_t *) calloc(nnz,sizeof(rapEntry_t));
  } else {
    newRAPEntries = (rapEntry_t *) calloc(1,sizeof(rapEntry_t));
  }

  //compress nonzeros
  nnz = 0;
  if (recvNtotal) newRAPEntries[nnz++] = recvRAPEntries[0];
  for (int i=1;i<recvNtotal;i++) {
    if ((recvRAPEntries[i].I!=recvRAPEntries[i-1].I)||
          (recvRAPEntries[i].J!=recvRAPEntries[i-1].J)) {
      newRAPEntries[nnz++] = recvRAPEntries[i];
    } else {
      newRAPEntries[nnz-1].coef += recvRAPEntries[i].coef;
    }
  }

  int numAggs = globalAggStarts[rank+1]-globalAggStarts[rank]; //local number of aggregates

  csr *RAP = (csr*) calloc(1,sizeof(csr));

  RAP->Nrows = numAggs;
  RAP->Ncols = numAggs;

  RAP->NlocalCols = numAggs;

  RAP->diagRowStarts = (int *) calloc(numAggs+1, sizeof(int));
  RAP->offdRowStarts = (int *) calloc(numAggs+1, sizeof(int));

  for (int n=0;n<nnz;n++) {
    int row = newRAPEntries[n].I - globalAggOffset;
    if ((newRAPEntries[n].J > globalAggStarts[rank]-1)&&
          (newRAPEntries[n].J < globalAggStarts[rank+1])) {
      RAP->diagRowStarts[row+1]++;
    } else {
      RAP->offdRowStarts[row+1]++;
    }
  }

  // cumulative sum
  for(int i=0; i<numAggs; i++) {
    RAP->diagRowStarts[i+1] += RAP->diagRowStarts[i];
    RAP->offdRowStarts[i+1] += RAP->offdRowStarts[i];
  }
  RAP->diagNNZ = RAP->diagRowStarts[numAggs];
  RAP->offdNNZ = RAP->offdRowStarts[numAggs];

  int *diagCols;
  dfloat *diagCoefs;
  if (RAP->diagNNZ) {
    RAP->diagCols  = (int *)   calloc(RAP->diagNNZ, sizeof(int));
    RAP->diagCoefs = (dfloat *) calloc(RAP->diagNNZ, sizeof(dfloat));
    diagCols  = (int *)   calloc(RAP->diagNNZ, sizeof(int));
    diagCoefs = (dfloat *) calloc(RAP->diagNNZ, sizeof(dfloat));
  }
  if (RAP->offdNNZ) {
    RAP->offdCols  = (int *)   calloc(RAP->offdNNZ,sizeof(int));
    RAP->offdCoefs = (dfloat *) calloc(RAP->offdNNZ, sizeof(dfloat));
  }

  int diagCnt =0;
  int offdCnt =0;
  for (int n=0;n<nnz;n++) {
    if ((newRAPEntries[n].J > globalAggStarts[rank]-1)&&
          (newRAPEntries[n].J < globalAggStarts[rank+1])) {
      diagCols[diagCnt]  = newRAPEntries[n].J - globalAggOffset;
      diagCoefs[diagCnt] = newRAPEntries[n].coef;
      diagCnt++;
    } else {
      RAP->offdCols[offdCnt]  = newRAPEntries[n].J;
      RAP->offdCoefs[offdCnt] = newRAPEntries[n].coef;
      offdCnt++;
    }
  }

  //move diagonal entries first
  for (int i=0;i<RAP->Nrows;i++) {
    int start = RAP->diagRowStarts[i];
    int cnt = 1;
    for (int j=RAP->diagRowStarts[i]; j<RAP->diagRowStarts[i+1]; j++) {
      if (diagCols[j] == i) { //move diagonal to first entry
        RAP->diagCols[start] = diagCols[j];
        RAP->diagCoefs[start] = diagCoefs[j];
      } else {
        RAP->diagCols[start+cnt] = diagCols[j];
        RAP->diagCoefs[start+cnt] = diagCoefs[j];
        cnt++;
      }
    }
  }

  //record global indexing of columns
  RAP->colMap = (int *)   calloc(RAP->Ncols, sizeof(int));
  for (int i=0;i<RAP->Ncols;i++)
    RAP->colMap[i] = i + globalAggOffset;

  if (RAP->offdNNZ) {
    //we now need to reorder the x vector for the halo, and shift the column indices
    int *col = (int *) calloc(RAP->offdNNZ,sizeof(int));
    for (int n=0;n<RAP->offdNNZ;n++)
      col[n] = RAP->offdCols[n]; //copy non-local column global ids

    //sort by global index
    std::sort(col,col+RAP->offdNNZ);

    //count unique non-local column ids
    RAP->NHalo = 0;
    for (int n=1;n<RAP->offdNNZ;n++)
      if (col[n]!=col[n-1])  col[++RAP->NHalo] = col[n];
    RAP->NHalo++; //number of unique columns

    RAP->Ncols += RAP->NHalo;

    //save global column ids in colMap
    RAP->colMap = (int *) realloc(RAP->colMap,RAP->Ncols*sizeof(int));
    for (int n=0; n<RAP->NHalo; n++)
      RAP->colMap[n+RAP->NlocalCols] = col[n];

    //shift the column indices to local indexing
    for (int n=0;n<RAP->offdNNZ;n++) {
      int gcol = RAP->offdCols[n];
      for (int m=RAP->NlocalCols;m<RAP->Ncols;m++) {
        if (gcol == RAP->colMap[m])
          RAP->offdCols[n] = m;
      }
    }
    free(col);
  }
  csrHaloSetup(RAP,globalAggStarts);

  //clean up
  MPI_Type_free(&MPI_RAPENTRY_T);

  if (M) free(PEntries);
  free(sendCounts); free(recvCounts);
  free(sendOffsets); free(recvOffsets);
  if (RAP->diagNNZ) {
    free(diagCols);
    free(diagCoefs);
  }
  free(RAPEntries);
  free(newRAPEntries);
  free(recvRAPEntries);

  return RAP;
}


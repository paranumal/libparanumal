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

csr *strong_graph(csr *A, dfloat threshold);
bool customLess(int smax, dfloat rmax, hlong imax, int s, dfloat r, hlong i);
hlong *form_aggregates(agmgLevel *level, csr *C);
void find_aggregate_owners(agmgLevel *level, hlong* FineToCoarse, setupAide options);
csr *construct_interpolator(agmgLevel *level, hlong *FineToCoarse, dfloat **nullCoarseA);
csr *transpose(agmgLevel* level, csr *A, hlong *globalRowStarts, hlong *globalColStarts);
csr *galerkinProd(agmgLevel *level, csr *R, csr *A, csr *P);
void coarsenAgmgLevel(agmgLevel *level, csr **coarseA, csr **P, csr **R, dfloat **nullCoarseA, setupAide options);


void agmgSetup(parAlmond_t *parAlmond, csr *A, dfloat *nullA, hlong *globalRowStarts, setupAide options){

  int rank, size;
  rank = agmg::rank;
  size = agmg::size;

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
  int ChebyshevIterations=2; //default to degree 2
  if (options.compareArgs("PARALMOND SMOOTHER", "CHEBYSHEV")) {
    smoothType = CHEBYSHEV;
    options.getArgs("PARALMOND CHEBYSHEV DEGREE", ChebyshevIterations);
  } else { //default to DAMPED_JACOBI
    smoothType = DAMPED_JACOBI;
  }
  levels[lev]->ChebyshevIterations = ChebyshevIterations;

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
  levels[lev]->globalRowStarts = (hlong *) calloc(size+1,sizeof(hlong));
  for (int r=0;r<size+1;r++)
      levels[lev]->globalRowStarts[r] = globalRowStarts[r];

  hlong localSize = (hlong) levels[lev]->A->Nrows;
  hlong globalSize = 0;
  MPI_Allreduce(&localSize, &globalSize, 1, MPI_HLONG, MPI_SUM, agmg::comm);

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
    
    levels[lev+1]->ChebyshevIterations = ChebyshevIterations;

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

    const hlong localCoarseDim = (hlong) levels[lev+1]->A->Nrows;
    hlong globalCoarseSize;
    MPI_Allreduce(&localCoarseDim, &globalCoarseSize, 1, MPI_HLONG, MPI_SUM, agmg::comm);

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
    dlong N = levels[n]->Nrows;
    dlong M = levels[n]->Ncols;

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
  //buffer for innerproducts in kcycle
  dlong numBlocks = ((levels[0]->Nrows+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  parAlmond->rho  = (dfloat*) calloc(3*numBlocks,sizeof(dfloat));
  parAlmond->o_rho  = device.malloc(3*numBlocks*sizeof(dfloat), parAlmond->rho); 
}

void parAlmondReport(parAlmond_t *parAlmond) {

  int rank, size;
  rank = agmg::rank;
  size = agmg::size;

  if(rank==0) {
    printf("------------------ParAlmond Report-----------------------------------\n");
    printf("---------------------------------------------------------------------\n");
    printf("level| active ranks |   dimension   |  nnzs         |  nnz/row      |\n");
    printf("     |              | (min,max,avg) | (min,max,avg) | (min,max,avg) |\n");
    printf("---------------------------------------------------------------------\n");
  }

  for(int lev=0; lev<parAlmond->numLevels; lev++){

    dlong Nrows = parAlmond->levels[lev]->Nrows;
    hlong hNrows = (hlong) parAlmond->levels[lev]->Nrows;

    int active = (Nrows>0) ? 1:0;
    int totalActive=0;
    MPI_Allreduce(&active, &totalActive, 1, MPI_INT, MPI_SUM, agmg::comm);

    dlong minNrows=0, maxNrows=0;
    hlong totalNrows=0;
    dfloat avgNrows;
    MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_DLONG, MPI_MAX, agmg::comm);
    MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, agmg::comm);
    avgNrows = (dfloat) totalNrows/totalActive;

    if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
    MPI_Allreduce(&Nrows, &minNrows, 1, MPI_DLONG, MPI_MIN, agmg::comm);


    long long int nnz;
    if (parAlmond->levels[lev]->A)
      nnz = parAlmond->levels[lev]->A->diagNNZ+parAlmond->levels[lev]->A->offdNNZ;
    else
      nnz =0;
    long long int minNnz=0, maxNnz=0, totalNnz=0;
    dfloat avgNnz;
    MPI_Allreduce(&nnz, &maxNnz, 1, MPI_LONG_LONG_INT, MPI_MAX, agmg::comm);
    MPI_Allreduce(&nnz, &totalNnz, 1, MPI_LONG_LONG_INT, MPI_SUM, agmg::comm);
    avgNnz = (dfloat) totalNnz/totalActive;

    if (nnz==0) nnz = maxNnz; //set this so it's ignored for the global min
    MPI_Allreduce(&nnz, &minNnz, 1, MPI_LONG_LONG_INT, MPI_MIN, agmg::comm);

    Nrows = parAlmond->levels[lev]->Nrows;
    dfloat nnzPerRow = (Nrows==0) ? 0 : (dfloat) nnz/Nrows;
    dfloat minNnzPerRow=0, maxNnzPerRow=0, avgNnzPerRow=0;
    MPI_Allreduce(&nnzPerRow, &maxNnzPerRow, 1, MPI_DFLOAT, MPI_MAX, agmg::comm);
    MPI_Allreduce(&nnzPerRow, &avgNnzPerRow, 1, MPI_DFLOAT, MPI_SUM, agmg::comm);
    avgNnzPerRow /= totalActive;

    if (Nrows==0) nnzPerRow = maxNnzPerRow;
    MPI_Allreduce(&nnzPerRow, &minNnzPerRow, 1, MPI_DFLOAT, MPI_MIN, agmg::comm);

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
void coarsenAgmgLevel(agmgLevel *level, csr **coarseA, csr **P, csr **R, dfloat **nullCoarseA, setupAide options){

  // establish the graph of strong connections
  level->threshold = 0.5;

  csr *C = strong_graph(level->A, level->threshold);

  hlong *FineToCoarse = form_aggregates(level, C);

  find_aggregate_owners(level,FineToCoarse,options);

  *P = construct_interpolator(level, FineToCoarse, nullCoarseA);
  *R = transpose(level, *P, level->globalRowStarts, level->globalAggStarts);
  *coarseA = galerkinProd(level, *R, level->A, *P);
}

csr * strong_graph(csr *A, dfloat threshold){

  const dlong N = A->Nrows;
  const dlong M = A->Ncols;

  csr *C = (csr *) calloc(1, sizeof(csr));

  C->Nrows = N;
  C->Ncols = M;

  C->diagRowStarts = (dlong *) calloc(N+1,sizeof(dlong));
  C->offdRowStarts = (dlong *) calloc(N+1,sizeof(dlong));

  dfloat *maxOD;
  if (N) maxOD = (dfloat *) calloc(N,sizeof(dfloat));

  //store the diagonal of A for all needed columns
  dfloat *diagA = (dfloat *) calloc(M,sizeof(dfloat));
  for (dlong i=0;i<N;i++)
    diagA[i] = A->diagCoefs[A->diagRowStarts[i]];
  csrHaloExchange(A, sizeof(dfloat), diagA, A->sendBuffer, diagA+A->NlocalCols);

  #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    dfloat sign = (diagA[i] >= 0) ? 1:-1;
    dfloat Aii = fabs(diagA[i]);

    //find maxOD
    //local entries
    dlong Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];
    for(dlong jj= Jstart+1; jj<Jend; jj++){
      dlong col = A->diagCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->diagCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > maxOD[i]) maxOD[i] = OD;
    }
    //non-local entries
    Jstart = A->offdRowStarts[i], Jend = A->offdRowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      dlong col = A->offdCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->offdCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > maxOD[i]) maxOD[i] = OD;
    }

    int diag_strong_per_row = 1; // diagonal entry
    //local entries
    Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];
    for(dlong jj = Jstart+1; jj<Jend; jj++){
      dlong col = A->diagCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->diagCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD[i]) diag_strong_per_row++;
    }
    int offd_strong_per_row = 0;
    //non-local entries
    Jstart = A->offdRowStarts[i], Jend = A->offdRowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      dlong col = A->offdCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->offdCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD[i]) offd_strong_per_row++;
    }

    C->diagRowStarts[i+1] = diag_strong_per_row;
    C->offdRowStarts[i+1] = offd_strong_per_row;
  }

  // cumulative sum
  for(dlong i=1; i<N+1 ; i++) {
    C->diagRowStarts[i] += C->diagRowStarts[i-1];
    C->offdRowStarts[i] += C->offdRowStarts[i-1];
  }

  C->diagNNZ = C->diagRowStarts[N];
  C->offdNNZ = C->offdRowStarts[N];

  if (C->diagNNZ) C->diagCols = (dlong *) calloc(C->diagNNZ, sizeof(dlong));
  if (C->offdNNZ) C->offdCols = (dlong *) calloc(C->offdNNZ, sizeof(dlong));

  // fill in the columns for strong connections
  #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    dfloat sign = (diagA[i] >= 0) ? 1:-1;
    dfloat Aii = fabs(diagA[i]);

    dlong diagCounter = C->diagRowStarts[i];
    dlong offdCounter = C->offdRowStarts[i];

    //local entries
    C->diagCols[diagCounter++] = i;// diag entry
    dlong Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];
    for(dlong jj = Jstart+1; jj<Jend; jj++){
      dlong col = A->diagCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->diagCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD[i])
        C->diagCols[diagCounter++] = A->diagCols[jj];
    }
    Jstart = A->offdRowStarts[i], Jend = A->offdRowStarts[i+1];
    for(dlong jj = Jstart; jj<Jend; jj++){
      dlong col = A->offdCols[jj];
      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->offdCoefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD[i])
        C->offdCols[offdCounter++] = A->offdCols[jj];
    }
  }
  if(N) free(maxOD);

  return C;
}

bool customLess(int smax, dfloat rmax, hlong imax, int s, dfloat r, hlong i){

  if(s > smax) return true;
  if(smax > s) return false;

  if(r > rmax) return true;
  if(rmax > r) return false;

  if(i > imax) return true;
  if(i < imax) return false;

  return false;
}

hlong * form_aggregates(agmgLevel *level, csr *C){

  int rank, size;
  rank = agmg::rank;
  size = agmg::size;

  const dlong N   = C->Nrows;
  const dlong M   = C->Ncols;
  const dlong diagNNZ = C->diagNNZ;
  const dlong offdNNZ = C->offdNNZ;

  hlong *FineToCoarse = (hlong *) calloc(M, sizeof(hlong));
  for (dlong i =0;i<M;i++) FineToCoarse[i] = -1;

  dfloat *rands  = (dfloat *) calloc(M, sizeof(dfloat));
  int   *states = (int *)   calloc(M, sizeof(int));

  dfloat *Tr = (dfloat *) calloc(M, sizeof(dfloat));
  int    *Ts = (int *)    calloc(M, sizeof(int));
  hlong  *Ti = (hlong *)  calloc(M, sizeof(hlong));
  hlong  *Tc = (hlong *)  calloc(M, sizeof(hlong));

  csr *A = level->A;
  hlong *globalRowStarts = level->globalRowStarts;

  int    *intSendBuffer;
  hlong  *hlongSendBuffer;
  dfloat *dfloatSendBuffer;
  if (level->A->NsendTotal) {
    intSendBuffer = (int *) calloc(A->NsendTotal,sizeof(int));
    hlongSendBuffer = (hlong *) calloc(A->NsendTotal,sizeof(hlong));
    dfloatSendBuffer = (dfloat *) calloc(A->NsendTotal,sizeof(dfloat));
  }

  for(dlong i=0; i<N; i++)
    rands[i] = (dfloat) drand48();

  for(dlong i=0; i<N; i++)
    states[i] = 0;

  // add the number of non-zeros in each column
  //local non-zeros
  for(dlong i=0; i<diagNNZ; i++)
    rands[C->diagCols[i]] += 1.;

  int *nnzCnt, *recvNnzCnt;
  if (A->NHalo) nnzCnt = (int *) calloc(A->NHalo,sizeof(int));
  if (A->NsendTotal) recvNnzCnt = (int *) calloc(A->NsendTotal,sizeof(int));

  //count the non-local non-zeros
  for (dlong i=0;i<offdNNZ;i++)
    nnzCnt[C->offdCols[i]-A->NlocalCols]++;

  //do a reverse halo exchange
  int tag = 999;

  // initiate immediate send  and receives to each other process as needed
  dlong recvOffset = 0;
  dlong sendOffset = 0;
  int sendMessage = 0, recvMessage = 0;
  for(int r=0;r<size;++r){
    if (A->NsendTotal) {
      if(A->NsendPairs[r]) {
        MPI_Irecv(recvNnzCnt+sendOffset, A->NsendPairs[r], MPI_INT, r, tag,
            agmg::comm, (MPI_Request*)A->haloSendRequests+sendMessage);
        sendOffset += A->NsendPairs[r];
        ++sendMessage;
      }
    }
    if (A->NrecvTotal) {
      if(A->NrecvPairs[r]){
        MPI_Isend(nnzCnt+recvOffset, A->NrecvPairs[r], MPI_INT, r, tag,
            agmg::comm, (MPI_Request*)A->haloRecvRequests+recvMessage);
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
    dlong id = A->haloElementList[i];

    rands[id] += recvNnzCnt[i];
  }

  if (A->NHalo) free(nnzCnt);
  if (A->NsendTotal) free(recvNnzCnt);

  //share randomizer values
  csrHaloExchange(A, sizeof(dfloat), rands, dfloatSendBuffer, rands+A->NlocalCols);



  hlong done = 0;
  while(!done){
    // first neighbours
    #pragma omp parallel for
    for(dlong i=0; i<N; i++){

      int smax = states[i];
      dfloat rmax = rands[i];
      hlong imax = i + globalRowStarts[rank];

      if(smax != 1){
        //local entries
        for(dlong jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
          const dlong col = C->diagCols[jj];
          if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
            smax = states[col];
            rmax = rands[col];
            imax = col + globalRowStarts[rank];
          }
        }
        //nonlocal entries
        for(dlong jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
          const dlong col = C->offdCols[jj];
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
    csrHaloExchange(A, sizeof(hlong), Ti, hlongSendBuffer, Ti+A->NlocalCols);

    // second neighbours
    #pragma omp parallel for
    for(dlong i=0; i<N; i++){
      int    smax = Ts[i];
      dfloat rmax = Tr[i];
      hlong  imax = Ti[i];

      //local entries
      for(dlong jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
        const dlong col = C->diagCols[jj];
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }
      //nonlocal entries
      for(dlong jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
        const dlong col = C->offdCols[jj];
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
    hlong cnt = std::count(states, states+N, 0);
    MPI_Allreduce(&cnt,&done,1,MPI_HLONG, MPI_SUM,agmg::comm);
    done = (done == 0) ? 1 : 0;
  }

  dlong numAggs = 0;
  dlong *gNumAggs = (dlong *) calloc(size,sizeof(dlong));
  level->globalAggStarts = (hlong *) calloc(size+1,sizeof(hlong));
  // count the coarse nodes/aggregates
  for(dlong i=0; i<N; i++)
    if(states[i] == 1) numAggs++;

  MPI_Allgather(&numAggs,1,MPI_DLONG,gNumAggs,1,MPI_DLONG,agmg::comm);

  level->globalAggStarts[0] = 0;
  for (int r=0;r<size;r++)
    level->globalAggStarts[r+1] = level->globalAggStarts[r] + gNumAggs[r];

  numAggs = 0;
  // enumerate the coarse nodes/aggregates
  for(dlong i=0; i<N; i++)
    if(states[i] == 1)
      FineToCoarse[i] = level->globalAggStarts[rank] + numAggs++;

  //share the initial aggregate flags
  csrHaloExchange(A, sizeof(hlong), FineToCoarse, hlongSendBuffer, FineToCoarse+A->NlocalCols);

  // form the aggregates
  #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int   smax = states[i];
    dfloat rmax = rands[i];
    hlong  imax = i + globalRowStarts[rank];
    hlong  cmax = FineToCoarse[i];

    if(smax != 1){
      //local entries
      for(dlong jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
        const dlong col = C->diagCols[jj];
        if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
          smax = states[col];
          rmax = rands[col];
          imax = col + globalRowStarts[rank];
          cmax = FineToCoarse[col];
        }
      }
      //nonlocal entries
      for(dlong jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
        const dlong col = C->offdCols[jj];
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

  csrHaloExchange(A, sizeof(hlong), FineToCoarse, hlongSendBuffer, FineToCoarse+A->NlocalCols);
  csrHaloExchange(A, sizeof(dfloat), Tr, dfloatSendBuffer, Tr+A->NlocalCols);
  csrHaloExchange(A, sizeof(int), Ts, intSendBuffer, Ts+A->NlocalCols);
  csrHaloExchange(A, sizeof(hlong), Ti, hlongSendBuffer, Ti+A->NlocalCols);
  csrHaloExchange(A, sizeof(hlong), Tc, hlongSendBuffer, Tc+A->NlocalCols);

  // second neighbours
  #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int    smax = Ts[i];
    dfloat rmax = Tr[i];
    hlong  imax = Ti[i];
    hlong  cmax = Tc[i];

    //local entries
    for(dlong jj=C->diagRowStarts[i]+1;jj<C->diagRowStarts[i+1];jj++){
      const dlong col = C->diagCols[jj];
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }
    //nonlocal entries
    for(dlong jj=C->offdRowStarts[i];jj<C->offdRowStarts[i+1];jj++){
      const dlong col = C->offdCols[jj];
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

  csrHaloExchange(A, sizeof(hlong), FineToCoarse, hlongSendBuffer, FineToCoarse+A->NlocalCols);

  free(rands);
  free(states);
  free(Tr);
  free(Ts);
  free(Ti);
  free(Tc);
  if (level->A->NsendTotal) {
    free(intSendBuffer);
    free(hlongSendBuffer);
    free(dfloatSendBuffer);
  }

  //TODO maybe free C here?

  return FineToCoarse;
}

typedef struct {

  dlong fineId;
  hlong coarseId;
  hlong newCoarseId;

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

void find_aggregate_owners(agmgLevel *level, hlong* FineToCoarse, setupAide options) {
  // MPI info
  int rank, size;
  rank = agmg::rank;
  size = agmg::size;

  dlong N = level->A->Nrows;

  //Need to establish 'ownership' of aggregates
  
  //Keep the current partitioning for STRONGNODES. 
  // The rank that had the strong node for each aggregate owns the aggregate
  if (options.compareArgs("PARALMOND PARTITION", "STRONGNODES")) return;

  //populate aggregate array
  hlong gNumAggs = level->globalAggStarts[size]; //total number of aggregates
  
  parallelAggregate_t *sendAggs;
  if (N) 
    sendAggs = (parallelAggregate_t *) calloc(N,sizeof(parallelAggregate_t));
  else 
    sendAggs = (parallelAggregate_t *) calloc(1,sizeof(parallelAggregate_t));

  for (dlong i=0;i<N;i++) {
    sendAggs[i].fineId = i;
    sendAggs[i].originRank = rank;

    sendAggs[i].coarseId = FineToCoarse[i];

    //set a temporary owner. Evenly distibute aggregates amoungst ranks
    sendAggs[i].ownerRank = (int) (FineToCoarse[i]*size)/gNumAggs;
  }

  // Make the MPI_PARALLEL_AGGREGATE data type
  MPI_Datatype MPI_PARALLEL_AGGREGATE;
  MPI_Datatype dtype[5] = {MPI_DLONG, MPI_HLONG, MPI_HLONG, MPI_INT, MPI_INT};
  int blength[5] = {1, 1, 1, 1, 1};
  MPI_Aint addr[5], displ[5];
  MPI_Get_address ( &(sendAggs[0]            ), addr+0);
  MPI_Get_address ( &(sendAggs[0].coarseId   ), addr+1);
  MPI_Get_address ( &(sendAggs[0].newCoarseId), addr+2);
  MPI_Get_address ( &(sendAggs[0].originRank ), addr+3);
  MPI_Get_address ( &(sendAggs[0].ownerRank  ), addr+4);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  displ[4] = addr[4] - addr[0];
  MPI_Type_create_struct (5, blength, displ, dtype, &MPI_PARALLEL_AGGREGATE);
  MPI_Type_commit (&MPI_PARALLEL_AGGREGATE);

  //sort by owning rank for all_reduce
  qsort(sendAggs, N, sizeof(parallelAggregate_t), compareOwner);

  int *sendCounts = (int *) calloc(size,sizeof(int));
  int *recvCounts = (int *) calloc(size,sizeof(int));
  int *sendOffsets = (int *) calloc(size+1,sizeof(int));
  int *recvOffsets = (int *) calloc(size+1,sizeof(int));

  for(dlong i=0;i<N;++i)
    sendCounts[sendAggs[i].ownerRank]++;

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, agmg::comm);

  // find send and recv offsets for gather
  dlong recvNtotal = 0;
  for(int r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r];
  }
  parallelAggregate_t *recvAggs = (parallelAggregate_t *) calloc(recvNtotal,sizeof(parallelAggregate_t));

  MPI_Alltoallv(sendAggs, sendCounts, sendOffsets, MPI_PARALLEL_AGGREGATE,
                recvAggs, recvCounts, recvOffsets, MPI_PARALLEL_AGGREGATE,
                agmg::comm);

  //sort by coarse aggregate number, and then by original rank
  qsort(recvAggs, recvNtotal, sizeof(parallelAggregate_t), compareAgg);

  //count the number of unique aggregates here
  dlong NumUniqueAggs =0;
  if (recvNtotal) NumUniqueAggs++;
  for (dlong i=1;i<recvNtotal;i++)
    if(recvAggs[i].coarseId!=recvAggs[i-1].coarseId) NumUniqueAggs++;

  //get their locations in the array
  dlong *aggStarts;
  if (NumUniqueAggs)
    aggStarts = (dlong *) calloc(NumUniqueAggs+1,sizeof(dlong));
  dlong cnt = 1;
  for (dlong i=1;i<recvNtotal;i++)
    if(recvAggs[i].coarseId!=recvAggs[i-1].coarseId) aggStarts[cnt++] = i;
  aggStarts[NumUniqueAggs] = recvNtotal;


  if (options.compareArgs("PARALMOND PARTITION", "DISTRIBUTED")) { //rank that contributes most to the aggregate ownes it
    //use a random dfloat for each rank to break ties.
    dfloat rand = (dfloat) drand48();
    dfloat *gRands = (dfloat *) calloc(size,sizeof(dfloat));
    MPI_Allgather(&rand, 1, MPI_DFLOAT, gRands, 1, MPI_DFLOAT, agmg::comm);

    //determine the aggregates majority owner
    int *rankCounts = (int *) calloc(size,sizeof(int));
    for (dlong n=0;n<NumUniqueAggs;n++) {
      //populate randomizer
      for (int r=0;r<size;r++)
        rankCounts[r] = gRands[r];

      //count the number of contributions to the aggregate from the separate ranks
      for (dlong i=aggStarts[n];i<aggStarts[n+1];i++)
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
      for (dlong i=aggStarts[n];i<aggStarts[n+1];i++)
        recvAggs[i].ownerRank = ownerRank;
    }
    free(gRands); free(rankCounts);
  } else { //default SATURATE: always choose the lowest rank to own the aggregate
    for (dlong n=0;n<NumUniqueAggs;n++) {
      
      int minrank = size;

      //count the number of contributions to the aggregate from the separate ranks
      for (dlong i=aggStarts[n];i<aggStarts[n+1];i++){

        minrank = (recvAggs[i].originRank<minrank) ? recvAggs[i].originRank : minrank;
      }

      //set this aggregate's owner
      for (dlong i=aggStarts[n];i<aggStarts[n+1];i++)
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

  for(dlong i=0;i<recvNtotal;++i)
    newSendCounts[recvAggs[i].ownerRank]++;

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(newSendCounts, 1, MPI_INT, newRecvCounts, 1, MPI_INT, agmg::comm);

  // find send and recv offsets for gather
  dlong newRecvNtotal = 0;
  for(int r=0;r<size;++r){
    newSendOffsets[r+1] = newSendOffsets[r] + newSendCounts[r];
    newRecvOffsets[r+1] = newRecvOffsets[r] + newRecvCounts[r];
    newRecvNtotal += newRecvCounts[r];
  }
  parallelAggregate_t *newRecvAggs = (parallelAggregate_t *) calloc(newRecvNtotal,sizeof(parallelAggregate_t));

  MPI_Alltoallv(   recvAggs, newSendCounts, newSendOffsets, MPI_PARALLEL_AGGREGATE,
                newRecvAggs, newRecvCounts, newRecvOffsets, MPI_PARALLEL_AGGREGATE,
                agmg::comm);

  //sort by coarse aggregate number, and then by original rank
  qsort(newRecvAggs, newRecvNtotal, sizeof(parallelAggregate_t), compareAgg);

  //count the number of unique aggregates this rank owns
  dlong numAggs = 0;
  if (newRecvNtotal) numAggs++;
  for (dlong i=1;i<newRecvNtotal;i++)
    if(newRecvAggs[i].coarseId!=newRecvAggs[i-1].coarseId) numAggs++;

  //determine a global numbering of the aggregates
  dlong *lNumAggs = (dlong*) calloc(size,sizeof(dlong));
  MPI_Allgather(&numAggs, 1, MPI_DLONG, lNumAggs, 1, MPI_INT, agmg::comm);

  level->globalAggStarts[0] = 0;
  for (int r=0;r<size;r++)
    level->globalAggStarts[r+1] = level->globalAggStarts[r] + lNumAggs[r];

  //set the new global coarse index
  cnt = level->globalAggStarts[rank];
  if (newRecvNtotal) newRecvAggs[0].newCoarseId = cnt;
  for (dlong i=1;i<newRecvNtotal;i++) {
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

  for(dlong i=0;i<newRecvNtotal;++i)
    sendCounts[newRecvAggs[i].originRank]++;

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, agmg::comm);

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
                agmg::comm);

  //clean up
  MPI_Barrier(agmg::comm);
  MPI_Type_free(&MPI_PARALLEL_AGGREGATE);

  free(recvAggs);
  free(sendCounts);  free(recvCounts);
  free(sendOffsets); free(recvOffsets);
  free(newRecvAggs);
  free(newSendCounts);  free(newRecvCounts);
  free(newSendOffsets); free(newRecvOffsets);

  //record the new FineToCoarse map
  for (dlong i=0;i<N;i++)
    FineToCoarse[sendAggs[i].fineId] = sendAggs[i].newCoarseId;

  free(sendAggs);
}


csr *construct_interpolator(agmgLevel *level, hlong *FineToCoarse, dfloat **nullCoarseA){
  // MPI info
  int rank, size;
  rank = agmg::rank;
  size = agmg::size;

  const dlong N = level->A->Nrows;
  // const dlong M = level->A->Ncols;

  hlong *globalAggStarts = level->globalAggStarts;

  const hlong globalAggOffset = level->globalAggStarts[rank];
  const dlong NCoarse = (dlong) (globalAggStarts[rank+1]-globalAggStarts[rank]); //local num agg

  csr* P = (csr *) calloc(1, sizeof(csr));

  P->Nrows = N;
  P->Ncols = NCoarse;

  P->NlocalCols = NCoarse;
  P->NHalo = 0;

  P->diagRowStarts = (dlong *) calloc(N+1, sizeof(dlong));
  P->offdRowStarts = (dlong *) calloc(N+1, sizeof(dlong));

  // each row has exactly one nonzero per row
  P->diagNNZ =0;
  P->offdNNZ =0;
  for(dlong i=0; i<N; i++) {
    hlong col = FineToCoarse[i];
    if ((col>globalAggOffset-1)&&(col<globalAggOffset+NCoarse)) {
      P->diagNNZ++;
      P->diagRowStarts[i+1]++;
    } else {
      P->offdNNZ++;
      P->offdRowStarts[i+1]++;
    }
  }
  for(dlong i=0; i<N; i++) {
    P->diagRowStarts[i+1] += P->diagRowStarts[i];
    P->offdRowStarts[i+1] += P->offdRowStarts[i];
  }

  if (P->diagNNZ) {
    P->diagCols  = (dlong *)  calloc(P->diagNNZ, sizeof(dlong));
    P->diagCoefs = (dfloat *) calloc(P->diagNNZ, sizeof(dfloat));
  }
  hlong *offdCols;
  if (P->offdNNZ) {
    offdCols  = (hlong *)  calloc(P->offdNNZ, sizeof(hlong));
    P->offdCols  = (dlong *)  calloc(P->offdNNZ, sizeof(dlong));
    P->offdCoefs = (dfloat *) calloc(P->offdNNZ, sizeof(dfloat));
  }

  dlong diagCnt = 0;
  dlong offdCnt = 0;
  for(dlong i=0; i<N; i++) {
    hlong col = FineToCoarse[i];
    if ((col>globalAggStarts[rank]-1)&&(col<globalAggStarts[rank+1])) {
      P->diagCols[diagCnt] = (dlong) (col - globalAggOffset); //local index
      P->diagCoefs[diagCnt++] = level->A->null[i];
    } else {
      offdCols[offdCnt] = col;
      P->offdCoefs[offdCnt++] = level->A->null[i];
    }
  }

  //record global indexing of columns
  P->colMap = (hlong *)   calloc(P->Ncols, sizeof(hlong));
  for (dlong i=0;i<P->Ncols;i++)
    P->colMap[i] = i + globalAggOffset;

  if (P->offdNNZ) {
    //we now need to reorder the x vector for the halo, and shift the column indices
    hlong *col = (hlong *) calloc(P->offdNNZ,sizeof(hlong));
    for (dlong i=0;i<P->offdNNZ;i++)
      col[i] = offdCols[i]; //copy non-local column global ids

    //sort by global index
    std::sort(col,col+P->offdNNZ);

    //count unique non-local column ids
    P->NHalo = 0;
    for (dlong i=1;i<P->offdNNZ;i++)
      if (col[i]!=col[i-1])  col[++P->NHalo] = col[i];
    P->NHalo++; //number of unique columns

    P->Ncols += P->NHalo;

    //save global column ids in colMap
    P->colMap = (hlong *) realloc(P->colMap, P->Ncols*sizeof(hlong));
    for (dlong i=0; i<P->NHalo; i++)
      P->colMap[i+P->NlocalCols] = col[i];
    free(col);

    //shift the column indices to local indexing
    for (dlong i=0;i<P->offdNNZ;i++) {
      hlong gcol = offdCols[i];
      for (dlong m=P->NlocalCols;m<P->Ncols;m++) {
        if (gcol == P->colMap[m])
          P->offdCols[i] = m;
      }
    }
    free(offdCols);
  }

  csrHaloSetup(P,globalAggStarts);

  // normalize the columns of P
  *nullCoarseA = (dfloat *) calloc(P->Ncols,sizeof(dfloat));

  //add local nonzeros
  for(dlong i=0; i<P->diagNNZ; i++)
    (*nullCoarseA)[P->diagCols[i]] += P->diagCoefs[i] * P->diagCoefs[i];

  dfloat *nnzSum, *recvNnzSum;
  if (P->NHalo) nnzSum = (dfloat *) calloc(P->NHalo,sizeof(dfloat));
  if (P->NsendTotal) recvNnzSum = (dfloat *) calloc(P->NsendTotal,sizeof(dfloat));

  //add the non-local non-zeros
  for (dlong i=0;i<P->offdNNZ;i++)
    nnzSum[P->offdCols[i]-P->NlocalCols] += P->offdCoefs[i] * P->offdCoefs[i];

  //do a reverse halo exchange
  int tag = 999;

  // initiate immediate send  and receives to each other process as needed
  dlong recvOffset = 0;
  dlong sendOffset = 0;
  int sendMessage = 0, recvMessage = 0;
  for(int r=0;r<size;++r){
    if (P->NsendTotal) {
      if(P->NsendPairs[r]) {
        MPI_Irecv(recvNnzSum+sendOffset, P->NsendPairs[r], MPI_DFLOAT, r, tag,
            agmg::comm, (MPI_Request*)P->haloSendRequests+sendMessage);
        sendOffset += P->NsendPairs[r];
        ++sendMessage;
      }
    }
    if (P->NrecvTotal) {
      if(P->NrecvPairs[r]){
        MPI_Isend(nnzSum+recvOffset, P->NrecvPairs[r], MPI_DFLOAT, r, tag,
            agmg::comm, (MPI_Request*)P->haloRecvRequests+recvMessage);
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

  for(dlong i=0;i<P->NsendTotal;++i){
    // local index of outgoing element in halo exchange
    dlong id = P->haloElementList[i];

    (*nullCoarseA)[id] += recvNnzSum[i];
  }

  if (P->NHalo) free(nnzSum);

  for(dlong i=0; i<NCoarse; i++)
    (*nullCoarseA)[i] = sqrt((*nullCoarseA)[i]);

  csrHaloExchange(P, sizeof(dfloat), *nullCoarseA, P->sendBuffer, *nullCoarseA+P->NlocalCols);

  for(dlong i=0; i<P->diagNNZ; i++)
    P->diagCoefs[i] /= (*nullCoarseA)[P->diagCols[i]];
  for(dlong i=0; i<P->offdNNZ; i++)
    P->offdCoefs[i] /= (*nullCoarseA)[P->offdCols[i]];

  MPI_Barrier(agmg::comm);
  if (P->NsendTotal) free(recvNnzSum);

  return P;
}

typedef struct {

  hlong row;
  hlong col;
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
                hlong *globalRowStarts, hlong *globalColStarts){

  // MPI info
  int rank, size;
  rank = agmg::rank;
  size = agmg::size;

  csr *At = (csr *) calloc(1,sizeof(csr));

  At->Nrows = A->Ncols-A->NHalo;
  At->Ncols = A->Nrows;
  At->diagNNZ   = A->diagNNZ; //local entries remain local

  At->NlocalCols = At->Ncols;

  At->diagRowStarts = (dlong *)   calloc(At->Nrows+1, sizeof(dlong));
  At->offdRowStarts = (dlong *)   calloc(At->Nrows+1, sizeof(dlong));

  //start with local entries
  if (A->diagNNZ) {
    At->diagCols      = (dlong *)  calloc(At->diagNNZ, sizeof(dlong));
    At->diagCoefs     = (dfloat *) calloc(At->diagNNZ, sizeof(dfloat));
  }

  // count the num of nonzeros per row for transpose
  for(dlong i=0; i<A->diagNNZ; i++){
    dlong row = A->diagCols[i];
    At->diagRowStarts[row+1]++;
  }

  // cumulative sum for rows
  for(dlong i=1; i<=At->Nrows; i++)
    At->diagRowStarts[i] += At->diagRowStarts[i-1];

  int *counter = (int *) calloc(At->Nrows+1,sizeof(int));
  for (dlong i=0; i<At->Nrows+1; i++)
    counter[i] = At->diagRowStarts[i];

  for(dlong i=0; i<A->Nrows; i++){
    const dlong Jstart = A->diagRowStarts[i], Jend = A->diagRowStarts[i+1];

    for(dlong jj=Jstart; jj<Jend; jj++){
      dlong row = A->diagCols[jj];
      At->diagCols[counter[row]]  = i;
      At->diagCoefs[counter[row]] = A->diagCoefs[jj];

      counter[row]++;
    }
  }
  free(counter);

  //record global indexing of columns
  At->colMap = (hlong *)   calloc(At->Ncols, sizeof(hlong));
  for (dlong i=0;i<At->Ncols;i++)
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
  for(dlong i=0;i<A->Nrows;++i){
    for (dlong j=A->offdRowStarts[i];j<A->offdRowStarts[i+1];j++) {
      hlong col =  A->colMap[A->offdCols[j]]; //global ids
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

  MPI_Alltoall(Nsend, 1, MPI_INT, Nrecv, 1, MPI_INT, agmg::comm);

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
                      MPI_CHAR, r, tag, agmg::comm,
                      (MPI_Request*)A->haloSendRequests+recvMessage);
        recvOffset += Nrecv[r]*sizeof(nonzero_t);
        ++recvMessage;
      }
    }
    if (A->offdNNZ) {
      if(Nsend[r]){
        MPI_Isend(((char*)sendNonZeros)+sendOffset, Nsend[r]*sizeof(nonzero_t),
                      MPI_CHAR, r, tag, agmg::comm,
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

    hlong *offdCols  = (hlong *)   calloc(At->offdNNZ,sizeof(hlong));
    At->offdCols  = (dlong *)   calloc(At->offdNNZ,sizeof(dlong));
    At->offdCoefs = (dfloat *) calloc(At->offdNNZ, sizeof(dfloat));

    //find row starts
    for(dlong n=0;n<At->offdNNZ;++n) {
      dlong row = (dlong) (recvNonZeros[n].row - globalColStarts[rank]);
      At->offdRowStarts[row+1]++;
    }
    //cumulative sum
    for (dlong i=0;i<At->Nrows;i++)
      At->offdRowStarts[i+1] += At->offdRowStarts[i];

    //fill cols and coefs
    for (dlong i=0; i<At->Nrows; i++) {
      for (dlong j=At->offdRowStarts[i]; j<At->offdRowStarts[i+1]; j++) {
        offdCols[j]  = recvNonZeros[j].col;
        At->offdCoefs[j] = recvNonZeros[j].val;
      }
    }
    free(recvNonZeros);

    //we now need to reorder the x vector for the halo, and shift the column indices
    hlong *col = (hlong *) calloc(At->offdNNZ,sizeof(hlong));
    for (dlong n=0;n<At->offdNNZ;n++)
      col[n] = offdCols[n]; //copy non-local column global ids

    //sort by global index
    std::sort(col,col+At->offdNNZ);

    //count unique non-local column ids
    At->NHalo = 0;
    for (dlong n=1;n<At->offdNNZ;n++)
      if (col[n]!=col[n-1])  col[++At->NHalo] = col[n];
    At->NHalo++; //number of unique columns

    At->Ncols += At->NHalo;

    //save global column ids in colMap
    At->colMap = (hlong *) realloc(At->colMap,At->Ncols*sizeof(hlong));
    for (dlong n=0; n<At->NHalo; n++)
      At->colMap[n+At->NlocalCols] = col[n];
    free(col);

    //shift the column indices to local indexing
    for (dlong n=0;n<At->offdNNZ;n++) {
      hlong gcol = offdCols[n];
      for (dlong m=At->NlocalCols;m<At->Ncols;m++) {
        if (gcol == At->colMap[m])
          At->offdCols[n] = m;
      }
    }
    free(offdCols);
  }

  csrHaloSetup(At,globalRowStarts);

  return At;
}

typedef struct {

  hlong coarseId;
  dfloat coef;

} pEntry_t;

typedef struct {

  hlong I;
  hlong J;
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
  rank = agmg::rank;
  size = agmg::size;

  hlong *globalAggStarts = level->globalAggStarts;
  // hlong *globalRowStarts = level->globalRowStarts;

  hlong globalAggOffset = globalAggStarts[rank];

  //The galerkin product can be computed as
  // (RAP)_IJ = sum_{i in Agg_I} sum_{j in Agg_j} P_iI A_ij P_jJ
  // Since each row of P has only one entry, we can share the ncessary
  // P entries, form the products, and send them to their destination rank

  dlong N = A->Nrows;
  dlong M = A->Ncols;

  //printf("Level has %d rows, and is making %d aggregates\n", N, globalAggStarts[rank+1]-globalAggStarts[rank]);

  pEntry_t *PEntries;
  if (M) 
    PEntries = (pEntry_t *) calloc(M,sizeof(pEntry_t));
  else 
    PEntries = (pEntry_t *) calloc(1,sizeof(pEntry_t));

  //record the entries of P that this rank has
  dlong cnt =0;
  for (dlong i=0;i<N;i++) {
    for (dlong j=P->diagRowStarts[i];j<P->diagRowStarts[i+1];j++) {
      PEntries[cnt].coarseId = P->diagCols[j] + globalAggOffset; //global ID
      PEntries[cnt].coef = P->diagCoefs[j];
      cnt++;
    }
    for (dlong j=P->offdRowStarts[i];j<P->offdRowStarts[i+1];j++) {
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
  dlong totalNNZ = A->diagNNZ+A->offdNNZ;
  if (totalNNZ) 
    RAPEntries = (rapEntry_t *) calloc(totalNNZ,sizeof(rapEntry_t));
  else 
    RAPEntries = (rapEntry_t *) calloc(1,sizeof(rapEntry_t)); //MPI_AlltoAll doesnt like null pointers
  
  // Make the MPI_RAPENTRY_T data type
  MPI_Datatype MPI_RAPENTRY_T;
  MPI_Datatype dtype[3] = {MPI_HLONG, MPI_HLONG, MPI_DFLOAT};
  int blength[3] = {1, 1, 1};
  MPI_Aint addr[3], displ[3];
  MPI_Get_address ( &(RAPEntries[0]     ), addr+0);
  MPI_Get_address ( &(RAPEntries[0].J   ), addr+1);
  MPI_Get_address ( &(RAPEntries[0].coef), addr+2);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  MPI_Type_create_struct (3, blength, displ, dtype, &MPI_RAPENTRY_T);
  MPI_Type_commit (&MPI_RAPENTRY_T);

  //for the RAP products
  cnt =0;
  for (dlong i=0;i<N;i++) {
    for (dlong j=A->diagRowStarts[i];j<A->diagRowStarts[i+1];j++) {
      dlong col  = A->diagCols[j];
      dfloat coef = A->diagCoefs[j];

      RAPEntries[cnt].I = PEntries[i].coarseId;
      RAPEntries[cnt].J = PEntries[col].coarseId;
      RAPEntries[cnt].coef = coef*PEntries[i].coef*PEntries[col].coef;
      cnt++;
    }
  }
  for (dlong i=0;i<N;i++) {
    for (dlong j=A->offdRowStarts[i];j<A->offdRowStarts[i+1];j++) {
      dlong col  = A->offdCols[j];
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

  for(dlong i=0;i<totalNNZ;++i) {
    hlong id = RAPEntries[i].I;
    for (int r=0;r<size;r++) {
      if (globalAggStarts[r]-1<id && id < globalAggStarts[r+1])
        sendCounts[r]++;
    }
  }

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, agmg::comm);

  // find send and recv offsets for gather
  dlong recvNtotal = 0;
  for(int r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r];
  }
  rapEntry_t *recvRAPEntries;
  if (recvNtotal) 
    recvRAPEntries = (rapEntry_t *) calloc(recvNtotal,sizeof(rapEntry_t));
  else 
    recvRAPEntries = (rapEntry_t *) calloc(1,sizeof(rapEntry_t));//MPI_AlltoAll doesnt like null pointers
  
  MPI_Alltoallv(    RAPEntries, sendCounts, sendOffsets, MPI_RAPENTRY_T,
                recvRAPEntries, recvCounts, recvOffsets, MPI_RAPENTRY_T,
                agmg::comm);

  //sort entries by the coarse row and col
  if (recvNtotal) qsort(recvRAPEntries, recvNtotal, sizeof(rapEntry_t), compareRAPEntries);

  //count total number of nonzeros;
  dlong nnz =0;
  if (recvNtotal) nnz++;
  for (dlong i=1;i<recvNtotal;i++)
    if ((recvRAPEntries[i].I!=recvRAPEntries[i-1].I)||
          (recvRAPEntries[i].J!=recvRAPEntries[i-1].J)) nnz++;

  rapEntry_t *newRAPEntries;
  if (nnz)
    newRAPEntries = (rapEntry_t *) calloc(nnz,sizeof(rapEntry_t));
  else 
    newRAPEntries = (rapEntry_t *) calloc(1,sizeof(rapEntry_t));
  
  //compress nonzeros
  nnz = 0;
  if (recvNtotal) newRAPEntries[nnz++] = recvRAPEntries[0];
  for (dlong i=1;i<recvNtotal;i++) {
    if ((recvRAPEntries[i].I!=recvRAPEntries[i-1].I)||
          (recvRAPEntries[i].J!=recvRAPEntries[i-1].J)) {
      newRAPEntries[nnz++] = recvRAPEntries[i];
    } else {
      newRAPEntries[nnz-1].coef += recvRAPEntries[i].coef;
    }
  }

  dlong numAggs = (dlong) (globalAggStarts[rank+1]-globalAggStarts[rank]); //local number of aggregates

  csr *RAP = (csr*) calloc(1,sizeof(csr));

  RAP->Nrows = numAggs;
  RAP->Ncols = numAggs;

  RAP->NlocalCols = numAggs;

  RAP->diagRowStarts = (dlong *) calloc(numAggs+1, sizeof(dlong));
  RAP->offdRowStarts = (dlong *) calloc(numAggs+1, sizeof(dlong));

  for (dlong n=0;n<nnz;n++) {
    dlong row = (dlong) (newRAPEntries[n].I - globalAggOffset);
    if ((newRAPEntries[n].J > globalAggStarts[rank]-1)&&
          (newRAPEntries[n].J < globalAggStarts[rank+1])) {
      RAP->diagRowStarts[row+1]++;
    } else {
      RAP->offdRowStarts[row+1]++;
    }
  }

  // cumulative sum
  for(dlong i=0; i<numAggs; i++) {
    RAP->diagRowStarts[i+1] += RAP->diagRowStarts[i];
    RAP->offdRowStarts[i+1] += RAP->offdRowStarts[i];
  }
  RAP->diagNNZ = RAP->diagRowStarts[numAggs];
  RAP->offdNNZ = RAP->offdRowStarts[numAggs];

  dlong *diagCols;
  dfloat *diagCoefs;
  if (RAP->diagNNZ) {
    RAP->diagCols  = (dlong *)   calloc(RAP->diagNNZ, sizeof(dlong));
    RAP->diagCoefs = (dfloat *) calloc(RAP->diagNNZ, sizeof(dfloat));
    diagCols  = (dlong *)   calloc(RAP->diagNNZ, sizeof(dlong));
    diagCoefs = (dfloat *) calloc(RAP->diagNNZ, sizeof(dfloat));
  }
  hlong *offdCols;
  if (RAP->offdNNZ) {
    offdCols  = (hlong *)   calloc(RAP->offdNNZ,sizeof(hlong));
    RAP->offdCols  = (dlong *)   calloc(RAP->offdNNZ,sizeof(dlong));
    RAP->offdCoefs = (dfloat *) calloc(RAP->offdNNZ, sizeof(dfloat));
  }

  dlong diagCnt =0;
  dlong offdCnt =0;
  for (dlong n=0;n<nnz;n++) {
    if ((newRAPEntries[n].J > globalAggStarts[rank]-1)&&
          (newRAPEntries[n].J < globalAggStarts[rank+1])) {
      diagCols[diagCnt]  = (dlong) (newRAPEntries[n].J - globalAggOffset);
      diagCoefs[diagCnt] = newRAPEntries[n].coef;
      diagCnt++;
    } else {
      offdCols[offdCnt]  = newRAPEntries[n].J;
      RAP->offdCoefs[offdCnt] = newRAPEntries[n].coef;
      offdCnt++;
    }
  }

  //move diagonal entries first
  for (dlong i=0;i<RAP->Nrows;i++) {
    dlong start = RAP->diagRowStarts[i];
    int cnt = 1;
    for (dlong j=RAP->diagRowStarts[i]; j<RAP->diagRowStarts[i+1]; j++) {
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
  RAP->colMap = (hlong *)   calloc(RAP->Ncols, sizeof(hlong));
  for (dlong i=0;i<RAP->Ncols;i++)
    RAP->colMap[i] = i + globalAggOffset;

  if (RAP->offdNNZ) {
    //we now need to reorder the x vector for the halo, and shift the column indices
    hlong *col = (hlong *) calloc(RAP->offdNNZ,sizeof(hlong));
    for (dlong n=0;n<RAP->offdNNZ;n++)
      col[n] = offdCols[n]; //copy non-local column global ids

    //sort by global index
    std::sort(col,col+RAP->offdNNZ);

    //count unique non-local column ids
    RAP->NHalo = 0;
    for (dlong n=1;n<RAP->offdNNZ;n++)
      if (col[n]!=col[n-1])  col[++RAP->NHalo] = col[n];
    RAP->NHalo++; //number of unique columns

    RAP->Ncols += RAP->NHalo;

    //save global column ids in colMap
    RAP->colMap = (hlong *) realloc(RAP->colMap,RAP->Ncols*sizeof(hlong));
    for (dlong n=0; n<RAP->NHalo; n++)
      RAP->colMap[n+RAP->NlocalCols] = col[n];

    //shift the column indices to local indexing
    for (dlong n=0;n<RAP->offdNNZ;n++) {
      hlong gcol = offdCols[n];
      for (dlong m=RAP->NlocalCols;m<RAP->Ncols;m++) {
        if (gcol == RAP->colMap[m])
          RAP->offdCols[n] = m;
      }
    }
    free(col);
    free(offdCols);
  }
  csrHaloSetup(RAP,globalAggStarts);

  //clean up
  MPI_Barrier(agmg::comm);
  MPI_Type_free(&MPI_RAPENTRY_T);

  free(PEntries);
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


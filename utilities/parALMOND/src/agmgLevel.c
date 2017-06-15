#include "parAlmond.h"

//TODO do we really need these functions? just replace their calls
void restrict(agmgLevel *level, dfloat *r, dfloat *Rr){
  axpy(level->R, 1.0, r, 0.0, Rr);
}

void restrict(parAlmond_t *parAlmond, agmgLevel *level, occa::memory o_r, occa::memory o_Rr){
  axpy(parAlmond, level->deviceR, 1.0, o_r, 0.0, o_Rr);
}

void interpolate(agmgLevel *level, dfloat *x, dfloat *Px){
  axpy(level->P, 1.0, x, 1.0, Px);
}

void interpolate(parAlmond_t *parAlmond, agmgLevel *level, occa::memory o_x, occa::memory o_Px){
  axpy(parAlmond, level->dcsrP, 1.0, o_x, 0.0, o_Px);
}

void allocate(agmgLevel *level){
  if(level->A->Nrows){
    iint m = level->A->Nrows;
    iint n = level->A->Ncols;

    level->x    = (dfloat *) calloc(n, sizeof(dfloat));
    level->rhs  = (dfloat *) calloc(m, sizeof(dfloat));
    level->res  = (dfloat *) calloc(n, sizeof(dfloat));
  }
}

csr* distribute(csr *A, iint *globalRowStarts, iint *globalColStarts) {
  iint rank,size;
  MPI_Comm_size(MPI_COMM_WORLD, &size );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  csr *localA = (csr *) calloc(1,sizeof(csr));

  iint rowStart = globalRowStarts[rank];
  iint rowEnd   = globalRowStarts[rank+1];
  iint colStart = globalColStarts[rank];
  iint colEnd   = globalColStarts[rank+1];

  localA->Nrows = rowEnd-rowStart;
  localA->numLocalIds = colEnd-colStart;
  localA->nnz   = A->rowStarts[rowEnd]-A->rowStarts[rowStart];
  localA->rowStarts = (iint *) calloc(rowEnd-rowStart+1,sizeof(iint));
  localA->cols      = (iint *) calloc(localA->nnz,sizeof(iint));
  localA->coefs     = (dfloat *) calloc(localA->nnz,sizeof(dfloat));

  //copy the slice of the global matrix
  for (iint n=0;n<localA->Nrows+1;n++)
    localA->rowStarts[n] = A->rowStarts[rowStart+n]-A->rowStarts[rowStart];

  for (iint n=0;n<localA->nnz;n++) {
    localA->cols[n]  = A->cols[A->rowStarts[rowStart]+n];
    localA->coefs[n] = A->coefs[A->rowStarts[rowStart]+n];
  }

  //we now need to reorder the x vector for the halo, and shift the column indices
  if (localA->nnz) {
    iint *col = (iint *) calloc(localA->nnz,sizeof(iint)); //list of global column ids for every nonzero
    for (iint n=0;n<localA->nnz;n++)
      col[n] = localA->cols[n]; //grab global ids
    
    //sort by global index
    std::sort(col,col+localA->nnz);

    //compress
    iint numcol = 0; 
    for (iint n=1;n<localA->nnz;n++)
      if (col[n]!=col[numcol]) 
        col[++numcol] = col[n];
    numcol++; //number of unique columns

    iint NHalo = 0;
    for (iint n=0;n<numcol;n++) 
      if ( (col[n] < colStart) || (col[n] > colEnd-1)) NHalo++; //number of halo columns

    localA->NHalo = NHalo;
    localA->Ncols = localA->numLocalIds + NHalo;
    
    localA->colMap = (iint*) calloc(localA->Ncols,sizeof(iint));
    
    for (iint n=0; n < localA->numLocalIds;n++)
      localA->colMap[n] = n + colStart;

    iint cnt = localA->numLocalIds;
    for (iint n=0; n<numcol;n++) 
      if ((col[n]<colStart) || (col[n]>colEnd-1)) 
        localA->colMap[cnt++] = col[n];
    
    free(col);

    //shift the column indices to local indexing
    for (iint n=0;n<localA->nnz;n++) {
      iint gcol = localA->cols[n];
      if ((gcol > colStart-1) && (gcol < colEnd)) {//local column
        localA->cols[n] -= colStart;
      } else {
        for (iint m = localA->numLocalIds;m<localA->Ncols;m++)
          if (gcol == localA->colMap[m])
            localA->cols[n] = m;
      }
    }
  }

  csrHaloSetup(localA, globalColStarts);

  return localA;
}



void setup_smoother(agmgLevel *level, SmoothType s){

  level->stype = s;

  if(s == JACOBI){
    return;
  }

  if(s == DAMPED_JACOBI){
    // estimate rho(invD * A)
    dfloat rho=0;

    dfloat *invD;
    if(level->A->Nrows)	{
      invD = (dfloat *) calloc(level->A->Nrows, sizeof(dfloat));
      for (iint i=0;i<level->A->Nrows;i++)
        invD[i] = 1.0/level->A->coefs[level->A->rowStarts[i]];

      rho = rhoDinvA(level->A, invD);

      free(invD);
    }

    level->smoother_params = (dfloat *) calloc(1,sizeof(dfloat));

    level->smoother_params[0] = (4./3.)/rho;
    return;
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

dfloat rhoDinvA(csr *A, dfloat *invD){

  const iint m = A->Nrows;
  const iint n = A->Ncols;

  int k = 10; 

  if(k > m)
    k = m;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double *H = new double [k*k];
  for(int i=0; i<k*k; i++)
    H[i] = 0.;

  // allocate memory for basis
  dfloat **V = (dfloat **) calloc(k+1, sizeof(dfloat *));
  dfloat *Vx = (dfloat *) calloc(n, sizeof(dfloat));

  for(int i=0; i<=k; i++)
    V[i] = (dfloat *) calloc(m, sizeof(dfloat));

  // generate a random vector for initial basis vector
  randomize(m, Vx);

  dfloat norm_vo = norm(m, Vx);
  scaleVector(m, Vx, 1./norm_vo);

  for (iint i=0;i<m;i++)
    V[0][i] = Vx[i];

  for(int j=0; j<k; j++){

    for (iint i=0;i<m;i++)
      Vx[i] = V[j][i];

    // v[j+1] = invD*(A*v[j])
    axpy(A, 1.0, Vx, 0., V[j+1]);

    dotStar(m, invD, V[j+1]);

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij = innerProd(m, V[i], V[j+1]);

      // v[j+1] = v[j+1] - hij*v[i]
      vectorAdd(m,-hij, V[i], 1.0, V[j+1]);

      H[i + j*k] = (double) hij;
    }

    if(j+1 < k){
      H[j+1+ j*k] = (double) norm(m,V[j+1]);

      scaleVector(m,V[j+1], 1./H[j+1 + j*k]);
    }
  }

  double *WR = new double[k];
  double *WI = new double[k];

  eig(k, H, WR, WI);

  double rho = 0.;

  for(int i=0; i<k; i++){
    double rho_i  = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);

    if(rho < rho_i) {
      rho = rho_i;
    }
  }

  delete [] H;
  delete [] WR;
  delete [] WI;

  // free memory
  for(int i=0; i<=k; i++){
    free(V[i]);
  }

  return rho;
}


void smooth(agmgLevel *level, dfloat *rhs, dfloat *x, bool x_is_zero){
  if(level->stype == JACOBI){
    smoothJacobi(level->A, rhs, x, x_is_zero);
    return;
  }

  if(level->stype == DAMPED_JACOBI){
    smoothDampedJacobi(level->A, rhs, x, level->smoother_params[0], x_is_zero);
    return;
  }
}


void smooth(parAlmond_t *parAlmond, agmgLevel *level, occa::memory o_rhs, occa::memory o_x, bool x_is_zero){

  if(level->stype == JACOBI){
    smoothJacobi(parAlmond, level->deviceA, o_rhs, o_x, x_is_zero);
    return;
  }

  if(level->stype == DAMPED_JACOBI){
    smoothDampedJacobi(parAlmond, level->deviceA, o_rhs, o_x, level->smoother_params[0], x_is_zero);
    return;
  }
}

void matFreeSmooth(parAlmond_t *parAlmond, agmgLevel *level, occa::memory &o_r, occa::memory &o_x, bool x_is_zero) {
  if(level->stype == JACOBI){
    matFreeSmoothJacobi(parAlmond, level->deviceA, o_r, o_x, x_is_zero);
    return;
  }

  if(level->stype == DAMPED_JACOBI){
    matFreeSmoothDampedJacobi(parAlmond, level->deviceA, o_r, o_x, level->smoother_params[0], x_is_zero);
    return;
  }
}

csr * strong_graph(csr *A, dfloat threshold){

  const iint m = A->Nrows;

  csr *C = (csr *) calloc(1, sizeof(csr));

  C->Nrows = m;
  C->Ncols = m;

  C->rowStarts = (iint *) calloc(m+1,sizeof(iint));

  for(iint i=0; i<m; i++){
    dfloat sign = (A->coefs[A->rowStarts[i]] >= 0) ? 1:-1;

    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    dfloat maxOD = 0.;

    dfloat Aii = fabs(A->coefs[A->rowStarts[i]]);
    iint jj = Jstart+1;
    for(; jj<Jend; jj++){
      iint col = A->cols[jj];
      dfloat Ajj = fabs(A->coefs[A->rowStarts[col]]);
      dfloat OD = -sign*A->coefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > maxOD) maxOD = OD;
    }

    jj = Jstart+1;
    iint strong_per_row = 1; // diagonal entry
    for(; jj<Jend; jj++){
      iint col = A->cols[jj];
      dfloat Ajj = fabs(A->coefs[A->rowStarts[col]]);
      dfloat OD = -sign*A->coefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD) strong_per_row++;
    }

    C->rowStarts[i+1] = strong_per_row;
  }

  // cumulative sum
  for(iint i=1; i<=m ; i++)
    C->rowStarts[i] += C->rowStarts[i-1];

  C->nnz = C->rowStarts[m];

  C->cols = (iint *) calloc(C->nnz, sizeof(iint));

  // fill in the columns for strong connections
  for(iint i=0; i<m; i++){

    dfloat sign = (A->coefs[A->rowStarts[i]] >= 0) ? 1:-1;

    dfloat Aii = fabs(A->coefs[A->rowStarts[i]]);

    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    dfloat maxOD = 0.;

    iint jj = Jstart+1;
    for(; jj<Jend; jj++){
      iint col = A->cols[jj];
      dfloat Ajj = fabs(A->coefs[A->rowStarts[col]]);
      dfloat OD = -sign*A->coefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > maxOD) maxOD = OD;
    }

    iint counter = C->rowStarts[i];

    // diag entry
    C->cols[counter++] = i;

    jj = Jstart+1;
    for(; jj<Jend; jj++){
      iint col = A->cols[jj];
      dfloat Ajj = fabs(A->coefs[A->rowStarts[col]]);
      dfloat OD = -sign*A->coefs[jj]/(sqrt(Aii)*sqrt(Ajj));
      if(OD > threshold*maxOD)
        C->cols[counter++] = A->cols[jj];
    }
  }

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

  const iint m   = C->Nrows;
  const iint nnz = C->nnz;

  iint *FineToCoarse = (iint *) calloc(m, sizeof(iint));
  for (iint i =0;i<m;i++) FineToCoarse[i] = -1;

  dfloat *rands  = (dfloat *) calloc(m, sizeof(dfloat)); 
  iint   *states = (iint *)   calloc(m, sizeof(iint));

  dfloat *Tr = (dfloat *) calloc(m, sizeof(dfloat));
  iint   *Ts = (iint *)   calloc(m, sizeof(iint)); 
  iint   *Ti = (iint *)   calloc(m, sizeof(iint)); 

  for(iint i=0; i<m; i++){
    rands[i] = (dfloat) drand48();
    states[i] = 0;
  }

  // add the number of strong connections
  for(iint i=0; i<nnz; i++)
    rands[C->cols[i]] += 1.;

  //-->Halo exchange for rands, fill halo zone of states with 0

  bool done = false;
  while(!done){
    // first neighbours
    for(iint i=0; i<m; i++){

      iint smax = states[i];
      dfloat rmax = rands[i];
      iint imax = i; //-->globalIndex

      if(smax != 1){
        iint jj=C->rowStarts[i], Jend=C->rowStarts[i+1];

        jj++;
        for(;jj<Jend;jj++){
          const iint col = C->cols[jj];
          if(customLess(smax, rmax, imax, states[col], rands[col], col)){//col-->globalIndex[col];
            smax = states[col];
            rmax = rands[col];
            imax = col; //-->globalIndex[col];
          }
        }
      }
      Ts[i] = smax;
      Tr[i] = rmax;
      Ti[i] = imax;
    }

    //-->Halo exchange for Ts,Tr,Ti

    // second neighbours
    for(iint i=0; i<m; i++){
      iint   smax = Ts[i];
      dfloat rmax = Tr[i];
      iint   imax = Ti[i];

      iint jj=C->rowStarts[i], Jend=C->rowStarts[i+1];

      jj++;
      for(;jj<Jend;jj++){
        const iint col = C->cols[jj];
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }

      // if I am the strongest among all the 1 and 2 ring neighbours
      // I am an MIS node
      if(states[i] == 0 && imax == i) //-->globalIndex[i]
        states[i] = 1;

      // if there is an MIS node within distance 2, I am removed
      if(states[i] == 0 && smax == 1)
        states[i] = -1;
    }

    //-->halo exchange for states

    // if number of undecided nodes = 0, algorithm terminates  
    done = (std::count(states, states+m, 0) == 0);
    //-->MPI All_reduce 'done'
  }

  level->numAggregates = 0;
  // enumerate the coarse nodes/aggregates
  for(iint i=0; i<m; i++){
    if(states[i] == 1)
      FineToCoarse[i] = level->numAggregates++;
  }
  //--> MPI All_reduce numAggregates, and make globalIndexing

  //-->Halo exchange for FineToCoarse

  // form the aggregates
  for(iint i=0; i<m; i++){
    iint   smax = states[i];
    dfloat rmax = rands[i];
    iint   imax = i; //-->globalIndex[i]

    if(smax != 1){
      iint jj=C->rowStarts[i], Jend=C->rowStarts[i+1];

      jj++;
      for(;jj<Jend;jj++){
        const iint col = C->cols[jj];
        if(customLess(smax, rmax, imax, states[col], rands[col], col)){ //-->globalIndex
          smax = states[col];
          rmax = rands[col];
          imax = col; //-->globalIndex
        }
      }
    }
    Ts[i] = smax;
    Tr[i] = rmax;
    Ti[i] = imax;

    //-->col = local index where globalIndex[col] = imax

    if(states[i] == -1 && smax == 1 && FineToCoarse[imax] > -1) //->FineToCoarse[col] 
      FineToCoarse[i] = FineToCoarse[imax];//->FineToCoarse[col]
  }

  //-->Halo exchange for FineToCoarse
  //-->Halo exchange for Ts,Tr,Ti

  // second neighbours
  for(iint i=0; i<m; i++){
    iint smax   = Ts[i];
    dfloat rmax = Tr[i];
    iint imax   = Ti[i];

    iint jj=C->rowStarts[i], Jend=C->rowStarts[i+1];

    jj++;
    for(;jj<Jend;jj++){
      const iint col = C->cols[jj];
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){ //-->globalIndex
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
      }
    }

    //-->col = localIndex where globalIndex[col] = imax

    if(states[i] == -1 && smax == 1 && FineToCoarse[imax] > -1)//->FineToCoarse[col] 
      FineToCoarse[i] = FineToCoarse[imax];//->FineToCoarse[col] 
  }

  free(rands); 
  free(states);
  free(Tr);
  free(Ts);
  free(Ti);

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
  iint numOwnedAggs = 0;
  if (newRecvNtotal) numOwnedAggs++;
  for (iint i=1;i<newRecvNtotal;i++) 
    if(newRecvAggs[i].coarseId!=newRecvAggs[i-1].coarseId) numOwnedAggs++; 
  
  //determine a global numbering of the aggregates
  iint *globalNumAggs = (iint *) calloc(size,sizeof(iint));
  MPI_Allgather(&numOwnedAggs, 1, MPI_IINT, globalNumAggs, 1, MPI_IINT, MPI_COMM_WORLD);
  
  iint *NumAggsOffset = (iint *) calloc(size+1,sizeof(iint));
  for (iint r=0;r<size;r++)
    NumAggsOffset[r+1] += globalNumAggs[r];
  
  //set the new global coarse index
  iint cnt = NumAggsOffset[rank];
  if (newRecvNtotal) newRecvAggs[0].newCoarseId = cnt;
  for (iint i=1;i<newRecvNtotal;i++) {
    if(newRecvAggs[i].coarseId!=newRecvAggs[i-1].coarseId) cnt++;

    newRecvAggs[i].newCoarseId = cnt; 
  }
  free(globalNumAggs); free(NumAggsOffset);


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

  //TODO probably need to record the owners as well for the halo setup
}


void construct_interpolator(agmgLevel *level, iint *FineToCoarse, dfloat **nullCoarseA){

  const iint m = level->A->Nrows;
  const iint n = level->numAggregates; //local num agg

  level->P = (csr *) calloc(1, sizeof(csr));

  level->P->Nrows = m;
  level->P->Ncols = n;
  level->P->nnz = m;

  level->P->rowStarts = (iint *) calloc(m+1, sizeof(iint));

  // each row has exactly one nonzero per row
  for(iint i=0; i<=m; i++)
    level->P->rowStarts[i] = i;

  level->P->cols = (iint *) calloc(m, sizeof(iint));
  level->P->coefs = (dfloat *) calloc(m, sizeof(dfloat));
  for (iint i=0; i<m; i++) {
    level->P->cols[i] = FineToCoarse[i];
    level->P->coefs[i] = level->nullA[i];
  }

  // normalize the columns of P
  *nullCoarseA = (dfloat *) calloc(n,sizeof(dfloat));
   
  for(iint i=0; i<m; i++)
    (*nullCoarseA)[level->P->cols[i]] += level->P->coefs[i] * level->P->coefs[i];

  for(iint i=0; i<n; i++)
    (*nullCoarseA)[i] = sqrt((*nullCoarseA)[i]);

  for(iint i=0; i<m; i++){
    level->P->coefs[i] /= (*nullCoarseA)[FineToCoarse[i]];
  }

  level->R = transpose(level->P);
}

csr *galerkinProd(agmgLevel *level){

  iint numAgg = level->numAggregates;

  if(numAgg == 0){
    for (iint i=0; i<level->P->nnz; i++) //find max column index
      numAgg = (numAgg<level->P->cols[i]) ? level->P->cols[i] : numAgg;
    numAgg++;
  }
  

  csr *RAP = (csr*) calloc(1,sizeof(csr));

  RAP->Nrows = numAgg;
  RAP->Ncols = numAgg;

  RAP->rowStarts = (iint *) calloc(numAgg+1, sizeof(iint));

  dfloat *dummyCoefs = (dfloat*) calloc(level->A->nnz,sizeof(dfloat));
  for (iint i=0; i<level->A->nnz;i++) //copy A coefs
    dummyCoefs[i] = level->A->coefs[i];

  long base = numAgg+1;
  struct key_value_pair1 *pair = (struct key_value_pair1 *) calloc(level->A->nnz,sizeof(struct key_value_pair1));

  for(iint i=0; i<level->A->Nrows; i++){
    for(iint jj=level->A->rowStarts[i]; jj<level->A->rowStarts[i+1]; jj++){
      iint j = level->A->cols[jj];

      dummyCoefs[jj] *= (level->P->coefs[i] * level->P->coefs[j]);

      iint I = level->P->cols[i];
      iint J = level->P->cols[j];

      pair[jj].value = jj;

      if (I == J)
        pair[jj].key = ((long)I) * base;
      else
        pair[jj].key = ((long)I) * base + ((long) (J+1));

    }
  }

  qsort(pair, level->A->nnz, sizeof(struct key_value_pair1), compare_key1);

  iint nnz = 1;
  for(iint i=1; i<level->A->nnz; i++)
    if(pair[i].key > pair[i-1].key) nnz++;

  RAP->nnz = nnz;

  RAP->cols  = (iint *) calloc(nnz, sizeof(iint));
  RAP->coefs = (dfloat *) calloc(nnz, sizeof(dfloat));

  iint count = 0;
  RAP->coefs[count] = dummyCoefs[pair[0].value];
  RAP->cols[count] = 0; // This assumes that first entry is diagonal
  RAP->rowStarts[1]++;
  for(iint i=1; i<level->A->nnz; i++){
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
  for(iint i=1; i<=numAgg; i++)
    RAP->rowStarts[i] += RAP->rowStarts[i-1];

  //MPI comms
  RAP->NrecvTotal =0;
  RAP->NsendTotal =0;

  return RAP;
}

csr * transpose(csr *A){

  csr *At = (csr *) calloc(1,sizeof(csr));

  At->Nrows = A->Ncols;
  At->Ncols = A->Nrows;
  At->nnz   = A->nnz;

  At->rowStarts = (iint *) calloc(At->Nrows+1, sizeof(iint));
  At->cols      = (iint *) calloc(At->nnz, sizeof(iint));
  At->coefs     = (dfloat *) calloc(At->nnz, sizeof(dfloat));

  // count the no of nonzeros per row for transpose
  for(iint i=0; i<A->nnz; i++){
    iint row = A->cols[i];
    At->rowStarts[row+1]++;
  }

  // cumulative sum for rows
  for(iint i=1; i<=At->Nrows; i++)
    At->rowStarts[i] += At->rowStarts[i-1];

  iint counter[At->Nrows+1];
  for (iint i=0; i<At->Nrows+1; i++)
    counter[i] = At->rowStarts[i];

  if(A->Nrows != A->Ncols){
    for(iint i=0; i<A->Nrows; i++){
      const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];
      
      for(iint jj=Jstart; jj<Jend; jj++){
        iint row = A->cols[jj];
        At->cols[counter[row]] = i;
        At->coefs[counter[row]] = A->coefs[jj];

        counter[row]++;
      }
    }
  } else{
    // fill in diagonals first
    for(iint i=0; i<A->Nrows; i++){
      At->cols[counter[i]] = i;

      At->coefs[counter[i]] = A->coefs[A->rowStarts[i]];
      counter[i]++;
    }

    // fill the remaining ones
    for(iint i=0; i<A->Nrows; i++){
      const iint Jstart = A->rowStarts[i]+1, Jend = A->rowStarts[i+1];

      for(iint jj=Jstart; jj<Jend; jj++){
        iint row = A->cols[jj];

        At->cols[counter[row]] = i;
        At->coefs[counter[row]] = A->coefs[jj];

        counter[row]++;
      }
    }
  }

  return At;
}

void coarsen(agmgLevel *level, csr **coarseA, dfloat **nullCoarseA){

  const iint m = level->A->Nrows;
  const iint n = level->A->Ncols;

  if(m != n){
    printf("Matrix must be square..exiting \n");
    exit(-1);
  }

  // establish the graph of strong connections
  level->threshold = 0.5;
  
  csr *C = strong_graph(level->A, level->threshold);

  iint *FineToCoarse = form_aggregates(level, C);

  construct_interpolator(level, FineToCoarse, nullCoarseA);

  *coarseA = galerkinProd(level);
}

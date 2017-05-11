#include "parAlmond.h"

//TODO do we really need these functions? just replace their calls
void restrict(agmgLevel *level, dfloat *r, dfloat *Rr){
  axpy(level->R, 1.0, r, 0.0, Rr);
}

void restrict(almond_t *almond, agmgLevel *level, occa::memory o_r, occa::memory o_Rr){
  axpy(almond, level->deviceR, 1.0, o_r, 0.0, o_Rr);
}

void interpolate(agmgLevel *level, dfloat *x, dfloat *Px){
  axpy(level->P, 1.0, x, 1.0, Px);
}

void interpolate(almond_t *almond, agmgLevel *level, occa::memory o_x, occa::memory o_Px){
  axpy(almond, level->dcsrP, 1.0, o_x, 0.0, o_Px);
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


void smooth(almond_t *almond, agmgLevel *level, occa::memory o_rhs, occa::memory o_x, bool x_is_zero){

  if(level->stype == JACOBI){
    smoothJacobi(almond, level->deviceA, o_rhs, o_x, x_is_zero);
    return;
  }

  if(level->stype == DAMPED_JACOBI){
    smoothDampedJacobi(almond, level->deviceA, o_rhs, o_x, level->smoother_params[0], x_is_zero);
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

    iint jj = Jstart+1;
    for(; jj<Jend; jj++){
      dfloat OD = -sign*A->coefs[jj];
      if(OD > maxOD) maxOD = OD;
    }

    jj = Jstart+1;
    iint strong_per_row = 1; // diagonal entry
    for(; jj<Jend; jj++){
      dfloat OD = -sign*A->coefs[jj];
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

    const iint Jstart = A->rowStarts[i], Jend = A->rowStarts[i+1];

    dfloat maxOD = 0.;

    iint jj = Jstart+1;
    for(; jj<Jend; jj++){
      dfloat OD = -sign*A->coefs[jj];
      if(OD > maxOD) maxOD = OD;
    }

    iint counter = C->rowStarts[i];

    // diag entry
    C->cols[counter++] = i;

    jj = Jstart+1;
    for(; jj<Jend; jj++){
      dfloat OD = -sign*A->coefs[jj];
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

  dfloat *Tr     = (dfloat *) calloc(m, sizeof(dfloat)); 
  dfloat *rands  = (dfloat *) calloc(m, sizeof(dfloat)); 
  dfloat *Tr_hat = (dfloat *) calloc(m, sizeof(dfloat));
  iint *Ts     = (iint *) calloc(m, sizeof(iint)); 
  iint *Ts_hat = (iint *) calloc(m, sizeof(iint)); 
  iint *Ti     = (iint *) calloc(m, sizeof(iint)); 
  iint *Ti_hat = (iint *) calloc(m, sizeof(iint)); 
  iint *states = (iint *) calloc(m, sizeof(iint));

  for(iint i=0; i<m; i++){
    rands[i] = (dfloat) drand48();
    Ti[i] = Ti_hat[i] = i;
    states[i] = 0;
  }

  // count the number of
  for(iint i=0; i<nnz; i++)
    rands[C->cols[i]] += 1.;

  for(iint i=0; i<m; i++) {
    Tr[i] = rands[i];
    Tr_hat[i] = rands[i];
    Ts_hat[i] = 0;
  }

  bool done = false;

  while(!done){
    // first neighbours
    for(iint i=0; i<m; i++){

      iint smax = states[i];

      dfloat rmax = rands[i];
      iint imax = i;

      if(smax != 1){
        iint jj=C->rowStarts[i], Jend=C->rowStarts[i+1];

        jj++;

        for(;jj<Jend;jj++){
          const iint col = C->cols[jj];
          
          // TODO : sixth argument needs attention for parallel implementation
          if(customLess(smax, rmax, imax, states[col], rands[col], col)){
            smax = states[col];
            rmax = rands[col];
            imax = col;
          }
        }
      }

      Ts_hat[i] = smax;
      Tr_hat[i] = rmax;
      Ti_hat[i] = imax;
    }


    // second neighbours
    for(iint i=0; i<m; i++){

      iint smax = Ts_hat[i];
      dfloat rmax = Tr_hat[i];
      iint imax = Ti_hat[i];

      iint jj=C->rowStarts[i], Jend=C->rowStarts[i+1];

      jj++;

      for(;jj<Jend;jj++){
        const iint col = C->cols[jj];
        if(customLess(smax, rmax, imax, Ts_hat[col], Tr_hat[col], Ti_hat[col])){
          smax = Ts_hat[col];
          rmax = Tr_hat[col];
          imax = Ti_hat[col];
        }
      }

      // if I am the strongest among all the 1 and 2 ring neighbours
      // I am an MIS node
      if(states[i] == 0 && imax == i)
        states[i] = 1;

      // if there is an MIS node within distance 2, I am removed
      if(states[i] == 0 && smax == 1)
        states[i] = -1;
    }

    // if number of undecided nodes = 0, algorithm terminates  
    done = (std::count(states, states+m, 0) == 0);
  }

  level->numAggregates = 0;
  // enumerate the coarse nodes/aggregates
  for(iint i=0; i<m; i++){
    if(states[i] == 1)
      FineToCoarse[i] = level->numAggregates++;
  }

  // TODO: cumulative scan for MPI to get global enumeration of aggregations

  // form the aggregates
  for(iint i=0; i<m; i++){
    iint smax = states[i];
    dfloat rmax = rands[i];
    iint imax = i;

    if(smax != 1){
      iint jj=C->rowStarts[i], Jend=C->rowStarts[i+1];

      jj++;

      for(;jj<Jend;jj++){
        const iint col = C->cols[jj];
        // TODO : sixth argument needs attention for parallel implementation
        if(customLess(smax, rmax, imax, states[col], rands[col], col)){
          smax = states[col];
          rmax = rands[col];
          imax = col;
        }
      }
    }

    Ts_hat[i] = smax;
    Tr_hat[i] = rmax;
    Ti_hat[i] = imax;

    if(states[i] == -1 && smax == 1 && FineToCoarse[imax] > -1)
      FineToCoarse[i] = FineToCoarse[imax];
  }


  // second neighbours
  for(iint i=0; i<m; i++){
    iint smax = Ts_hat[i];
    dfloat rmax = Tr_hat[i];
    iint imax = Ti_hat[i];

    iint jj=C->rowStarts[i], Jend=C->rowStarts[i+1];

    jj++;

    for(;jj<Jend;jj++){
      const iint col = C->cols[jj];
      if(customLess(smax, rmax, imax, Ts_hat[col], Tr_hat[col], Ti_hat[col])){
        smax = Ts_hat[col];
        rmax = Tr_hat[col];
        imax = Ti_hat[col];
      }
    }

    if(states[i] == -1 && smax == 1 && FineToCoarse[imax] > -1)
      FineToCoarse[i] = FineToCoarse[imax];
  }

  free(Tr    ); 
  free(rands ); 
  free(Tr_hat);
  free(Ts    );
  free(Ts_hat);
  free(Ti    );
  free(Ti_hat);
  free(states);

  return FineToCoarse;
}


void construct_interpolator(agmgLevel *level, iint *FineToCoarse, dfloat **nullCoarseA){

  const iint m = level->A->Nrows;
  const iint n = level->numAggregates;

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

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
        invD[i] = 1.0/level->A->diagCoefs[level->A->diagRowStarts[i]];

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





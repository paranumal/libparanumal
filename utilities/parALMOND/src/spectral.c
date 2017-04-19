#include "parAlmond.h"

extern "C"{
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
	double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
}


static void eig(const int Nrows, double *A, double *WR,
	  double *WI){

  int NB	= 256;
  char JOBVL	= 'V';
  char JOBVR	= 'V';
  int     N	= Nrows;
  int   LDA	= Nrows;
  int  LWORK	= (NB+2)*N;

  double *WORK	= new double[LWORK];
  double *VL	= new double[Nrows*Nrows];
  double *VR	= new double[Nrows*Nrows];

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
    VL, &LDA, VR, &LDA, WORK, &LWORK, &INFO);


  assert(INFO == 0);

  delete [] VL;
  delete [] VR;
  delete [] WORK;
}

dfloat rhoDinvA(csr *A, dfloat *invD, int l){

  const iint m = A->Nrows;
  const iint n = A->Ncols;

  int k = l;

  if(k > m)
    k = m;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double *H = new double [k*k];
  for(int i=0; i<k*k; i++)
    H[i] = 0.;

  // allocate memory for basis
  dfloat **V = (dfloat **) calloc(k+1, sizeof(dfloat *))
  dfloat *Vx = (dfloat *) calloc(n, sizeof(dfloat));

  for(int i=0; i<=k; i++)
    V[i] = (dfloat *) calloc(m, sizeof(dfloat));

  // generate a random vector for initial basis vector
  randomize(Vx);

  dfloat norm_vo = norm(Vx);
  scaleVector(Vx, 1./norm_vo);

  for (iint i=0;i<m;i++)
    V[0][i] = Vx[i];

  for(int j=0; j<k; j++){

    for (iint i=0;i<m;i++)
      Vx[i] = V[j][i];

    // v[j+1] = invD*(A*v[j])
    axpy(A, 1.0, Vx, 0., V[j+1]);

    dotStar(invD, V[j+1]);

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
    	// H(i,j) = v[i]'*A*v[j]
    	dfloat hij = innerProd(*V[i], *V[j+1]);

    	// v[j+1] = v[j+1] - hij*v[i]
    	vectorAdd(-hij, V[i], 1.0, V[j+1]);

    	H[i + j*k] = (double) hij;
    }

    if(j+1 < k){
    	H[j+1+ j*k] = (double) norm(V[j+1]);

    	scaleVector(V[j+1], 1./H[j+1 + j*k]);
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

#include "parAlmond.h"


void pcg(parAlmond_t *parAlmond,
   dfloat *b,
   dfloat *x,
   iint maxIt,
   dfloat tol){

  csr *A = parAlmond->levels[0]->A;

  const iint m = A->Nrows;
  const iint n = A->Ncols;

  iint size, rank;
  dfloat localRes, globalRes;
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  parAlmond->ktype = PCG;

  // initial residual
  dfloat localNB = innerProd(m, b, b);
  dfloat globalNB;
  MPI_Allreduce(&localNB,&globalNB,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

  globalNB = sqrt(globalNB);

  // initial guess
  dfloat *x0 = (dfloat *) calloc(n, sizeof(dfloat));
  dfloat *r0 = (dfloat *) calloc(m, sizeof(dfloat));

  // initial residue  r0 = b - A*x0;
  zeqaxpy(parAlmond->levels[0]->A,-1.0, x0, 1.0, b, r0);

  dfloat *r = (dfloat *) calloc(m, sizeof(dfloat));
  for (iint i=0; i<m; i++) {
    r[i] = r0[i];
    x[i] = x0[i];
  }
  
  dfloat *di = (dfloat *) calloc(m, sizeof(dfloat));
  dfloat *Adi = (dfloat *) calloc(m, sizeof(dfloat));
  dfloat *di_previous = (dfloat *) calloc(m, sizeof(dfloat));

  dfloat rho, alpha, beta, rho_previous;
  dfloat rhoLocal, rhoGlobal, diAdiLocal, diAdiGlobal;

  iint flag = 1;

  dfloat *resvec = (dfloat *) calloc(m+1, sizeof(dfloat));
  resvec[0] = globalNB;

  for (iint i=1; i<=maxIt; i++){

    // apply preconditioner (make sure that precond is symmetric)
    // di = M^{-1}*r
    for (iint k=0;k<m;k++) 
      parAlmond->levels[0]->rhs[k] = r[k];
    kcycle(parAlmond,0);
    for (iint k=0;k<m;k++)
      di[k] = parAlmond->levels[0]->x[k];

    // rho = di'*r
    rhoLocal = innerProd(m, di, r);
    MPI_Allreduce(&rhoLocal,&rhoGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

    // TODO : try flexible conjugate gradient
    if(i > 1){
      beta = rhoGlobal/rho_previous;

      // di = di + beta*di_previous
      vectorAdd(m, beta, di_previous, 1.0, di);
    }

    //   Adi = A*di;
    axpy(A, 1.0, di, 0.0, Adi);
    diAdiLocal = innerProd(m, di, Adi);
    MPI_Allreduce(&diAdiLocal,&diAdiGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
    
    alpha =  rhoGlobal/ diAdiGlobal;

    // update solution
    //    x = x + alpha * di;
    vectorAdd(m, alpha, di, 1.0, x);

    // update residue
    // r = r - alpha * Adi;
    vectorAdd(m, -alpha, Adi, 1.0, r);

    localRes = innerProd(m, r, r);
    MPI_Allreduce(&localRes,&globalRes,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
    globalRes = sqrt(globalRes);

    for (iint k=0;k<m;k++)
      di_previous[k] = di[k];

    rho_previous = rhoGlobal;

    resvec[i] = globalRes;

    //printf("iter = %d, globalRes = %g\n",i, globalRes);

    if(globalRes < tol){ //*globalNB){
      flag = 0;
      break;
    }
  }
}

//TODO need to link the MPI processes in this solve
void pcg(parAlmond_t *parAlmond,
   occa::memory o_b,
   occa::memory o_x,
   iint maxIt,
   dfloat tol){

  hyb* A = parAlmond->levels[0]->deviceA;

  const iint n = A->Nrows;
  const iint m = A->Ncols;

  parAlmond->ktype = PCG;

  // initial residual
  dfloat nb = innerProd(parAlmond, m, o_b, o_b);
  nb = sqrt(nb);

  occa::memory o_x0, o_r0, o_r, o_di, o_Adi, o_di_previous;

  dfloat *dummy = (dfloat *) calloc(m, sizeof(dfloat));

  // initial guess
  o_x0 = parAlmond->device.malloc(m*sizeof(dfloat), dummy);
  o_r0 = parAlmond->device.malloc(m*sizeof(dfloat), dummy);
  o_r = parAlmond->device.malloc(m*sizeof(dfloat), dummy);
  o_di = parAlmond->device.malloc(m*sizeof(dfloat), dummy);
  o_Adi = parAlmond->device.malloc(m*sizeof(dfloat), dummy);
  o_di_previous = parAlmond->device.malloc(m*sizeof(dfloat), dummy);

  // initial residue  r0 = b - A*x0;
  zeqaxpy(parAlmond, A, -1.0, o_x0, 1.0, o_b, o_r0);

  //    r = r0;
  o_r.copyFrom(o_r0);

  //    x = x0;
  o_x.copyFrom(o_x0);

  dfloat rho, alpha, beta, rho_previous;

  iint flag = 1;


  dfloat *resvec = (dfloat *) calloc(maxIt+1, sizeof(dfloat));
  resvec[0] = nb;

  for (iint i=1; i<=maxIt; i++){
    // apply preconditioner (make sure that precond is symmetric)
    // di = M^{-1}*r
    parAlmond->levels[0]->o_rhs.copyFrom(o_r);
    device_kcycle(parAlmond,0);
    o_di.copyFrom(parAlmond->levels[0]->o_x);

    // rho = di'*r
    rho = innerProd(parAlmond, m, o_di, o_r);

    // TODO : try flexible conjugate gradient

    if(i > 1){
      beta = rho/rho_previous;

      // di = di + beta*di_previous
      vectorAdd(parAlmond, m, beta, o_di_previous, 1.0, o_di);
    }

    //   Adi = A*di;
    axpy(parAlmond, A, 1.0, o_di, 0.0, o_Adi);

    alpha =  rho/ innerProd(parAlmond, m, o_di, o_Adi);

    // update solution
    //    x = x + alpha * di;
    vectorAdd(parAlmond, m, alpha, o_di, 1.0, o_x);

    // update residue
    // r = r - alpha * Adi;
    vectorAdd(parAlmond, m, -alpha, o_Adi, 1.0, o_r);

    resvec[i] = innerProd(parAlmond, m, o_r, o_r);
    resvec[i] = sqrt(resvec[i]);

    //      di_previous = di;
    o_di_previous.copyFrom(o_di);
    rho_previous = rho;

    if(resvec[i] < tol*nb){
      flag = 0;
      break;
    }
  }
}



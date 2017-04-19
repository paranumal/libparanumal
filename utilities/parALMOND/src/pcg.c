#include "parAlmond.h"


void pcg(almond_t *almond,
   csr *A,
   dfloat *b,
   dfloat *x,
   iint maxIt,
   dfloat tol){

  const iint m = A->Nrows;
  const iint n = A->Ncols;

  iint size, rank;
  dfloat localRes, globalRes;
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  almond->ktype = PCG;

  // initial residual
  dfloat localNB = norm(m,b);
  dfloat globalNB;
  MPI_Allreduce(&localNB,&globalNB,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

  // initial guess
  dfloat *x0 = (dfloat *) calloc(n, sizeof(dfloat));
  dfloat *r0 = (dfloat *) calloc(m, sizeof(dfloat));

  // initial residue  r0 = b - A*x0;
  zeqaxpy(A, -1.0, x0, 1.0, b, r0);

  dfloat *r = (dfloat *) calloc(m, sizeof(dfloat));
  for (iint i=0; i<n; i++) {
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
    solve(almond, r, di);

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

    alpha =  rho/ diAdiGlobal;

    // update solution
    //    x = x + alpha * di;
    vectorAdd(m, alpha, di, 1.0, x);

    // update residue
    // r = r - alpha * Adi;
    vectorAdd(m, -alpha, Adi, 1.0, r);

    resvec[i] = norm(m, r);

    di_previous = di;
    rho_previous = rhoGlobal;

    localRes = resvec[i];

    MPI_Allreduce(&localRes,&globalRes,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

    if(globalRes < tol*globalNB){
      flag = 0;
      break;
    }
  }

  // r = b - A*x
  zeqaxpy(A, -1.0, x, 1.0, b, r);

  dfloat relres = norm(m, r)/globalNB;
}

//TODO need to link the MPI processes in this solve
void pcg(almond_t *almond,
   hyb *A,
   occa::memory o_b,
   occa::memory o_x,
   iint maxIt,
   dfloat tol){

  const iint m = A->Nrows;
  const iint n = A->Ncols;

  almond->ktype = PCG;

  // initial residual
  dfloat nb = innerProd(almond, m, o_b, o_b);
  nb = sqrt(nb);

  occa::memory o_x0, o_r0, o_r, o_di, o_Adi, o_di_previous;

  dfloat *dummy = (dfloat *) calloc(m, sizeof(dfloat));

  // initial guess
  o_x0 = almond->device.malloc(m*sizeof(dfloat), dummy);
  o_r0 = almond->device.malloc(m*sizeof(dfloat), dummy);
  o_r = almond->device.malloc(m*sizeof(dfloat), dummy);
  o_di = almond->device.malloc(m*sizeof(dfloat), dummy);
  o_Adi = almond->device.malloc(m*sizeof(dfloat), dummy);
  o_di_previous = almond->device.malloc(m*sizeof(dfloat), dummy);

  // initial residue  r0 = b - A*x0;
  zeqaxpy(almond, A, -1.0, o_x0, 1.0, o_b, o_r0);

  //    r = r0;
  copyVector(almond, m, o_r0, o_r);

  //    x = x0;
  copyVector(almond, m, o_x0, o_x);


  dfloat rho, alpha, beta, rho_previous;

  iint flag = 1;


  dfloat *resvec = (dfloat *) calloc(maxIt+1, sizeof(dfloat));
  resvec[0] = nb;

  for (iint i=1; i<=maxIt; i++){
    // apply preconditioner (make sure that precond is symmetric)
    // di = M^{-1}*r
    solve(almond, o_r, o_di);

    // rho = di'*r
    rho = innerProd(almond, m, o_di, o_r);

    // TODO : try flexible conjugate gradient

    if(i > 1){
      beta = rho/rho_previous;

      // di = di + beta*di_previous
      vectorAdd(almond, m, beta, o_di_previous, 1.0, o_di);
    }

    //   Adi = A*di;
    axpy(almond, A, 1.0, o_di, 0.0, o_Adi);

    alpha =  rho/ innerProd(almond, m, o_di, o_Adi);

    // update solution
    //    x = x + alpha * di;
    vectorAdd(almond, m, alpha, o_di, 1.0, o_x);

    // update residue
    // r = r - alpha * Adi;
    vectorAdd(almond, m, -alpha, o_Adi, 1.0, o_r);

    resvec[i] = innerProd(almond, m, o_r, o_r);
    resvec[i] = sqrt(resvec[i]);

    //      di_previous = di;
    copyVector(almond, m, o_di, o_di_previous);
    rho_previous = rho;

    if(resvec[i] < tol*nb){
      flag = 0;
      break;
    }
  }

  // r = b - A*x
  zeqaxpy(almond, A, -1.0, o_x, 1.0, o_b, o_r);

  dfloat relres = sqrt(innerProd(almond, m, o_r, o_r))/nb;
}



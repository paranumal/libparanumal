#include "agmg.h"


void pcg(parAlmond_t *parAlmond,
   iint maxIt,
   dfloat tol){

  csr *A = parAlmond->levels[0]->A;

  const iint m = A->Nrows;
  const iint n = A->Ncols;

  parAlmond->ktype = PCG;

  // use parAlmond's buffers
  dfloat *r = parAlmond->levels[0]->rhs;
  dfloat *z = parAlmond->levels[0]->x;

  // initial residual
  dfloat rdotr0Local = innerProd(m, r, r);
  dfloat rdotr0 = 0;
  MPI_Allreduce(&rdotr0Local,&rdotr0,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

  dfloat *x, *p, *Ap;

  x  = (dfloat *) calloc(m,sizeof(dfloat));
  Ap = (dfloat *) calloc(m,sizeof(dfloat));
  p  = (dfloat *) calloc(m,sizeof(dfloat));

  //    x = 0;
  setVector(m, x, 0.0);

  //sanity check
  if (rdotr0<=(tol*tol)) {
    for (iint i=0;i<m;i++)
      parAlmond->levels[0]->x[i] = x[i];

    free(x); free(p); free(Ap);
    return;
  }

  // Precondition, z = M^{-1}*r
  kcycle(parAlmond,0);
  for (iint i=0;i<m;i++)
    p[i] = z[i];

  dfloat rdotz0Local = innerProd(m, r, z);
  dfloat rdotz0 = 0;
  MPI_Allreduce(&rdotz0Local,&rdotz0,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  dfloat alpha, beta, pAp;

  iint Niter = 0;
  while(rdotr0>(tol*tol)){
    //   Ap = A*p;
    axpy(A, 1.0, p, 0.0, Ap,parAlmond->nullSpace,parAlmond->nullSpacePenalty);

    dfloat pApLocal = innerProd(m, p, Ap);
    pAp = 0;
    MPI_Allreduce(&pApLocal,&pAp,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

    alpha = rdotz0/pAp;

    // update solution
    //    x = x + alpha * p;
    vectorAdd(m, alpha, p, 1.0, x);

    // update residual
    // r = r - alpha * Ap;
    vectorAdd(m, -alpha, Ap, 1.0, r);


    dfloat rdotr1Local = innerProd(m, r, r);
    rdotr1 = 0;
    MPI_Allreduce(&rdotr1Local,&rdotr1,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

    if(rdotr1 < tol*tol) {
      rdotr0 = rdotr1;
      break;
    }

    // Precondition, z = M^{-1}*r
    kcycle(parAlmond,0);

    dfloat rdotz1Local = innerProd(m, r, z);
    rdotz1 = 0;
    MPI_Allreduce(&rdotz1Local,&rdotz1,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

  #if 1
    // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
    dfloat zdotApLocal = innerProd(m, z, Ap);
    dfloat zdotAp = 0;
    MPI_Allreduce(&zdotApLocal,&zdotAp,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
    beta = -alpha*zdotAp/rdotz0;
  #else
    beta = rdotz1/rdotz0;
  #endif

    // p = z + beta*p
    vectorAdd(m, 1.0, z, beta, p);

    // switch rdotz0 <= rdotz1
    rdotz0 = rdotz1;

    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;

    Niter++;

    printf("Almond PCG iter %d, res = %g\n", Niter, sqrt(rdotr0));

    if(Niter==maxIt) break;
  }

  //copy result back to parAlmond's x storage
  for (iint i=0;i<m;i++)
    parAlmond->levels[0]->x[i] = x[i];

  free(x); free(p); free(Ap);
}

void device_pcg(parAlmond_t *parAlmond, iint maxIt, dfloat tol){

  hyb* A = parAlmond->levels[0]->deviceA;

  const iint m = A->Nrows;
  const iint n = A->Ncols;

  parAlmond->ktype = PCG;

  // use parAlmond's buffers
  occa::memory &o_r = parAlmond->levels[0]->o_rhs;
  occa::memory &o_z = parAlmond->levels[0]->o_x;

  // initial residual
  dfloat rdotr0Local = innerProd(parAlmond, m, o_r, o_r);
  dfloat rdotr0 = 0;
  MPI_Allreduce(&rdotr0Local,&rdotr0,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

  occa::memory o_x, o_p, o_Ap;

  o_x  = parAlmond->device.malloc(n*sizeof(dfloat),parAlmond->levels[0]->x);
  o_Ap = parAlmond->device.malloc(n*sizeof(dfloat),parAlmond->levels[0]->x);
  o_p  = parAlmond->device.malloc(n*sizeof(dfloat),parAlmond->levels[0]->x);

  //    x = 0;
  setVector(parAlmond, m, o_x, 0.0);

  //sanity check
  if (rdotr0<=(tol*tol)) {
    parAlmond->levels[0]->o_x.copyFrom(o_x);
    printf("Almond PCG iter %d, res = %g\n", 0, sqrt(rdotr0));
    o_x.free(); o_p.free(); o_Ap.free();
    return;
  }

  // Precondition, z = M^{-1}*r
  device_kcycle(parAlmond,0);
  o_p.copyFrom(o_z);

  dfloat rdotz0Local = innerProd(parAlmond, m, o_r, o_z);
  dfloat rdotz0 = 0;
  MPI_Allreduce(&rdotz0Local,&rdotz0,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  dfloat alpha, beta, pAp;

  iint Niter = 0;
  while(rdotr0>(tol*tol)){
    //   Ap = A*p;
    axpy(parAlmond, A, 1.0, o_p, 0.0, o_Ap,parAlmond->nullSpace,parAlmond->nullSpacePenalty);

    dfloat pApLocal = innerProd(parAlmond, m, o_p, o_Ap);
    pAp = 0;
    MPI_Allreduce(&pApLocal,&pAp,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

    alpha = rdotz0/pAp;

    // update solution
    //    x = x + alpha * p;
    vectorAdd(parAlmond, m, alpha, o_p, 1.0, o_x);

    // update residual
    // r = r - alpha * Ap;
    vectorAdd(parAlmond, m, -alpha, o_Ap, 1.0, o_r);


    dfloat rdotr1Local = innerProd(parAlmond, m, o_r, o_r);
    rdotr1 = 0.;
    MPI_Allreduce(&rdotr1Local,&rdotr1,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

    if(rdotr1 < tol*tol) {
      rdotr0 = rdotr1;
      break;
    }

    // Precondition, z = M^{-1}*r
    device_kcycle(parAlmond,0);

    dfloat rdotz1Local = innerProd(parAlmond, m, o_r, o_z);
    rdotz1 = 0;
    MPI_Allreduce(&rdotz1Local,&rdotz1,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

  #if 1
    // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
    dfloat zdotApLocal = innerProd(parAlmond, m, o_z, o_Ap);
    dfloat zdotAp = 0;
    MPI_Allreduce(&zdotApLocal,&zdotAp,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
    beta = -alpha*zdotAp/rdotz0;
  #else
    beta = rdotz1/rdotz0;
  #endif

    // p = z + beta*p
    vectorAdd(parAlmond, m, 1.0, o_z, beta, o_p);

    // switch rdotz0 <= rdotz1
    rdotz0 = rdotz1;

    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;

    Niter++;

    //printf("Almond PCG iter %d, res = %g\n", Niter, sqrt(rdotr0));

    if(Niter==maxIt) break;
  }

  //copy result back to parAlmond's x storage
  parAlmond->levels[0]->o_x.copyFrom(o_x);

  printf("Almond PCG iter %d, res = %g\n", Niter, sqrt(rdotr0));

  o_x.free(); o_p.free(); o_Ap.free();
}



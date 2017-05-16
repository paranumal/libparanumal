#include "parAlmond.h"

parAlmond_t * agmgSetup(csr *A, dfloat *nullA, iint *globalRowStarts, const char* options){
  iint rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  parAlmond_t *parAlmond = (parAlmond_t *) calloc(1,sizeof(parAlmond_t));

  const iint coarseSize = 10;

  double seed = MPI_Wtime();
  double gSeed;
  MPI_Allreduce(&seed, &gSeed, 1, MPI_LONG, MPI_BXOR, MPI_COMM_WORLD);
  srand48(gSeed);

  agmgLevel **levels = (agmgLevel **) calloc(MAX_LEVELS,sizeof(agmgLevel *));

  levels[0] = (agmgLevel *) calloc(1,sizeof(agmgLevel));

  //copy A matrix
  levels[0]->A = (csr *) calloc(1,sizeof(csr));

  levels[0]->A->Nrows = A->Nrows;
  levels[0]->A->Ncols = A->Ncols;
  levels[0]->A->nnz = A->nnz;

  levels[0]->A->rowStarts = (iint *) calloc(A->Nrows+1,sizeof(iint));
  levels[0]->A->cols      = (iint *) calloc(A->nnz, sizeof(iint));
  levels[0]->A->coefs     = (dfloat *) calloc(A->nnz, sizeof(dfloat));
  for (iint i=0; i<A->Nrows+1; i++)
    levels[0]->A->rowStarts[i] = A->rowStarts[i];

  for (iint i=0; i<A->nnz; i++) {
    levels[0]->A->cols[i] = A->cols[i];   
    levels[0]->A->coefs[i] = A->coefs[i];   
  }
  levels[0]->A->NsendTotal = A->NsendTotal;
  levels[0]->A->NrecvTotal = A->NrecvTotal;
  levels[0]->A->NHalo      = A->NHalo;

  //set up level size
  levels[0]->nullA = nullA;
  levels[0]->Nrows = A->Nrows;
  levels[0]->Ncols = A->Ncols;

  if (globalRowStarts) {
    levels[0]->globalRowStarts = (iint *) calloc(size+1,sizeof(iint));
    for (iint r=0;r<size+1;r++) 
      levels[0]->globalRowStarts[r] = globalRowStarts[r];
  }

  int numLevels = 1;
  int lev =0;

  bool done = false;
  while(!done){
    const iint dim = levels[lev]->A->Nrows;
    csr *coarseA = (csr *) calloc(1,sizeof(csr));
    dfloat *nullCoarseA;

    coarsen(levels[lev], &coarseA, &nullCoarseA); 

    const iint coarseDim = coarseA->Nrows;

    // allocate vectors required
    allocate(levels[lev]);

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
      allocate(levels[lev+1]);
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

      levels[n]->x    = (dfloat *) calloc(Nmax,sizeof(dfloat));
      levels[n]->rhs  = (dfloat *) calloc(M,sizeof(dfloat));
      levels[n]->res  = (dfloat *) calloc(Nmax,sizeof(dfloat));
    }
    levels[numLevels-1]->A = distribute(levels[numLevels-1]->A,
                                  levels[numLevels-1]->globalRowStarts,
                                  levels[numLevels-1]->globalRowStarts);

    iint M    = levels[numLevels-1]->A->Nrows;
    iint Nmax = levels[numLevels-1]->A->Ncols;

    if (numLevels>1) Nmax = levels[numLevels-2]->P->Ncols > Nmax ? levels[numLevels-2]->P->Ncols : Nmax;

    levels[numLevels-1]->Nrows = M;
    levels[numLevels-1]->Ncols = Nmax;

    levels[numLevels-1]->x    = (dfloat *) calloc(Nmax,sizeof(dfloat));
    levels[numLevels-1]->rhs  = (dfloat *) calloc(M,sizeof(dfloat));
    levels[numLevels-1]->res  = (dfloat *) calloc(Nmax,sizeof(dfloat));
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
    parAlmond->levels[i]->deviceA = newHYB(parAlmond, parAlmond->levels[i]->A);
    //parAlmond->levels[i]->deviceA = newDCSR(parAlmond, parAlmond->levels[i]->A);
    if (i < parAlmond->numLevels-1) {
      parAlmond->levels[i]->dcsrP   = newDCSR(parAlmond, parAlmond->levels[i]->P);
      parAlmond->levels[i]->deviceR = newHYB(parAlmond, parAlmond->levels[i]->R);
      //parAlmond->levels[i]->deviceR = newDCSR(parAlmond, parAlmond->levels[i]->R);
    }

    iint N = parAlmond->levels[i]->Ncols;
    iint M = parAlmond->levels[i]->Nrows;

    parAlmond->levels[i]->o_x   = parAlmond->device.malloc(N*sizeof(dfloat), parAlmond->levels[i]->x);
    parAlmond->levels[i]->o_rhs = parAlmond->device.malloc(M*sizeof(dfloat), parAlmond->levels[i]->rhs);
    parAlmond->levels[i]->o_res = parAlmond->device.malloc(N*sizeof(dfloat), parAlmond->levels[i]->res);

    if(i > 0){
      parAlmond->levels[i]->o_ckp1 = parAlmond->device.malloc(N*sizeof(dfloat), parAlmond->levels[i]->x);
      parAlmond->levels[i]->o_dkp1 = parAlmond->device.malloc(N*sizeof(dfloat), parAlmond->levels[i]->x);
      parAlmond->levels[i]->o_vkp1 = parAlmond->device.malloc(M*sizeof(dfloat), parAlmond->levels[i]->x);
      parAlmond->levels[i]->o_wkp1 = parAlmond->device.malloc(M*sizeof(dfloat), parAlmond->levels[i]->x);
      parAlmond->levels[i]->o_rkp1 = parAlmond->device.malloc(M*sizeof(dfloat), parAlmond->levels[i]->x);
    }
  }

  //if using matrix-free A action, free unnecessary buffers
  if (strstr(parAlmond->options,"MATRIXFREE")) {
    parAlmond->levels[0]->deviceA->E->o_cols.free();
    parAlmond->levels[0]->deviceA->E->o_coefs.free();
    parAlmond->levels[0]->deviceA->C->o_rows.free();
    parAlmond->levels[0]->deviceA->C->o_cols.free();
    parAlmond->levels[0]->deviceA->C->o_coefs.free();
    parAlmond->levels[0]->deviceA->C->o_temp_rows.free();
    parAlmond->levels[0]->deviceA->C->o_temp_Ax.free();
    if(parAlmond->levels[0]->deviceA->NsendTotal) {
      parAlmond->levels[0]->deviceA->o_haloElementList.free();
      parAlmond->levels[0]->deviceA->o_haloBuffer.free();
    }
  }
}

void kcycle(parAlmond_t *parAlmond, int k){

  iint m = parAlmond->levels[k]->Nrows;
  iint n = parAlmond->levels[k]->Ncols;
  iint mCoarse = parAlmond->levels[k+1]->Nrows;
  iint nCoarse = parAlmond->levels[k+1]->Ncols;

  char name[BUFSIZ];
  sprintf(name, "host kcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  scaleVector(m, parAlmond->levels[k]->x, 0.0);

  smooth(parAlmond->levels[k], parAlmond->levels[k]->rhs, parAlmond->levels[k]->x, true);

    // res = - A*x + rhs (i.e., rhs - A*x)
  zeqaxpy(parAlmond->levels[k]->A, -1.0, parAlmond->levels[k]->x, 1.0, 
          parAlmond->levels[k]->rhs, parAlmond->levels[k]->res);

  // restrict the residual to next level
  restrict(parAlmond->levels[k], parAlmond->levels[k]->res, parAlmond->levels[k+1]->rhs);

  dfloat *ckp1 = (dfloat *) calloc(nCoarse,sizeof(dfloat)); 
  dfloat *vkp1 = (dfloat *) calloc(mCoarse,sizeof(dfloat)); 
  dfloat *dkp1 = (dfloat *) calloc(nCoarse,sizeof(dfloat)); 
  dfloat *wkp1 = (dfloat *) calloc(mCoarse,sizeof(dfloat)); 
  dfloat *rkp1 = (dfloat *) calloc(mCoarse,sizeof(dfloat));

  if(k+1 < parAlmond->numLevels - 1){

    // first inner krylov iteration
    kcycle(parAlmond, k+1);

    for (iint i=0;i<mCoarse; i++)
      ckp1[i] = parAlmond->levels[k+1]->x[i];

    // v = A*c
    axpy(parAlmond->levels[k+1]->A, 1.0, ckp1, 0.0, vkp1);

    for (iint i=0;i<mCoarse; i++)
      rkp1[i] = parAlmond->levels[k+1]->rhs[i];

    dfloat rho1, alpha1, norm_rkp1, norm_rktilde_p;
    dfloat rho1Global, alpha1Global, norm_rkp1Global, norm_rktilde_pGlobal;

    if(parAlmond->ktype == PCG){
      rho1   = innerProd(mCoarse, vkp1, ckp1);
      alpha1 = innerProd(mCoarse, rkp1, ckp1);
      MPI_Allreduce(&rho1,&rho1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&alpha1,&alpha1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
    }

    if(parAlmond->ktype == GMRES){
      rho1   = innerProd(mCoarse,vkp1, vkp1);
      alpha1 = innerProd(mCoarse,vkp1, rkp1);
      MPI_Allreduce(&rho1,&rho1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&alpha1,&alpha1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
    }

    norm_rkp1 = innerProd(mCoarse,rkp1,rkp1);
    MPI_Allreduce(&norm_rkp1,&norm_rkp1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
    norm_rkp1Global = sqrt(norm_rkp1Global);

    // rkp1 = rkp1 - (alpha1/rho1)*vkp1
    vectorAdd(mCoarse, -alpha1Global/rho1Global, vkp1, 1.0, rkp1);

    norm_rktilde_p = innerProd(mCoarse,rkp1,rkp1);
    MPI_Allreduce(&norm_rktilde_p,&norm_rktilde_pGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
    norm_rktilde_pGlobal = sqrt(norm_rktilde_pGlobal);

    for (iint i=0;i<mCoarse; i++)
      parAlmond->levels[k+1]->rhs[i] = rkp1[i];

    dfloat t = 0.2;

    if(norm_rktilde_pGlobal < t*norm_rkp1Global){
      //      parAlmond->levels[k+1]->x = (alpha1/rho1)*ckp1
      dfloat adivrho = alpha1Global/rho1Global;
      vectorAdd(mCoarse, adivrho, ckp1, 0., parAlmond->levels[k+1]->x);
    } else{
    
      kcycle(parAlmond, k+1);

      dfloat *dkp1 = (dfloat *) calloc(nCoarse,sizeof(dfloat)); 
      dfloat *wkp1 = (dfloat *) calloc(mCoarse,sizeof(dfloat));

      for (iint i=0;i<mCoarse; i++)
        dkp1[i] = parAlmond->levels[k+1]->x[i];

      // w = A*d
      axpy(parAlmond->levels[k+1]->A, 1.0, dkp1, 0., wkp1);

      dfloat gamma, beta, alpha2;
      dfloat gammaGlobal, betaGlobal, alpha2Global;

      if(parAlmond->ktype == PCG){
        gamma  = innerProd(mCoarse, vkp1, dkp1);
        beta   = innerProd(mCoarse, wkp1, dkp1);
        alpha2 = innerProd(mCoarse, rkp1, dkp1);
        MPI_Allreduce(&gamma,&gammaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&beta,&betaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&alpha2,&alpha2Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      }
      if(parAlmond->ktype == GMRES){
        gamma  = innerProd(mCoarse, wkp1, vkp1);
        beta   = innerProd(mCoarse, wkp1, wkp1);
        alpha2 = innerProd(mCoarse, wkp1, rkp1);
        MPI_Allreduce(&gamma,&gammaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&beta,&betaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&alpha2,&alpha2Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      }

      if(fabs(rho1Global) > (dfloat) 1e-20){

        dfloat rho2 = betaGlobal - gammaGlobal*gammaGlobal/rho1Global;

        if(fabs(rho2) > (dfloat) 1e-20){

          // parAlmond->levels[k+1]->x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ckp1 + (alpha2/rho2)*dkp1
          dfloat a = alpha1Global/rho1Global - gammaGlobal*alpha2Global/(rho1Global*rho2);
          dfloat b = alpha2Global/rho2;

          vectorAdd(mCoarse, a, ckp1, 0.0, parAlmond->levels[k+1]->x);
          vectorAdd(mCoarse, b, dkp1, 1.0, parAlmond->levels[k+1]->x);
        }
      }
    }
  } else {
    if (parAlmond->Acoarse != NULL) {
      //use coarse sovler 
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      for (iint n=0;n<mCoarse;n++) 
        parAlmond->rhsCoarse[n+parAlmond->coarseOffset] = parAlmond->levels[k+1]->rhs[n];

      xxtSolve(parAlmond->xCoarse, parAlmond->Acoarse, parAlmond->rhsCoarse); 

      for (iint n=0;n<mCoarse;n++) 
        parAlmond->levels[k+1]->x[n] = parAlmond->xCoarse[n+parAlmond->coarseOffset];
    } else {
      scaleVector(mCoarse, parAlmond->levels[k+1]->x, 0.);
      smooth(parAlmond->levels[k+1], parAlmond->levels[k+1]->rhs, parAlmond->levels[k+1]->x, true);
    }
  }


  interpolate(parAlmond->levels[k], parAlmond->levels[k+1]->x, parAlmond->levels[k]->x);

  smooth(parAlmond->levels[k], parAlmond->levels[k]->rhs, parAlmond->levels[k]->x,false);

  occaTimerToc(parAlmond->device,name);
}


void device_kcycle(parAlmond_t *parAlmond, int k){

  iint m = parAlmond->levels[k]->Nrows;
  iint n = parAlmond->levels[k]->Ncols;
  iint mCoarse = parAlmond->levels[k+1]->Nrows;
  iint nCoarse = parAlmond->levels[k+1]->Ncols;

  char name[BUFSIZ];
  sprintf(name, "device kcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  scaleVector(parAlmond, m, parAlmond->levels[k]->o_x, 0.0);

  //use matrix free action if its been given
  if ((k==0)&&strstr(parAlmond->options,"MATRIXFREE")) {
    matFreeSmooth(parAlmond, parAlmond->levels[k], parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_x, true);

    matFreeZeqAXPY(parAlmond, parAlmond->levels[k]->deviceA, -1.0, parAlmond->levels[k]->o_x,  1.0,
             parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_res);
  } else {
    smooth(parAlmond, parAlmond->levels[k], parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_x, true);
  
    // res = - A*x + rhs (i.e., rhs - A*x)
    zeqaxpy(parAlmond, parAlmond->levels[k]->deviceA, -1.0, parAlmond->levels[k]->o_x,  1.0,
             parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_res);
  }

  // restrict the residual to next level
  restrict(parAlmond, parAlmond->levels[k], parAlmond->levels[k]->o_res, parAlmond->levels[k+1]->o_rhs);

  if(k+1 < parAlmond->numLevels - 1){
    if(k>2) {
      device_vcycle(parAlmond,k+1);
      //device_kcycle(parAlmond, k+1);
    } else{
      // first inner krylov iteration
      device_kcycle(parAlmond,k+1);

      //ckp1 = parAlmond->levels[k+1]->x;
      copyVector(parAlmond, mCoarse, parAlmond->levels[k+1]->o_x, parAlmond->levels[k+1]->o_ckp1);

      // v = A*c
      axpy(parAlmond, parAlmond->levels[k+1]->deviceA, 1.0, parAlmond->levels[k+1]->o_ckp1, 0.0,
              parAlmond->levels[k+1]->o_vkp1);

      // rkp1 = parAlmond->levels[k+1]->rhs;
      copyVector(parAlmond, mCoarse, parAlmond->levels[k+1]->o_rhs, parAlmond->levels[k+1]->o_rkp1);

      dfloat rho1Local, alpha1Local, norm_rkp1Local, norm_rktilde_pLocal;
      dfloat rho1Global, alpha1Global, norm_rkp1Global, norm_rktilde_pGlobal;

      if(parAlmond->ktype == PCG){
        rho1Local   = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_ckp1, parAlmond->levels[k+1]->o_vkp1);
        alpha1Local = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_ckp1, parAlmond->levels[k+1]->o_rkp1);
        MPI_Allreduce(&rho1Local,&rho1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&alpha1Local,&alpha1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      }

      if(parAlmond->ktype == GMRES){
        rho1Local   = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_vkp1, parAlmond->levels[k+1]->o_vkp1);
        alpha1Local = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_vkp1, parAlmond->levels[k+1]->o_rkp1);
        MPI_Allreduce(&rho1Local,&rho1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&alpha1Local,&alpha1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      }
   
      norm_rkp1Local = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_rkp1,parAlmond->levels[k+1]->o_rkp1);
      MPI_Allreduce(&norm_rkp1Local,&norm_rkp1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      norm_rkp1Global = sqrt(norm_rkp1Global);

      // rkp1 = rkp1 - (alpha1/rho1)*vkp1
      vectorAdd(parAlmond, mCoarse, -alpha1Global/rho1Global, parAlmond->levels[k+1]->o_vkp1, 1.0,
          parAlmond->levels[k+1]->o_rkp1 );

      norm_rktilde_pLocal = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_rkp1,parAlmond->levels[k+1]->o_rkp1);
      MPI_Allreduce(&norm_rktilde_pLocal,&norm_rktilde_pGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      norm_rktilde_pGlobal = sqrt(norm_rktilde_pGlobal);

      //      parAlmond->levels[k+1]->rhs = rkp1;
      copyVector(parAlmond, mCoarse, parAlmond->levels[k+1]->o_rkp1, parAlmond->levels[k+1]->o_rhs);

      dfloat t = 0.2;

      if(norm_rktilde_pGlobal < t*norm_rkp1Global){
        //      parAlmond->levels[k+1]->x = (alpha1/rho1)*ckp1
        vectorAdd(parAlmond, mCoarse, alpha1Global/rho1Global, parAlmond->levels[k+1]->o_ckp1, 0., parAlmond->levels[k+1]->o_x);
      } else{

        device_kcycle(parAlmond,k+1);

        //  dkp1 = parAlmond->levels[k+1]->x;
        copyVector(parAlmond, mCoarse, parAlmond->levels[k+1]->o_x, parAlmond->levels[k+1]->o_dkp1);

        // w = A*d
        axpy(parAlmond, parAlmond->levels[k+1]->deviceA, 1.0, parAlmond->levels[k+1]->o_dkp1, 0.,
                parAlmond->levels[k+1]->o_wkp1);

        dfloat gammaLocal, betaLocal, alpha2Local;
        dfloat gammaGlobal, betaGlobal, alpha2Global;

        if(parAlmond->ktype == PCG){
          gammaLocal  = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_dkp1, parAlmond->levels[k+1]->o_vkp1);
          betaLocal   = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_dkp1, parAlmond->levels[k+1]->o_wkp1);
          alpha2Local = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_dkp1, parAlmond->levels[k+1]->o_rkp1);
          MPI_Allreduce(&gammaLocal,&gammaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
          MPI_Allreduce(&betaLocal,&betaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
          MPI_Allreduce(&alpha2Local,&alpha2Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        }
        if(parAlmond->ktype == GMRES){
          gammaLocal  = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_wkp1, parAlmond->levels[k+1]->o_vkp1);
          betaLocal   = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_wkp1, parAlmond->levels[k+1]->o_wkp1);
          alpha2Local = innerProd(parAlmond, mCoarse, parAlmond->levels[k+1]->o_wkp1, parAlmond->levels[k+1]->o_rkp1);
          MPI_Allreduce(&gammaLocal,&gammaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
          MPI_Allreduce(&betaLocal,&betaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
          MPI_Allreduce(&alpha2Local,&alpha2Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        }

        if(fabs(rho1Global) > (dfloat) 1e-20){

          dfloat rho2 = betaGlobal - gammaGlobal*gammaGlobal/rho1Global;

          if(fabs(rho2) > (dfloat) 1e-20){
            // parAlmond->levels[k+1]->x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ckp1 + (alpha2/rho2)*dkp1

            dfloat a = alpha1Global/rho1Global - gammaGlobal*alpha2Global/(rho1Global*rho2);
            dfloat b = alpha2Global/rho2;

            vectorAdd(parAlmond, mCoarse, a, parAlmond->levels[k+1]->o_ckp1, 0.0, parAlmond->levels[k+1]->o_x);
            vectorAdd(parAlmond, mCoarse, b, parAlmond->levels[k+1]->o_dkp1, 1.0, parAlmond->levels[k+1]->o_x);
          }
        }
      }
    }
  } else{
    if (parAlmond->Acoarse != NULL) {
      //use coarse sovler 
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      parAlmond->levels[k+1]->o_rhs.copyTo(parAlmond->rhsCoarse+parAlmond->coarseOffset);
      xxtSolve(parAlmond->xCoarse, parAlmond->Acoarse, parAlmond->rhsCoarse); 
      parAlmond->levels[k+1]->o_x.copyFrom(parAlmond->xCoarse+parAlmond->coarseOffset,mCoarse*sizeof(dfloat));  
    } else {
      scaleVector(parAlmond, mCoarse, parAlmond->levels[k+1]->o_x, 0.);
      smooth(parAlmond, parAlmond->levels[k+1], parAlmond->levels[k+1]->o_rhs, parAlmond->levels[k+1]->o_x, true);
    }
  }

  interpolate(parAlmond, parAlmond->levels[k], parAlmond->levels[k+1]->o_x, parAlmond->levels[k]->o_x);

  //use matrix free action if its been given
  if ((k==0)&&strstr(parAlmond->options,"MATRIXFREE")) {
    matFreeSmooth(parAlmond, parAlmond->levels[k], parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_x, false);
  } else {
    smooth(parAlmond, parAlmond->levels[k], parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_x, false);
  }

  occaTimerToc(parAlmond->device,name);
}



void vcycle(parAlmond_t *parAlmond, int k) {

  const iint m = parAlmond->levels[k]->Nrows;
  const iint mCoarse = parAlmond->levels[k+1]->Nrows;

  char name[BUFSIZ];
  sprintf(name, "host vcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  scaleVector(m, parAlmond->levels[k]->x,  0.0);

  smooth(parAlmond->levels[k], parAlmond->levels[k]->rhs, parAlmond->levels[k]->x, true);

  // res = rhs - A*x
  zeqaxpy(parAlmond->levels[k]->A, -1.0, parAlmond->levels[k]->x, 1.0, parAlmond->levels[k]->rhs,
     parAlmond->levels[k]->res);

  // restrict the residual to next level
  restrict(parAlmond->levels[k], parAlmond->levels[k]->res, parAlmond->levels[k+1]->rhs);

  if(k+1 < parAlmond->numLevels - 1){
    vcycle(parAlmond,k+1);
  } else{
    if (parAlmond->Acoarse != NULL) {
      //use coarse sovler 
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      for (iint n=0;n<mCoarse;n++) 
        parAlmond->rhsCoarse[n+parAlmond->coarseOffset] = parAlmond->levels[k+1]->rhs[n];

      xxtSolve(parAlmond->xCoarse, parAlmond->Acoarse, parAlmond->rhsCoarse); 

      for (iint n=0;n<mCoarse;n++) 
        parAlmond->levels[k+1]->x[n] = parAlmond->xCoarse[n+parAlmond->coarseOffset];
    } else {
      scaleVector(mCoarse, parAlmond->levels[k+1]->x, 0.);
      smooth(parAlmond->levels[k+1], parAlmond->levels[k+1]->rhs, parAlmond->levels[k+1]->x,true);
    }
  }

  interpolate(parAlmond->levels[k], parAlmond->levels[k+1]->x, parAlmond->levels[k]->x);
  smooth(parAlmond->levels[k], parAlmond->levels[k]->rhs, parAlmond->levels[k]->x,false);

  occaTimerToc(parAlmond->device,name);
}


void device_vcycle(parAlmond_t *parAlmond, int k){

#define GPU_CPU_SWITCH_SIZE 1024 //TODO move this the the almond struct?

  const iint m = parAlmond->levels[k]->Nrows;
  const iint mCoarse = parAlmond->levels[k+1]->Nrows;

  // switch to cpu if the problem size is too small for gpu
  if(m < GPU_CPU_SWITCH_SIZE){
    parAlmond->levels[k]->o_rhs.copyTo(&(parAlmond->levels[k]->rhs[0]), m*sizeof(dfloat));
    vcycle(parAlmond, k);
    parAlmond->levels[k]->o_x.copyFrom(&(parAlmond->levels[k]->x[0]), m*sizeof(dfloat));
    return;
  }

  char name[BUFSIZ];
  sprintf(name, "device vcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  scaleVector(parAlmond, m, parAlmond->levels[k]->o_x, 0.0);

  if ((k==0)&&strstr(parAlmond->options,"MATRIXFREE")){
    matFreeSmooth(parAlmond, parAlmond->levels[k], parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_x, true);
    matFreeZeqAXPY(parAlmond, parAlmond->levels[k]->deviceA,-1.0, parAlmond->levels[k]->o_x,  1.0,
             parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_res);
  } else {
    smooth(parAlmond, parAlmond->levels[k], parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_x, true);
    // res = rhs - A*x
    zeqaxpy(parAlmond, parAlmond->levels[k]->deviceA, -1.0, parAlmond->levels[k]->o_x,  1.0,
             parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_res);
  }

  // restrict the residual to next level
  restrict(parAlmond, parAlmond->levels[k], parAlmond->levels[k]->o_res, parAlmond->levels[k+1]->o_rhs);


  if(k+1 < parAlmond->numLevels - 1){
    device_vcycle(parAlmond, k+1);
  }else{
    if (parAlmond->Acoarse != NULL) {
      //use coarse sovler 
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      parAlmond->levels[k+1]->o_rhs.copyTo(parAlmond->rhsCoarse+parAlmond->coarseOffset);
      xxtSolve(parAlmond->xCoarse, parAlmond->Acoarse, parAlmond->rhsCoarse); 
      parAlmond->levels[k+1]->o_x.copyFrom(parAlmond->xCoarse+parAlmond->coarseOffset,mCoarse*sizeof(dfloat));  
    } else {
      scaleVector(parAlmond, mCoarse, parAlmond->levels[k+1]->o_x, 0.);
      smooth(parAlmond, parAlmond->levels[k+1], parAlmond->levels[k+1]->o_rhs, parAlmond->levels[k+1]->o_x, true);
    }
  }

  interpolate(parAlmond, parAlmond->levels[k], parAlmond->levels[k+1]->o_x, parAlmond->levels[k]->o_x);

  if ((k==0)&&strstr(parAlmond->options,"MATRIXFREE")){
    matFreeSmooth(parAlmond, parAlmond->levels[k], parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_x, false);
  } else {
    smooth(parAlmond, parAlmond->levels[k], parAlmond->levels[k]->o_rhs, parAlmond->levels[k]->o_x,false);
  }

  occaTimerToc(parAlmond->device,name);
}

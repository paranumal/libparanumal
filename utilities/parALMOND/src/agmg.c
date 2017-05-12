#include "parAlmond.h"

almond_t * setup(csr *A, dfloat *nullA, iint *globalRowStarts){
  iint rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  almond_t *almond = (almond_t *) calloc(1,sizeof(almond_t));

  const iint coarseSize = 10;

  double seed = MPI_Wtime();
  double gSeed;
  MPI_Allreduce(&seed, &gSeed, 1, MPI_LONG, MPI_BXOR, MPI_COMM_WORLD);
  srand48(gSeed);

  agmgLevel **levels = (agmgLevel **) calloc(MAX_LEVELS,sizeof(agmgLevel *));

  levels[0] = (agmgLevel *) calloc(1,sizeof(agmgLevel));
  levels[0]->A = A;           //TODO maybe these should be copies, not just pointer assignments
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

  almond->ktype = PCG;


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

  almond->levels = levels;
  almond->numLevels = numLevels;

  return almond;
}


void sync_setup_on_device(almond_t *almond, occa::device dev){
  //set occa device pointer
  almond->device = dev;
  buildAlmondKernels(almond);

  for(int i=0; i<almond->numLevels; i++){
    almond->levels[i]->deviceA = newHYB(almond, almond->levels[i]->A);
    //almond->levels[i]->deviceA = newDCSR(almond, almond->levels[i]->A);
    if (i < almond->numLevels-1) {
      almond->levels[i]->dcsrP   = newDCSR(almond, almond->levels[i]->P);
      almond->levels[i]->deviceR = newHYB(almond, almond->levels[i]->R);
      //almond->levels[i]->deviceR = newDCSR(almond, almond->levels[i]->R);
    }

    iint N = almond->levels[i]->Ncols;
    iint M = almond->levels[i]->Nrows;

    almond->levels[i]->o_x   = almond->device.malloc(N*sizeof(dfloat), almond->levels[i]->x);
    almond->levels[i]->o_rhs = almond->device.malloc(M*sizeof(dfloat), almond->levels[i]->rhs);
    almond->levels[i]->o_res = almond->device.malloc(N*sizeof(dfloat), almond->levels[i]->res);

    if(i > 0){
      almond->levels[i]->o_ckp1 = almond->device.malloc(N*sizeof(dfloat), almond->levels[i]->x);
      almond->levels[i]->o_dkp1 = almond->device.malloc(N*sizeof(dfloat), almond->levels[i]->x);
      almond->levels[i]->o_vkp1 = almond->device.malloc(M*sizeof(dfloat), almond->levels[i]->x);
      almond->levels[i]->o_wkp1 = almond->device.malloc(M*sizeof(dfloat), almond->levels[i]->x);
      almond->levels[i]->o_rkp1 = almond->device.malloc(M*sizeof(dfloat), almond->levels[i]->x);
    }
  }
}

void solve(almond_t *almond, dfloat *rhs, dfloat *x){
  //copy rhs and zero x
  for (iint n=0;n<almond->levels[0]->A->Nrows;n++) {
    almond->levels[0]->rhs[n] = rhs[n];
    almond->levels[0]->x[n] = 0.0;
  }

  kcycle(almond, 0);
  //vcycle(almond, 0);

  //copy back
  for (iint n=0;n<almond->levels[0]->A->Nrows;n++)
    x[n] = almond->levels[0]->x[n];
}

void solve(almond_t *almond, occa::memory o_rhs, occa::memory o_x){
  const iint N = almond->levels[0]->deviceA->Nrows;

  //copy rhs and zero x
  copyVector(almond, N, o_rhs, almond->levels[0]->o_rhs);
  scaleVector(almond, N, almond->levels[0]->o_x, 0.);

  device_kcycle(almond, 0);
  //device_vcycle(almond, 0);

  //copy back
  copyVector(almond, N, almond->levels[0]->o_x, o_x);
}

void kcycle(almond_t *almond, int k){

  iint m = almond->levels[k]->Nrows;
  iint n = almond->levels[k]->Ncols;
  iint mCoarse = almond->levels[k+1]->Nrows;
  iint nCoarse = almond->levels[k+1]->Ncols;

  char name[BUFSIZ];
  sprintf(name, "host kcycle level %d", k);
  occaTimerTic(almond->device,name);

  // if its not first level zero it out
  if(k>0)
    scaleVector(m, almond->levels[k]->x, 0.0);


  smooth(almond->levels[k], almond->levels[k]->rhs, almond->levels[k]->x, true);

    // res = - A*x + rhs (i.e., rhs - A*x)
  zeqaxpy(almond->levels[k]->A, -1.0, almond->levels[k]->x, 1.0, 
          almond->levels[k]->rhs, almond->levels[k]->res);

  // restrict the residual to next level
  restrict(almond->levels[k], almond->levels[k]->res, almond->levels[k+1]->rhs);

  dfloat *ckp1 = (dfloat *) calloc(nCoarse,sizeof(dfloat)); 
  dfloat *vkp1 = (dfloat *) calloc(mCoarse,sizeof(dfloat)); 
  dfloat *dkp1 = (dfloat *) calloc(nCoarse,sizeof(dfloat)); 
  dfloat *wkp1 = (dfloat *) calloc(mCoarse,sizeof(dfloat)); 
  dfloat *rkp1 = (dfloat *) calloc(mCoarse,sizeof(dfloat));

  if(k+1 < almond->numLevels - 1){

    // first inner krylov iteration
    kcycle(almond, k+1);

    for (iint i=0;i<mCoarse; i++)
      ckp1[i] = almond->levels[k+1]->x[i];

    // v = A*c
    axpy(almond->levels[k+1]->A, 1.0, ckp1, 0.0, vkp1);

    for (iint i=0;i<mCoarse; i++)
      rkp1[i] = almond->levels[k+1]->rhs[i];

    dfloat rho1, alpha1, norm_rkp1, norm_rktilde_p;
    dfloat rho1Global, alpha1Global, norm_rkp1Global, norm_rktilde_pGlobal;

    if(almond->ktype == PCG){
      rho1   = innerProd(mCoarse, vkp1, ckp1);
      alpha1 = innerProd(mCoarse, rkp1, ckp1);
      MPI_Allreduce(&rho1,&rho1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&alpha1,&alpha1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
    }

    if(almond->ktype == GMRES){
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
      almond->levels[k+1]->rhs[i] = rkp1[i];

    dfloat t = 0.2;

    if(norm_rktilde_pGlobal < t*norm_rkp1Global){
      //      almond->levels[k+1]->x = (alpha1/rho1)*ckp1
      dfloat adivrho = alpha1Global/rho1Global;
      vectorAdd(mCoarse, adivrho, ckp1, 0., almond->levels[k+1]->x);
    } else{
    
      kcycle(almond, k+1);

      dfloat *dkp1 = (dfloat *) calloc(nCoarse,sizeof(dfloat)); 
      dfloat *wkp1 = (dfloat *) calloc(mCoarse,sizeof(dfloat));

      for (iint i=0;i<mCoarse; i++)
        dkp1[i] = almond->levels[k+1]->x[i];

      // w = A*d
      axpy(almond->levels[k+1]->A, 1.0, dkp1, 0., wkp1);

      dfloat gamma, beta, alpha2;
      dfloat gammaGlobal, betaGlobal, alpha2Global;

      if(almond->ktype == PCG){
        gamma  = innerProd(mCoarse, vkp1, dkp1);
        beta   = innerProd(mCoarse, wkp1, dkp1);
        alpha2 = innerProd(mCoarse, rkp1, dkp1);
        MPI_Allreduce(&gamma,&gammaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&beta,&betaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&alpha2,&alpha2Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      }
      if(almond->ktype == GMRES){
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

          // almond->levels[k+1]->x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ckp1 + (alpha2/rho2)*dkp1
          dfloat a = alpha1Global/rho1Global - gammaGlobal*alpha2Global/(rho1Global*rho2);
          dfloat b = alpha2Global/rho2;

          vectorAdd(mCoarse, a, ckp1, 0.0, almond->levels[k+1]->x);
          vectorAdd(mCoarse, b, dkp1, 1.0, almond->levels[k+1]->x);
        }
      }
    }
  } else {
    scaleVector(mCoarse, almond->levels[k+1]->x, 0.);
    if (almond->coarseSolve != NULL) {
      //use coarse sovler 
      dfloat *xCoarse   = (dfloat*) calloc(almond->coarseTotal,sizeof(dfloat));
      dfloat *rhsCoarse = (dfloat*) calloc(almond->coarseTotal,sizeof(dfloat));

      for (iint n=0;n<mCoarse;n++) 
        rhsCoarse[n+almond->coarseOffset] = almond->levels[k+1]->rhs[n];

      almond->coarseSolve((void *) xCoarse, almond->ACoarse, (void *) rhsCoarse); 

      for (iint n=0;n<mCoarse;n++) 
        almond->levels[k+1]->x[n] = xCoarse[n+almond->coarseOffset];
    } else {
      smooth(almond->levels[k+1], almond->levels[k+1]->rhs, almond->levels[k+1]->x, true);
    }
  }


  interpolate(almond->levels[k], almond->levels[k+1]->x, almond->levels[k]->x);

  smooth(almond->levels[k], almond->levels[k]->rhs, almond->levels[k]->x,false);

  occaTimerToc(almond->device,name);
}


void device_kcycle(almond_t *almond, int k){

  iint m = almond->levels[k]->Nrows;
  iint n = almond->levels[k]->Ncols;
  iint mCoarse = almond->levels[k+1]->Nrows;
  iint nCoarse = almond->levels[k+1]->Ncols;

  char name[BUFSIZ];
  sprintf(name, "device kcycle level %d", k);
  occaTimerTic(almond->device,name);

  // if its not first level zero it out
  if(k>0)
    scaleVector(almond, m, almond->levels[k]->o_x, 0.0);

  //use matrix free action if its been given
  //if ((k==0)&&almond->matFreeAX) {
  //  matFreeSmooth(almond, almond->levels[k]->o_rhs, almond->levels[k]->o_x, true);
  //  matFreeZeqAXPY(almond, -1.0, almond->levels[k]->o_x,  1.0,
  //           almond->levels[k]->o_rhs, almond->levels[k]->o_res);
  //} else {
    smooth(almond, almond->levels[k], almond->levels[k]->o_rhs, almond->levels[k]->o_x, true);

    // res = - A*x + rhs (i.e., rhs - A*x)
    zeqaxpy(almond, almond->levels[k]->deviceA, -1.0, almond->levels[k]->o_x,  1.0,
             almond->levels[k]->o_rhs, almond->levels[k]->o_res);
  //}

  // restrict the residual to next level
  restrict(almond, almond->levels[k], almond->levels[k]->o_res, almond->levels[k+1]->o_rhs);

  if(k+1 < almond->numLevels - 1){
    if(k>2) {
      device_vcycle(almond,k+1);
      //device_kcycle(almond, k+1);
    } else{
      // first inner krylov iteration
      device_kcycle(almond,k+1);

      //ckp1 = almond->levels[k+1]->x;
      copyVector(almond, mCoarse, almond->levels[k+1]->o_x, almond->levels[k+1]->o_ckp1);

      // v = A*c
      axpy(almond, almond->levels[k+1]->deviceA, 1.0, almond->levels[k+1]->o_ckp1, 0.0,
              almond->levels[k+1]->o_vkp1);

      // rkp1 = almond->levels[k+1]->rhs;
      copyVector(almond, mCoarse, almond->levels[k+1]->o_rhs, almond->levels[k+1]->o_rkp1);

      dfloat rho1Local, alpha1Local, norm_rkp1Local, norm_rktilde_pLocal;
      dfloat rho1Global, alpha1Global, norm_rkp1Global, norm_rktilde_pGlobal;

      if(almond->ktype == PCG){
        rho1Local   = innerProd(almond, mCoarse, almond->levels[k+1]->o_ckp1, almond->levels[k+1]->o_vkp1);
        alpha1Local = innerProd(almond, mCoarse, almond->levels[k+1]->o_ckp1, almond->levels[k+1]->o_rkp1);
        MPI_Allreduce(&rho1Local,&rho1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&alpha1Local,&alpha1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      }

      if(almond->ktype == GMRES){
        rho1Local   = innerProd(almond, mCoarse, almond->levels[k+1]->o_vkp1, almond->levels[k+1]->o_vkp1);
        alpha1Local = innerProd(almond, mCoarse, almond->levels[k+1]->o_vkp1, almond->levels[k+1]->o_rkp1);
        MPI_Allreduce(&rho1Local,&rho1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&alpha1Local,&alpha1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      }
   
      norm_rkp1Local = innerProd(almond, mCoarse, almond->levels[k+1]->o_rkp1,almond->levels[k+1]->o_rkp1);
      MPI_Allreduce(&norm_rkp1Local,&norm_rkp1Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      norm_rkp1Global = sqrt(norm_rkp1Global);

      // rkp1 = rkp1 - (alpha1/rho1)*vkp1
      vectorAdd(almond, mCoarse, -alpha1Global/rho1Global, almond->levels[k+1]->o_vkp1, 1.0,
          almond->levels[k+1]->o_rkp1 );

      norm_rktilde_pLocal = innerProd(almond, mCoarse, almond->levels[k+1]->o_rkp1,almond->levels[k+1]->o_rkp1);
      MPI_Allreduce(&norm_rktilde_pLocal,&norm_rktilde_pGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      norm_rktilde_pGlobal = sqrt(norm_rktilde_pGlobal);

      //      almond->levels[k+1]->rhs = rkp1;
      copyVector(almond, mCoarse, almond->levels[k+1]->o_rkp1, almond->levels[k+1]->o_rhs);

      dfloat t = 0.2;

      if(norm_rktilde_pGlobal < t*norm_rkp1Global){
        //      almond->levels[k+1]->x = (alpha1/rho1)*ckp1
        vectorAdd(almond, mCoarse, alpha1Global/rho1Global, almond->levels[k+1]->o_ckp1, 0., almond->levels[k+1]->o_x);
      } else{

        device_kcycle(almond,k+1);

        //  dkp1 = almond->levels[k+1]->x;
        copyVector(almond, mCoarse, almond->levels[k+1]->o_x, almond->levels[k+1]->o_dkp1);

        // w = A*d
        axpy(almond, almond->levels[k+1]->deviceA, 1.0, almond->levels[k+1]->o_dkp1, 0.,
                almond->levels[k+1]->o_wkp1);

        dfloat gammaLocal, betaLocal, alpha2Local;
        dfloat gammaGlobal, betaGlobal, alpha2Global;

        if(almond->ktype == PCG){
          gammaLocal  = innerProd(almond, mCoarse, almond->levels[k+1]->o_dkp1, almond->levels[k+1]->o_vkp1);
          betaLocal   = innerProd(almond, mCoarse, almond->levels[k+1]->o_dkp1, almond->levels[k+1]->o_wkp1);
          alpha2Local = innerProd(almond, mCoarse, almond->levels[k+1]->o_dkp1, almond->levels[k+1]->o_rkp1);
          MPI_Allreduce(&gammaLocal,&gammaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
          MPI_Allreduce(&betaLocal,&betaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
          MPI_Allreduce(&alpha2Local,&alpha2Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        }
        if(almond->ktype == GMRES){
          gammaLocal  = innerProd(almond, mCoarse, almond->levels[k+1]->o_wkp1, almond->levels[k+1]->o_vkp1);
          betaLocal   = innerProd(almond, mCoarse, almond->levels[k+1]->o_wkp1, almond->levels[k+1]->o_wkp1);
          alpha2Local = innerProd(almond, mCoarse, almond->levels[k+1]->o_wkp1, almond->levels[k+1]->o_rkp1);
          MPI_Allreduce(&gammaLocal,&gammaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
          MPI_Allreduce(&betaLocal,&betaGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
          MPI_Allreduce(&alpha2Local,&alpha2Global,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
        }

        if(fabs(rho1Global) > (dfloat) 1e-20){

          dfloat rho2 = betaGlobal - gammaGlobal*gammaGlobal/rho1Global;

          if(fabs(rho2) > (dfloat) 1e-20){
            // almond->levels[k+1]->x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ckp1 + (alpha2/rho2)*dkp1

            dfloat a = alpha1Global/rho1Global - gammaGlobal*alpha2Global/(rho1Global*rho2);
            dfloat b = alpha2Global/rho2;

            vectorAdd(almond, mCoarse, a, almond->levels[k+1]->o_ckp1, 0.0, almond->levels[k+1]->o_x);
            vectorAdd(almond, mCoarse, b, almond->levels[k+1]->o_dkp1, 1.0, almond->levels[k+1]->o_x);
          }
        }
      }
    }
  } else{
    scaleVector(almond, mCoarse, almond->levels[k+1]->o_x, 0.);
    if (almond->coarseSolve != NULL) {
      //use direct sovler passed as input
      dfloat *xCoarse   = (dfloat*) calloc(almond->coarseTotal,sizeof(dfloat));
      dfloat *rhsCoarse = (dfloat*) calloc(almond->coarseTotal,sizeof(dfloat));

      almond->levels[k+1]->o_rhs.copyTo(rhsCoarse+almond->coarseOffset);

      almond->coarseSolve((void *) xCoarse, almond->ACoarse, (void *) rhsCoarse); 

      almond->levels[k+1]->o_x.copyFrom(xCoarse+almond->coarseOffset,mCoarse*sizeof(dfloat));  
    } else {
      smooth(almond, almond->levels[k+1], almond->levels[k+1]->o_rhs, almond->levels[k+1]->o_x, true);
    }
  }

  interpolate(almond, almond->levels[k], almond->levels[k+1]->o_x, almond->levels[k]->o_x);

  //use matrix free action if its been given
  //if ((k==0)&&almond->matFreeAX) {
  //  matFreeSmooth(almond, almond->levels[k]->o_rhs, almond->levels[k]->o_x, false);
  //} else {
    smooth(almond, almond->levels[k], almond->levels[k]->o_rhs, almond->levels[k]->o_x, false);
  //}

  occaTimerToc(almond->device,name);
}



void vcycle(almond_t *almond, int k) {

  const iint m = almond->levels[k]->Nrows;
  const iint mCoarse = almond->levels[k+1]->Nrows;

  char name[BUFSIZ];
  sprintf(name, "host vcycle level %d", k);
  occaTimerTic(almond->device,name);

  // if its not first level zero it out
  if(k>0)
    scaleVector(m, almond->levels[k]->x,  0.0);

  smooth(almond->levels[k], almond->levels[k]->rhs, almond->levels[k]->x, true);

  // res = rhs - A*x
  zeqaxpy(almond->levels[k]->A, -1.0, almond->levels[k]->x, 1.0, almond->levels[k]->rhs,
     almond->levels[k]->res);

  // restrict the residual to next level
  restrict(almond->levels[k], almond->levels[k]->res, almond->levels[k+1]->rhs);

  if(k+1 < almond->numLevels - 1){
    vcycle(almond,k+1);
  } else{
    scaleVector(mCoarse, almond->levels[k+1]->x, 0.);
    if (almond->coarseSolve != NULL) {
      //use direct sovler passed as input
      dfloat *xCoarse   = (dfloat*) calloc(almond->coarseTotal,sizeof(dfloat));
      dfloat *rhsCoarse = (dfloat*) calloc(almond->coarseTotal,sizeof(dfloat));

      for (iint n=0;n<mCoarse;n++)
        rhsCoarse[n+almond->coarseOffset] = almond->levels[k+1]->rhs[n];

      almond->coarseSolve((void *) xCoarse, almond->ACoarse, (void *) rhsCoarse); 

      for (iint n=0;n<mCoarse;n++) 
        almond->levels[k+1]->x[n] = xCoarse[n+almond->coarseOffset];
    } else {
      smooth(almond->levels[k+1], almond->levels[k+1]->rhs, almond->levels[k+1]->x,true);
    }
  }

  interpolate(almond->levels[k], almond->levels[k+1]->x, almond->levels[k]->x);
  smooth(almond->levels[k], almond->levels[k]->rhs, almond->levels[k]->x,false);

  occaTimerToc(almond->device,name);
}


void device_vcycle(almond_t *almond, int k){

#define GPU_CPU_SWITCH_SIZE 0 //TODO move this the the almond struct?

  const iint m = almond->levels[k]->Nrows;
  const iint mCoarse = almond->levels[k+1]->Nrows;

  // switch to cpu if the problem size is too small for gpu
  if(m < GPU_CPU_SWITCH_SIZE){
    almond->levels[k]->o_rhs.copyTo(&(almond->levels[k]->rhs[0]), m*sizeof(dfloat));
    vcycle(almond, k);
    almond->levels[k]->o_x.copyFrom(&(almond->levels[k]->x[0]), m*sizeof(dfloat));
    return;
  }

  char name[BUFSIZ];
  sprintf(name, "device vcycle level %d", k);
  occaTimerTic(almond->device,name);

  // if its not first level zero it out
  if(k>0)
    scaleVector(almond, m, almond->levels[k]->o_x, 0.0);

  //if ((k==0)&&almond->matFreeAX) {
  //  matFreeSmooth(almond, almond->levels[k]->o_rhs, almond->levels[k]->o_x, true);
  //  matFreeZeqAXPY(almond, -1.0, almond->levels[k]->o_x,  1.0,
  //           almond->levels[k]->o_rhs, almond->levels[k]->o_res);
  //} else {
    smooth(almond, almond->levels[k], almond->levels[k]->o_rhs, almond->levels[k]->o_x, true);
    // res = rhs - A*x
    zeqaxpy(almond, almond->levels[k]->deviceA, -1.0, almond->levels[k]->o_x,  1.0,
             almond->levels[k]->o_rhs, almond->levels[k]->o_res);
  //}

  // restrict the residual to next level
  restrict(almond, almond->levels[k], almond->levels[k]->o_res, almond->levels[k+1]->o_rhs);


  if(k+1 < almond->numLevels - 1){
    device_vcycle(almond, k+1);
  }else{
    scaleVector(almond, mCoarse, almond->levels[k+1]->o_x, 0.);
    if (almond->coarseSolve != NULL) {
      //use direct sovler passed as input
      dfloat *xCoarse   = (dfloat*) calloc(almond->coarseTotal,sizeof(dfloat));
      dfloat *rhsCoarse = (dfloat*) calloc(almond->coarseTotal,sizeof(dfloat));

      almond->levels[k+1]->o_rhs.copyTo(rhsCoarse+almond->coarseOffset);

      almond->coarseSolve((void *) xCoarse, almond->ACoarse, (void *) rhsCoarse); 

      almond->levels[k+1]->o_x.copyFrom(xCoarse+almond->coarseOffset,mCoarse*sizeof(dfloat));  
    } else {
      smooth(almond, almond->levels[k+1], almond->levels[k+1]->o_rhs, almond->levels[k+1]->o_x, true);
    }
  }

  interpolate(almond, almond->levels[k], almond->levels[k+1]->o_x, almond->levels[k]->o_x);

  //if ((k==0)&&almond->matFreeAX) {
  //  matFreeSmooth(almond, almond->levels[k]->o_rhs, almond->levels[k]->o_x, false);
  //} else {
    smooth(almond, almond->levels[k], almond->levels[k]->o_rhs, almond->levels[k]->o_x,false);
  //}

  occaTimerToc(almond->device,name);
}

#include "parAlmond.h"

void kcycle(parAlmond_t *parAlmond, int k){

  iint m = parAlmond->levels[k]->Nrows;
  iint n = parAlmond->levels[k]->Ncols;
  iint mCoarse = parAlmond->levels[k+1]->Nrows;
  iint nCoarse = parAlmond->levels[k+1]->Ncols;

  char name[BUFSIZ];
  sprintf(name, "host kcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  setVector(m, parAlmond->levels[k]->x, 0.0);

  smooth(parAlmond->levels[k], parAlmond->levels[k]->rhs, parAlmond->levels[k]->x, true);

    // res = - A*x + rhs (i.e., rhs - A*x)
  zeqaxpy(parAlmond->levels[k]->A, -1.0, parAlmond->levels[k]->x, 1.0,
          parAlmond->levels[k]->rhs, parAlmond->levels[k]->res);

  // restrict the residual to next level
  restrict(parAlmond->levels[k], parAlmond->levels[k]->res, parAlmond->levels[k+1]->rhs);

  if(k+1 < parAlmond->numLevels - 1){
    dfloat *ckp1 = parAlmond->levels[k+1]->ckp1;
    dfloat *vkp1 = parAlmond->levels[k+1]->vkp1;
    dfloat *wkp1 = parAlmond->levels[k+1]->wkp1;
    dfloat *dkp1 = parAlmond->levels[k+1]->x;
    dfloat *rkp1 = parAlmond->levels[k+1]->rhs;

    // first inner krylov iteration
    kcycle(parAlmond, k+1);

    //ckp1 = x
    memcpy(ckp1,parAlmond->levels[k+1]->x,mCoarse*sizeof(dfloat));

    // v = A*c
    axpy(parAlmond->levels[k+1]->A, 1.0, ckp1, 0.0, vkp1);

    dfloat rhoLocal[3], rhoGlobal[3];

    dfloat rho1, alpha1, norm_rkp1;
    dfloat norm_rktilde_p, norm_rktilde_pGlobal;

    if(parAlmond->ktype == PCG)
      kcycleCombinedOp1(mCoarse, rhoLocal, ckp1, rkp1, vkp1);

    if(parAlmond->ktype == GMRES)
      kcycleCombinedOp1(mCoarse, rhoLocal, vkp1, rkp1, vkp1);

    MPI_Allreduce(rhoLocal,rhoGlobal,3,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

    alpha1 = rhoGlobal[0];
    rho1   = rhoGlobal[1];
    norm_rkp1 = sqrt(rhoGlobal[2]);

    // rkp1 = rkp1 - (alpha1/rho1)*vkp1
    norm_rktilde_p = vectorAddInnerProd(mCoarse, -alpha1/rho1, vkp1, 1.0, rkp1);
    MPI_Allreduce(&norm_rktilde_p,&norm_rktilde_pGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
    norm_rktilde_pGlobal = sqrt(norm_rktilde_pGlobal);

    dfloat t = 0.2;

    if(norm_rktilde_pGlobal < t*norm_rkp1){
      // x = (alpha1/rho1)*x
      scaleVector(mCoarse, parAlmond->levels[k+1]->x, alpha1/rho1);
    } else{

      kcycle(parAlmond, k+1);

      // w = A*d
      axpy(parAlmond->levels[k+1]->A, 1.0, dkp1, 0., wkp1);

      dfloat gamma, beta, alpha2;

      if(parAlmond->ktype == PCG)
        kcycleCombinedOp2(mCoarse,rhoLocal,dkp1,vkp1,wkp1,rkp1);

      if(parAlmond->ktype == GMRES)
        kcycleCombinedOp2(mCoarse,rhoLocal,wkp1,vkp1,wkp1,rkp1);

      MPI_Allreduce(rhoLocal,rhoGlobal,3,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

      gamma  = rhoGlobal[0];
      beta   = rhoGlobal[1];
      alpha2 = rhoGlobal[2];

      if(fabs(rho1) > (dfloat) 1e-20){

        dfloat rho2 = beta - gamma*gamma/rho1;

        if(fabs(rho2) > (dfloat) 1e-20){
          // parAlmond->levels[k+1]->x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ckp1 + (alpha2/rho2)*dkp1
          dfloat a = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
          dfloat b = alpha2/rho2;

          vectorAdd(mCoarse, a, ckp1, b, parAlmond->levels[k+1]->x);
        }
      }
    }
  } else {
    if (parAlmond->ExactSolve != NULL) {
      //use coarse sovler
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      for (iint n=0;n<mCoarse;n++)
        parAlmond->rhsCoarse[n+parAlmond->coarseOffset] = parAlmond->levels[k+1]->rhs[n];

      xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);

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
  setVector(parAlmond, m, parAlmond->levels[k]->o_x, 0.0);

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

  if(mCoarse)
    parAlmond->levels[k+1]->o_rhs.copyTo(parAlmond->levels[k+1]->rhs);

  if(k+1 < parAlmond->numLevels - 1){
    if(k>2) {
      device_vcycle(parAlmond,k+1);
      //device_kcycle(parAlmond, k+1);
    } else{
      // first inner krylov iteration
      device_kcycle(parAlmond,k+1);

      //ckp1 = parAlmond->levels[k+1]->x;
      parAlmond->levels[k+1]->o_ckp1.copyFrom(parAlmond->levels[k+1]->o_x);

      // v = A*c
      axpy(parAlmond, parAlmond->levels[k+1]->deviceA, 1.0, parAlmond->levels[k+1]->o_ckp1, 0.0,
              parAlmond->levels[k+1]->o_vkp1);

      dfloat rhoLocal[3], rhoGlobal[3];

      dfloat rho1, alpha1, norm_rkp1;
      dfloat norm_rktilde_pLocal, norm_rktilde_pGlobal;

      if(parAlmond->ktype == PCG)
        kcycleCombinedOp1(parAlmond, mCoarse, rhoLocal,
                          parAlmond->levels[k+1]->o_ckp1,
                          parAlmond->levels[k+1]->o_rhs,
                          parAlmond->levels[k+1]->o_vkp1);

      if(parAlmond->ktype == GMRES)
        kcycleCombinedOp1(parAlmond, mCoarse, rhoLocal,
                          parAlmond->levels[k+1]->o_vkp1,
                          parAlmond->levels[k+1]->o_rhs,
                          parAlmond->levels[k+1]->o_vkp1);

      MPI_Allreduce(rhoLocal,rhoGlobal,3,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

      alpha1 = rhoGlobal[0];
      rho1   = rhoGlobal[1];
      norm_rkp1 = sqrt(rhoGlobal[2]);

      // rkp1 = rkp1 - (alpha1/rho1)*vkp1
      norm_rktilde_pLocal = vectorAddInnerProd(parAlmond, mCoarse, -alpha1/rho1,
                                                parAlmond->levels[k+1]->o_vkp1, 1.0,
                                                parAlmond->levels[k+1]->o_rhs);
      MPI_Allreduce(&norm_rktilde_pLocal,&norm_rktilde_pGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      norm_rktilde_pGlobal = sqrt(norm_rktilde_pGlobal);

      dfloat t = 0.2;
      if(norm_rktilde_pGlobal < t*norm_rkp1){
        //      parAlmond->levels[k+1]->x = (alpha1/rho1)*x
        scaleVector(parAlmond,mCoarse, parAlmond->levels[k+1]->o_x, alpha1/rho1);
      } else{
        device_kcycle(parAlmond,k+1);

        // w = A*d
        axpy(parAlmond, parAlmond->levels[k+1]->deviceA, 1.0, parAlmond->levels[k+1]->o_x, 0.,
                parAlmond->levels[k+1]->o_wkp1);

        dfloat gamma, beta, alpha2;

        if(parAlmond->ktype == PCG)
          kcycleCombinedOp2(parAlmond,mCoarse,rhoLocal,
                            parAlmond->levels[k+1]->o_x,
                            parAlmond->levels[k+1]->o_vkp1,
                            parAlmond->levels[k+1]->o_wkp1,
                            parAlmond->levels[k+1]->o_rhs);

        if(parAlmond->ktype == GMRES)
          kcycleCombinedOp2(parAlmond,mCoarse,rhoLocal,
                            parAlmond->levels[k+1]->o_wkp1,
                            parAlmond->levels[k+1]->o_vkp1,
                            parAlmond->levels[k+1]->o_wkp1,
                            parAlmond->levels[k+1]->o_rhs);

        MPI_Allreduce(rhoLocal,rhoGlobal,3,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

        gamma  = rhoGlobal[0];
        beta   = rhoGlobal[1];
        alpha2 = rhoGlobal[2];

        if(fabs(rho1) > (dfloat) 1e-20){

          dfloat rho2 = beta - gamma*gamma/rho1;

          if(fabs(rho2) > (dfloat) 1e-20){
            // parAlmond->levels[k+1]->x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ckp1 + (alpha2/rho2)*dkp1
            dfloat a = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
            dfloat b = alpha2/rho2;

            vectorAdd(parAlmond, mCoarse, a, parAlmond->levels[k+1]->o_ckp1,
                                          b, parAlmond->levels[k+1]->o_x);
          }
        }
      }
    }
  } else{
    if (parAlmond->ExactSolve != NULL) {
      //use coarse sovler
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      parAlmond->levels[k+1]->o_rhs.copyTo(parAlmond->rhsCoarse+parAlmond->coarseOffset);
      xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);
      parAlmond->levels[k+1]->o_x.copyFrom(parAlmond->xCoarse+parAlmond->coarseOffset,mCoarse*sizeof(dfloat));
    } else {
      setVector(parAlmond, mCoarse, parAlmond->levels[k+1]->o_x, 0.);
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
  setVector(m, parAlmond->levels[k]->x,  0.0);

  smooth(parAlmond->levels[k], parAlmond->levels[k]->rhs, parAlmond->levels[k]->x, true);

  // res = rhs - A*x
  zeqaxpy(parAlmond->levels[k]->A, -1.0, parAlmond->levels[k]->x, 1.0, parAlmond->levels[k]->rhs,
     parAlmond->levels[k]->res);

  // restrict the residual to next level
  restrict(parAlmond->levels[k], parAlmond->levels[k]->res, parAlmond->levels[k+1]->rhs);

  if(k+1 < parAlmond->numLevels - 1){
    vcycle(parAlmond,k+1);
  } else{
    if (parAlmond->ExactSolve != NULL) {
      //use coarse sovler
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      for (iint n=0;n<mCoarse;n++)
        parAlmond->rhsCoarse[n+parAlmond->coarseOffset] = parAlmond->levels[k+1]->rhs[n];

      xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);

      for (iint n=0;n<mCoarse;n++)
        parAlmond->levels[k+1]->x[n] = parAlmond->xCoarse[n+parAlmond->coarseOffset];
    } else {
      // scaleVector(mCoarse, parAlmond->levels[k+1]->x, 0.);
      smooth(parAlmond->levels[k+1], parAlmond->levels[k+1]->rhs, parAlmond->levels[k+1]->x,true);
    }
  }

  interpolate(parAlmond->levels[k], parAlmond->levels[k+1]->x, parAlmond->levels[k]->x);
  smooth(parAlmond->levels[k], parAlmond->levels[k]->rhs, parAlmond->levels[k]->x,false);

  occaTimerToc(parAlmond->device,name);
}


void device_vcycle(parAlmond_t *parAlmond, int k){

#define GPU_CPU_SWITCH_SIZE 0 //TODO move this the the almond struct?

  const iint m = parAlmond->levels[k]->Nrows;
  const iint mCoarse = parAlmond->levels[k+1]->Nrows;

  // switch to cpu if the problem size is too small for gpu
  if(m < GPU_CPU_SWITCH_SIZE){
    parAlmond->levels[k]->o_rhs.copyTo(parAlmond->levels[k]->rhs, m*sizeof(dfloat));
    vcycle(parAlmond, k);
    parAlmond->levels[k]->o_x.copyFrom(parAlmond->levels[k]->x, m*sizeof(dfloat));
    return;
  }

  char name[BUFSIZ];
  sprintf(name, "device vcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  setVector(parAlmond, m, parAlmond->levels[k]->o_x, 0.0);

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
    if (parAlmond->ExactSolve != NULL) {
      //use coarse sovler
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      parAlmond->levels[k+1]->o_rhs.copyTo(parAlmond->rhsCoarse+parAlmond->coarseOffset);
      xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);
      parAlmond->levels[k+1]->o_x.copyFrom(parAlmond->xCoarse+parAlmond->coarseOffset,mCoarse*sizeof(dfloat));
    } else {
      //      setVector(parAlmond, mCoarse, parAlmond->levels[k+1]->o_x, 0.);
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

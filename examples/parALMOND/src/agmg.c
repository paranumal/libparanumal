#include "agmg.h"

void kcycle(parAlmond_t *parAlmond, int k){

  agmgLevel **levels = parAlmond->levels;

  iint m = levels[k]->Nrows;
  iint n = levels[k]->Ncols;
  iint mCoarse = levels[k+1]->Nrows;
  iint nCoarse = levels[k+1]->Ncols;

  char name[BUFSIZ];
  sprintf(name, "host kcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  //setVector(m, levels[k]->x, 0.0);

  levels[k]->smooth(levels[k]->smoothArgs, levels[k]->rhs, levels[k]->x, true);

  // res = r - A*x
  levels[k]->Ax(levels[k]->AxArgs,levels[k]->x,levels[k]->res);
  vectorAdd(m, 1.0, levels[k]->rhs, -1.0, levels[k]->res);

  // coarsen the residual to next level
  levels[k+1]->coarsen(levels[k+1]->coarsenArgs, levels[k]->res, levels[k+1]->rhs);

  if(k+1 < parAlmond->numLevels - 1){
    if(k>2) {
      vcycle(parAlmond,k+1);
      //kcycle(parAlmond, k+1);
    } else{
      dfloat *ckp1 = levels[k+1]->ckp1;
      dfloat *vkp1 = levels[k+1]->vkp1;
      dfloat *wkp1 = levels[k+1]->wkp1;
      dfloat *dkp1 = levels[k+1]->x;
      dfloat *rkp1 = levels[k+1]->rhs;

      // first inner krylov iteration
      kcycle(parAlmond, k+1);

      //ckp1 = x
      memcpy(ckp1,levels[k+1]->x,mCoarse*sizeof(dfloat));

      // v = A*c
      levels[k+1]->Ax(levels[k+1]->AxArgs,ckp1,vkp1);

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
        scaleVector(mCoarse, levels[k+1]->x, alpha1/rho1);
      } else{

        kcycle(parAlmond, k+1);

        // w = A*d
        levels[k+1]->Ax(levels[k+1]->AxArgs,dkp1,wkp1);

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
            // levels[k+1]->x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ckp1 + (alpha2/rho2)*dkp1
            dfloat a = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
            dfloat b = alpha2/rho2;

            vectorAdd(mCoarse, a, ckp1, b, levels[k+1]->x);
          }
        }
      }
    }
  } else {
    if (parAlmond->ExactSolve != NULL) {
      //use coarse sovler
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      for (iint n=0;n<mCoarse;n++)
        parAlmond->rhsCoarse[n+parAlmond->coarseOffset] = levels[k+1]->rhs[n];

      xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);

      for (iint n=0;n<mCoarse;n++)
        levels[k+1]->x[n] = parAlmond->xCoarse[n+parAlmond->coarseOffset];
    } else {
      levels[k+1]->smooth(levels[k+1]->smoothArgs, levels[k+1]->rhs, levels[k+1]->x, true);
    }
  }

  levels[k+1]->prolongate(levels[k+1]->prolongateArgs, levels[k+1]->x, levels[k]->x);

  levels[k]->smooth(levels[k]->smoothArgs, levels[k]->rhs, levels[k]->x, false);

  occaTimerToc(parAlmond->device,name);
}


void device_kcycle(parAlmond_t *parAlmond, int k){

  agmgLevel **levels = parAlmond->levels;

  iint m = levels[k]->Nrows;
  iint n = levels[k]->Ncols;
  iint mCoarse = levels[k+1]->Nrows;
  iint nCoarse = levels[k+1]->Ncols;

  char name[BUFSIZ];
  sprintf(name, "device kcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  //setVector(parAlmond, m, levels[k]->o_x, 0.0);

  levels[k]->device_smooth(levels[k]->smoothArgs, levels[k]->o_rhs, levels[k]->o_x, true);

  // res = rhs - A*x
  levels[k]->device_Ax(levels[k]->AxArgs,levels[k]->o_x,levels[k]->o_res);
  vectorAdd(parAlmond, m, 1.0, levels[k]->o_rhs, -1.0, levels[k]->o_res);

  // coarsen the residual to next level
  levels[k+1]->device_coarsen(levels[k+1]->coarsenArgs, levels[k]->o_res, levels[k+1]->o_rhs);


  if(k+1 < parAlmond->numLevels - 1){
    if(k>2) {
      device_vcycle(parAlmond,k+1);
      //device_kcycle(parAlmond, k+1);
    } else{
      // first inner krylov iteration
      device_kcycle(parAlmond,k+1);

      //ckp1 = levels[k+1]->x;
      if (mCoarse)
        levels[k+1]->o_ckp1.copyFrom(levels[k+1]->o_x);

      // v = A*c
      levels[k+1]->device_Ax(levels[k+1]->AxArgs,levels[k+1]->o_ckp1,levels[k+1]->o_vkp1);

      dfloat rhoLocal[3], rhoGlobal[3];
      dfloat rho1, alpha1, norm_rkp1;
      dfloat norm_rktilde_pLocal, norm_rktilde_pGlobal;

      if(parAlmond->ktype == PCG)
        kcycleCombinedOp1(parAlmond, mCoarse, rhoLocal,
                          levels[k+1]->o_ckp1,
                          levels[k+1]->o_rhs,
                          levels[k+1]->o_vkp1);

      if(parAlmond->ktype == GMRES)
        kcycleCombinedOp1(parAlmond, mCoarse, rhoLocal,
                          levels[k+1]->o_vkp1,
                          levels[k+1]->o_rhs,
                          levels[k+1]->o_vkp1);

      MPI_Allreduce(rhoLocal,rhoGlobal,3,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

      alpha1 = rhoGlobal[0];
      rho1   = rhoGlobal[1];
      norm_rkp1 = sqrt(rhoGlobal[2]);

      // rkp1 = rkp1 - (alpha1/rho1)*vkp1
      norm_rktilde_pLocal = vectorAddInnerProd(parAlmond, mCoarse, -alpha1/rho1,
                                                levels[k+1]->o_vkp1, 1.0,
                                                levels[k+1]->o_rhs);
      MPI_Allreduce(&norm_rktilde_pLocal,&norm_rktilde_pGlobal,1,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);
      norm_rktilde_pGlobal = sqrt(norm_rktilde_pGlobal);

      dfloat t = 0.2;
      if(norm_rktilde_pGlobal < t*norm_rkp1){
        //      levels[k+1]->x = (alpha1/rho1)*x
        scaleVector(parAlmond,mCoarse, levels[k+1]->o_x, alpha1/rho1);
      } else{
        device_kcycle(parAlmond,k+1);

        // w = A*x
        levels[k+1]->device_Ax(levels[k+1]->AxArgs,levels[k+1]->o_x,levels[k+1]->o_wkp1);

        dfloat gamma, beta, alpha2;

        if(parAlmond->ktype == PCG)
          kcycleCombinedOp2(parAlmond,mCoarse,rhoLocal,
                            levels[k+1]->o_x,
                            levels[k+1]->o_vkp1,
                            levels[k+1]->o_wkp1,
                            levels[k+1]->o_rhs);

        if(parAlmond->ktype == GMRES)
          kcycleCombinedOp2(parAlmond,mCoarse,rhoLocal,
                            levels[k+1]->o_wkp1,
                            levels[k+1]->o_vkp1,
                            levels[k+1]->o_wkp1,
                            levels[k+1]->o_rhs);

        MPI_Allreduce(rhoLocal,rhoGlobal,3,MPI_DFLOAT,MPI_SUM,MPI_COMM_WORLD);

        gamma  = rhoGlobal[0];
        beta   = rhoGlobal[1];
        alpha2 = rhoGlobal[2];

        if(fabs(rho1) > (dfloat) 1e-20){

          dfloat rho2 = beta - gamma*gamma/rho1;

          if(fabs(rho2) > (dfloat) 1e-20){
            // levels[k+1]->x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ckp1 + (alpha2/rho2)*dkp1
            dfloat a = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
            dfloat b = alpha2/rho2;

            vectorAdd(parAlmond, mCoarse, a, levels[k+1]->o_ckp1,
                                          b, levels[k+1]->o_x);
          }
        }
      }
    }
  } else{
    if (parAlmond->ExactSolve != NULL) {
      //use coarse sovler
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      levels[k+1]->o_rhs.copyTo(parAlmond->rhsCoarse+parAlmond->coarseOffset);
      xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);
      levels[k+1]->o_x.copyFrom(parAlmond->xCoarse+parAlmond->coarseOffset,mCoarse*sizeof(dfloat));
    } else {
      levels[k+1]->device_smooth(levels[k+1]->smoothArgs, levels[k+1]->o_rhs, levels[k+1]->o_x, true);
    }
  }

  levels[k+1]->device_prolongate(levels[k+1]->prolongateArgs, levels[k+1]->o_x, levels[k]->o_x);

  levels[k]->device_smooth(levels[k]->smoothArgs, levels[k]->o_rhs, levels[k]->o_x, false);

  occaTimerToc(parAlmond->device,name);
}



void vcycle(parAlmond_t *parAlmond, int k) {

  agmgLevel **levels = parAlmond->levels;

  const iint m = levels[k]->Nrows;
  const iint mCoarse = levels[k+1]->Nrows;

  char name[BUFSIZ];
  sprintf(name, "host vcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  //setVector(m, levels[k]->x,  0.0);

  levels[k]->smooth(levels[k]->smoothArgs, levels[k]->rhs, levels[k]->x, true);

  // res = rhs - A*x
  levels[k]->Ax(levels[k]->AxArgs,levels[k]->x,levels[k]->res);
  vectorAdd(m, 1.0, levels[k]->rhs, -1.0, levels[k]->res);

  // coarsen the residual to next level
  levels[k+1]->coarsen(levels[k+1]->coarsenArgs, levels[k]->res, levels[k+1]->rhs);

  if(k+1 < parAlmond->numLevels - 1){
    vcycle(parAlmond,k+1);
  } else{
    if (parAlmond->ExactSolve != NULL) {
      //use coarse sovler
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      for (iint n=0;n<mCoarse;n++)
        parAlmond->rhsCoarse[n+parAlmond->coarseOffset] = levels[k+1]->rhs[n];

      xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);

      for (iint n=0;n<mCoarse;n++)
        levels[k+1]->x[n] = parAlmond->xCoarse[n+parAlmond->coarseOffset];
    } else {
      levels[k+1]->smooth(levels[k+1]->smoothArgs, levels[k+1]->rhs, levels[k+1]->x,true);
    }
  }

  levels[k+1]->prolongate(levels[k+1]->prolongateArgs, levels[k+1]->x, levels[k]->x);
  levels[k]->smooth(levels[k]->smoothArgs, levels[k]->rhs, levels[k]->x,false);

  occaTimerToc(parAlmond->device,name);
}


void device_vcycle(parAlmond_t *parAlmond, int k){

  agmgLevel **levels = parAlmond->levels;

  const iint m = levels[k]->Nrows;
  const iint mCoarse = levels[k+1]->Nrows;

  // switch to cpu if the problem size is too small for gpu
  if(m < GPU_CPU_SWITCH_SIZE){
    levels[k]->o_rhs.copyTo(levels[k]->rhs, m*sizeof(dfloat));
    vcycle(parAlmond, k);
    levels[k]->o_x.copyFrom(levels[k]->x, m*sizeof(dfloat));
    return;
  }

  char name[BUFSIZ];
  sprintf(name, "device vcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  //setVector(parAlmond, m, levels[k]->o_x, 0.0);

  levels[k]->device_smooth(levels[k]->smoothArgs, levels[k]->o_rhs, levels[k]->o_x, true);

  // res = rhs - A*x
  levels[k]->device_Ax(levels[k]->AxArgs,levels[k]->o_x,levels[k]->o_res);
  vectorAdd(parAlmond, m, 1.0, levels[k]->o_rhs, -1.0, levels[k]->o_res);

  // coarsen the residual to next level
  levels[k+1]->device_coarsen(levels[k+1]->coarsenArgs, levels[k]->o_res, levels[k+1]->o_rhs);

  if(k+1 < parAlmond->numLevels - 1){
    device_vcycle(parAlmond, k+1);
  }else{
    if (parAlmond->ExactSolve != NULL) {
      //use coarse sovler
      for (iint n=0;n<parAlmond->coarseTotal;n++)
        parAlmond->rhsCoarse[n] =0.;

      levels[k+1]->o_rhs.copyTo(parAlmond->rhsCoarse+parAlmond->coarseOffset);
      xxtSolve(parAlmond->xCoarse, parAlmond->ExactSolve, parAlmond->rhsCoarse);
      levels[k+1]->o_x.copyFrom(parAlmond->xCoarse+parAlmond->coarseOffset,mCoarse*sizeof(dfloat));
    } else {
      levels[k+1]->device_smooth(levels[k+1]->smoothArgs, levels[k+1]->o_rhs, levels[k+1]->o_x, true);
    }
  }

  levels[k+1]->device_prolongate(levels[k+1]->prolongateArgs, levels[k+1]->o_x, levels[k]->o_x);

  levels[k]->device_smooth(levels[k]->smoothArgs, levels[k]->o_rhs, levels[k]->o_x,false);

  occaTimerToc(parAlmond->device,name);
}

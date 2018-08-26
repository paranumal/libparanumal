/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "agmg.h"

namespace agmg{

  int rank;
  int size;
  MPI_Comm comm;

}

void kcycle(parAlmond_t *parAlmond, int k){

  agmgLevel **levels = parAlmond->levels;

  dlong m = levels[k]->Nrows;
  // dlong n = levels[k]->Ncols;

  //check for base level
  if(k==parAlmond->numLevels-1) {
    if (parAlmond->invCoarseA != NULL) {
      //use exact sovler
      exactCoarseSolve(parAlmond, m, levels[k]->rhs, levels[k]->x);
    } else {
      levels[k]->smooth(levels[k]->smoothArgs, levels[k]->rhs, levels[k]->x, true);
    }
    return;
  }

  char name[BUFSIZ];
  sprintf(name, "host kcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  dlong mCoarse = levels[k+1]->Nrows;
  // dlong nCoarse = levels[k+1]->Ncols;

  // zero out x
  //setVector(m, levels[k]->x, 0.0);

  levels[k]->smooth(levels[k]->smoothArgs, levels[k]->rhs, levels[k]->x, true);

  // res = r - A*x
  levels[k]->Ax(levels[k]->AxArgs,levels[k]->x,levels[k]->res);
  vectorAdd(m, 1.0, levels[k]->rhs, -1.0, levels[k]->res);

  // coarsen the residual to next level, checking if the residual needs to be gathered after
  if (levels[k+1]->gatherLevel==true) {
    levels[k+1]->coarsen(levels[k+1]->coarsenArgs, levels[k]->res, levels[k+1]->Srhs);
    levels[k+1]->gather (levels[k+1]->gatherArgs,  levels[k+1]->Srhs, levels[k+1]->rhs);
  } else {
    levels[k+1]->coarsen(levels[k+1]->coarsenArgs, levels[k]->res, levels[k+1]->rhs);
  }

  if(k>2) {
    vcycle(parAlmond,k+1);
    //kcycle(parAlmond, k+1);
  } else{
    dfloat *ckp1 = levels[k+1]->ckp1;
    dfloat *vkp1 = levels[k+1]->vkp1;
    dfloat *wkp1 = levels[k+1]->wkp1;
    dfloat *dkp1 = levels[k+1]->x;
    dfloat *rkp1 = levels[k+1]->rhs;
    dfloat *w = levels[k+1]->weight;
    bool weighted = levels[k+1]->weightedInnerProds;

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
      kcycleCombinedOp1(mCoarse, rhoLocal, ckp1, rkp1, vkp1, w, weighted);

    if(parAlmond->ktype == GMRES)
      kcycleCombinedOp1(mCoarse, rhoLocal, vkp1, rkp1, vkp1, w, weighted);

    MPI_Allreduce(rhoLocal,rhoGlobal,3,MPI_DFLOAT,MPI_SUM,agmg::comm);

    alpha1 = rhoGlobal[0];
    rho1   = rhoGlobal[1];
    norm_rkp1 = sqrt(rhoGlobal[2]);

    // rkp1 = rkp1 - (alpha1/rho1)*vkp1
    norm_rktilde_p = vectorAddInnerProd(mCoarse, -alpha1/rho1, vkp1, 1.0, rkp1, w, weighted);
    MPI_Allreduce(&norm_rktilde_p,&norm_rktilde_pGlobal,1,MPI_DFLOAT,MPI_SUM,agmg::comm);
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
        kcycleCombinedOp2(mCoarse,rhoLocal,dkp1,vkp1,wkp1,rkp1, w, weighted);

      if(parAlmond->ktype == GMRES)
        kcycleCombinedOp2(mCoarse,rhoLocal,wkp1,vkp1,wkp1,rkp1, w, weighted);

      MPI_Allreduce(rhoLocal,rhoGlobal,3,MPI_DFLOAT,MPI_SUM,agmg::comm);

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

  if (levels[k+1]->gatherLevel==true) {
    levels[k+1]->scatter(levels[k+1]->scatterArgs,  levels[k+1]->x, levels[k+1]->Sx);
    levels[k+1]->prolongate(levels[k+1]->prolongateArgs, levels[k+1]->Sx, levels[k]->x);
  } else {
    levels[k+1]->prolongate(levels[k+1]->prolongateArgs, levels[k+1]->x, levels[k]->x);
  }

  levels[k]->smooth(levels[k]->smoothArgs, levels[k]->rhs, levels[k]->x, false);

  occaTimerToc(parAlmond->device,name);
}


void device_kcycle(parAlmond_t *parAlmond, int k){

  agmgLevel **levels = parAlmond->levels;

  dlong m = levels[k]->Nrows;
  // dlong n = levels[k]->Ncols;

  if(m < GPU_CPU_SWITCH_SIZE){
    levels[k]->o_rhs.copyTo(levels[k]->rhs, m*sizeof(dfloat));
    kcycle(parAlmond, k);
    levels[k]->o_x.copyFrom(levels[k]->x, m*sizeof(dfloat));
    return;
  }

  //check for base level
  if(k==parAlmond->numLevels-1) {
    if (parAlmond->invCoarseA != NULL) {
      //use exact sovler
      device_exactCoarseSolve(parAlmond, m, levels[k]->o_rhs, levels[k]->o_x);
    } else {
      levels[k]->device_smooth(levels[k]->smoothArgs, levels[k]->o_rhs, levels[k]->o_x, true);
    }
    return;
  }

  dlong mCoarse = levels[k+1]->Nrows;
  // dlong nCoarse = levels[k+1]->Ncols;

  char name[BUFSIZ];
  sprintf(name, "device kcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // zero out x
  //setVector(parAlmond, m, levels[k]->o_x, 0.0);

  levels[k]->device_smooth(levels[k]->smoothArgs, levels[k]->o_rhs, levels[k]->o_x, true);

  // res = rhs - A*x
  levels[k]->device_Ax(levels[k]->AxArgs,levels[k]->o_x,levels[k]->o_res);
  vectorAdd(parAlmond, m, 1.0, levels[k]->o_rhs, -1.0, levels[k]->o_res);

  // coarsen the residual to next level, checking if the residual needs to be gathered after
  if (levels[k+1]->gatherLevel==true) {
    levels[k+1]->device_coarsen(levels[k+1]->coarsenArgs, levels[k]->o_res, levels[k+1]->o_Srhs);
    levels[k+1]->device_gather (levels[k+1]->gatherArgs,  levels[k+1]->o_Srhs, levels[k+1]->o_rhs);
  } else {
    levels[k+1]->device_coarsen(levels[k+1]->coarsenArgs, levels[k]->o_res, levels[k+1]->o_rhs);
  }

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

    // kcycleCombinedOp1(parAlmond,N,aDotbc,a,b,c,w,bool) 
    //    returns aDotbc[0] = a.b, aDotbc[1] = a.c, aDotbc[2] = b.b
    //       or aDotbc[0] = w.a.b, aDotbc[1] = w.a.c, aDotbc[2] = w.b.b
    if(parAlmond->ktype == PCG)
      kcycleCombinedOp1(parAlmond, mCoarse, rhoLocal,
                        levels[k+1]->o_ckp1,
                        levels[k+1]->o_rhs,
                        levels[k+1]->o_vkp1,
                        levels[k+1]->o_weight,
                        levels[k+1]->weightedInnerProds);

    if(parAlmond->ktype == GMRES)
      kcycleCombinedOp1(parAlmond, mCoarse, rhoLocal,
                        levels[k+1]->o_vkp1,
                        levels[k+1]->o_rhs,
                        levels[k+1]->o_vkp1,
                        levels[k+1]->o_weight,
                        levels[k+1]->weightedInnerProds);

    MPI_Allreduce(rhoLocal,rhoGlobal,3,MPI_DFLOAT,MPI_SUM,agmg::comm);

    alpha1 = rhoGlobal[0];
    rho1   = rhoGlobal[1];
    norm_rkp1 = sqrt(rhoGlobal[2]);

    // rkp1 = rkp1 - (alpha1/rho1)*vkp1
    norm_rktilde_pLocal = vectorAddInnerProd(parAlmond, mCoarse, -alpha1/rho1,
                                              levels[k+1]->o_vkp1, 1.0,
                                              levels[k+1]->o_rhs,
                                              levels[k+1]->o_weight,
                                              levels[k+1]->weightedInnerProds);
    MPI_Allreduce(&norm_rktilde_pLocal,&norm_rktilde_pGlobal,1,MPI_DFLOAT,MPI_SUM,agmg::comm);
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

      // kcycleCombinedOp2(parAlmond,N,aDotbc,a,b,c,d,w,bool) 
      //   returns aDotbcd[0] = a.b, aDotbcd[1] = a.c, aDotbcd[2] = a.d,
      //      or aDotbcd[0] = w.a.b, aDotbcd[1] = w.a.c, aDotbcd[2] = w.a.d,
      if(parAlmond->ktype == PCG)
        kcycleCombinedOp2(parAlmond,mCoarse,rhoLocal,
                          levels[k+1]->o_x,
                          levels[k+1]->o_vkp1,
                          levels[k+1]->o_wkp1,
                          levels[k+1]->o_rhs,
                          levels[k+1]->o_weight,
                          levels[k+1]->weightedInnerProds);

      if(parAlmond->ktype == GMRES)
        kcycleCombinedOp2(parAlmond,mCoarse,rhoLocal,
                          levels[k+1]->o_wkp1,
                          levels[k+1]->o_vkp1,
                          levels[k+1]->o_wkp1,
                          levels[k+1]->o_rhs,
                          levels[k+1]->o_weight,
                          levels[k+1]->weightedInnerProds);

      MPI_Allreduce(rhoLocal,rhoGlobal,3,MPI_DFLOAT,MPI_SUM,agmg::comm);

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

  if (levels[k+1]->gatherLevel==true) {
    levels[k+1]->device_scatter   (levels[k+1]->scatterArgs,  levels[k+1]->o_x, levels[k+1]->o_Sx);
    levels[k+1]->device_prolongate(levels[k+1]->prolongateArgs, levels[k+1]->o_Sx, levels[k]->o_x);
  } else {
    levels[k+1]->device_prolongate(levels[k+1]->prolongateArgs, levels[k+1]->o_x, levels[k]->o_x);
  }

  levels[k]->device_smooth(levels[k]->smoothArgs, levels[k]->o_rhs, levels[k]->o_x, false);

  occaTimerToc(parAlmond->device,name);
}



void vcycle(parAlmond_t *parAlmond, int k) {

  agmgLevel **levels = parAlmond->levels;

  const dlong m = levels[k]->Nrows;

  //check for base level
  if(k==parAlmond->numLevels-1) {
    if (parAlmond->invCoarseA != NULL) {
      //use exact sovler
      exactCoarseSolve(parAlmond, m, levels[k]->rhs, levels[k]->x);
    } else {
      levels[k]->smooth(levels[k]->smoothArgs, levels[k]->rhs, levels[k]->x, true);
    }
    return;
  }

  char name[BUFSIZ];
  sprintf(name, "host vcycle level %d", k);
  occaTimerTic(parAlmond->device,name);

  // const int mCoarse = levels[k+1]->Nrows;

  // zero out x
  //setVector(m, levels[k]->x,  0.0);

  levels[k]->smooth(levels[k]->smoothArgs, levels[k]->rhs, levels[k]->x, true);

  // res = rhs - A*x
  levels[k]->Ax(levels[k]->AxArgs,levels[k]->x,levels[k]->res);
  vectorAdd(m, 1.0, levels[k]->rhs, -1.0, levels[k]->res);

  // coarsen the residual to next level, checking if the residual needs to be gathered after
  if (levels[k+1]->gatherLevel==true) {
    levels[k+1]->coarsen(levels[k+1]->coarsenArgs, levels[k]->res, levels[k+1]->Srhs);
    levels[k+1]->gather (levels[k+1]->gatherArgs,  levels[k+1]->Srhs, levels[k+1]->rhs);
  } else {
    levels[k+1]->coarsen(levels[k+1]->coarsenArgs, levels[k]->res, levels[k+1]->rhs);
  }

  vcycle(parAlmond,k+1);

  if (levels[k+1]->gatherLevel==true) {
    levels[k+1]->scatter(levels[k+1]->scatterArgs,  levels[k+1]->x, levels[k+1]->Sx);
    levels[k+1]->prolongate(levels[k+1]->prolongateArgs, levels[k+1]->Sx, levels[k]->x);
  } else {
    levels[k+1]->prolongate(levels[k+1]->prolongateArgs, levels[k+1]->x, levels[k]->x);
  }

  levels[k]->smooth(levels[k]->smoothArgs, levels[k]->rhs, levels[k]->x,false);

  occaTimerToc(parAlmond->device,name);
}


void device_vcycle(parAlmond_t *parAlmond, int k){

  agmgLevel **levels = parAlmond->levels;

  const dlong m = levels[k]->Nrows;
  // const dlong mCoarse = levels[k+1]->Nrows;

  // switch to cpu if the problem size is too small for gpu
  if(m < GPU_CPU_SWITCH_SIZE){
    levels[k]->o_rhs.copyTo(levels[k]->rhs, m*sizeof(dfloat));
    vcycle(parAlmond, k);
    levels[k]->o_x.copyFrom(levels[k]->x, m*sizeof(dfloat));
    return;
  }

  //check for base level
  if (k==parAlmond->numLevels-1) {
    if (parAlmond->invCoarseA != NULL) {
      //use exact sovler
      device_exactCoarseSolve(parAlmond, m, levels[k]->o_rhs, levels[k]->o_x);
    } else {
      levels[k]->device_smooth(levels[k]->smoothArgs, levels[k]->o_rhs, levels[k]->o_x, true);
    }
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

  // coarsen the residual to next level, checking if the residual needs to be gathered after
  if (levels[k+1]->gatherLevel==true) {
    levels[k+1]->device_coarsen(levels[k+1]->coarsenArgs, levels[k]->o_res, levels[k+1]->o_Srhs);
    levels[k+1]->device_gather (levels[k+1]->gatherArgs,  levels[k+1]->o_Srhs, levels[k+1]->o_rhs);
  } else {
    levels[k+1]->device_coarsen(levels[k+1]->coarsenArgs, levels[k]->o_res, levels[k+1]->o_rhs);
  }

  device_vcycle(parAlmond, k+1);

  if (levels[k+1]->gatherLevel==true) {
    levels[k+1]->device_scatter   (levels[k+1]->scatterArgs,  levels[k+1]->o_x, levels[k+1]->o_Sx);
    levels[k+1]->device_prolongate(levels[k+1]->prolongateArgs, levels[k+1]->o_Sx, levels[k]->o_x);
  } else {
    levels[k+1]->device_prolongate(levels[k+1]->prolongateArgs, levels[k+1]->o_x, levels[k]->o_x);
  }

  levels[k]->device_smooth(levels[k]->smoothArgs, levels[k]->o_rhs, levels[k]->o_x,false);

  occaTimerToc(parAlmond->device,name);
}

/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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
#include "cds.h"

void extbdfCoefficents(cds_t *cds, int order);

void cdsRunEXTBDF(cds_t *cds){

  mesh_t *mesh = cds->mesh;
  
  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"CDS");

  // Write Initial Data
   if(cds->outputStep) cdsReport(cds, cds->startTime, 0.0);

   for(int tstep=0;tstep<cds->NtimeSteps;++tstep){
     // for(int tstep=0;tstep<5;++tstep){

    if(tstep<1) 
      extbdfCoefficents(cds,tstep+1);
    else if(tstep<2 && cds->temporalOrder>=2) 
      extbdfCoefficents(cds,tstep+1);
    else if(tstep<3 && cds->temporalOrder>=3) 
      extbdfCoefficents(cds,tstep+1);

    dfloat time = cds->startTime + tstep*cds->dt;

#if 1  
    hlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
    if(cds->Nsubsteps) {
      cdsSubCycle(cds, time, cds->Nstages, cds->o_U, cds->o_S,  cds->o_NS);
    } else {
      cdsAdvection(cds, time, cds->o_U, cds->o_S, cds->o_NS);
    }    

    cdsHelmholtzRhs  (cds, time+cds->dt, cds->Nstages, cds->o_rhsS);

    cdsHelmholtzSolve(cds, time+cds->dt, cds->Nstages, cds->o_rhsS, cds->o_rkS);
   
    // cycle history
    for (int s=cds->Nstages;s>1;s--) {
      cds->o_S.copyFrom(cds->o_S, cds->Ntotal*cds->NSfields*sizeof(dfloat), 
                            			(s-1)*cds->Ntotal*cds->NSfields*sizeof(dfloat), 
                            			(s-2)*cds->Ntotal*cds->NSfields*sizeof(dfloat));

      cds->o_NS.copyFrom(cds->o_NS, cds->Ntotal*cds->NSfields*sizeof(dfloat), 
                            			 (s-1)*cds->Ntotal*cds->NSfields*sizeof(dfloat), 
                            			 (s-2)*cds->Ntotal*cds->NSfields*sizeof(dfloat));
    }

    //copy updated scalar
    cds->o_S.copyFrom(cds->o_rkS, cds->NSfields*cds->Ntotal*sizeof(dfloat)); 
#else 
    cdsSolveStep(cds, time, cds->dt, cds->o_U, cds->o_S);
#endif
 
    occaTimerTic(mesh->device,"Report");

    if(cds->outputStep){
      if(((tstep+1)%(cds->outputStep))==0){
         printf("\rtstep = %d, solver iteration: S - %3d \n", tstep+1, cds->Niter);
         cdsReport(cds, time+cds->dt, tstep+1);
      }
    }
    
    //  if (cds->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, cds->NiterU, cds->NiterV, cds->NiterP); fflush(stdout);
    // if (cds->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", tstep+1, cds->NiterU, cds->NiterV, cds->NiterW, cds->NiterP); fflush(stdout);
    
    occaTimerToc(mesh->device,"Report");
  }
  occaTimerToc(mesh->device,"CDS");


  dfloat finalTime = cds->NtimeSteps*cds->dt; printf("Writing Final Data...\n");

  if(cds->outputStep) cdsReport(cds, finalTime,cds->NtimeSteps);
  
  if(mesh->rank==0) occa::printTimer();
}


void extbdfCoefficents(cds_t *cds, int order) {

  if(order==1) {
     //advection, first order in time, increment
    cds->g0 =  1.0f; 
    dfloat extbdfB[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfA[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};
    
    memcpy(cds->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(cds->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(cds->extbdfC, extbdfC, 3*sizeof(dfloat));

    cds->o_extbdfB.copyFrom(extbdfB);
    cds->o_extbdfA.copyFrom(extbdfA);
    cds->o_extbdfC.copyFrom(extbdfC);

    cds->ExplicitOrder = 1;    
    
    cds->lambda = cds->g0 / (cds->dt * cds->alf);
    cds->ig0 = 1.0/cds->g0; 
  } else if(order==2) {
    //advection, second order in time, increment
    cds->g0 =  1.5f;
    dfloat extbdfB[3] = {2.0f,-0.5f, 0.0f};
    dfloat extbdfA[3] = {2.0f,-1.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};

    memcpy(cds->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(cds->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(cds->extbdfC, extbdfC, 3*sizeof(dfloat));

    cds->o_extbdfB.copyFrom(extbdfB);
    cds->o_extbdfA.copyFrom(extbdfA);
    cds->o_extbdfC.copyFrom(extbdfC);

    cds->ExplicitOrder=2;

    cds->lambda = cds->g0 / (cds->dt * cds->alf);
    cds->ig0 = 1.0/cds->g0; 
  } else if(order==3) {
    //advection, third order in time, increment
    cds->g0 =  11.f/6.f;
    dfloat extbdfB[3] = {3.0f,-1.5f, 1.0f/3.0f};
    dfloat extbdfA[3] = {3.0f,-3.0f, 1.0f};
    dfloat extbdfC[3] = {2.0f,-1.0f, 0.0f};
    
    memcpy(cds->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(cds->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(cds->extbdfC, extbdfC, 3*sizeof(dfloat));

    cds->o_extbdfB.copyFrom(extbdfB);
    cds->o_extbdfA.copyFrom(extbdfA);
    cds->o_extbdfC.copyFrom(extbdfC);

    cds->ExplicitOrder=3;

    cds->lambda = cds->g0 / (cds->dt * cds->alf); // alf = k/rho*cp thermal diffusivity
    cds->ig0 = 1.0/cds->g0; 
  }
}

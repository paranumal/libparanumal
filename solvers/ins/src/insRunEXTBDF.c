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

#include "ins.h"

void extbdfCoefficents(ins_t *ins, int order);

void insRunEXTBDF(ins_t *ins){

  mesh_t *mesh = ins->mesh;
  
  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"INS");

  int NstokesSteps = 0;
  dfloat oldDt = ins->dt;
  ins->dt *= 100;

  if (mesh->rank==0) printf("Number of Timesteps: %d\n", ins->NtimeSteps);
  for(int tstep=0;tstep<NstokesSteps;++tstep){
    if(tstep<1) 
      extbdfCoefficents(ins,tstep+1);
    else if(tstep<2 && ins->temporalOrder>=2) 
      extbdfCoefficents(ins,tstep+1);
    else if(tstep<3 && ins->temporalOrder>=3) 
      extbdfCoefficents(ins,tstep+1);

    insGradient (ins, 0, ins->o_P, ins->o_GP);

    insVelocityRhs  (ins, 0, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
    insVelocitySolve(ins, 0, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);

    insPressureRhs  (ins, 0, ins->Nstages);
    insPressureSolve(ins, 0, ins->Nstages); 

    insPressureUpdate(ins, 0, ins->Nstages, ins->o_rkP);
    insGradient(ins, 0, ins->o_rkP, ins->o_rkGP);

    //cycle history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_U.copyFrom(ins->o_U, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      ins->o_P.copyFrom(ins->o_P, ins->Ntotal*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*sizeof(dfloat));
    }

    //copy updated pressure
    ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat)); 

    //update velocity
    insVelocityUpdate(ins, 0, ins->Nstages, ins->o_rkGP, ins->o_rkU);

    //copy updated pressure
    ins->o_U.copyFrom(ins->o_rkU, ins->NVfields*ins->Ntotal*sizeof(dfloat)); 

    //cycle rhs history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_GP.copyFrom(ins->o_GP, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    }

    if (mesh->rank==0) printf("\rSstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
  }

  if (mesh->rank==0) printf("\n");

  ins->dt = oldDt;
  // Write Initial Data
  // if(ins->outputStep) insReport(ins, ins->startTime, 0);

  // for(int tstep=0;tstep<ins->NtimeSteps;++tstep){
  for(int tstep=0;tstep<1;++tstep){

    // if(ins->restartedFromFile){
      // if(tstep=0 && ins->temporalOrder>=2) 
      //   extbdfCoefficents(ins,2);
      // else if(tstep=0 && ins->temporalOrder>=3) 
      //   extbdfCoefficents(ins,3);
    // }else{
      if(tstep<1) 
        extbdfCoefficents(ins,tstep+1);
      else if(tstep<2 && ins->temporalOrder>=2) 
        extbdfCoefficents(ins,tstep+1);
      else if(tstep<3 && ins->temporalOrder>=3) 
        extbdfCoefficents(ins,tstep+1);
    // }

    dfloat time = ins->startTime + tstep*ins->dt;

    if(ins->Nsubsteps) {
      insSubCycle(ins, time, ins->Nstages, ins->o_U, ins->o_NU);
    } else {
      insAdvection(ins, time, ins->o_U, ins->o_NU);
    } 
    insGradient (ins, time, ins->o_P, ins->o_GP);

    insVelocityRhs  (ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
    insVelocitySolve(ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);

    // insPressureRhs  (ins, time+ins->dt, ins->Nstages);
    // insPressureSolve(ins, time+ins->dt, ins->Nstages); 

    // insPressureUpdate(ins, time+ins->dt, ins->Nstages, ins->o_rkP);
    // insGradient(ins, time+ins->dt, ins->o_rkP, ins->o_rkGP);

    // //cycle history
    // for (int s=ins->Nstages;s>1;s--) {
    //   ins->o_U.copyFrom(ins->o_U, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
    //                               (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
    //                               (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    //   ins->o_P.copyFrom(ins->o_P, ins->Ntotal*sizeof(dfloat), 
    //                               (s-1)*ins->Ntotal*sizeof(dfloat), 
    //                               (s-2)*ins->Ntotal*sizeof(dfloat));
    // }

    // //copy updated pressure
    // ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat)); 

    // //update velocity
    // insVelocityUpdate(ins, time+ins->dt, ins->Nstages, ins->o_rkGP, ins->o_rkU);

    // //copy updated pressure
    // ins->o_U.copyFrom(ins->o_rkU, ins->NVfields*ins->Ntotal*sizeof(dfloat)); 

    // //cycle rhs history
    // for (int s=ins->Nstages;s>1;s--) {
    //   ins->o_NU.copyFrom(ins->o_NU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
    //                               (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
    //                               (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    //   ins->o_GP.copyFrom(ins->o_GP, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
    //                               (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
    //                               (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    // }

    // occaTimerTic(mesh->device,"Report");

    // if(ins->outputStep){
    //   if(((tstep+1)%(ins->outputStep))==0){
    //     if (ins->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP);
    //     if (ins->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);
    //     insReport(ins, time+ins->dt, tstep+1);

    //     // Write a restart file
    //     if(ins->writeRestartFile){
    //       if(mesh->rank==0) printf("\nWriting Binary Restart File....");
    //         insRestartWrite(ins, ins->options, time+ins->dt);
    //       if(mesh->rank==0) printf("done\n");
    //     }

    //     // // Update Time-Step Size
    //     // if(ins->dtAdaptStep){
    //     //   if(((ins->tstep)%(ins->dtAdaptStep))==0){
    //     //     if(rank==0) printf("\n Adapting time Step Size to ");
    //     //       insComputeDt(ins, ins->time);
    //     //     if(rank==0) printf("%.4e\n", ins->dt);
    //     //      // Interpolate history for the new time step size
    //     //       insInterpolateHistory(ins, ins->dtold, ins->dt);

    //     //   }
    //     // } 
    //   }
    // }

    // if (ins->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
    // if (ins->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); fflush(stdout);
    
    // occaTimerToc(mesh->device,"Report");
  }
  occaTimerToc(mesh->device,"INS");


  dfloat finalTime = ins->NtimeSteps*ins->dt;
  printf("\n");

  // if(ins->outputStep) insReport(ins, finalTime,ins->NtimeSteps);
  
  if(mesh->rank==0) occa::printTimer();
}


void extbdfCoefficents(ins_t *ins, int order) {

  if(order==1) {
     //advection, first order in time, increment
    ins->g0 =  1.0f; 
    dfloat extbdfB[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfA[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};
    
    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder = 1;    
    
    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0; 
  } else if(order==2) {
    //advection, second order in time, increment
    ins->g0 =  1.5f;
    dfloat extbdfB[3] = {2.0f,-0.5f, 0.0f};
    dfloat extbdfA[3] = {2.0f,-1.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};

    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder=2;

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0; 
  } else if(order==3) {
    //advection, third order in time, increment
    ins->g0 =  11.f/6.f;
    dfloat extbdfB[3] = {3.0f,-1.5f, 1.0f/3.0f};
    dfloat extbdfA[3] = {3.0f,-3.0f, 1.0f};
    dfloat extbdfC[3] = {2.0f,-1.0f, 0.0f};
    
    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder=3;

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0; 
  }
}

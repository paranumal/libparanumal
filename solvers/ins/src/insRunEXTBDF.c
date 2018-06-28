#include "ins.h"

void extbdfCoefficents(ins_t *ins, int order);

void insRunEXTBDF(ins_t *ins){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh_t *mesh = ins->mesh;
  
  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"INS");

  int NstokesSteps = 0;
  dfloat oldDt = ins->dt;
  ins->dt *= 100;

  if (rank ==0) printf("Number of Timesteps: %d\n", ins->NtimeSteps);
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
      ins->o_GP.copyFrom(ins->o_NU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    }

    if (rank==0) printf("\rSstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
  }

  if (rank==0) printf("\n");

  ins->dt = oldDt;
  // Write Initial Data
  if(ins->outputStep) insReport(ins, 0.0, 0);

  for(int tstep=0;tstep<ins->NtimeSteps;++tstep){
    if(tstep<1) 
      extbdfCoefficents(ins,tstep+1);
    else if(tstep<2 && ins->temporalOrder>=2) 
      extbdfCoefficents(ins,tstep+1);
    else if(tstep<3 && ins->temporalOrder>=3) 
      extbdfCoefficents(ins,tstep+1);
    
    dfloat time = ins->startTime + tstep*ins->dt;

    if(ins->Nsubsteps) {
      insSubCycle(ins, time, ins->Nstages, ins->o_U, ins->o_NU);
    } else {
      insAdvection(ins, time, ins->o_U, ins->o_NU);
    } 
    insGradient (ins, time, ins->o_P, ins->o_GP);

    insVelocityRhs  (ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
    insVelocitySolve(ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);

    insPressureRhs  (ins, time+ins->dt, ins->Nstages);
    insPressureSolve(ins, time+ins->dt, ins->Nstages); 

    insPressureUpdate(ins, time+ins->dt, ins->Nstages, ins->o_rkP);
    insGradient(ins, time+ins->dt, ins->o_rkP, ins->o_rkGP);

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
    insVelocityUpdate(ins, time+ins->dt, ins->Nstages, ins->o_rkGP, ins->o_rkU);

    //copy updated pressure
    ins->o_U.copyFrom(ins->o_rkU, ins->NVfields*ins->Ntotal*sizeof(dfloat)); 

    //cycle rhs history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_NU.copyFrom(ins->o_NU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      // ins->o_GP.copyFrom(ins->o_GP, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
      //                             (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
      //                             (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      ins->o_GP.copyFrom(ins->o_NU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    }

    occaTimerTic(mesh->device,"Report");

    if(ins->outputStep){
      if(((tstep+1)%(ins->outputStep))==0){
        if (ins->dim==2 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP);
        if (ins->dim==3 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);
        insReport(ins, time+ins->dt, tstep+1);

        // Write a restart file
        if(ins->writeRestartFile){
          if(rank==0) printf("\nWriting Binary Restart File....");
            insRestartWrite(ins, ins->options, time+ins->dt);
          if(rank==0) printf("done\n");
        } 
      }
    }

    if (ins->dim==2 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
    if (ins->dim==3 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); fflush(stdout);
    
    occaTimerToc(mesh->device,"Report");
  }
  occaTimerToc(mesh->device,"INS");


  dfloat finalTime = ins->NtimeSteps*ins->dt;
  printf("\n");

  if(ins->outputStep) insReport(ins, finalTime,ins->NtimeSteps);
  
  if(rank==0) occa::printTimer();
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
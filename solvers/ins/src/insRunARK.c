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

void insRunARK(ins_t *ins){

  mesh_t *mesh = ins->mesh;
  
  //PID Control
  dfloat safe  = 0.8;   //safety factor

  //error control parameters
  dfloat beta    = 0.05; 
  dfloat factor1 = 0.2;
  dfloat factor2 = 5.0;
  dfloat exp1       = 0.25 - 0.75*beta;
  dfloat invfactor1 = 1.0/factor1;
  dfloat invfactor2 = 1.0/factor2;
  dfloat facold     = 1E-4;


  int profile_application = 1; 
  timer profiler(profile_application); profiler.initTimer(mesh->device);
  occaTimerTic(mesh->device,"INS");

  // occa::initTimer(mesh->device);
  // occaTimerTic(mesh->device,"INS");
#if 0
  // 
  int NstokesSteps = 100;
  dfloat oldDt     = ins->dt;
  ins->dt         *= 10;

  ins->time = 0.0; 

  if (mesh->rank ==0) printf("Running Initial Stokes Solve:\n");
  
  for(int tstep=0;tstep<NstokesSteps;++tstep){

    insDiffusion(ins, ins->time, ins->o_U, ins->o_LU);
    insGradient (ins, ins->time, ins->o_P, ins->o_GP);
    for(int stage=1;stage<=ins->Nstages;++stage){
      // intermediate stage time
      dfloat stageTime = ins->time + ins->rkC[stage]*ins->dt;
      insVelocityRhs  (ins, stageTime, stage, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
      insVelocitySolve(ins, stageTime, stage, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);

      insPressureRhs  (ins, stageTime, stage);
      insPressureSolve(ins, stageTime, stage);      

      insPressureUpdate(ins, stageTime, stage, ins->o_rkP);
      insGradient(ins, stageTime, ins->o_rkP, ins->o_rkGP);

      insVelocityUpdate(ins, stageTime, stage, ins->o_rkGP, ins->o_rkU);
      
      //compute and save NU and LU
      // insAdvection(ins, stageTime, ins->o_rkU, ins->o_rkNU);
      insDiffusion(ins, stageTime, ins->o_rkU, ins->o_rkLU); 


      ins->o_NU.copyFrom(ins->o_rkNU, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
      ins->o_LU.copyFrom(ins->o_rkLU, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);

      if (stage==ins->Nstages) break; //final stage
      ins->o_U.copyFrom(ins->o_rkU, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
      ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat), stage*ins->Ntotal*sizeof(dfloat), 0);
      ins->o_GP.copyFrom(ins->o_rkGP, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
    } 

    //accept the step and proceed
    ins->o_U.copyFrom(ins->o_rkU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
    ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat), 0);
    ins->tstep++;
    ins->time += ins->dt;

    if (ins->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", ins->tstep, ins->NiterU, ins->NiterV, ins->NiterP);
  }

  printf("done \n");

  insReport(ins, 0.0, 0);
#endif

  ins->tstep = 0;
  int done = 0;
  ins->time = ins->startTime;

  // Compute dt according the local CFL condition
  if(ins->dtAdaptStep) insComputeDt(ins, ins->time); 
  // Write Initial Data
  if(ins->outputStep) insReport(ins, 0.0, 0);
  // Write Initial Force Data (assumes U = o_U )
  if(ins->outputForceStep) insForces(ins, ins->time); 

  while (!done) {

    if (ins->dt<ins->dtMIN){
      printf("ERROR: Time step became too small at time step=%d\n", ins->tstep);
      exit (-1);
    }
    if (isnan(ins->dt)) {
      printf("ERROR: Solution became unstable at time step=%d\n", ins->tstep);
      exit (-1);
    }


    //check for final timestep
    if (ins->time+ins->dt > ins->finalTime){
      ins->dt = ins->finalTime-ins->time;
      done = 1;
    }

    insAdvection(ins, ins->time, ins->o_U, ins->o_NU);
    insDiffusion(ins, ins->time, ins->o_U, ins->o_LU);
    insGradient (ins, ins->time, ins->o_P, ins->o_GP);

    for(int stage=1;stage<=ins->Nstages;++stage){

      // intermediate stage time
      dfloat stageTime = ins->time + ins->rkC[stage]*ins->dt;

      insVelocityRhs  (ins, stageTime, stage, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
      insVelocitySolve(ins, stageTime, stage, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);

      insPressureRhs  (ins, stageTime, stage);
      insPressureSolve(ins, stageTime, stage);      

      insPressureUpdate(ins, stageTime, stage, ins->o_rkP);
      insGradient(ins, stageTime, ins->o_rkP, ins->o_rkGP);

      insVelocityUpdate(ins, stageTime, stage, ins->o_rkGP, ins->o_rkU);
      
      //compute and save NU and LU
      insAdvection(ins, stageTime, ins->o_rkU, ins->o_rkNU);
      insDiffusion(ins, stageTime, ins->o_rkU, ins->o_rkLU); 


      ins->o_NU.copyFrom(ins->o_rkNU, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
      ins->o_LU.copyFrom(ins->o_rkLU, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);

      if (stage==ins->Nstages) break; //final stage
      ins->o_U.copyFrom(ins->o_rkU, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
      ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat), stage*ins->Ntotal*sizeof(dfloat), 0);
      ins->o_GP.copyFrom(ins->o_rkGP, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
    } 

    if (ins->embeddedRKFlag==0) {//check if an embedded rk method is being used
      //accept the step and proceed
      ins->o_U.copyFrom(ins->o_rkU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
      ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat), 0);
      ins->tstep++;
      ins->time += ins->dt;

      occaTimerTic(mesh->device,"Report");
      if(ins->outputStep){
        if(((ins->tstep)%(ins->outputStep))==0){
          if (ins->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", ins->tstep+1, ins->NiterU, ins->NiterV, ins->NiterP);
          if (ins->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", ins->tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);
          insReport(ins, ins->time, ins->tstep);
        }
      }
      
      if(ins->outputForceStep){
        if(((ins->tstep)%(ins->outputForceStep))==0){
          ins->o_U.copyTo(ins->U);
          ins->o_P.copyTo(ins->P);
          insForces(ins, ins->time);
        }
      }

      // Update Time-Step Size
      if(ins->dtAdaptStep){
        if(((ins->tstep)%(ins->dtAdaptStep))==0){
          if(mesh->rank==0) printf("\n Adapting time Step Size to ");
          insComputeDt(ins, ins->time);
          if(mesh->rank==0) printf("%.4e\n", ins->dt);
        }
      }
    }



      
    
    if (ins->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", ins->tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
    if (ins->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", ins->tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); fflush(stdout);
    occaTimerToc(mesh->device,"Report");


/*
    dlong Nlocal = mesh->Nelements*mesh->Np*ins->NVfields;
    ins->errorEstimateKernel(Nlocal, 
                            ins->ATOL,
                            ins->RTOL,
                            ins->o_U,
                            ins->o_rkU,
                            ins->o_NU,
                            ins->o_LU,
                            ins->o_erkerr,
                            ins->o_irkerr,
                            ins->o_errtmp);

    ins->o_errtmp.copyTo(ins->errtmp);
    dfloat localerr = 0;
    dfloat err = 0;
    for(int n=0;n<ins->Nblock;++n){
      localerr += ins->errtmp[n];
    }
    MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

    err = sqrt(err/(ins->totalElements*mesh->Np*ins->NVfields));

    dfloat fac1 = pow(err,exp1);
    dfloat fac  = fac1/pow(facold,beta);

    fac = mymax(invfactor2, mymin(invfactor1,fac/safe));
    dfloat dtnew = ins->dt/fac;

    if(err<1.0){
      ins->o_q.copyFrom(ins->o_rkq);
      ins->o_pmlqx.copyFrom(ins->o_rkqx);
      ins->o_pmlqy.copyFrom(ins->o_rkqy);

      facold = mymax(err,1E-4);
      ins->time += ins->dt;

      ins->ins->tstep++;
    }else{
      ins->rtstep++; 
      dtnew = ins->dt/(mymax(invfactor1,fac1/safe));
      done =0;
    }

    ins->dt = dtnew;
    ins->atstep++;

    printf("\rTime = %.4e (%d). Average Dt = %.4e, Rejection rate = %.2g   ", time, ins->tstep, time/(dfloat)ins->tstep, Nregect/(dfloat) ins->tstep); fflush(stdout);
  */
  }
  occaTimerToc(mesh->device,"INS");


  dfloat finalTime = ins->NtimeSteps*ins->dt;
  printf("\n");
  insReport(ins, finalTime, ins->NtimeSteps);
  
  if(mesh->rank==0) profiler.printTimer();
}

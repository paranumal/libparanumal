#include "ins.h"

void insRunARK(ins_t *ins){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

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



  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"INS");

  // int NstokesSteps = 0;
  // dfloat oldDt = ins->dt;
  // ins->dt *= 100;

  // if (rank ==0) printf("Running Initial Stokes Solve:\n");
  // for(int tstep=0;tstep<NstokesSteps;++tstep){
    

  //   if (rank==0) printf("\rSstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
  // }

  // if (rank==0) printf("\n");

  // ins->dt = oldDt;


  // Write Initial Data
  insReport(ins, 0);


  ins->tstep = 0;
  int done = 0;
  dfloat time = ins->startTime;
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

    insAdvection(ins, time, ins->o_U, ins->o_NU);
    insDiffusion(ins, time, ins->o_U, ins->o_LU);
    insGrad(ins, time, ins->o_P, ins->o_GP);

    for(int stage=1;stage<ins->NrkStages;++stage){

      // intermediate stage time
      dfloat stageTime = time + ins->rkC[stage]*ins->dt;

      insVelocityRhs(ins, time);
      insVelocitySolve(ins, time);

      insPressureRhs(ins, time);
      insPressureSolve(ins, time);      

      insPressureUpdate(ins, ins->o_rkP);
      insGrad(ins, time, ins->o_rkP, ins->o_rkGP);

      insVelocityUpdate(ins, ins->o_rkU);

      insAdvection(ins, time, ins->o_rkU, ins->o_rkNU);
      insDiffusion(ins, time, ins->o_rkU, ins->o_rkLU); 
    } 

    occaTimerTic(mesh->device,"Report");
    if(((tstep+1)%(ins->outputStep))==0){
      if (ins->dim==2 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP);
      if (ins->dim==3 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);
      insReport(ins, tstep+1);
    }

    if (ins->dim==2 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
    if (ins->dim==3 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); fflush(stdout);
    
    occaTimerToc(mesh->device,"Report");

    dlong Ntotal = mesh->Nelements*mesh->Np*ins->NVfields;
    ins->errorEstimateKernel(Ntotal, 
                            ins->ATOL,
                            ins->RTOL,
                            ins->o_q,
                            ins->o_rkq,
                            ins->o_rkerr,
                            ins->o_errtmp);

    ins->o_errtmp.copyTo(ins->errtmp);
    dfloat localerr = 0;
    dfloat err = 0;
    for(int n=0;n<ins->Nblock;++n){
      localerr += ins->errtmp[n];
    }
    MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

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

      ins->tstep++;
    }else{
      ins->rtstep++; 
      dtnew = ins->dt/(mymax(invfactor1,fac1/safe));
      done =0;
    }

    ins->dt = dtnew;
    ins->atstep++;

    printf("\rTime = %.4e (%d). Average Dt = %.4e, Rejection rate = %.2g   ", time, tstep, time/(dfloat)tstep, Nregect/(dfloat) tstep); fflush(stdout);
  }
  occaTimerToc(mesh->device,"INS");


  dfloat finalTime = ins->NtimeSteps*ins->dt;
  printf("\n");
  insReport(ins, ins->NtimeSteps);
  
  if(rank==0) occa::printTimer();
}

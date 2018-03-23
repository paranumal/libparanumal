#include "insTet3D.h"

void insRunTet3D(ins_t *ins, char *options){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh3D *mesh = ins->mesh;
  
  // Write Initial Data
  insReportTet3D(ins, 0, options);

  occa::initTimer(mesh->device);

  //ins->NtimeSteps = 271000; 
  ins->NtimeSteps =10;

  // double tic_tot = 0.f, toc_tot = 0.f; 
  // double tic_adv = 0.f, toc_adv = 0.f;
  // double tic_pre = 0.f, toc_pre = 0.f;
  // double tic_vel = 0.f, toc_vel = 0.f;
  // double tic_upd = 0.f, toc_upd = 0.f;

  for(int tstep=0;tstep<ins->NtimeSteps;++tstep){
    if(tstep<1){
      //advection, first order in time, increment
      ins->b0 =  1.f,  ins->a0 =  1.0f, ins->c0 = 1.0f;  // 2
      ins->b1 =  0.f,  ins->a1 =  0.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.f; 
      ins->ExplicitOrder = 1;    
      
      ins->lambda = ins->g0 / (ins->dt * ins->nu);
      ins->idt = 1.0/ins->dt; 
      ins->ig0 = 1.0/ins->g0; 
    } else if(tstep<2) {
      //advection, second order in time, increment
      ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 1.0f;  // 2
      ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.5f;
      ins->ExplicitOrder=2;

      ins->lambda = ins->g0 / (ins->dt * ins->nu);
      ins->idt = 1.0/ins->dt; 
      ins->ig0 = 1.0/ins->g0; 
    // } else {
    //   //advection, third order in time, increment
    //   ins->b0 =  3.f,       ins->a0  =  3.0f, ins->c0 = 2.0f;
    //   ins->b1 = -1.5f,      ins->a1  = -3.0f, ins->c1 = -1.0f;
    //   ins->b2 =  1.f/3.f,   ins->a2  =  1.0f, ins->c2 =  0.0f;
    //   ins->g0 =  11.f/6.f;
    //   ins->ExplicitOrder=3;

    //   ins->lambda = ins->g0 / (ins->dt * ins->nu);
    //   ins->idt = 1.0/ins->dt; 
    //   ins->ig0 = 1.0/ins->g0; 
    }

    // mesh->device.finish();
    // MPI_Barrier(MPI_COMM_WORLD);
    // tic_tot = MPI_Wtime(); 
    // tic_adv = MPI_Wtime(); 
    if(strstr(options,"SUBCYCLING")) {
      insAdvectionSubCycleStepTet3D(ins, tstep, options);
    } else {
      insAdvectionStepTet3D(ins, tstep, options);
    }
    // mesh->device.finish();
    // MPI_Barrier(MPI_COMM_WORLD);
    // toc_adv = MPI_Wtime(); 

    // tic_vel = MPI_Wtime(); 
    insHelmholtzStepTet3D(ins, tstep, options);
    // mesh->device.finish();
    // MPI_Barrier(MPI_COMM_WORLD);
    // toc_vel = MPI_Wtime(); 

    // tic_pre = MPI_Wtime(); 
    insPoissonStepTet3D(  ins, tstep, options);
    // mesh->device.finish();
    // MPI_Barrier(MPI_COMM_WORLD);
    // toc_pre = MPI_Wtime(); 

    // tic_upd = MPI_Wtime(); 
    insUpdateStepTet3D(   ins, tstep, options);
    // mesh->device.finish();
    // MPI_Barrier(MPI_COMM_WORLD);
    // toc_upd = MPI_Wtime(); 
    // toc_tot = MPI_Wtime(); 

    if(((tstep+1)%(ins->errorStep))==0){
      if (rank==0) printf("\rtstep = %d, time = %3.2E, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", tstep+1, (tstep+1)*ins->dt, ins->NiterU, ins->NiterV, ins->NiterW,  ins->NiterP);
      insReportTet3D(ins, tstep+1,options);
    }

    if (rank==0) printf("\rtstep = %d, time = %3.2E, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", tstep+1, (tstep+1)*ins->dt, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); fflush(stdout);
    //if (rank==0) printf("\ntotaltime = %3.2E, advectiontime = %3.2E, velocitytime = %3.2E, pressuretime = %3.2E, updatetime = %3.2E \n", toc_tot- tic_tot, toc_adv- tic_adv, toc_vel- tic_vel, toc_pre- tic_pre, toc_upd- tic_upd );
  }
    
  // For Final Time
  dfloat finaltime = (ins->NtimeSteps)*ins->dt;
  printf("\n");
  insReportTet3D(ins, ins->NtimeSteps,options);
  //insErrorNormsTet3D(ins, ins->finalTime, options);
}




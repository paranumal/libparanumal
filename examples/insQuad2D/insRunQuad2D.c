#include "insQuad2D.h"

void insRunQuad2D(ins_t *ins, char *options){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = ins->mesh;
  
  // Write Initial Data
  insReportQuad2D(ins, 0, options);
  
  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"INS");
  
  // if(ins->Nsubsteps)
  // ins->NtimeSteps = 160/ins->Nsubsteps;
  // else
  ins->NtimeSteps=32000;

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

    // if(strstr(options,"ALGEBRAIC")){
    if(strstr(options,"SUBCYCLING")) {
        occaTimerTic(mesh->device,"AdvectionSubStep");
        insAdvectionSubCycleStepQuad2D(ins, tstep, options);
        occaTimerToc(mesh->device,"AdvectionSubStep");
    } else {
      occaTimerTic(mesh->device,"AdvectionStep");
      insAdvectionStepQuad2D(ins, tstep, options);
      occaTimerToc(mesh->device,"AdvectionStep");
    }

    occaTimerTic(mesh->device,"HelmholtzStep");
    insHelmholtzStepQuad2D(ins, tstep, options); 
    occaTimerToc(mesh->device,"HelmholtzStep");

    occaTimerTic(mesh->device,"PoissonStep");
    insPoissonStepQuad2D(ins, tstep, options);
    occaTimerToc(mesh->device,"PoissonStep");

    occaTimerTic(mesh->device,"UpdateStep");
    insUpdateStepQuad2D(ins, tstep, options);
    occaTimerToc(mesh->device,"UpdateStep"); 
    
    occaTimerTic(mesh->device,"Report");
    if(strstr(options, "VTU")){
      if(((tstep+1)%(ins->errorStep))==0){
        if (rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP);
        insReportQuad2D(ins, tstep+1,options);
      }
    }

    if (rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
    
    occaTimerToc(mesh->device,"Report");
  }

  occaTimerToc(mesh->device,"INS");

  dfloat finalTime = ins->NtimeSteps*ins->dt;
  printf("\n");
  insReportQuad2D(ins, ins->NtimeSteps,options);
  if(rank==0) occa::printTimer();
}




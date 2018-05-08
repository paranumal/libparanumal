#include "ins.h"

void insRun(ins_t *ins){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh_t *mesh = ins->mesh;
  
  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"INS");

  int NstokesSteps = 0;
  dfloat oldDt = ins->dt;
  ins->dt *= 100;

  if (rank ==0) printf("Running Initial Stokes Solve:\n");
  for(int tstep=0;tstep<NstokesSteps;++tstep){
    if(tstep<1){
       //advection, first order in time, increment
      ins->b0 =  1.f,  ins->a0 =  0.0f, ins->c0 = 1.0f;  // 2
      ins->b1 =  0.f,  ins->a1 =  0.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.f; 
      ins->ExplicitOrder = 1;    
      
      ins->lambda = ins->g0 / (ins->dt * ins->nu);
      ins->idt = 1.0/ins->dt; 
      ins->ig0 = 1.0/ins->g0; 
    } else if(tstep<2) {
      //advection, second order in time, increment
      ins->b0 =  2.f,  ins->a0 =  0.0f, ins->c0 = 1.0f;  // 2
      ins->b1 = -0.5f, ins->a1 =  0.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.5f;
      ins->ExplicitOrder=2;

      ins->lambda = ins->g0 / (ins->dt * ins->nu);
      ins->idt = 1.0/ins->dt; 
      ins->ig0 = 1.0/ins->g0; 
    } else {
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

    if(ins->Nsubsteps) {
      occaTimerTic(mesh->device,"AdvectionSubStep");
      insAdvectionSubCycleStep(ins, 0.0);
      occaTimerToc(mesh->device,"AdvectionSubStep");
    } else {
      occaTimerTic(mesh->device,"AdvectionStep");
      insAdvectionStep(ins, 0.0);
      occaTimerToc(mesh->device,"AdvectionStep");
    }

    occaTimerTic(mesh->device,"HelmholtzStep");
    insHelmholtzStep(ins, 0.0); 
    occaTimerToc(mesh->device,"HelmholtzStep");

    occaTimerTic(mesh->device,"PoissonStep");
    insPoissonStep(ins, 0.0);
    occaTimerToc(mesh->device,"PoissonStep");

    occaTimerTic(mesh->device,"UpdateStep");
    insUpdateStep(ins, 0.0);
    occaTimerToc(mesh->device,"UpdateStep"); 

    if (rank==0) printf("\rSstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
  }

  if (rank==0) printf("\n");

  ins->dt = oldDt;
  // Write Initial Data
  insReport(ins, 0);

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
    
    dfloat time = tstep*ins->dt;

    // if(strstr(options,"ALGEBRAIC")){
    if(ins->Nsubsteps) {
      occaTimerTic(mesh->device,"AdvectionSubStep");
      insAdvectionSubCycleStep(ins, time);
      occaTimerToc(mesh->device,"AdvectionSubStep");
    } else {
      occaTimerTic(mesh->device,"AdvectionStep");
      insAdvectionStep(ins, time);
      occaTimerToc(mesh->device,"AdvectionStep");
    } 

    occaTimerTic(mesh->device,"HelmholtzStep");
    insHelmholtzStep(ins, time+ins->dt); 
    occaTimerToc(mesh->device,"HelmholtzStep");

    occaTimerTic(mesh->device,"PoissonStep");
    insPoissonStep(ins, time+ins->dt);
    occaTimerToc(mesh->device,"PoissonStep");

    occaTimerTic(mesh->device,"UpdateStep");
    insUpdateStep(ins, time+ins->dt);
    occaTimerToc(mesh->device,"UpdateStep"); 

    occaTimerTic(mesh->device,"Report");
    if(((tstep+1)%(ins->outputStep))==0){
      if (ins->dim==2 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP);
      if (ins->dim==3 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);
      insReport(ins, tstep+1);
    }

    if (ins->dim==2 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
    if (ins->dim==3 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); fflush(stdout);
    
    occaTimerToc(mesh->device,"Report");
  }
  occaTimerToc(mesh->device,"INS");


  dfloat finalTime = ins->NtimeSteps*ins->dt;
  printf("\n");
  insReport(ins, ins->NtimeSteps);
  
  if(rank==0) occa::printTimer();
}

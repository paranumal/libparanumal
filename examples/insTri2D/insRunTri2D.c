#include "insTri2D.h"

void insRunTri2D(ins_t *ins, char *options){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = ins->mesh;
  
  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"INS");


  ins->adv_time  = 0.f;
  ins->velx_time = 0.f;
  ins->vely_time = 0.f;
  ins->pr_time   = 0.f;
  ins->ins_time  = 0.f;

  // if(ins->Nsubsteps)
  // ins->NtimeSteps = 160/ins->Nsubsteps;
  // else
  //ins->NtimeSteps=32000;
  
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

    // if(strstr(options,"ALGEBRAIC")){
    if(strstr(options,"SUBCYCLING")) {
      occaTimerTic(mesh->device,"AdvectionSubStep");
      insAdvectionSubCycleStepTri2D(ins, 0, options);
      occaTimerToc(mesh->device,"AdvectionSubStep");
    } else {
      occaTimerTic(mesh->device,"AdvectionStep");
      insAdvectionStepTri2D(ins, 0, options);
      occaTimerToc(mesh->device,"AdvectionStep");
    }

    occaTimerTic(mesh->device,"HelmholtzStep");
    insHelmholtzStepTri2D(ins, 0, options); 
    occaTimerToc(mesh->device,"HelmholtzStep");

    occaTimerTic(mesh->device,"PoissonStep");
    insPoissonStepTri2D(ins, 0, options);
    occaTimerToc(mesh->device,"PoissonStep");

    occaTimerTic(mesh->device,"UpdateStep");
    insUpdateStepTri2D(ins, 0, options);
    occaTimerToc(mesh->device,"UpdateStep"); 

    if (rank==0) printf("\rSstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
  }

  if (rank==0) printf("\n");

  ins->dt = oldDt;
  // Write Initial Data
  insReportTri2D(ins, 0, options);


  dfloat ins_start = MPI_Wtime();

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
    
    dfloat adv_start = MPI_Wtime();
    // if(strstr(options,"ALGEBRAIC")){
    if(strstr(options,"SUBCYCLING")) {
      occaTimerTic(mesh->device,"AdvectionSubStep");
      insAdvectionSubCycleStepTri2D(ins, tstep, options);
      occaTimerToc(mesh->device,"AdvectionSubStep");
    } else {
      occaTimerTic(mesh->device,"AdvectionStep");
      insAdvectionStepTri2D(ins, tstep, options);
      occaTimerToc(mesh->device,"AdvectionStep");
    } 
    dfloat adv_end = MPI_Wtime();
    ins->adv_time    += (adv_end - adv_start); 




    occaTimerTic(mesh->device,"HelmholtzStep");
    insHelmholtzStepTri2D(ins, tstep, options); 
    occaTimerToc(mesh->device,"HelmholtzStep");

    occaTimerTic(mesh->device,"PoissonStep");
    insPoissonStepTri2D(ins, tstep, options);
    occaTimerToc(mesh->device,"PoissonStep");

    occaTimerTic(mesh->device,"UpdateStep");
    insUpdateStepTri2D(ins, tstep, options);
    occaTimerToc(mesh->device,"UpdateStep"); 

    occaTimerTic(mesh->device,"Report");
    if(strstr(options, "VTU")){
      if(((tstep+1)%(ins->errorStep))==0){
        if (rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP);
        insReportTri2D(ins, tstep+1,options);
      }
    }

    if(((tstep+1)%(ins->errorStep))==0){
      int Ns = 0;

      if(strstr(options,"SUBCYCLING"))
        Ns = 0;
      else
      Ns = ins->Nsubsteps;

       char fname[BUFSIZ];
       sprintf(fname, "IterationNumbers_%04d_%04d.dat",mesh->N, Ns);

       FILE *fp; fp = fopen(fname, "a");
       fprintf(fp, "%.4e %2d %2d %2d\n", (tstep+1)*ins->dt, ins->NiterU, ins->NiterV, ins->NiterP);
       fclose(fp);
      }


    if (rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
    
    occaTimerToc(mesh->device,"Report");
  }
  
  dfloat ins_end = MPI_Wtime();
  ins->ins_time = ins_end -ins_start; 

  occaTimerToc(mesh->device,"INS");


  char fname[BUFSIZ];
  sprintf(fname, "TimerResult.dat");
  

  int Ns = 0;
   if(strstr(options,"SUBCYCLING"))
      Ns = 0;
   else
     Ns = ins->Nsubsteps;
  FILE *fp; fp = fopen(fname, "a");
  fprintf(fp, "%2d %2d %.4e %.4e %.4e %.4e %.4e\n", mesh->N,Ns,ins->adv_time, ins->velx_time, 
                                   ins->vely_time, ins->pr_time,ins->ins_time);
  fclose(fp);





  dfloat finalTime = ins->NtimeSteps*ins->dt;
  printf("\n");
  insReportTri2D(ins, ins->NtimeSteps,options);
  //insErrorNorms2D(ins, finalTime, options);
  
  if(rank==0) occa::printTimer();
}

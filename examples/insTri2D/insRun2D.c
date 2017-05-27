#include "ins2D.h"

void insRun2D(ins_t *ins, char *options, char *velSolverOptions, char *prSolverOptions){

  mesh2D *mesh = ins->mesh; 
  occa::kernelInfo kernelInfo = ins->kernelInfo;
  
    
  // Set-up Pressure Incriment Solver
  solver_t *prsolver = ellipticSolveSetupTri2D(mesh,1.0f, kernelInfo, prSolverOptions); 
  ins->prsolver = prsolver; 

















  

  // Allocate MPI send buffer
  iint helmholtzHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NTfields)*sizeof(dfloat);
  dfloat *helmholtzSendBuffer = (dfloat*) malloc(helmholtzHaloBytes);
  dfloat *helmholtzRecvBuffer = (dfloat*) malloc(helmholtzHaloBytes);
  //
  iint poissonHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat);
  dfloat *poissonSendBuffer = (dfloat*) malloc(poissonHaloBytes);
  dfloat *poissonRecvBuffer = (dfloat*) malloc(poissonHaloBytes);

  // No need to do like this, just for consistency 
  iint updateHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  dfloat *updateSendBuffer = (dfloat*) malloc(updateHaloBytes);
  dfloat *updateRecvBuffer = (dfloat*) malloc(updateHaloBytes);

  occa::initTimer(mesh->device);

 
  //First Order Coefficients
  ins->a0 = 1.0, ins->b0 = 1.f;
  ins->a1 = 0.f, ins->b1 = 0.f;
  ins->b2 = 0.f, ins->b2 = 0.f;
  ins->g0 = 1.f; 
  

  // // Set-up First Order Solver
  // ins->lamda = ins->g0 / (ins->dt * ins->nu);
  // solver_t *velsolver = ellipticSolveSetupTri2D(mesh, ins->lamda, kernelInfo, velSolverOptions); 
  // ins->velsolver = velsolver; 

  // One First Order Step 
  //insHelmholtzStep2D(ins, 0, helmholtzHaloBytes, helmholtzSendBuffer, helmholtzRecvBuffer, options, velSolverOptions);
  //insPoissonStep2D(ins, 0, poissonHaloBytes, poissonSendBuffer, poissonRecvBuffer, options, prSolverOptions);
  //insUpdateStep2D(ins, 0, updateHaloBytes, updateSendBuffer, updateRecvBuffer, options);
  //
  printf("Finished first order step\n");
  insReport2D(ins, ins->NtimeSteps,options);



  // // Switch to second order
  // ins->a0 =  2.0,  ins->b0  = 2.0f;
  // ins->a1 = -0.5f, ins->b1 = -1.0f;
  // ins->b2 =  0.f,  ins->b2  = 0.f;
  // ins->g0 =  1.5f; 

  // // Set-up First Order Solver
  // ins->lamda = ins->g0 / (ins->dt * ins->nu);
  // velsolver = ellipticSolveSetupTri2D(mesh, ins->lamda, kernelInfo, velSolverOptions); 
  // ins->velsolver = velsolver; 

//   insHelmholtzStep2D(ins, 1, helmholtzHaloBytes, helmholtzSendBuffer, helmholtzRecvBuffer, options, velSolverOptions);
//   insPoissonStep2D(ins, 1, poissonHaloBytes, poissonSendBuffer, poissonRecvBuffer, options,prSolverOptions);
//   insUpdateStep2D(ins, 1, updateHaloBytes, updateSendBuffer, updateRecvBuffer, options);
//   printf("Finished second order step\n");

// insReport2D(ins, ins->NtimeSteps,options);
  

   // Switch to third order
  // ins->a0 =  3.f,       ins->b0  =  3.0f;
  // ins->a1 = -1.5f,      ins->b1  = -3.0f;
  // ins->b2 =  1.f/3.f,   ins->b2  =  1.0f;
  // ins->g0 =  11.f/6.f; 

  // // Set-up First Order Solver
  // ins->lamda = ins->g0 / (ins->dt * ins->nu);
  // velsolver = ellipticSolveSetupTri2D(mesh, ins->lamda, kernelInfo, velSolverOptions); 
  // ins->velsolver = velsolver; 
  // printf("Finished third order solver-setup\n");


  // for(iint tstep=2;tstep<ins->NtimeSteps;++tstep){
  //      //
  //     insHelmholtzStep2D(ins, tstep, helmholtzHaloBytes, helmholtzSendBuffer, helmholtzRecvBuffer, options);
  //     insPoissonStep2D(ins, tstep, poissonHaloBytes, poissonSendBuffer, poissonRecvBuffer, options);
  //     insUpdateStep2D(ins, tstep, updateHaloBytes, updateSendBuffer, updateRecvBuffer, options);

  //    if(strstr(options, "REPORT")){
  //     if((tstep%ins->errorStep)==0){
  //       insReport2D(ins, tstep,options);
  //      }
  //    }
  // }
  
  // // For Final Time
  // insReport2D(ins, ins->NtimeSteps,options);

  occa::printTimer();

  // Deallocate Halo MPI storage
  free(helmholtzSendBuffer);
  free(helmholtzRecvBuffer);
  //
  free(poissonSendBuffer);
  free(poissonRecvBuffer);
  //
  free(updateSendBuffer);
  free(updateRecvBuffer);
  //

}




#include "ins2D.h"

void insRun2D(solver_t *ins, char *options){

  mesh2D *mesh = ins->mesh; 

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

   
  insHelmholtzStep2D(ins, 0, helmholtzHaloBytes, helmholtzSendBuffer, helmholtzRecvBuffer, options);
  insPoissonStep2D(ins, 0, poissonHaloBytes, poissonSendBuffer, poissonRecvBuffer, options);
  insUpdateStep2D(ins, 0, updateHaloBytes, updateSendBuffer, updateRecvBuffer, options);

  // Switch to second order
  //ins->a0 = 2.0, ins->b0 = 2.0,  ins->a1 = -0.5, ins->b1 =-1.0, ins->g0 = 1.5; 


  occa::initTimer(mesh->device);

  for(iint tstep=2;tstep<ins->NtimeSteps;++tstep){
       //
       //insHelmholtzStep2D(ins, tstep, helmholtzHaloBytes, helmholtzSendBuffer, helmholtzRecvBuffer, options);
       
       //insPressureStep2D(ins, tstep, haloBytes, sendBuffer, recvBuffer,options);

      // insUpdateStep2D(ins, tstep, updateHaloBytes, updateSendBuffer, updateRecvBuffer, options);










     if(strstr(options, "REPORT")){
      if((tstep%ins->errorStep)==0){
        insReport2D(ins, tstep,options);
      }
    }
    }
  
  // // For Final Time
  insReport2D(ins, ins->NtimeSteps,options);

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




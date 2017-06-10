#include "ins2D.h"

void insRun2D(ins_t *ins, char *options){

  mesh2D *mesh = ins->mesh; 
  // Write Initial Data
  insReport2D(ins, 0, options);
  // Allocate MPI buffer for velocity step solver!! May Change Later!!!!!!
  iint tHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NTfields)*sizeof(dfloat);
  dfloat *tSendBuffer = (dfloat*) malloc(tHaloBytes);
  dfloat *tRecvBuffer = (dfloat*) malloc(tHaloBytes);
  
  iint vHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat);
  dfloat *vSendBuffer = (dfloat*) malloc(vHaloBytes);
  dfloat *vRecvBuffer = (dfloat*) malloc(vHaloBytes);

  // No need to do like this, just for consistency 
  iint pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  dfloat *pSendBuffer = (dfloat*) malloc(pHaloBytes);
  dfloat *pRecvBuffer = (dfloat*) malloc(pHaloBytes);

  occa::initTimer(mesh->device);
  
  for(iint tstep=0;tstep<10;++tstep){

    switch(tstep){
    case 0:
      ins->b0 = 1.f, ins->a0 = 1.f;
      ins->b1 = 0.f, ins->a1 = 0.f;
      ins->b2 = 0.f, ins->a2 = 0.f;
      ins->g0 = 1.f; 
      break;

    case 1:
      ins->b0 =  2.f,  ins->a0  = 2.0f;
      ins->b1 = -0.5f, ins->a1 = -1.0f;
      ins->b2 =  0.f,  ins->a2  = 0.f;
      ins->g0 =  1.5f;
      break;

    case 2:
      ins->b0 =  3.f,       ins->a0  =  3.0f;
      ins->b1 = -1.5f,      ins->a1  = -3.0f;
      ins->b2 =  1.f/3.f,   ins->a2  =  1.0f;
      ins->g0 =  11.f/6.f; 
      break;
    }
    
    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    
    insAdvectionStep2D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
    insHelmholtzStep2D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
    insPoissonStep2D(  ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    insUpdateStep2D(   ins, tstep, pHaloBytes,pSendBuffer,pRecvBuffer, options);
    
    if(strstr(options, "REPORT")){
      if((tstep%ins->errorStep)==0){
	      insReport2D(ins, tstep,options);
      }
    }
  }
  
  // For Final Time
  insReport2D(ins, ins->NtimeSteps,options);

  // Deallocate Halo MPI storage
  free(tSendBuffer);
  free(tRecvBuffer);
  free(vSendBuffer);
  free(vRecvBuffer);
  free(pSendBuffer);
  free(pRecvBuffer);
}




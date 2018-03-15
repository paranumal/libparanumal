#include "ins2D.h"

void insRunQuad2D(ins_t *ins, char *options){
 mesh2D *mesh = ins->mesh;
  // Write Initial Data
  //insReport2D(ins, 0, options);
  // Allocate MPI buffer for velocity step solver!! May Change Later!!!!!!
  int tHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NTfields)*sizeof(dfloat);
  dfloat *tSendBuffer = (dfloat*) malloc(tHaloBytes);
  dfloat *tRecvBuffer = (dfloat*) malloc(tHaloBytes);

  int vHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat);
  dfloat *vSendBuffer = (dfloat*) malloc(vHaloBytes);
  dfloat *vRecvBuffer = (dfloat*) malloc(vHaloBytes);

  // No need to do like this, just for consistency
  int pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  dfloat *pSendBuffer = (dfloat*) malloc(pHaloBytes);
  dfloat *pRecvBuffer = (dfloat*) malloc(pHaloBytes);

  // Set subscycling
  int subcycling =0;
  if(strstr(options,"SUBCYCLING")){ subcycling = 1; }

  occa::initTimer(mesh->device);

  for(int tstep=0;tstep<ins->NtimeSteps;++tstep){
   if(tstep<1){
      //advection, first order in time, no increment
     ins->b0 =  1.f,  ins->a0 =  1.0f, ins->c0 = 1.0f;  // 2
     ins->b1 =  0.f,  ins->a1 =  0.0f, ins->c1 = 0.0f; // -1
     ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
     ins->g0 =  1.f;      
   }
   else 
   //  //if(tstep<2)
    {
    // advection, second order in time, no increment
     ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 1.0f;  // 2
     ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
     ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
     ins->g0 =  1.5f;
     }
  
    ins->lambda = ins->g0 / (ins->dt * ins->nu);



    #if 0

    switch(subcycling){
      case 1:
        insAdvectionSubCycleStep2D(ins, tstep,tSendBuffer,tRecvBuffer,vSendBuffer,vRecvBuffer, options);
      break;

      case 0:
       insAdvectionStep2D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
      break;
    }
   
    insHelmholtzStep2D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
    insPoissonStep2D(  ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    insUpdateStep2D(   ins, tstep, pHaloBytes,pSendBuffer,pRecvBuffer, options);
   
    printf("tstep = %d\n", tstep);
    if(strstr(options, "REPORT")){
      if(((tstep+1)%ins->errorStep)==0){
        insReport2D(ins, tstep+1,options);
      }
    }
    #endif
  }


  

  #if 0
   // For Final Time
  insReport2D(ins, ins->NtimeSteps+1,options);
  insErrorNorms2D(ins, ins->finalTime, options);
  #endif
  
  // Deallocate Halo MPI storage
  free(tSendBuffer);
  free(tRecvBuffer);
  free(vSendBuffer);
  free(vRecvBuffer);
  free(pSendBuffer);
  free(pRecvBuffer);
}




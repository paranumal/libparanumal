#include "ins2D.h"

void insRun2D(solver_t *ins, char *options){

  mesh2D *mesh = ins->mesh; 

  // Allocate MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*ins->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  occa::initTimer(mesh->device);

  for(iint tstep=0;tstep<ins->NtimeSteps;++tstep){
      
       insAdvectionStep2D(ins, tstep, haloBytes, sendBuffer, recvBuffer,options);

     if(strstr(options, "REPORT")){
      if((tstep%mesh->errorStep)==0){
        insReport2D(ins, tstep,options);
      }
    }
    }
  
  // // For Final Time
  insReport2D(ins, ins->NtimeSteps,options);

  occa::printTimer();

  // Deallocate Halo MPI storage
  free(recvBuffer);
  free(sendBuffer);
}




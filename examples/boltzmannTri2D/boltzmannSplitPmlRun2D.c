#include "boltzmann2D.h"

void boltzmannSplitPmlRun2D(mesh2D *mesh, char *options){

  // Allocate MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  occa::initTimer(mesh->device);


   if(strstr(options, "LSERK")){

    for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

      boltzmannSplitPmlLserkStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer, options);

      if((tstep%mesh->errorStep)==0){
        boltzmannReport2D(mesh, tstep,options);
      }
    }
  }


  if(strstr(options, "LSIMEX")){

    for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

      boltzmannSplitPmlLsimexStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer,options);

      if((tstep%mesh->errorStep)==0){
        boltzmannReport2D(mesh, tstep,options);
      }
    }
  }


  if(strstr(options, "SARK3")){

    for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

      boltzmannSplitPmlSark3Step2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer,options);

      if((tstep%mesh->errorStep)==0){
        boltzmannReport2D(mesh, tstep,options);
      }
    }
  }

   if(strstr(options, "SAAB3")){

    for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

      boltzmannSplitPmlSaab3Step2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer,options);

      if((tstep%mesh->errorStep)==0){
        boltzmannReport2D(mesh, tstep,options);
      }
    }
  }
 
  
  boltzmannReport2D(mesh, mesh->NtimeSteps,options);

  occa::printTimer();

  // Deallocate Halo MPI storage
  free(recvBuffer);
  free(sendBuffer);
}




#include "boltzmann2D.h"

void boltzmannRun2D(mesh2D *mesh, char *options){

  // Allocate MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  occa::initTimer(mesh->device);


   

    for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){
      if(strstr(options, "LSERK")){
      boltzmannLserkStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer, options);
      }

      if(strstr(options, "LSIMEX")){

        if(strstr(options, "UNSPLITPML")){
          boltzmannUnsplitPmlLsimexStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer,options);
          }
        else{
          boltzmannSplitPmlLsimexStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer,options);

        }
      }
      
      if(strstr(options, "SARK3")){
       boltzmannSark3Step2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer,options);
      }

      if(strstr(options, "SAAB3")){

       boltzmannSaab3Step2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer,options);
      }

     if(strstr(options, "REPORT")){
      if((tstep%mesh->errorStep)==0){
        boltzmannReport2D(mesh, tstep,options);
      }
     }
    }
  
  // For Final Time
  boltzmannReport2D(mesh, mesh->NtimeSteps,options);

  occa::printTimer();

  // Deallocate Halo MPI storage
  free(recvBuffer);
  free(sendBuffer);
}




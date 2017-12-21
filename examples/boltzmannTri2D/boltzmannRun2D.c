#include "boltzmann2D.h"

void boltzmannRun2D(mesh2D *mesh, char *options){

   // For initial data
  boltzmannReport2D(mesh, 0 ,options);

  occa::initTimer(mesh->device);

      // Allocate MPI send buffer for single rate integrators
    iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
    dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
    dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
    // VOLUME KERNELS
    mesh->device.finish();
    occa::tic("Boltzmann Solver");

  // for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){
    for(iint tstep=0;tstep<100;++tstep){
      if(strstr(options, "LSERK")){
      boltzmannLserkStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer, options);
      }

      if(strstr(options, "LSIMEX")){
      boltzmannLsimexStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer,options);
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

  mesh->device.finish();
  occa::toc("Boltzmann Solver");
  
  // For Final Time
  boltzmannReport2D(mesh, mesh->NtimeSteps,options);

  occa::printTimer();

  // Deallocate Halo MPI storage
  free(recvBuffer);
  free(sendBuffer);
}




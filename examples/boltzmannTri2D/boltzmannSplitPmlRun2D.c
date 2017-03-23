#include "boltzmann2D.h"

void boltzmannSplitPmlRun2D(mesh2D *mesh){

  // Allocate MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  occa::initTimer(mesh->device);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    // perform a 5-stage LSERK4 time step
    #if TIME_DISC==LSERK
      boltzmannSplitPmlLserkStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer);
    #endif

    #if TIME_DISC==LSIMEX
      boltzmannSplitPmlLsimexStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer);
    #endif  
      
    #if TIME_DISC==MRAB
      boltzmannSplitPmlMrabStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer);
    #endif
   // Perform semi-analytic integration for pml damping term
    #if TIME_DISC==SAAB
      boltzmannSplitPmlSaabStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer);
    #endif
    // Perform semi-analytic integration for pml damping term
    #if TIME_DISC==SARK
      boltzmannSplitPmlSarkStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer);
    #endif

    // output statistics if this is an output step
    if((tstep%mesh->errorStep)==0){
      boltzmannReport2D(mesh, tstep);
    }
  }
  
  boltzmannReport2D(mesh, mesh->NtimeSteps);

  
  occa::printTimer();

  // Deallocate Halo MPI storage
  free(recvBuffer);
  free(sendBuffer);
}




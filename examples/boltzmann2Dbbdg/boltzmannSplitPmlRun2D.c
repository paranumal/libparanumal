#include "boltzmann2D.h"

void boltzmannSplitPmlRun2D(mesh2D *mesh){

  // Allocate MPI send buffer
  int haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  occa::initTimer(mesh->device);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

    // perform a 5-stage LSERK4 time step
    boltzmannSplitPmlLserkStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer);

    // output statistics if this is an output step
    if((tstep%mesh->errorStep)==0){
      boltzmannReport2D(mesh, tstep);
    }
  }

  occa::printTimer();

  // Deallocate Halo MPI storage
  free(recvBuffer);
  free(sendBuffer);
}




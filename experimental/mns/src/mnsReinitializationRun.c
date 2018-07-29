/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "mns.h"

void mnsReinitializationRun(mns_t *mns){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh_t *mesh = mns->mesh; 

  // MPI send buffer
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  int haloBytes;

  haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);

  if (haloBytes) {
  occa::memory o_sendBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
  occa::memory o_recvBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
  sendBuffer = (dfloat*) o_sendBufferPinned.getMappedPointer();
  recvBuffer = (dfloat*) o_recvBufferPinned.getMappedPointer();
  }

  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device, "REINITIALIZATION_RUN");

  
  // Compute Gradient of Phi0 = Phi and Compute SgnPhi
  mnsComputeSignumTerm(mns, haloBytes, sendBuffer, recvBuffer);

  for(int tstep=0;tstep<mns->NtimeSteps;++tstep){

    if((tstep%mns->outputStep)==0){
      occaTimerTic(mesh->device, "REINITIALIZATION_OUTPUT");
      dfloat time = mns->dt*tstep; 
      mnsReport(mns, time, tstep);
      occaTimerToc(mesh->device, "REINITIALIZATION_OUTPUT");
    }

    occaTimerTic(mesh->device, "REINITIALIZATION_STEP");  
    mnsReinitializationStep(mns, tstep, haloBytes, sendBuffer, recvBuffer);
    occaTimerToc(mesh->device, "REINITIALIZATION_STEP");  
  }


  printf("writing Final data\n");  
  occaTimerTic(mesh->device, "LEVEL_SET_OUTPUT");   
    mnsReport(mns, mns->finalTime, mns->NtimeSteps);
  occaTimerToc(mesh->device, "LEVEL_SET_OUTPUT");



  occaTimerToc(mesh->device, "REINITIALIZATION_RUN");

  occa::printTimer();
}

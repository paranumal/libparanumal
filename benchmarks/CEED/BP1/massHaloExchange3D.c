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

#include "massHex3D.h"

void massStartHaloExchange3D(solver_t *solver, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer){

  mesh3D *mesh = solver->mesh;
  
  // count size of halo for this process
  int haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  int haloOffset = mesh->Nelements*mesh->Np*sizeof(dfloat);
  
  // extract halo on DEVICE
  if(haloBytes){

    // make sure compute device is ready to perform halo extract
    //    mesh->device.finish();

    // switch to data stream
    //    mesh->device.setStream(solver->dataStream);

    // extract halo on data stream
    mesh->haloExtractKernel(mesh->totalHaloPairs, mesh->Np, mesh->o_haloElementList,
			    o_q, mesh->o_haloBuffer);

    // queue up async copy of halo on data stream
    mesh->o_haloBuffer.copyTo(sendBuffer);
    
    //    mesh->device.setStream(solver->defaultStream);
  }
}

void massInterimHaloExchange3D(solver_t *solver, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer){

  mesh3D *mesh = solver->mesh;

  // count size of halo for this process
  int haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  int haloOffset = mesh->Nelements*mesh->Np*sizeof(dfloat);
  
  // extract halo on DEVICE
  if(haloBytes){
    
    // copy extracted halo to HOST
    //    mesh->device.setStream(solver->dataStream);

    // make sure async copy finished
    //    mesh->device.finish(); 
    
    // start halo exchange HOST<>HOST
    meshHaloExchangeStart(mesh,
			  mesh->Np*sizeof(dfloat),
			  sendBuffer,
			  recvBuffer);
    
    //    mesh->device.setStream(solver->defaultStream);

  }
}
    

void massEndHaloExchange3D(solver_t *solver, occa::memory &o_q, dfloat *recvBuffer){

  mesh3D *mesh = solver->mesh;
  
  // count size of halo for this process
  int haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  int haloOffset = mesh->Nelements*mesh->Np*sizeof(dfloat);
  
  // extract halo on DEVICE
  if(haloBytes){
    // finalize recv on HOST
    meshHaloExchangeFinish(mesh);
    
    // copy into halo zone of o_r  HOST>DEVICE
    //    mesh->device.setStream(solver->dataStream);
    o_q.copyFrom(recvBuffer, haloBytes, haloOffset);
    
    //    mesh->device.setStream(solver->defaultStream);
  }
}

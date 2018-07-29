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

#include "ellipticHex3D.h"

void ellipticParallelGatherScatter(mesh3D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op){

  // use gather map for gather and scatter
  occaTimerTic(mesh->device,"meshParallelGatherScatter3D");
  meshParallelGatherScatter(mesh, ogs, o_q, o_gsq, type, op);
  occaTimerToc(mesh->device,"meshParallelGatherScatter3D");
}


void ellipticHaloGatherScatter(solver_t *solver, 
			       ogs_t *halo, 
			       occa::memory &o_v,
			       const char *type,
			       const char *op){


  mesh3D *mesh = solver->mesh;

  if(halo->Ngather){
    // rough way to make sure 
    mesh->device.finish();
    
    // set stream to halo stream
    //    mesh->device.setStream(solver->dataStream);
    mesh->device.finish();
    
    // gather halo nodes on DEVICE
    mesh->gatherKernel(halo->Ngather, halo->o_gatherOffsets, halo->o_gatherLocalIds, o_v, halo->o_gatherTmp);

    mesh->device.finish();
    
    // copy partially gathered halo data from DEVICE to HOST
    halo->o_gatherTmp.asyncCopyTo(halo->gatherTmp);

    mesh->device.finish();
    
    // wait for async copy
    occa::streamTag tag = mesh->device.tagStream();
    mesh->device.waitFor(tag);
    
    // gather across MPI processes then scatter back
    gsParallelGatherScatter(halo->gatherGsh, halo->gatherTmp, dfloatString, op); // danger on hardwired type

    mesh->device.finish();
    
    // copy totally gather halo data back from HOST to DEVICE
    halo->o_gatherTmp.asyncCopyFrom(halo->gatherTmp); 

    mesh->device.finish();
    
    tag = mesh->device.tagStream();
    mesh->device.waitFor(tag);
    
    // do scatter back to local nodes
    mesh->scatterKernel(halo->Ngather, halo->o_gatherOffsets, halo->o_gatherLocalIds, halo->o_gatherTmp, o_v);
    mesh->device.finish();
    
    // revert back to default stream
    mesh->device.setStream(solver->defaultStream);
  }
  
}


void ellipticNonHaloGatherScatter(solver_t *solver, 
				  ogs_t *nonHalo, 
				  occa::memory &o_v,
				  const char *type,
				  const char *op){


  mesh3D *mesh = solver->mesh;

  // set stream to default stream
  //  mesh->device.setStream(solver->defaultStream);
  
  // unified gather-scatter operation on DEVICE for non-halo nodes
  mesh->gatherScatterKernel(nonHalo->Ngather, nonHalo->o_gatherOffsets, nonHalo->o_gatherLocalIds, o_v);
  
}

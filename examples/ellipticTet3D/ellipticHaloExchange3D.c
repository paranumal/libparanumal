#include "ellipticTet3D.h"

void ellipticStartHaloExchange3D(solver_t *solver, occa::memory &o_q, int Nentries, dfloat *sendBuffer, dfloat *recvBuffer){

  mesh3D *mesh = solver->mesh;
  
  // count size of halo for this process
  int haloBytes = mesh->totalHaloPairs*Nentries*sizeof(dfloat);
  int haloOffset = mesh->Nelements*Nentries*sizeof(dfloat);
  
  // extract halo on DEVICE
  if(haloBytes){

    // make sure compute device is ready to perform halo extract
    mesh->device.finish();

    // switch to data stream
    mesh->device.setStream(solver->dataStream);

    // extract halo on data stream
    mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList,
          o_q, mesh->o_haloBuffer);

    // queue up async copy of halo on data stream
    mesh->o_haloBuffer.asyncCopyTo(sendBuffer);

    mesh->device.setStream(solver->defaultStream);
  }
}

void ellipticInterimHaloExchange3D(solver_t *solver, occa::memory &o_q, int Nentries, dfloat *sendBuffer, dfloat *recvBuffer){

  mesh3D *mesh = solver->mesh;

  // count size of halo for this process
  int haloBytes = mesh->totalHaloPairs*Nentries*sizeof(dfloat);
  int haloOffset = mesh->Nelements*Nentries*sizeof(dfloat);
  
  // extract halo on DEVICE
  if(haloBytes){
    
    // copy extracted halo to HOST
    mesh->device.setStream(solver->dataStream);

    // make sure async copy finished
    mesh->device.finish(); 
    
    // start halo exchange HOST<>HOST
    meshHaloExchangeStart(mesh,
        Nentries*sizeof(dfloat),
        sendBuffer,
        recvBuffer);
    
    mesh->device.setStream(solver->defaultStream);

  }
}
    

void ellipticEndHaloExchange3D(solver_t *solver, occa::memory &o_q, int Nentries, dfloat *recvBuffer){

  mesh3D *mesh = solver->mesh;
  
  // count size of halo for this process
  int haloBytes = mesh->totalHaloPairs*Nentries*sizeof(dfloat);
  int haloOffset = mesh->Nelements*Nentries*sizeof(dfloat);
  
  // extract halo on DEVICE
  if(haloBytes){
    // finalize recv on HOST
    meshHaloExchangeFinish(mesh);
    
    // copy into halo zone of o_r  HOST>DEVICE
    mesh->device.setStream(solver->dataStream);
    o_q.asyncCopyFrom(recvBuffer, haloBytes, haloOffset);
    mesh->device.finish();
    
    mesh->device.setStream(solver->defaultStream);
    mesh->device.finish();
  }
}


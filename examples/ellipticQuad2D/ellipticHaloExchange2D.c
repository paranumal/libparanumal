#include "ellipticQuad2D.h"

void ellipticStartHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer){

  // count size of halo for this process
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  iint haloOffset = mesh->Nelements*mesh->Np*sizeof(dfloat);
  
  // extract halo on DEVICE
  if(haloBytes){
    
    // WARNING: uses dfloats
    mesh->haloExtractKernel(mesh->totalHaloPairs, mesh->Np, mesh->o_haloElementList,
			    o_q, mesh->o_haloBuffer);
    
    // copy extracted halo to HOST 
    mesh->o_haloBuffer.copyTo(sendBuffer);
    
    // start halo exchange HOST<>HOST
    meshHaloExchangeStart(mesh,
			  mesh->Np*sizeof(dfloat),
			  sendBuffer,
			  recvBuffer);
  }
}
    

void ellipticEndHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *recvBuffer){
  
  // count size of halo for this process
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  iint haloOffset = mesh->Nelements*mesh->Np*sizeof(dfloat);
  
  // extract halo on DEVICE
  if(haloBytes){
    // finalize recv on HOST
    meshHaloExchangeFinish(mesh);
    
    // copy into halo zone of o_r  HOST>DEVICE
    o_q.copyFrom(recvBuffer, haloBytes, haloOffset);
  }
}
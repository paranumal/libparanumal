#include "mns.h"

// complete a time step using LSERK4
void mnsComputeSignumTerm(mns_t *mns, int haloBytes, dfloat * sendBuffer, dfloat *recvBuffer){

  mesh_t *mesh = mns->mesh; 

  dfloat time = 0.0; 

  if(mesh->totalHaloPairs>0){
    #if ASYNC 
      mesh->device.setStream(dataStream);
    #endif

    mesh->haloExtractKernel(mesh->totalHaloPairs,
                            mesh->Np,
                            mesh->o_haloElementList,
                            mns->o_Phi,
                            mesh->o_haloBuffer);

    // copy extracted halo to HOST
    mesh->o_haloBuffer.asyncCopyTo(sendBuffer);

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }


  // Compute Gradient volume Kernel
  occaTimerTic(mesh->device,"GradientVolume");
  // Compute Volume Contribution
  mns->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_Dmatrices,
                            mns->fieldOffset,
                            mns->o_Phi,
                            mns->o_GPhi);
  occaTimerToc(mesh->device,"GradientVolume");



  if(mesh->totalHaloPairs>0){

    #if ASYNC 
      mesh->device.setStream(dataStream);
    #endif

    //make sure the async copy is finished
    mesh->device.finish();
    // start halo exchange
    meshHaloExchangeStart(mesh,
                          mesh->Np*sizeof(dfloat),
                          sendBuffer,
                          recvBuffer);
    // wait for halo data to arrive
    meshHaloExchangeFinish(mesh);
    // copy halo data to DEVICE
    size_t offset = mesh->Np*mesh->Nelements*sizeof(dfloat); // offset for halo data
    mns->o_Phi.asyncCopyFrom(recvBuffer, haloBytes, offset);
    mesh->device.finish();        

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }


  occaTimerTic(mesh->device,"GradientSurface");
  // Compute Surface Conribution
  mns->gradientSurfaceKernel(mesh->Nelements,
                            mesh->o_sgeo,
                            mesh->o_LIFTT,
                            mesh->o_vmapM,
                            mesh->o_vmapP,
                            mesh->o_EToB,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            time,
                            mns->fieldOffset,
                            mns->o_Phi,
                            mns->o_GPhi);
  occaTimerToc(mesh->device,"GradientSurface");


 
  occaTimerTic(mesh->device,"RegularizedSignum");
  mns->regularizedSignumKernel(mesh->Nelements,
                              mns->hmin,
                              mns->fieldOffset,
                              mns->o_Phi,
                              mns->o_GPhi,
                              mns->o_SPhi);

  occaTimerToc(mesh->device,"RegularizedSignum");



}
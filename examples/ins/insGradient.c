#include "ins.h"

// complete a time step using LSERK4
void insGradient(ins_t *ins, dfloat time, occa::memory o_P, occa::memory o_GP){

  mesh_t *mesh = ins->mesh;
  
  if (ins->pOptions.compareArgs("DISCRETIZATION","IPDG")) {
    if(mesh->totalHaloPairs>0){
      ins->pressureHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 o_P,
                                 o_pHaloBuffer);

      // copy extracted halo to HOST
      o_pHaloBuffer.copyTo(ins->pSendBuffer);

      // start halo exchange
      meshHaloExchangeStart(mesh,
                           mesh->Np*sizeof(dfloat),
                           ins->pSendBuffer,
                           ins->pRecvBuffer);
    }
  }

  occaTimerTic(mesh->device,"GradientVolume");
  // Compute Volume Contribution
  ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_Dmatrices,
                            ins->fieldOffset,
                            o_P,
                            o_GP);
  occaTimerToc(mesh->device,"GradientVolume");

  // COMPLETE HALO EXCHANGE
  if (ins->pOptions.compareArgs("DISCRETIZATION","IPDG")) {
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);

      ins->o_pHaloBuffer.copyFrom(ins->pRecvBuffer);

      ins->pressureHaloScatterKernel(mesh->Nelements,
                                      mesh->totalHaloPairs,
                                      o_P,
                                      ins->o_pHaloBuffer);
    }

    occaTimerTic(mesh->device,"GradientSurface");
    // Compute Surface Conribution
    ins->gradientSurfaceKernel(mesh->Nelements,
                               mesh->o_sgeo,
                               mesh->o_LIFTT,
                               mesh->o_vmapM,
                               mesh->o_vmapP,
                               mesh->o_EToB,
                               mesh->o_x,
                               mesh->o_y,
                               mesh->o_z,
                               time,
                               ins->fieldOffset,
                               o_P,
                               o_GP);
    occaTimerToc(mesh->device,"GradientSurface");
  }
}

#include "ins.h"

// Compute divergence of the velocity field using physical boundary data at t = time. 
void insDivergence(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_DU){

  mesh_t *mesh = ins->mesh;

  //if (ins->vOptions.compareArgs("DISCRETIZATION","IPDG")) {
    if(mesh->totalHaloPairs>0){
      ins->velocityHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 ins->fieldOffset,
                                 o_U,
                                 ins->o_vHaloBuffer);

      // copy extracted halo to HOST 
      ins->o_vHaloBuffer.copyTo(ins->vSendBuffer);           
    
      // start halo exchange
      meshHaloExchangeStart(mesh,
                           mesh->Np*(ins->NVfields)*sizeof(dfloat),
                           ins->vSendBuffer,
                           ins->vRecvBuffer);
    }
  //}

  // computes div u^(n+1) volume term
  occaTimerTic(mesh->device,"DivergenceVolume");
  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_Dmatrices,
                             ins->fieldOffset,
                             o_U,
                             o_DU);
  occaTimerToc(mesh->device,"DivergenceVolume");

  //if (ins->vOptions.compareArgs("DISCRETIZATION","IPDG")) {
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);

      ins->o_vHaloBuffer.copyFrom(ins->vRecvBuffer); 

      ins->velocityHaloScatterKernel(mesh->Nelements,
                                    mesh->totalHaloPairs,
                                    ins->fieldOffset,
                                    o_U,
                                    ins->o_vHaloBuffer);
    }

    //computes div u^(n+1) surface term
    occaTimerTic(mesh->device,"DivergenceSurface");
    ins->divergenceSurfaceKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_LIFTT,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                time,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                ins->fieldOffset,
                                o_U,
                                o_DU);
    occaTimerToc(mesh->device,"DivergenceSurface");
  //}
}
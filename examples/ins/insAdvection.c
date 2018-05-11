#include "ins.h"

// compute NU = N(U)
void insAdvection(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_NU){

  mesh_t *mesh = ins->mesh;
  
  //Exctract Halo On Device, all fields
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

  // Compute Volume Contribution
  occaTimerTic(mesh->device,"AdvectionVolume");
  if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    ins->advectionCubatureVolumeKernel(mesh->Nelements,
                                       mesh->o_vgeo,
                                       mesh->o_cubvgeo,
                                       mesh->o_cubDWmatrices,
                                       mesh->o_cubInterpT,
                                       mesh->o_cubProjectT,
                                       ins->fieldOffset,
                                       o_U,
                                       ins->o_cU,
                                       o_NU);
  } else {
    ins->advectionVolumeKernel(mesh->Nelements,
                               mesh->o_vgeo,
                               mesh->o_Dmatrices,
                               ins->fieldOffset,
                               o_U,
                               o_NU);
  }
  occaTimerToc(mesh->device,"AdvectionVolume");

  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){
    meshHaloExchangeFinish(mesh);

    ins->o_vHaloBuffer.copyFrom(ins->vRecvBuffer); 

    ins->velocityHaloScatterKernel(mesh->Nelements,
                                  mesh->totalHaloPairs,
                                  ins->fieldOffset,
                                  o_U,
                                  ins->o_vHaloBuffer);
  }

  occaTimerTic(mesh->device,"AdvectionSurface");
  if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    ins->advectionCubatureSurfaceKernel(mesh->Nelements,
                                        mesh->o_vgeo,
                                        mesh->o_sgeo,
                                        mesh->o_cubsgeo,
                                        mesh->o_intInterpT,
                                        mesh->o_intLIFTT,
                                        mesh->o_cubInterpT,
                                        mesh->o_cubProjectT,
                                        mesh->o_vmapM,
                                        mesh->o_vmapP,
                                        mesh->o_EToB,
                                        time,
                                        mesh->o_intx,
                                        mesh->o_inty,
                                        mesh->o_intz,
                                        ins->fieldOffset,
                                        o_U,
                                        o_NU);
  } else {
    ins->advectionSurfaceKernel(mesh->Nelements,
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
                                o_NU);
  }
  occaTimerToc(mesh->device,"AdvectionSurface");
}

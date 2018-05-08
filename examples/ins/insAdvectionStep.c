#include "ins.h"

// complete a time step using LSERK4
void insAdvectionStep(ins_t *ins, dfloat time){

  mesh_t *mesh = ins->mesh;
  
  // field offset
  dlong  offset  = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;  

  //Exctract Halo On Device, all fields
  if(mesh->totalHaloPairs>0){
    ins->totalHaloExtractKernel(mesh->Nelements,
                                mesh->totalHaloPairs,
                                mesh->o_haloElementList,
                                offset,
                                ins->o_U,
                                ins->o_P,
                                ins->o_tHaloBuffer);

    // copy extracted halo to HOST
    ins->o_tHaloBuffer.copyTo(ins->tSendBuffer);

    // start halo exchange
    meshHaloExchangeStart(mesh,
                          mesh->Np*(ins->NTfields)*sizeof(dfloat),
                          ins->tSendBuffer,
                          ins->tRecvBuffer);
  }

  occaTimerTic(mesh->device,"AdvectionVolume");
  // Compute Volume Contribution
  if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    ins->advectionCubatureVolumeKernel(mesh->Nelements,
                                       mesh->o_vgeo,
                                       mesh->o_cubDWmatrices,
                                       mesh->o_cubInterpT,
                                       mesh->o_cubProjectT,
                                       offset,
                                       ins->o_U,
                                       ins->o_NU);
  } else {
    ins->advectionVolumeKernel(mesh->Nelements,
                               mesh->o_vgeo,
                               mesh->o_Dmatrices,
                               offset,
                               ins->o_U,
                               ins->o_NU);
  }

  occaTimerToc(mesh->device,"AdvectionVolume");


  occaTimerTic(mesh->device,"GradientVolume");
  // Compute Volume Contribution
  ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_Dmatrices,
                            offset,
                            ins->o_P,
                            ins->o_gradP);
  occaTimerToc(mesh->device,"GradientVolume");

  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){

    meshHaloExchangeFinish(mesh);

    ins->o_tHaloBuffer.copyFrom(ins->tRecvBuffer);
    ins->totalHaloScatterKernel(mesh->Nelements,
                                mesh->totalHaloPairs,
                                mesh->o_haloElementList,
                                offset,
                                ins->o_U,
                                ins->o_P,
                                ins->o_tHaloBuffer);
  }

  occaTimerTic(mesh->device,"AdvectionSurface");
  if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    ins->advectionCubatureSurfaceKernel(mesh->Nelements,
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
                                        offset,
                                        ins->o_U,
                                        ins->o_NU);
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
                                offset,
                                ins->o_U,
                                ins->o_NU);
  }

  occaTimerToc(mesh->device,"AdvectionSurface");
  // Solve pressure gradient for time^(n+1) grad(p^(n+1))
  time += ins->dt;
  

  if (ins->pOptions.compareArgs("DISCRETIZATION","IPDG")) {
    const int solverid = 0; // Pressure Solve, use BCs for pressure

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
                               offset,
                               solverid, // pressure BCs
                               ins->o_P,
                               ins->o_gradP);
    occaTimerToc(mesh->device,"GradientSurface");
  }
}

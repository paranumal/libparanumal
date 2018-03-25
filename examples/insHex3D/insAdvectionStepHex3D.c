#include "insHex3D.h"

// complete a time step using LSERK4
void insAdvectionStepHex3D(ins_t *ins, int tstep, char *options){

  mesh3D *mesh = ins->mesh;
  dfloat t = (tstep+0)*ins->dt;  // to compute N(U^n) set t=tn 
  
  // field ioffset at this step
  dlong  offset  = mesh->Nelements+mesh->totalHaloPairs;  
  dlong ioffset  = ins->index*offset;

  //Exctract Halo On Device, all fields
  if(mesh->totalHaloPairs>0){
    ins->totalHaloExtractKernel(mesh->Nelements,
                                mesh->totalHaloPairs,
                                mesh->o_haloElementList,
                                ioffset,
                                ins->o_U,
                                ins->o_V,
                                ins->o_W,
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
  if(strstr(options, "CUBATURE")){
    ins->advectionCubatureVolumeKernel(mesh->Nelements,
                                       mesh->o_vgeo,
                                       mesh->o_cubvgeo,
                                       mesh->o_cubDWT,
                                       mesh->o_cubInterpT,
                                       mesh->o_cubProjectT,
                                       ioffset,
                                       ins->o_U,
                                       ins->o_V,
                                       ins->o_W,
                                       ins->o_NU,
                                       ins->o_NV,
                                       ins->o_NW);
  } else {
    ins->advectionVolumeKernel(mesh->Nelements,
                               mesh->o_vgeo,
                               mesh->o_D,
                               ioffset,
                               ins->o_U,
                               ins->o_V,
                               ins->o_W,
                               ins->o_NU,
                               ins->o_NV,
                               ins->o_NW);
  }

  occaTimerToc(mesh->device,"AdvectionVolume");


  occaTimerTic(mesh->device,"GradientVolume");
  // Compute Volume Contribution
  ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_D,
                            ioffset,
                            ins->o_P,
                            ins->o_Px,
                            ins->o_Py,
                            ins->o_Pz);
  occaTimerToc(mesh->device,"GradientVolume");

  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){

    meshHaloExchangeFinish(mesh);

    ins->o_tHaloBuffer.copyFrom(ins->tRecvBuffer);
    ins->totalHaloScatterKernel(mesh->Nelements,
                                mesh->totalHaloPairs,
                                mesh->o_haloElementList,
                                ioffset,
                                ins->o_U,
                                ins->o_V,
                                ins->o_W,
                                ins->o_P,
                                ins->o_tHaloBuffer);
  }

  occaTimerTic(mesh->device,"AdvectionSurface");
  if(strstr(options, "CUBATURE")){
    ins->advectionCubatureSurfaceKernel(mesh->Nelements,
                                        mesh->o_vgeo,
                                        mesh->o_cubsgeo,
                                        mesh->o_vmapM,
                                        mesh->o_vmapP,
                                        mesh->o_EToB,
                                        mesh->o_cubInterpT,
                                        mesh->o_cubProjectT,
                                        t,
                                        mesh->o_intx,
                                        mesh->o_inty,
                                        mesh->o_intz,
                                        ioffset,
                                        ins->o_U,
                                        ins->o_V,
                                        ins->o_W,
                                        ins->o_NU,
                                        ins->o_NV,
                                        ins->o_NW);
  } else {
    ins->advectionSurfaceKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                ioffset,
                                ins->o_U,
                                ins->o_V,
                                ins->o_W,
                                ins->o_NU,
                                ins->o_NV,
                                ins->o_NW);
  }

  occaTimerToc(mesh->device,"AdvectionSurface");
  // Solve pressure gradient for t^(n+1) grad(p^(n+1))
   t += ins->dt;
  

  if (strstr(ins->pSolverOptions,"IPDG")) {
    const int solverid = 0; // Pressure Solve, use BCs for pressure

    occaTimerTic(mesh->device,"GradientSurface");
    // Compute Surface Conribution
    ins->gradientSurfaceKernel(mesh->Nelements,
                               mesh->o_sgeo,
                               mesh->o_vmapM,
                               mesh->o_vmapP,
                               mesh->o_EToB,
                               mesh->o_x,
                               mesh->o_y,
                               mesh->o_z,
                               t,
                               ins->dt,
                               ins->c0,
                               ins->c1,
                               ins->c2,
                               ins->index,
                               offset,
                               solverid, // pressure BCs
                               ins->o_PI, //not used
                               ins->o_P,
                               ins->o_Px,
                               ins->o_Py,
                               ins->o_Pz);
    occaTimerToc(mesh->device,"GradientSurface");
  }
}

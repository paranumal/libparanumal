#include "ins2D.h"

// complete a time step using LSERK4
void insAdvectionStepSS2D(ins_t *ins, iint tstep,  iint haloBytes,
      dfloat * sendBuffer, dfloat * recvBuffer,
      char   * options){

  mesh2D *mesh = ins->mesh;
  dfloat t = (tstep+0)*ins->dt;
  
  const iint Ntotal = (mesh->Nelements+mesh->totalHaloPairs); 
  iint offset       = ins->index*Ntotal;

  //Exctract Halo On Device
  if(mesh->totalHaloPairs>0){
    ins->velocityHaloExtractKernel(mesh->Nelements,
                                  mesh->totalHaloPairs,
                                  mesh->o_haloElementList,
                                  offset,
                                  ins->o_U,
                                  ins->o_V,
                                  ins->o_vHaloBuffer);

    // copy extracted halo to HOST
    ins->o_vHaloBuffer.copyTo(sendBuffer);

    // start halo exchange
    meshHaloExchangeStart(mesh,
                          mesh->Np*(ins->NVfields)*sizeof(dfloat),
                          sendBuffer,
                          recvBuffer);
  }

  // Compute Volume Contribution
  if(strstr(options, "CUBATURE")){
    ins->advectionCubatureVolumeKernel(mesh->Nelements,
                                      mesh->o_vgeo,
                                      mesh->o_cubDrWT,
                                      mesh->o_cubDsWT,
                                      mesh->o_cubInterpT,
                                      offset,
                                      ins->o_U,
                                      ins->o_V,
                                      ins->o_NU,
                                      ins->o_NV);
  } else {
    ins->advectionVolumeKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_DrT,
                              mesh->o_DsT,
                              offset,
                              ins->o_U,
                              ins->o_V,
                              ins->o_NU,
                              ins->o_NV);
  }

  
  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){

    meshHaloExchangeFinish(mesh);

    ins->o_vHaloBuffer.copyFrom(recvBuffer);

    ins->velocityHaloScatterKernel(mesh->Nelements,
                                  mesh->totalHaloPairs,
                                  mesh->o_haloElementList,
                                  offset,
                                  ins->o_U,
                                  ins->o_V,
                                  ins->o_vHaloBuffer);
  }


  if(strstr(options, "CUBATURE")){
    ins->advectionCubatureSurfaceKernel(mesh->Nelements,
                                        mesh->o_sgeo,
                                        mesh->o_intInterpT,
                                        mesh->o_intLIFTT,
                                        mesh->o_vmapM,
                                        mesh->o_vmapP,
                                        mesh->o_EToB,
                                        t,
                                        mesh->o_intx,
                                        mesh->o_inty,
                                        offset,
                                        ins->o_U,
                                        ins->o_V,
                                        ins->o_NU,
                                        ins->o_NV);
  } else {
    ins->advectionSurfaceKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_LIFTT,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                offset,
                                ins->o_U,
                                ins->o_V,
                                ins->o_NU,
                                ins->o_NV);
  }



  ins->advectionUpdateKernel(mesh->Nelements,
                            ins->index,
                            Ntotal,
                            ins->dt,
                            ins->g0,
                            ins->a0,
                            ins->a1,
                            ins->a2,
                            ins->b0,
                            ins->b1,
                            ins->b2,
                            ins->o_U,
                            ins->o_V,
                            ins->o_NU,
                            ins->o_NV,
                            ins->o_Ut,
                            ins->o_Vt);
}

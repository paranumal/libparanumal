#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannSAABStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
                        dfloat * sendBuffer, dfloat *recvBuffer, char * options){



  iint mrab_order = 2;  // Third order

  if(tstep==0)
    mrab_order = 0; // first order
  else if(tstep==1)
    mrab_order = 1; // second order

  // intermediate stage time
  dfloat t = mesh->startTime + (tstep+1)*mesh->dt;

  // COMPUTE RAMP FUNCTION 
  dfloat ramp, drampdt;
  boltzmannRampFunction2D(t, &ramp, &drampdt);


  if(mesh->totalHaloPairs>0){
    // extract halo on DEVICE
    #if ASYNC 
      mesh->device.setStream(dataStream);
    #endif

    iint Nentries = mesh->Np*mesh->Nfields;
    mesh->haloExtractKernel(mesh->totalHaloPairs,
                            Nentries,
                            mesh->o_haloElementList,
                            mesh->o_q,
                            mesh->o_haloBuffer);

    // copy extracted halo to HOST
    mesh->o_haloBuffer.asyncCopyTo(sendBuffer);

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }


  if (mesh->nonPmlNelements)
    mesh->volumeKernel(mesh->nonPmlNelements,
                      mesh->o_nonPmlElementIds,
                      ramp, 
                      drampdt,
                      mesh->Nrhs,
                      mesh->shiftIndex,
                      mesh->o_vgeo,
                      mesh->o_DrT,
                      mesh->o_DsT,
                      mesh->o_q,
                      mesh->o_rhsq);

  if (mesh->pmlNelements)
    mesh->pmlVolumeKernel(mesh->pmlNelements,
                          mesh->o_pmlElementIds,
                          mesh->o_pmlIds,
                          ramp, 
                          drampdt,
                          mesh->Nrhs,
                          mesh->shiftIndex,
                          mesh->o_vgeo,
                          mesh->o_pmlSigmaX,
                          mesh->o_pmlSigmaY,
                          mesh->o_DrT,
                          mesh->o_DsT,
                          mesh->o_q,
                          mesh->o_pmlqx,
                          mesh->o_pmlqy,
                          mesh->o_rhsq,
                          mesh->o_pmlrhsqx,
                          mesh->o_pmlrhsqy);


  if(strstr(options, "CUBATURE")){ 

    if (mesh->nonPmlNelements)
      mesh->relaxationKernel(mesh->nonPmlNelements,
                            mesh->o_nonPmlElementIds,
                            mesh->Nrhs,
                            mesh->shiftIndex,
                            mesh->o_cubInterpT,
                            mesh->o_cubProjectT,
                            mesh->o_q,
                            mesh->o_rhsq);  

    if (mesh->pmlNelements)  
      mesh->pmlRelaxationKernel(mesh->pmlNelements,
                                mesh->o_pmlElementIds,
                                mesh->o_pmlIds,
                                mesh->Nrhs,
                                mesh->shiftIndex,
                                mesh->o_cubInterpT,
                                mesh->o_cubProjectT,
                                mesh->o_pmlSigmaX,
                                mesh->o_pmlSigmaY,
                                mesh->o_q,
                                mesh->o_pmlqx,
                                mesh->o_pmlqy,
                                mesh->o_rhsq,
                                mesh->o_pmlrhsqx,
                                mesh->o_pmlrhsqy);

      // mesh->pmlRelaxationKernel(mesh->pmlNelements,
      //                           mesh->o_pmlElementIds,
      //                           mesh->o_pmlIds,
      //                           mesh->Nrhs,
      //                           mesh->shiftIndex,
      //                           mesh->o_cubInterpT,
      //                           mesh->o_cubProjectT,
      //                           mesh->o_q,
      //                           mesh->o_rhsq);
  }



  if(mesh->totalHaloPairs>0){
    #if ASYNC 
      mesh->device.setStream(dataStream);
    #endif

    //make sure the async copy is finished
    mesh->device.finish();

    // start halo exchange
    meshHaloExchangeStart(mesh,
                          mesh->Nfields*mesh->Np*sizeof(dfloat),
                          sendBuffer,
                          recvBuffer);

    // wait for halo data to arrive
    meshHaloExchangeFinish(mesh);


    // copy halo data to DEVICE
    size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
    mesh->o_q.asyncCopyFrom(recvBuffer, haloBytes, offset);
    mesh->device.finish();        

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }


  if (mesh->nonPmlNelements)
    mesh->surfaceKernel(mesh->nonPmlNelements,
                        mesh->o_nonPmlElementIds,
                        t,
                        ramp,
                        mesh->Nrhs,
                        mesh->shiftIndex,
                        mesh->o_sgeo,
                        mesh->o_LIFTT,
                        mesh->o_vmapM,
                        mesh->o_vmapP,
                        mesh->o_EToB,
                        mesh->o_x,
                        mesh->o_y,
                        mesh->o_q,
                        mesh->o_rhsq);

  if (mesh->pmlNelements)
    mesh->pmlSurfaceKernel(mesh->pmlNelements,
                          mesh->o_pmlElementIds,
                          mesh->o_pmlIds,
                          t,
                          ramp,
                          mesh->Nrhs,
                          mesh->shiftIndex,
                          mesh->o_sgeo,
                          mesh->o_LIFTT,
                          mesh->o_vmapM,
                          mesh->o_vmapP,
                          mesh->o_EToB,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_q,
                          mesh->o_rhsq,
                          mesh->o_pmlrhsqx,
                          mesh->o_pmlrhsqy);

  const iint id = mrab_order*3;

  if (mesh->nonPmlNelements)
    mesh->updateKernel(mesh->nonPmlNelements,
                      mesh->o_nonPmlElementIds,
                      mesh->MRSAAB_C[0], //
                      mesh->MRAB_A[id+0], //
                      mesh->MRAB_A[id+1],
                      mesh->MRAB_A[id+2], //
                      mesh->MRSAAB_A[id+0], //
                      mesh->MRSAAB_A[id+1],
                      mesh->MRSAAB_A[id+2], 
                      mesh->shiftIndex,
                      mesh->o_rhsq,
                      mesh->o_q);

  if (mesh->pmlNelements)
    mesh->pmlUpdateKernel(mesh->pmlNelements,
                          mesh->o_pmlElementIds,
                          mesh->o_pmlIds,
                          mesh->MRSAAB_C[0], //
                          mesh->MRAB_A[id+0], //
                          mesh->MRAB_A[id+1],
                          mesh->MRAB_A[id+2], //
                          mesh->MRSAAB_A[id+0], //
                          mesh->MRSAAB_A[id+1],
                          mesh->MRSAAB_A[id+2], 
                          mesh->shiftIndex,
                          mesh->o_rhsq,
                          mesh->o_pmlrhsqx,
                          mesh->o_pmlrhsqy,
                          mesh->o_q,
                          mesh->o_pmlqx,
                          mesh->o_pmlqy);

  //rotate index
  mesh->shiftIndex = (mesh->shiftIndex+1)%3;
}

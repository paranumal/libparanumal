#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannSRABStep2D(bns_t *bns, int tstep, int haloBytes,
                        dfloat * sendBuffer, dfloat *recvBuffer, setupAide &options){


  int mrab_order = 2;  // Third order
  dfloat a1, a2, a3;
  if (tstep==0) {
    mrab_order = 0;  
    a1 = 1.0*bns->dt;
    a2 = 0.0;
    a3 = 0.0;
  } else if (tstep ==1) {
    mrab_order = 1;  
    a1 =  3.0/2.0*bns->dt;
    a2 = -1.0/2.0*bns->dt;
    a3 =  0.0;
  } else {
  a1 =  23./12.0*bns->dt;
  a2 = -16./12.0*bns->dt;
  a3 =   5./12.0*bns->dt;
  }
  
  mesh2D *mesh = bns->mesh; 
  // intermediate stage time
  dfloat t = bns->startTime + (tstep+1)*bns->dt;

  // COMPUTE RAMP FUNCTION 
  dfloat ramp, drampdt;
  boltzmannRampFunction2D(t, &ramp, &drampdt);


  if(mesh->totalHaloPairs>0){
    // extract halo on DEVICE
    #if ASYNC 
      mesh->device.setStream(dataStream);
    #endif

    int Nentries = mesh->Np*bns->Nfields;
    mesh->haloExtractKernel(mesh->totalHaloPairs,
                            Nentries,
                            mesh->o_haloElementList,
                            bns->o_q,
                            mesh->o_haloBuffer);

    // copy extracted halo to HOST
    mesh->o_haloBuffer.asyncCopyTo(sendBuffer);

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }


  if (mesh->nonPmlNelements)  
      bns->volumeKernel(mesh->nonPmlNelements,
                        mesh->o_nonPmlElementIds,
                        ramp, 
                        drampdt,
                        bns->Nrhs,
                        bns->shiftIndex,
                        mesh->o_vgeo,
                        mesh->o_DrT,
                        mesh->o_DsT,
                        bns->o_q,
                        bns->o_rhsq);

  if (mesh->pmlNelements)
     bns->pmlVolumeKernel(mesh->pmlNelements,
                          mesh->o_pmlElementIds,
                          mesh->o_pmlIds,
                          ramp, 
                          drampdt,
                          bns->Nrhs,
                          bns->shiftIndex,
                          mesh->o_vgeo,
                          bns->o_pmlSigmaX,
                          bns->o_pmlSigmaY,
                          mesh->o_DrT,
                          mesh->o_DsT,
                          bns->o_q,
                          bns->o_pmlqx,
                          bns->o_pmlqy,
                          bns->o_rhsq,
                          bns->o_pmlrhsqx,
                          bns->o_pmlrhsqy);


  if(options.compareArgs("RELAXATION TYPE","CUBATURE")){ 

    if (mesh->nonPmlNelements)
        bns->relaxationKernel(mesh->nonPmlNelements,
                            mesh->o_nonPmlElementIds,
                            bns->Nrhs,
                            bns->shiftIndex,
                            mesh->o_cubInterpT,
                            mesh->o_cubProjectT,
                            bns->o_q,
                            bns->o_rhsq);  

    if (mesh->pmlNelements)    
      bns->pmlRelaxationKernel(mesh->pmlNelements,
                                mesh->o_pmlElementIds,
                                mesh->o_pmlIds,
                                bns->Nrhs,
                                bns->shiftIndex,
                                mesh->o_cubInterpT,
                                mesh->o_cubProjectT,
                                bns->o_pmlSigmaX,
                                bns->o_pmlSigmaY,
                                bns->o_q,
                                bns->o_pmlqx,
                                bns->o_pmlqy,
                                bns->o_rhsq,
                                bns->o_pmlrhsqx,
                                bns->o_pmlrhsqy);
  }



  if(mesh->totalHaloPairs>0){
    #if ASYNC 
      mesh->device.setStream(dataStream);
    #endif

    //make sure the async copy is finished
    mesh->device.finish();

    // start halo exchange
    meshHaloExchangeStart(mesh,
                          bns->Nfields*mesh->Np*sizeof(dfloat),
                          sendBuffer,
                          recvBuffer);

    // wait for halo data to arrive
    meshHaloExchangeFinish(mesh);


    // copy halo data to DEVICE
    size_t offset = mesh->Np*bns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
    bns->o_q.asyncCopyFrom(recvBuffer, haloBytes, offset);
    mesh->device.finish();        

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }


// printf("Starting Surface step: mesh Nrhs: %d  %d \n", bns->Nrhs, bns->shiftIndex);

  if (mesh->nonPmlNelements)
    bns->surfaceKernel(mesh->nonPmlNelements,
                        mesh->o_nonPmlElementIds,
                        t,
                        ramp,
                        bns->Nrhs,
                        bns->shiftIndex,
                        mesh->o_sgeo,
                        mesh->o_LIFTT,
                        mesh->o_vmapM,
                        mesh->o_vmapP,
                        mesh->o_EToB,
                        mesh->o_x,
                        mesh->o_y,
                        bns->o_q,
                        bns->o_rhsq);

  if (mesh->pmlNelements)
    bns->pmlSurfaceKernel(mesh->pmlNelements,
                          mesh->o_pmlElementIds,
                          mesh->o_pmlIds,
                          t,
                          ramp,
                          bns->Nrhs,
                          bns->shiftIndex,
                          mesh->o_sgeo,
                          mesh->o_LIFTT,
                          mesh->o_vmapM,
                          mesh->o_vmapP,
                          mesh->o_EToB,
                          mesh->o_x,
                          mesh->o_y,
                          bns->o_q,
                          bns->o_rhsq,
                          bns->o_pmlrhsqx,
                          bns->o_pmlrhsqy);

  if (mesh->nonPmlNelements)
    bns->updateKernel(mesh->nonPmlNelements,
                      mesh->o_nonPmlElementIds,
                      bns->dt,
                      a1,a2,a3,
                      bns->shiftIndex,
                      bns->o_rhsq,
                      bns->o_q);

  if (mesh->pmlNelements)
    bns->pmlUpdateKernel(mesh->pmlNelements,
                          mesh->o_pmlElementIds,
                          mesh->o_pmlIds,
                          bns->dt,
                          a1,a2,a3,
                          bns->shiftIndex,
                          bns->o_rhsq,
                          bns->o_pmlrhsqx,
                          bns->o_pmlrhsqy,
                          bns->o_q,
                          bns->o_pmlqx,
                          bns->o_pmlqy);

  //rotate index
  bns->shiftIndex = (bns->shiftIndex+1)%3;


// }


}

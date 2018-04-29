#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannSARKStep2D(bns_t *bns, int tstep, int haloBytes,
                  dfloat * sendBuffer, dfloat *recvBuffer,setupAide &options){


  mesh2D *mesh = bns->mesh; 
  const int shift_base  = 0; 
  bns->shiftIndex       = 0;

 
  for(int s=0; s<3; ++s){

    // Stage time
    dfloat t = bns->startTime+ tstep*bns->dt + bns->dt*bns->RK_C[s];
    //

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



    occaTimerTic(mesh->device, "VolumeKernel");

    if (mesh->nonPmlNelements){
      occaTimerTic(mesh->device,"NonPmlVolumeKernel"); 
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
      occaTimerToc(mesh->device,"NonPmlVolumeKernel"); 
  }

  if (mesh->pmlNelements){
    occaTimerTic(mesh->device,"PmlVolumeKernel"); 
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
    occaTimerToc(mesh->device,"PmlVolumeKernel"); 
  }

  occaTimerToc(mesh->device, "VolumeKernel");    



  if(options.compareArgs("RELAXATION TYPE","CUBATURE")){
    occaTimerTic(mesh->device, "RelaxationKernel");
          
    if (mesh->nonPmlNelements){
      occaTimerTic(mesh->device, "NonPmlRelaxationKernel");
      bns->relaxationKernel(mesh->nonPmlNelements,
                            mesh->o_nonPmlElementIds,
                            bns->Nrhs,
                            bns->shiftIndex,
                            mesh->o_cubInterpT,
                            mesh->o_cubProjectT,
                            bns->o_q,
                            bns->o_rhsq);
      occaTimerToc(mesh->device, "NonPmlRelaxationKernel");
      } 

    if (mesh->pmlNelements){
      occaTimerTic(mesh->device, "PmlRelaxationKernel"); 
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
      occaTimerToc(mesh->device, "PmlRelaxationKernel");
    }

    occaTimerToc(mesh->device, "RelaxationKernel");
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

  


  occaTimerTic(mesh->device,"SurfaceKernel");  
  if (mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NonPmlSurfaceKernel");
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
    occaTimerToc(mesh->device,"NonPmlSurfaceKernel");

  }

  if (mesh->pmlNelements){
    occaTimerTic(mesh->device,"PmlSurfaceKernel");
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
    occaTimerToc(mesh->device,"PmlSurfaceKernel");

  }

  occaTimerToc(mesh->device,"SurfaceKernel");


    //rotate index
    bns->shiftIndex = (bns->shiftIndex+1)%3;

  occaTimerTic(mesh->device,"UpdateKernel");
 if(s<2){ 
      // Stage update
      if (mesh->pmlNelements){
        occaTimerTic(mesh->device,"PmlUpdateKernel");
        bns->pmlUpdateStageKernel(mesh->pmlNelements,
                  mesh->o_pmlElementIds,
                  mesh->o_pmlIds,
                  bns->dt,
                  bns->SARK_C[s],
                  bns->RK_A[s+1][0], 
                  bns->SARK_A[s+1][0], 
                  bns->RK_A[s+1][1], 
                  bns->SARK_A[s+1][1], 
                  ramp,
                  shift_base,
                  bns->o_rhsq,
                  bns->o_pmlrhsqx,
                  bns->o_pmlrhsqy,
                  bns->o_qS,
                  bns->o_qSx,
                  bns->o_qSy,
                  bns->o_pmlqx,
                  bns->o_pmlqy,
                  bns->o_q);
        occaTimerToc(mesh->device,"PmlUpdateKernel");
      }
        
      if(mesh->nonPmlNelements){
        occaTimerTic(mesh->device,"NonPmlUpdateKernel");
        bns->updateStageKernel(mesh->nonPmlNelements,
                  mesh->o_nonPmlElementIds,
                  bns->dt,
                  bns->SARK_C[s],
                  bns->RK_A[s+1][0], 
                  bns->SARK_A[s+1][0], 
                  bns->RK_A[s+1][1], 
                  bns->SARK_A[s+1][1],
                  shift_base,  
                  bns->o_rhsq,
                  bns->o_qS,
                  bns->o_q);
        occaTimerToc(mesh->device,"NonPmlUpdateKernel");
      }

    }

    else{
     // Final update s==2
     if (mesh->pmlNelements){
      occaTimerTic(mesh->device,"PmlUpdateKernel");
        bns->pmlUpdateKernel(mesh->pmlNelements,
                  mesh->o_pmlElementIds,
                  mesh->o_pmlIds,
                  bns->dt,
                  bns->SARK_C[s],
                  bns->RK_B[0], 
                  bns->SARK_B[0], 
                  bns->RK_B[1], 
                  bns->SARK_B[1], 
                  bns->RK_B[2], 
                  bns->SARK_B[2], 
                  ramp,
                  shift_base,
                  bns->o_rhsq,
                  bns->o_pmlrhsqx,
                  bns->o_pmlrhsqy,
                  bns->o_qS,
                  bns->o_qSx,
                  bns->o_qSy,
                  bns->o_pmlqx,
                  bns->o_pmlqy,
                  bns->o_q);
      occaTimerToc(mesh->device,"PmlUpdateKernel");
      }
        
      if(mesh->nonPmlNelements){
        occaTimerTic(mesh->device,"NonPmlUpdateKernel");
        bns->updateKernel(mesh->nonPmlNelements,
                  mesh->o_nonPmlElementIds,
                  bns->dt,
                  bns->SARK_C[s],
                  bns->RK_B[0], 
                  bns->SARK_B[0], 
                  bns->RK_B[1], 
                  bns->SARK_B[1], 
                  bns->RK_B[2], 
                  bns->SARK_B[2], 
                  shift_base,  
                  bns->o_rhsq,
                  bns->o_qS,
                  bns->o_q);
        occaTimerToc(mesh->device,"NonPmlUpdateKernel");
      }
    }

  occaTimerToc(mesh->device,"UpdateKernel");
        
  }   
}

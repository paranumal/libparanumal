#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannSARKStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
                  dfloat * sendBuffer, dfloat *recvBuffer,char * options){



  const iint shift_base      = 0; 
  mesh->shiftIndex = 0;

 
  for(iint s=0; s<3; ++s){

    // Stage time
    dfloat t = mesh->startTime+ tstep*mesh->dt + mesh->dt*mesh->RK_C[s];
    //

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


    occaTimerTic(mesh->device, "VolumeKernel");

    if (mesh->nonPmlNelements){
      occaTimerTic(mesh->device,"NonPmlVolumeKernel"); 
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
      occaTimerToc(mesh->device,"NonPmlVolumeKernel"); 
  }

  if (mesh->pmlNelements){
    occaTimerTic(mesh->device,"PmlVolumeKernel"); 
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
    occaTimerToc(mesh->device,"PmlVolumeKernel"); 
  }

  occaTimerToc(mesh->device, "VolumeKernel");    



  if(strstr(options, "CUBATURE")){
    occaTimerTic(mesh->device, "RelaxationKernel");
          
    if (mesh->nonPmlNelements){
      occaTimerTic(mesh->device, "NonPmlRelaxationKernel");
      mesh->relaxationKernel(mesh->nonPmlNelements,
                            mesh->o_nonPmlElementIds,
                            mesh->Nrhs,
                            mesh->shiftIndex,
                            mesh->o_cubInterpT,
                            mesh->o_cubProjectT,
                            mesh->o_q,
                            mesh->o_rhsq);
      occaTimerToc(mesh->device, "NonPmlRelaxationKernel");
      } 

    if (mesh->pmlNelements){
      occaTimerTic(mesh->device, "PmlRelaxationKernel"); 
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

  


  occaTimerTic(mesh->device,"SurfaceKernel");  
  if (mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NonPmlSurfaceKernel");
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
    occaTimerToc(mesh->device,"NonPmlSurfaceKernel");

  }

  if (mesh->pmlNelements){
    occaTimerTic(mesh->device,"PmlSurfaceKernel");
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
    occaTimerToc(mesh->device,"PmlSurfaceKernel");

  }

  occaTimerToc(mesh->device,"SurfaceKernel");


    //rotate index
    mesh->shiftIndex = (mesh->shiftIndex+1)%3;

  occaTimerTic(mesh->device,"UpdateKernel");
 if(s<2){ 
      // Stage update
      if (mesh->pmlNelements){
        occaTimerTic(mesh->device,"PmlUpdateKernel");
        mesh->pmlUpdateStageKernel(mesh->pmlNelements,
                  mesh->o_pmlElementIds,
                  mesh->o_pmlIds,
                  mesh->dt,
                  mesh->SARK_C[s],
                  mesh->RK_A[s+1][0], 
                  mesh->SARK_A[s+1][0], 
                  mesh->RK_A[s+1][1], 
                  mesh->SARK_A[s+1][1], 
                  ramp,
                  shift_base,
                  mesh->o_rhsq,
                  mesh->o_pmlrhsqx,
                  mesh->o_pmlrhsqy,
                  mesh->o_qS,
                  mesh->o_qSx,
                  mesh->o_qSy,
                  mesh->o_pmlqx,
                  mesh->o_pmlqy,
                  mesh->o_q);
        occaTimerToc(mesh->device,"PmlUpdateKernel");
      }
        
      if(mesh->nonPmlNelements){
        occaTimerTic(mesh->device,"NonPmlUpdateKernel");
        mesh->updateStageKernel(mesh->nonPmlNelements,
                  mesh->o_nonPmlElementIds,
                  mesh->dt,
                  mesh->SARK_C[s],
                  mesh->RK_A[s+1][0], 
                  mesh->SARK_A[s+1][0], 
                  mesh->RK_A[s+1][1], 
                  mesh->SARK_A[s+1][1],
                  shift_base,  
                  mesh->o_rhsq,
                  mesh->o_qS,
                  mesh->o_q);
        occaTimerToc(mesh->device,"NonPmlUpdateKernel");
      }

    }

    else{
     // Final update s==2
     if (mesh->pmlNelements){
      occaTimerTic(mesh->device,"PmlUpdateKernel");
        mesh->pmlUpdateKernel(mesh->pmlNelements,
                  mesh->o_pmlElementIds,
                  mesh->o_pmlIds,
                  mesh->dt,
                  mesh->SARK_C[s],
                  mesh->RK_B[0], 
                  mesh->SARK_B[0], 
                  mesh->RK_B[1], 
                  mesh->SARK_B[1], 
                  mesh->RK_B[2], 
                  mesh->SARK_B[2], 
                  ramp,
                  shift_base,
                  mesh->o_rhsq,
                  mesh->o_pmlrhsqx,
                  mesh->o_pmlrhsqy,
                  mesh->o_qS,
                  mesh->o_qSx,
                  mesh->o_qSy,
                  mesh->o_pmlqx,
                  mesh->o_pmlqy,
                  mesh->o_q);
      occaTimerToc(mesh->device,"PmlUpdateKernel");
      }
        
      if(mesh->nonPmlNelements){
        occaTimerTic(mesh->device,"NonPmlUpdateKernel");
        mesh->updateKernel(mesh->nonPmlNelements,
                  mesh->o_nonPmlElementIds,
                  mesh->dt,
                  mesh->SARK_C[s],
                  mesh->RK_B[0], 
                  mesh->SARK_B[0], 
                  mesh->RK_B[1], 
                  mesh->SARK_B[1], 
                  mesh->RK_B[2], 
                  mesh->SARK_B[2], 
                  shift_base,  
                  mesh->o_rhsq,
                  mesh->o_qS,
                  mesh->o_q);
        occaTimerToc(mesh->device,"NonPmlUpdateKernel");
      }
    }

  occaTimerToc(mesh->device,"UpdateKernel");
        
  }   
}

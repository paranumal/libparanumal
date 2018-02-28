#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannLSERKStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char * options){


   mesh->shiftIndex = 0; 	

// LSERK4 stages
for(iint rk=0;rk<mesh->Nrk;++rk){

  // intermediate stage time
  dfloat t = mesh->startTime + tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

  if(mesh->totalHaloPairs>0){
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

    // COMPUTE RAMP FUNCTION 
    dfloat ramp, drampdt;
    boltzmannRampFunction2D(t, &ramp, &drampdt);
      

    occaTimerTic(mesh->device, "VolumeKernel");    
    // compute volume contribution to DG boltzmann RHS
    if(mesh->pmlNelements){	
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

    // compute volume contribution to DG boltzmann RHS added d/dt (ramp(qbar)) to RHS
    if(mesh->nonPmlNelements){
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
  occaTimerToc(mesh->device, "VolumeKernel");    
    

	if(strstr(options, "CUBATURE")){ 
    occaTimerTic(mesh->device, "RelaxationKernel");
		if(mesh->pmlNelements){
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

		// compute relaxation terms using cubature
		if(mesh->nonPmlNelements){
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
		 // VOLUME KERNELS
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



  // SURFACE KERNELS
  occaTimerTic(mesh->device,"SurfaceKernel");

  if(mesh->pmlNelements){
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

  if(mesh->nonPmlNelements){
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
  occaTimerToc(mesh->device,"SurfaceKernel");

    
  // ramp function for flow at next RK stage
  dfloat tupdate = tstep*mesh->dt + mesh->dt*mesh->rkc[rk+1];
  dfloat rampUpdate, drampdtUpdate;
  boltzmannRampFunction2D(tupdate, &rampUpdate, &drampdtUpdate);

  //UPDATE
  occaTimerTic(mesh->device,"UpdateKernel");


  //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
  if (mesh->pmlNelements){   
    occaTimerTic(mesh->device,"PmlUpdateKernel");
    mesh->pmlUpdateKernel(mesh->pmlNelements,
                          mesh->o_pmlElementIds,
                          mesh->o_pmlIds,
                          mesh->dt,
                          mesh->rka[rk],
                          mesh->rkb[rk],
                          rampUpdate,
                          mesh->o_rhsq,
                          mesh->o_pmlrhsqx,
                          mesh->o_pmlrhsqy,
                          mesh->o_resq,
                          mesh->o_pmlresqx,
                          mesh->o_pmlresqy,
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
                      mesh->rka[rk],
                      mesh->rkb[rk],
                      mesh->o_rhsq,
                      mesh->o_resq,
                      mesh->o_q);
    occaTimerToc(mesh->device,"NonPmlUpdateKernel");
  }

 occaTimerToc(mesh->device,"UpdateKernel");
    
}


}

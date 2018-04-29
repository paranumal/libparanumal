#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannLSERKStep2D(bns_t *bns, int tstep, int haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, setupAide &options){


 bns->shiftIndex = 0; 
 mesh2D *mesh = bns->mesh; 	

// LSERK4 stages
for(int rk=0;rk<mesh->Nrk;++rk){

  // intermediate stage time
  dfloat t = bns->startTime + tstep*bns->dt + bns->dt*mesh->rkc[rk];

  if(mesh->totalHaloPairs>0){
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

    // COMPUTE RAMP FUNCTION 
    dfloat ramp, drampdt;
    boltzmannRampFunction2D(t, &ramp, &drampdt);
      

    occaTimerTic(mesh->device, "VolumeKernel");    
    // compute volume contribution to DG boltzmann RHS
    if(mesh->pmlNelements){	
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

    // compute volume contribution to DG boltzmann RHS added d/dt (ramp(qbar)) to RHS
    if(mesh->nonPmlNelements){
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
  occaTimerToc(mesh->device, "VolumeKernel");    
    

	if(options.compareArgs("RELAXATION TYPE","CUBATURE")){ 
    occaTimerTic(mesh->device, "RelaxationKernel");
		if(mesh->pmlNelements){
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

		// compute relaxation terms using cubature
		if(mesh->nonPmlNelements){
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



  // SURFACE KERNELS
  occaTimerTic(mesh->device,"SurfaceKernel");

  if(mesh->pmlNelements){
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

  if(mesh->nonPmlNelements){
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
  occaTimerToc(mesh->device,"SurfaceKernel");

    
  // ramp function for flow at next RK stage
  dfloat tupdate = tstep*bns->dt + bns->dt*mesh->rkc[rk+1];
  dfloat rampUpdate, drampdtUpdate;
  boltzmannRampFunction2D(tupdate, &rampUpdate, &drampdtUpdate);

  //UPDATE
  occaTimerTic(mesh->device,"UpdateKernel");


  //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
  if (mesh->pmlNelements){   
    occaTimerTic(mesh->device,"PmlUpdateKernel");
    bns->pmlUpdateKernel(mesh->pmlNelements,
                          mesh->o_pmlElementIds,
                          mesh->o_pmlIds,
                          bns->dt,
                          mesh->rka[rk],
                          mesh->rkb[rk],
                          rampUpdate,
                          bns->o_rhsq,
                          bns->o_pmlrhsqx,
                          bns->o_pmlrhsqy,
                          bns->o_resq,
                          bns->o_pmlresqx,
                          bns->o_pmlresqy,
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
                      mesh->rka[rk],
                      mesh->rkb[rk],
                      bns->o_rhsq,
                      bns->o_resq,
                      bns->o_q);
    occaTimerToc(mesh->device,"NonPmlUpdateKernel");
  }

 occaTimerToc(mesh->device,"UpdateKernel");
    
}


}

#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannSAADRKStep2D(bns_t *bns, dfloat time, int haloBytes,
				                   dfloat * sendBuffer, dfloat *recvBuffer, char * options){


 bns->shiftIndex = 0; 

 mesh2D *mesh = bns->mesh; 	

// LSERK4 stages
for(int rk=0;rk<bns->NrkStages;++rk){

  // intermediate stage time
  dfloat currentTime = time + bns->rkC[rk]*bns->dt;
  
   // printf("%.10e %.10e %.10e %.10e %.10e %.10e\n", expca, expmca, expma, bns->dt, bns->rkC[rk], bns->tauInv);


  
  occaTimerTic(mesh->device, "RKStageKernel");  
  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device, "NonPmlRKStageKernel");  
    bns->updateStageKernel(mesh->nonPmlNelements,
                           mesh->o_nonPmlElementIds,
                           mesh->Nelements,
                           rk,
                           bns->dt,
                           bns->o_sarkC,
                           bns->o_rkA,
                           bns->o_sarkA,
                           bns->o_q,
                           bns->o_rkrhsq,
                           bns->o_rkq);
    occaTimerToc(mesh->device, "NonPmlRKStageKernel");  
  }
  
  if(mesh->pmlNelements){
  occaTimerTic(mesh->device, "PmlRKStageKernel");  
  bns->pmlUpdateStageKernel(mesh->pmlNelements,
                        mesh->o_pmlElementIds,
                        mesh->o_pmlIds,
                        mesh->Nelements,
                        rk,
                        bns->dt,
                        bns->o_sarkC,
                        bns->o_rkA,
                        bns->o_sarkA,
                        bns->o_q,
                        bns->o_pmlqx,
                        bns->o_pmlqy,
                        bns->o_rkrhsq,
                        bns->o_rkrhsqx,
                        bns->o_rkrhsqy,
                        bns->o_rkq,
                        bns->o_rkqx,
                        bns->o_rkqy);
  occaTimerToc(mesh->device, "PmlRKStageKernel");
  }

  occaTimerToc(mesh->device, "RKStageKernel");  


  if(mesh->totalHaloPairs>0){
    #if ASYNC 
      mesh->device.setStream(dataStream);
    #endif

    int Nentries = mesh->Np*bns->Nfields;
    mesh->haloExtractKernel(mesh->totalHaloPairs,
                Nentries,
                mesh->o_haloElementList,
                bns->o_rkq,
                mesh->o_haloBuffer);

    // copy extracted halo to HOST
    mesh->o_haloBuffer.asyncCopyTo(sendBuffer);

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }

    // COMPUTE RAMP FUNCTION 
    dfloat ramp, drampdt;
    boltzmannRampFunction2D(currentTime, &ramp, &drampdt);
      

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
                        bns->o_rkq,
                        bns->o_rkqx,
                        bns->o_rkqy,
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
                          bns->o_rkq,
                          bns->o_rhsq);
      occaTimerToc(mesh->device,"NonPmlVolumeKernel");
  }
  occaTimerToc(mesh->device, "VolumeKernel");    
    

	if(strstr(options, "CUBATURE")){ 
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
                                bns->o_rkq,
                                bns->o_rkqx,
                                bns->o_rkqy,
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
                            bns->o_rkq,
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
    bns->o_rkq.asyncCopyFrom(recvBuffer, haloBytes, offset);
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
                          currentTime,
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
                          bns->o_rkq,
                          bns->o_rhsq,
                          bns->o_pmlrhsqx,
                          bns->o_pmlrhsqy);
    occaTimerToc(mesh->device,"PmlSurfaceKernel");
  }

  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NonPmlSurfaceKernel");
    bns->surfaceKernel(mesh->nonPmlNelements,
                        mesh->o_nonPmlElementIds,
                        currentTime,
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
                        bns->o_rkq,
                        bns->o_rhsq);
    occaTimerToc(mesh->device,"NonPmlSurfaceKernel");
  }
  occaTimerToc(mesh->device,"SurfaceKernel");

    
 //  //UPDATE
  occaTimerTic(mesh->device,"UpdateKernel");


  //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
  if (mesh->pmlNelements){   
    occaTimerTic(mesh->device,"PmlUpdateKernel");
    bns->pmlUpdateKernel(mesh->pmlNelements,
                        mesh->o_pmlElementIds,
                        mesh->o_pmlIds,
                        mesh->Nelements,
                        rk,
                        bns->dt,
                        bns->o_sarkC,
                        bns->o_rkA, 
                        bns->o_rkE,
                        bns->o_sarkA, 
                        bns->o_sarkE,
                        bns->o_q,
                        bns->o_pmlqx,
                        bns->o_pmlqy,
                        bns->o_rhsq,
                        bns->o_pmlrhsqx,
                        bns->o_pmlrhsqy,
                        bns->o_rkrhsq,
                        bns->o_rkrhsqx,
                        bns->o_rkrhsqy,
                        bns->o_rkq,
                        bns->o_rkqx,
                        bns->o_rkqy,
                        bns->o_rkerr);
    occaTimerToc(mesh->device,"PmlUpdateKernel");

  }

  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NonPmlUpdateKernel");
    bns->updateKernel(mesh->nonPmlNelements,
                      mesh->o_nonPmlElementIds,
                      mesh->Nelements,
                      rk,
                      bns->dt,
                      bns->o_sarkC,
                      bns->o_rkA, 
                      bns->o_rkE,
                      bns->o_sarkA, 
                      bns->o_sarkE,
                      bns->o_q,
                      bns->o_rhsq,
                      bns->o_rkrhsq,
                      bns->o_rkq,
                      bns->o_rkerr);
    occaTimerToc(mesh->device,"NonPmlUpdateKernel");
  }

  occaTimerToc(mesh->device,"UpdateKernel");
    
}


}

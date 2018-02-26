#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannLSIMEXStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				   dfloat * sendBuffer, dfloat *recvBuffer, char * options){


mesh->shiftIndex =0;
for(iint k=0;k<mesh->Nimex;++k){

  // intermediate stage time
  dfloat t = mesh->startTime+ tstep*mesh->dt + mesh->dt*mesh->LSIMEX_C[k];

  dfloat ramp, drampdt;
  boltzmannRampFunction2D(t, &ramp, &drampdt);


  // RESIDUAL UPDATE, i.e. Y = Q+ (a(k,k-1)-b(k-1))_ex *Y + (a(k,k-1)-b(k-1))_im *Z
  occaTimerTic(mesh->device,"ResidualUpdateKernel");
  // compute volume contribution to DG boltzmann RHS
  if(mesh->pmlNelements){

    occaTimerTic(mesh->device,"PMLResidualUpdateKernel");
    //	printf("pmlNel = %d\n", mesh->pmlNelements);
    mesh->pmlResidualUpdateKernel(mesh->pmlNelements,
                                  mesh->o_pmlElementIds,
                                  mesh->o_pmlIds,
                                  mesh->dt,
                                  ramp,
                                  mesh->LSIMEX_ABi[k],
                                  mesh->LSIMEX_ABe[k],
                                  mesh->o_q,
                                  mesh->o_pmlqx,
                                  mesh->o_pmlqy,
                                  mesh->o_qZ,
                                  mesh->o_qY,
                                  mesh->o_qYx,
                                  mesh->o_qYy);
    occaTimerToc(mesh->device,"PMLResidualUpdateKernel");
  }
  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NONPMLResidualUpdateKernel");
    mesh->residualUpdateKernel(mesh->nonPmlNelements,
                              mesh->o_nonPmlElementIds,
                              mesh->dt,
                              mesh->LSIMEX_ABi[k],
                              mesh->LSIMEX_ABe[k],
                              mesh->o_q,
                              mesh->o_qZ,
                              mesh->o_qY);
    occaTimerToc(mesh->device,"NONPMLResidualUpdateKernel");
  }
  
  occaTimerToc(mesh->device,"ResidualUpdateKernel");

      

  occaTimerTic(mesh->device,"PMLImplicit");
  // Compute Implicit Part of Boltzmann, node based no communication
  if(mesh->pmlNelements){
    occaTimerTic(mesh->device,"PMLImplicitSolve");
    mesh->pmlImplicitSolveKernel(mesh->pmlNelements,
                                mesh->o_pmlElementIds,
                                mesh->dt,
                                mesh->LSIMEX_Ad[k],
                                mesh->o_cubInterpT,
                                mesh->o_cubProjectT,
                                mesh->o_qY,
                                mesh->o_qZ); 
    occaTimerToc(mesh->device,"PMLImplicitSolve");

    //
    occaTimerTic(mesh->device,"PMLImplicitUpdate");
    // No surface term for implicit part
    mesh->pmlImplicitUpdateKernel(mesh->pmlNelements,
                                  mesh->o_pmlElementIds,
                                  mesh->o_pmlIds,
                                  mesh->dt,
                                  ramp,
                                  mesh->LSIMEX_Ad[k],
                                  mesh->o_qY,
                                  mesh->o_qYx,
                                  mesh->o_qYy,
                                  mesh->o_qZ,
                                  mesh->o_pmlqx,
                                  mesh->o_pmlqy,
                                  mesh->o_qSx,
                                  mesh->o_qSy,
                                  mesh->o_qS,
                                  mesh->o_q);
    occaTimerToc(mesh->device,"PMLImplicitUpdate");
  }

  occaTimerToc(mesh->device,"PMLImplicit");

 //      // compute volume contribution to DG boltzmann RHS
      
  occaTimerTic(mesh->device,"NONPMLImplicit");
  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NONPMLImplicitSolve");
    mesh->implicitSolveKernel(mesh->nonPmlNelements,
                              mesh->o_nonPmlElementIds,
                              mesh->dt,
                              mesh->LSIMEX_Ad[k],
                              mesh->o_cubInterpT,
                              mesh->o_cubProjectT,
                              mesh->o_qY,
                              mesh->o_qZ); 
    occaTimerToc(mesh->device,"NONPMLImplicitSolve");

    occaTimerTic(mesh->device,"NONPMLImplicitUpdate");

    //No surface term for implicit part
    mesh->implicitUpdateKernel(mesh->nonPmlNelements,
                              mesh->o_nonPmlElementIds,
                              mesh->dt,
                              mesh->LSIMEX_Ad[k],
                              mesh->o_qZ,
                              mesh->o_qY,
                              mesh->o_q,
                              mesh->o_qS);	
    occaTimerToc(mesh->device,"NONPMLImplicitUpdate");
  }

  occaTimerToc(mesh->device,"NONPMLImplicit");


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

   
  // VOLUME KERNELS
  occaTimerTic(mesh->device,"VolumeKernel");

  // // compute volume contribution to DG boltzmann RHS
  if(mesh->pmlNelements){	
    occaTimerTic(mesh->device,"PMLVolumeKernel");
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
                          mesh->o_qY,
                          mesh->o_qYx,
                          mesh->o_qYy);
    occaTimerToc(mesh->device,"PMLVolumeKernel");	

  }

    // compute volume contribution to DG boltzmann RHS added d/dt (ramp(qbar)) to RHS
    if(mesh->nonPmlNelements){

      occaTimerTic(mesh->device,"NONPMLVolumeKernel");      	
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
                         mesh->o_qY);
      occaTimerToc(mesh->device,"NONPMLVolumeKernel");
	}
    
    
  occaTimerToc(mesh->device,"VolumeKernel");


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
    // SURFACE KERNELS
    occaTimerTic(mesh->device,"PMLSurfaceKernel"); 	
    // compute surface contribution to DG boltzmann RHS
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
                            mesh->o_qY,
                            mesh->o_qYx,
                            mesh->o_qYy);

    occaTimerToc(mesh->device,"PMLSurfaceKernel"); 	
  }

  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NONPMLSurfaceKernel"); 

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
                        mesh->o_qY);
    occaTimerToc(mesh->device,"NONPMLSurfaceKernel"); 
  }

  occaTimerToc(mesh->device,"SurfaceKernel");



  occaTimerTic(mesh->device,"UpdateKernel");

  //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
  if (mesh->pmlNelements){ 
    occaTimerTic(mesh->device,"PMLUpdateKernel");  

    mesh->pmlUpdateKernel(mesh->pmlNelements,
                          mesh->o_pmlElementIds,
                          mesh->o_pmlIds,
                          mesh->dt,
                          mesh->LSIMEX_B[k],
                          ramp,
                          mesh->o_qZ,
                          mesh->o_qY,
                          mesh->o_qYx,
                          mesh->o_qYy,
                          mesh->o_qS,
                          mesh->o_qSx,
                          mesh->o_qSy,
                          mesh->o_pmlqx,
                          mesh->o_pmlqy,
                          mesh->o_q);

    occaTimerToc(mesh->device,"PMLUpdateKernel");  
  }

  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NONPMLUpdateKernel");   	
    mesh->updateKernel(mesh->nonPmlNelements,
                      mesh->o_nonPmlElementIds,
                      mesh->dt,
                      mesh->LSIMEX_B[k],
                      mesh->o_qZ,
                      mesh->o_qY,
                      mesh->o_qS,
                      mesh->o_q);

    occaTimerToc(mesh->device,"NONPMLUpdateKernel");  
  }

  occaTimerToc(mesh->device,"UpdateKernel");      
    
}

}

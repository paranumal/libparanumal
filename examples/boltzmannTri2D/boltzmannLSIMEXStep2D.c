#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannLSIMEXStep2D(bns_t *bns, int tstep, int haloBytes,
				   dfloat * sendBuffer, dfloat *recvBuffer, setupAide &options){


bns->shiftIndex =0;
mesh2D *mesh = bns->mesh; 

for(int k=0;k<bns->Nimex;++k){

  // intermediate stage time
  dfloat t = bns->startTime+ tstep*bns->dt + bns->dt*bns->LSIMEX_C[k];

  dfloat ramp, drampdt;
  boltzmannRampFunction2D(t, &ramp, &drampdt);


  // RESIDUAL UPDATE, i.e. Y = Q+ (a(k,k-1)-b(k-1))_ex *Y + (a(k,k-1)-b(k-1))_im *Z
  occaTimerTic(mesh->device,"ResidualUpdateKernel");
  // compute volume contribution to DG boltzmann RHS
  if(mesh->pmlNelements){

    occaTimerTic(mesh->device,"PMLResidualUpdateKernel");
    //	printf("pmlNel = %d\n", mesh->pmlNelements);
    bns->pmlResidualUpdateKernel( mesh->pmlNelements,
                                  mesh->o_pmlElementIds,
                                  mesh->o_pmlIds,
                                  bns->dt,
                                  ramp,
                                  bns->LSIMEX_ABi[k],
                                  bns->LSIMEX_ABe[k],
                                  bns->o_q,
                                  bns->o_pmlqx,
                                  bns->o_pmlqy,
                                  bns->o_qZ,
                                  bns->o_rhsq,
                                  bns->o_pmlrhsqx,
                                  bns->o_pmlrhsqy,
                                  bns->o_qY,
                                  bns->o_qYx,
                                  bns->o_qYy);
    occaTimerToc(mesh->device,"PMLResidualUpdateKernel");
  }
  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NONPMLResidualUpdateKernel");
    bns->residualUpdateKernel(mesh->nonPmlNelements,
                              mesh->o_nonPmlElementIds,
                              bns->dt,
                              bns->LSIMEX_ABi[k],
                              bns->LSIMEX_ABe[k],
                              bns->o_q,
                              bns->o_qZ,
                              bns->o_rhsq,
                              bns->o_qY);
    occaTimerToc(mesh->device,"NONPMLResidualUpdateKernel");
  }
  
  occaTimerToc(mesh->device,"ResidualUpdateKernel");

      

  occaTimerTic(mesh->device,"PMLImplicit");
  // Compute Implicit Part of Boltzmann, node based no communication
  // Solves imlicit part and adds sigma volume terms also.
  if(mesh->pmlNelements){
    occaTimerTic(mesh->device,"PMLImplicitSolve");
    bns->pmlImplicitSolveKernel(mesh->pmlNelements,
                                mesh->o_pmlElementIds,
                                mesh->o_pmlIds,
                                bns->dt,
                                bns->LSIMEX_Ad[k],
                                mesh->o_cubInterpT,
                                mesh->o_cubProjectT,
                                bns->o_pmlSigmaX,
                                bns->o_pmlSigmaY,
                                bns->o_qY,
                                bns->o_qYx,
                                bns->o_qYy,
                                bns->o_rhsq,
                                bns->o_pmlrhsqx,
                                bns->o_pmlrhsqy,
                                bns->o_qZ); 
    occaTimerToc(mesh->device,"PMLImplicitSolve");

    //
    occaTimerTic(mesh->device,"PMLImplicitUpdate");
    // No surface term for implicit part
    bns->pmlImplicitUpdateKernel(mesh->pmlNelements,
                                  mesh->o_pmlElementIds,
                                  mesh->o_pmlIds,
                                  bns->dt,
                                  ramp,
                                  bns->LSIMEX_Ad[k],
                                  bns->o_qZ,
                                  bns->o_qY);
    occaTimerToc(mesh->device,"PMLImplicitUpdate");
  }

  occaTimerToc(mesh->device,"PMLImplicit");

   // compute volume contribution to DG boltzmann RHS   
  occaTimerTic(mesh->device,"NONPMLImplicit");
  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NONPMLImplicitSolve");
    bns->implicitSolveKernel(mesh->nonPmlNelements,
                              mesh->o_nonPmlElementIds,
                              bns->dt,
                              bns->LSIMEX_Ad[k],
                              mesh->o_cubInterpT,
                              mesh->o_cubProjectT,
                              bns->o_qY,
                              bns->o_qZ); 
    occaTimerToc(mesh->device,"NONPMLImplicitSolve");

    occaTimerTic(mesh->device,"NONPMLImplicitUpdate");

    //No surface term for implicit part
    bns->implicitUpdateKernel(mesh->nonPmlNelements,
                              mesh->o_nonPmlElementIds,
                              bns->dt,
                              bns->LSIMEX_Ad[k],
                              bns->o_qZ,
                              bns->o_qY);	
    occaTimerToc(mesh->device,"NONPMLImplicitUpdate");
  }

  occaTimerToc(mesh->device,"NONPMLImplicit");


  if(mesh->totalHaloPairs>0){
    // extract halo on DEVICE
    #if ASYNC 
      mesh->device.setStream(dataStream);
    #endif

    int Nentries = mesh->Np*bns->Nfields;
    mesh->haloExtractKernel(mesh->totalHaloPairs,
                            Nentries,
                            mesh->o_haloElementList,
                            bns->o_qY,
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
                          bns->o_qY,
                          bns->o_qYx,
                          bns->o_qYy,
                          bns->o_rhsq,
                          bns->o_pmlrhsqx,
                          bns->o_pmlrhsqy);
    occaTimerToc(mesh->device,"PMLVolumeKernel");	

  }

    // compute volume contribution to DG boltzmann RHS added d/dt (ramp(qbar)) to RHS
    if(mesh->nonPmlNelements){

      occaTimerTic(mesh->device,"NONPMLVolumeKernel");      	
      bns->volumeKernel(mesh->nonPmlNelements,
                         mesh->o_nonPmlElementIds,
                         ramp, 
                         drampdt,
                         bns->Nrhs,
                         bns->shiftIndex,
                         mesh->o_vgeo,
                         mesh->o_DrT,
                         mesh->o_DsT,
                         bns->o_qY,
                         bns->o_rhsq);
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
                          bns->Nfields*mesh->Np*sizeof(dfloat),
                          sendBuffer,
                          recvBuffer);

    // wait for halo data to arrive
    meshHaloExchangeFinish(mesh);


    // copy halo data to DEVICE
    size_t offset = bns->Nfields*mesh->Np*mesh->Nelements*sizeof(dfloat); // offset for halo data
    bns->o_qY.asyncCopyFrom(recvBuffer, haloBytes, offset);
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
                            bns->o_qY,
                            bns->o_rhsq,
                            bns->o_pmlrhsqx,
                            bns->o_pmlrhsqy);

    occaTimerToc(mesh->device,"PMLSurfaceKernel"); 	
  }

  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NONPMLSurfaceKernel"); 

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
                        bns->o_qY,
                        bns->o_rhsq);
    occaTimerToc(mesh->device,"NONPMLSurfaceKernel"); 
  }

  occaTimerToc(mesh->device,"SurfaceKernel");



  occaTimerTic(mesh->device,"UpdateKernel");

  //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
  if (mesh->pmlNelements){ 
    occaTimerTic(mesh->device,"PMLUpdateKernel");  

    bns->pmlUpdateKernel(mesh->pmlNelements,
                          mesh->o_pmlElementIds,
                          mesh->o_pmlIds,
                          bns->dt,
                          bns->LSIMEX_B[k],
                          ramp,
                          bns->o_qZ,
                          bns->o_qY,
                          bns->o_qYx,
                          bns->o_qYy,
                          bns->o_rhsq,
                          bns->o_pmlrhsqx,
                          bns->o_pmlrhsqy,
                          bns->o_pmlqx,
                          bns->o_pmlqy,
                          bns->o_q);

    occaTimerToc(mesh->device,"PMLUpdateKernel");  
  }

  if(mesh->nonPmlNelements){
    occaTimerTic(mesh->device,"NONPMLUpdateKernel");   	
    bns->updateKernel(mesh->nonPmlNelements,
                      mesh->o_nonPmlElementIds,
                      bns->dt,
                      bns->LSIMEX_B[k],
                      bns->o_qZ,
                      bns->o_qY,
                      bns->o_rhsq,
                      bns->o_q);

    occaTimerToc(mesh->device,"NONPMLUpdateKernel");  
  }

  occaTimerToc(mesh->device,"UpdateKernel");      
    
}

}

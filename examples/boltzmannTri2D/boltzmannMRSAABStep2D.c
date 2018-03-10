#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannMRSAABStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
                         dfloat * sendBuffer, dfloat *recvBuffer, char * options){

for (iint Ntick=0; Ntick < pow(2,mesh->MRABNlevels-1);Ntick++) {

  // intermediate stage time
  dfloat t = mesh->dt*(tstep*pow(2.0,mesh->MRABNlevels-1) + Ntick);

  iint mrab_order = 0; 

  if(tstep==0)
    mrab_order = 0; // first order
  else if(tstep==1)
    mrab_order = 1; // second order
  else
    mrab_order = 2; // third order 



  // COMPUTE RAMP FUNCTION 
  dfloat ramp, drampdt;
  boltzmannRampFunction2D(t, &ramp, &drampdt);





  iint lev;
  for (lev=0;lev<mesh->MRABNlevels;lev++)
    if (Ntick % (1<<lev) != 0) break; //find the max lev to compute rhs
  
  if(mesh->totalHaloPairs>0){
    // extract halo on DEVICE
    #if ASYNC 
      mesh->device.setStream(dataStream);
    #endif

    iint Nentries = mesh->Nfp*mesh->Nfields*mesh->Nfaces;
    mesh->haloExtractKernel(mesh->totalHaloPairs,
                            Nentries,
                            mesh->o_haloElementList,
                            mesh->o_fQM,
                            mesh->o_haloBuffer);

    // copy extracted halo to HOST
    mesh->o_haloBuffer.asyncCopyTo(sendBuffer);

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }


  occaTimerTic(mesh->device, "VolumeKernel");  

  for (iint l=0;l<lev;l++) {

    if (mesh->MRABNelements[l]){
      occaTimerTic(mesh->device, "NonPmlVolumeKernel"); 
      mesh->volumeKernel(mesh->MRABNelements[l],
                        mesh->o_MRABelementIds[l],
                        ramp, 
                        drampdt,
                        mesh->Nrhs,
                        mesh->MRABshiftIndex[l],
                        mesh->o_vgeo,
                        mesh->o_DrT,
                        mesh->o_DsT,
                        mesh->o_q,
                        mesh->o_rhsq);
      occaTimerToc(mesh->device, "NonPmlVolumeKernel"); 
    }

    if (mesh->MRABpmlNelements[l]){
      occaTimerTic(mesh->device, "PmlVolumeKernel");
      mesh->pmlVolumeKernel(mesh->MRABpmlNelements[l],
                        mesh->o_MRABpmlElementIds[l],
                        mesh->o_MRABpmlIds[l],
                        ramp, 
                        drampdt,
                        mesh->Nrhs,
                        mesh->MRABshiftIndex[l],
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
      occaTimerToc(mesh->device, "PmlVolumeKernel");
    }
  }

  occaTimerToc(mesh->device, "VolumeKernel");   

     

  if(strstr(options, "CUBATURE")){ 
    occaTimerTic(mesh->device, "RelaxationKernel");
    for (iint l=0;l<lev;l++) {
      if (mesh->MRABNelements[l]){
        occaTimerTic(mesh->device,"NonPmlRelaxationKernel");
        mesh->relaxationKernel(mesh->MRABNelements[l],
                      mesh->o_MRABelementIds[l],
                      mesh->Nrhs,
                      mesh->MRABshiftIndex[l],
                      mesh->o_cubInterpT,
                      mesh->o_cubProjectT,
                      mesh->o_q,
                      mesh->o_rhsq); 
        occaTimerToc(mesh->device,"NonPmlRelaxationKernel");
        } 

      if (mesh->MRABpmlNelements[l]){
        occaTimerTic(mesh->device,"PmlRelaxationKernel");
        mesh->pmlRelaxationKernel(mesh->MRABpmlNelements[l],
                      mesh->o_MRABpmlElementIds[l],
                      mesh->o_MRABpmlIds[l],
                      mesh->Nrhs,
                      mesh->MRABshiftIndex[l],
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
        occaTimerToc(mesh->device,"PmlRelaxationKernel");
      }
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
          mesh->Nfields*mesh->Nfp*mesh->Nfaces*sizeof(dfloat),
          sendBuffer,
          recvBuffer);

    // wait for halo data to arrive
    meshHaloExchangeFinish(mesh);

    // copy halo data to DEVICE
    size_t offset = mesh->Nfaces*mesh->Nfp*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
    mesh->o_fQM.asyncCopyFrom(recvBuffer, haloBytes, offset);
    mesh->device.finish();        

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }


  // SURFACE KERNELS for boltzmann Nodal DG
  for (iint l=0;l<lev;l++) {
    occaTimerTic(mesh->device,"SurfaceKernel");
    if (mesh->MRABNelements[l]){
      occaTimerTic(mesh->device,"NonPmlSurfaceKernel");
      mesh->surfaceKernel(mesh->MRABNelements[l],
                          mesh->o_MRABelementIds[l],
                          t,
                          ramp,
                          mesh->Nrhs,
                          mesh->MRABshiftIndex[l],
                          mesh->o_sgeo,
                          mesh->o_LIFTT,
                          mesh->o_vmapM,
                          mesh->o_mapP,
                          mesh->o_EToB,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_q,
                          mesh->o_fQM,
                          mesh->o_rhsq);
      occaTimerToc(mesh->device,"NonPmlSurfaceKernel");
    }

    if (mesh->MRABpmlNelements[l]){
      occaTimerTic(mesh->device,"PmlSurfaceKernel");
      mesh->pmlSurfaceKernel(mesh->MRABpmlNelements[l],
                          mesh->o_MRABpmlElementIds[l],
                          mesh->o_MRABpmlIds[l],
                          t,
                          ramp,
                          mesh->Nrhs,
                          mesh->MRABshiftIndex[l],
                          mesh->o_sgeo,
                          mesh->o_LIFTT,
                          mesh->o_vmapM,
                          mesh->o_mapP,
                          mesh->o_EToB,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_q,
                          mesh->o_fQM,
                          mesh->o_rhsq,
                          mesh->o_pmlrhsqx,
                          mesh->o_pmlrhsqy);
       occaTimerToc(mesh->device,"PmlSurfaceKernel");
    }

    occaTimerToc(mesh->device,"SurfaceKernel");
  }


  for (lev=0;lev<mesh->MRABNlevels;lev++)
    if ((Ntick+1) % (1<<lev) !=0) break; //find the max lev to update

      // printf(" UPDATE:\n  Ntick :%d  lev: %d", Ntick, lev);

  for (iint l=0; l<lev; l++) {

    const iint id = mrab_order*mesh->MRABNlevels*3 + l*3;
    occaTimerTic(mesh->device,"UpdateKernel");

    if (mesh->MRABNelements[l]){
      occaTimerTic(mesh->device,"NonPmlUpdateKernel");
      mesh->updateKernel(mesh->MRABNelements[l],
                        mesh->o_MRABelementIds[l],
                        mesh->MRSAAB_C[l], //
                        mesh->MRAB_A[id+0], //
                        mesh->MRAB_A[id+1],
                        mesh->MRAB_A[id+2], //
                        mesh->MRSAAB_A[id+0], //
                        mesh->MRSAAB_A[id+1],
                        mesh->MRSAAB_A[id+2], //
                        mesh->MRABshiftIndex[l],
                        mesh->o_vmapM,
                        mesh->o_rhsq,
                        mesh->o_fQM,
                        mesh->o_q);
      occaTimerToc(mesh->device,"NonPmlUpdateKernel");
    }

    if (mesh->MRABpmlNelements[l]){
      occaTimerTic(mesh->device,"PmlUpdateKernel");
      mesh->pmlUpdateKernel(mesh->MRABpmlNelements[l],
                            mesh->o_MRABpmlElementIds[l],
                            mesh->o_MRABpmlIds[l],
                            mesh->MRSAAB_C[l], //
                            mesh->MRAB_A[id+0], //
                            mesh->MRAB_A[id+1],
                            mesh->MRAB_A[id+2], //
                            mesh->MRSAAB_A[id+0], //
                            mesh->MRSAAB_A[id+1],
                            mesh->MRSAAB_A[id+2], //
                            mesh->MRABshiftIndex[l],
                            mesh->o_vmapM,
                            mesh->o_rhsq,
                            mesh->o_pmlrhsqx,
                            mesh->o_pmlrhsqy,
                            mesh->o_q,
                            mesh->o_pmlqx,
                            mesh->o_pmlqy,
                            mesh->o_fQM);
      occaTimerToc(mesh->device,"PmlUpdateKernel");
    }

    occaTimerToc(mesh->device,"UpdateKernel");

    //rotate index
    mesh->MRABshiftIndex[l] = (mesh->MRABshiftIndex[l]+1)%3;
  }


  if (lev<mesh->MRABNlevels) {    
    // const iint id = mrab_order*mesh->MRABNlevels*3 + (lev-1)*3; // !!!!!
    const iint id = mrab_order*mesh->MRABNlevels*3 + (lev)*3;
    occaTimerTic(mesh->device,"TraceUpdateKernel");

    if (mesh->MRABNhaloElements[lev]){
      occaTimerTic(mesh->device,"NonPmlTraceUpdateKernel");
      mesh->traceUpdateKernel(mesh->MRABNhaloElements[lev],
                            mesh->o_MRABhaloIds[lev],
                            mesh->MRSAAB_C[lev-1], //
                            mesh->MRAB_B[id+0], //
                            mesh->MRAB_B[id+1],
                            mesh->MRAB_B[id+2], //
                            mesh->MRSAAB_B[id+0], //
                            mesh->MRSAAB_B[id+1],
                            mesh->MRSAAB_B[id+2], 
                            mesh->MRABshiftIndex[lev],
                            mesh->o_vmapM,
                            mesh->o_rhsq,
                            mesh->o_fQM,
                            mesh->o_q);
      occaTimerToc(mesh->device,"NonPmlTraceUpdateKernel");
    }

    if (mesh->MRABpmlNhaloElements[lev]){
      occaTimerTic(mesh->device,"PmlTraceUpdateKernel");
      mesh->pmlTraceUpdateKernel(mesh->MRABpmlNhaloElements[lev],
                              mesh->o_MRABpmlHaloElementIds[lev],
                              mesh->o_MRABpmlHaloIds[lev],
                              mesh->MRSAAB_C[lev-1], //
                              mesh->MRAB_B[id+0], //
                              mesh->MRAB_B[id+1],
                              mesh->MRAB_B[id+2], //
                              mesh->MRSAAB_B[id+0], //
                              mesh->MRSAAB_B[id+1],
                              mesh->MRSAAB_B[id+2], 
                              mesh->MRABshiftIndex[lev],
                              mesh->o_vmapM,
                              mesh->o_rhsq,
                              mesh->o_pmlrhsqx,
                              mesh->o_pmlrhsqy,
                              mesh->o_q,
                              mesh->o_pmlqx,
                              mesh->o_pmlqy,
                              mesh->o_fQM);
      occaTimerToc(mesh->device,"PmlTraceUpdateKernel");
    }


   occaTimerToc(mesh->device,"TraceUpdateKernel");
  }



}

  
}

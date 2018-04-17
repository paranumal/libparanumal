#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannMRSAABStep2D(bns_t *bns, int tstep, int haloBytes,
                         dfloat * sendBuffer, dfloat *recvBuffer, char * options){


  mesh2D *mesh = bns->mesh; 

for (int Ntick=0; Ntick < pow(2,mesh->MRABNlevels-1);Ntick++) {

  // intermediate stage time
  dfloat t = bns->dt*(tstep*pow(2.0,mesh->MRABNlevels-1) + Ntick);

  int mrab_order = 0; 

  if(tstep==0)
    mrab_order = 0; // first order
  else if(tstep==1)
    mrab_order = 1; // second order
  else
    mrab_order = 2; // third order 



  // COMPUTE RAMP FUNCTION 
  dfloat ramp, drampdt;
  boltzmannRampFunction2D(t, &ramp, &drampdt);





  int lev;
  for (lev=0;lev<mesh->MRABNlevels;lev++)
    if (Ntick % (1<<lev) != 0) break; //find the max lev to compute rhs
  
  if(mesh->totalHaloPairs>0){
    // extract halo on DEVICE
    #if ASYNC 
      mesh->device.setStream(dataStream);
    #endif

    int Nentries = mesh->Nfp*mesh->Nfields*mesh->Nfaces;
    mesh->haloExtractKernel(mesh->totalHaloPairs,
                            Nentries,
                            mesh->o_haloElementList,
                            bns->o_fQM,
                            mesh->o_haloBuffer);

    // copy extracted halo to HOST
    mesh->o_haloBuffer.asyncCopyTo(sendBuffer);

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }


  occaTimerTic(mesh->device, "VolumeKernel");  

  for (int l=0;l<lev;l++) {

    if (mesh->MRABNelements[l]){
      occaTimerTic(mesh->device, "NonPmlVolumeKernel"); 
      bns->volumeKernel(mesh->MRABNelements[l],
                        mesh->o_MRABelementIds[l],
                        ramp, 
                        drampdt,
                        bns->Nrhs,
                        mesh->MRABshiftIndex[l],
                        mesh->o_vgeo,
                        mesh->o_DrT,
                        mesh->o_DsT,
                        bns->o_q,
                        bns->o_rhsq);
      occaTimerToc(mesh->device, "NonPmlVolumeKernel"); 
    }

    if (mesh->MRABpmlNelements[l]){
      occaTimerTic(mesh->device, "PmlVolumeKernel");
      bns->pmlVolumeKernel(mesh->MRABpmlNelements[l],
                          mesh->o_MRABpmlElementIds[l],
                          mesh->o_MRABpmlIds[l],
                          ramp, 
                          drampdt,
                          bns->Nrhs,
                          mesh->MRABshiftIndex[l],
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
      occaTimerToc(mesh->device, "PmlVolumeKernel");
    }
  }

  occaTimerToc(mesh->device, "VolumeKernel");   

     

  if(strstr(options, "CUBATURE")){ 
    occaTimerTic(mesh->device, "RelaxationKernel");
    for (int l=0;l<lev;l++) {
      if (mesh->MRABNelements[l]){
        occaTimerTic(mesh->device,"NonPmlRelaxationKernel");
        bns->relaxationKernel(mesh->MRABNelements[l],
                      mesh->o_MRABelementIds[l],
                      bns->Nrhs,
                      mesh->MRABshiftIndex[l],
                      mesh->o_cubInterpT,
                      mesh->o_cubProjectT,
                      bns->o_q,
                      bns->o_rhsq); 
        occaTimerToc(mesh->device,"NonPmlRelaxationKernel");
        } 

      if (mesh->MRABpmlNelements[l]){
        occaTimerTic(mesh->device,"PmlRelaxationKernel");
        bns->pmlRelaxationKernel(mesh->MRABpmlNelements[l],
                                mesh->o_MRABpmlElementIds[l],
                                mesh->o_MRABpmlIds[l],
                                bns->Nrhs,
                                mesh->MRABshiftIndex[l],
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
    bns->o_fQM.asyncCopyFrom(recvBuffer, haloBytes, offset);
    mesh->device.finish();        

    #if ASYNC 
      mesh->device.setStream(defaultStream);
    #endif
  }


  // SURFACE KERNELS for boltzmann Nodal DG
  for (int l=0;l<lev;l++) {
    occaTimerTic(mesh->device,"SurfaceKernel");
    if (mesh->MRABNelements[l]){
      occaTimerTic(mesh->device,"NonPmlSurfaceKernel");
      bns->surfaceKernel(mesh->MRABNelements[l],
                          mesh->o_MRABelementIds[l],
                          t,
                          ramp,
                          bns->Nrhs,
                          mesh->MRABshiftIndex[l],
                          mesh->o_sgeo,
                          mesh->o_LIFTT,
                          mesh->o_vmapM,
                          mesh->o_mapP,
                          mesh->o_EToB,
                          mesh->o_x,
                          mesh->o_y,
                          bns->o_q,
                          bns->o_fQM,
                          bns->o_rhsq);
      occaTimerToc(mesh->device,"NonPmlSurfaceKernel");
    }

    if (mesh->MRABpmlNelements[l]){
      occaTimerTic(mesh->device,"PmlSurfaceKernel");
      bns->pmlSurfaceKernel(mesh->MRABpmlNelements[l],
                          mesh->o_MRABpmlElementIds[l],
                          mesh->o_MRABpmlIds[l],
                          t,
                          ramp,
                          bns->Nrhs,
                          mesh->MRABshiftIndex[l],
                          mesh->o_sgeo,
                          mesh->o_LIFTT,
                          mesh->o_vmapM,
                          mesh->o_mapP,
                          mesh->o_EToB,
                          mesh->o_x,
                          mesh->o_y,
                          bns->o_q,
                          bns->o_fQM,
                          bns->o_rhsq,
                          bns->o_pmlrhsqx,
                          bns->o_pmlrhsqy);
       occaTimerToc(mesh->device,"PmlSurfaceKernel");
    }

    occaTimerToc(mesh->device,"SurfaceKernel");
  }


  for (lev=0;lev<mesh->MRABNlevels;lev++)
    if ((Ntick+1) % (1<<lev) !=0) break; //find the max lev to update

      // printf(" UPDATE:\n  Ntick :%d  lev: %d", Ntick, lev);

  for (int l=0; l<lev; l++) {

    const int id = mrab_order*mesh->MRABNlevels*3 + l*3;
    occaTimerTic(mesh->device,"UpdateKernel");

    if (mesh->MRABNelements[l]){
      occaTimerTic(mesh->device,"NonPmlUpdateKernel");
      bns->updateKernel(mesh->MRABNelements[l],
                        mesh->o_MRABelementIds[l],
                        bns->MRSAAB_C[l], //
                        bns->MRAB_A[id+0], //
                        bns->MRAB_A[id+1],
                        bns->MRAB_A[id+2], //
                        bns->MRSAAB_A[id+0], //
                        bns->MRSAAB_A[id+1],
                        bns->MRSAAB_A[id+2], //
                        mesh->MRABshiftIndex[l],
                        mesh->o_vmapM,
                        bns->o_rhsq,
                        bns->o_fQM,
                        bns->o_q);
      occaTimerToc(mesh->device,"NonPmlUpdateKernel");
    }

    if (mesh->MRABpmlNelements[l]){
      occaTimerTic(mesh->device,"PmlUpdateKernel");
      bns->pmlUpdateKernel(mesh->MRABpmlNelements[l],
                            mesh->o_MRABpmlElementIds[l],
                            mesh->o_MRABpmlIds[l],
                            bns->MRSAAB_C[l], //
                            bns->MRAB_A[id+0], //
                            bns->MRAB_A[id+1],
                            bns->MRAB_A[id+2], //
                            bns->MRSAAB_A[id+0], //
                            bns->MRSAAB_A[id+1],
                            bns->MRSAAB_A[id+2], //
                            mesh->MRABshiftIndex[l],
                            mesh->o_vmapM,
                            bns->o_rhsq,
                            bns->o_pmlrhsqx,
                            bns->o_pmlrhsqy,
                            bns->o_q,
                            bns->o_pmlqx,
                            bns->o_pmlqy,
                            bns->o_fQM);
      occaTimerToc(mesh->device,"PmlUpdateKernel");
    }

    occaTimerToc(mesh->device,"UpdateKernel");

    //rotate index
    mesh->MRABshiftIndex[l] = (mesh->MRABshiftIndex[l]+1)%3;
  }


  if (lev<mesh->MRABNlevels) {    
    // const int id = mrab_order*mesh->MRABNlevels*3 + (lev-1)*3; // !!!!!
    const int id = mrab_order*mesh->MRABNlevels*3 + (lev)*3;
    occaTimerTic(mesh->device,"TraceUpdateKernel");

    if (mesh->MRABNhaloElements[lev]){
      occaTimerTic(mesh->device,"NonPmlTraceUpdateKernel");
      bns->traceUpdateKernel(mesh->MRABNhaloElements[lev],
                            mesh->o_MRABhaloIds[lev],
                            bns->MRSAAB_C[lev-1], //
                            bns->MRAB_B[id+0], //
                            bns->MRAB_B[id+1],
                            bns->MRAB_B[id+2], //
                            bns->MRSAAB_B[id+0], //
                            bns->MRSAAB_B[id+1],
                            bns->MRSAAB_B[id+2], 
                            mesh->MRABshiftIndex[lev],
                            mesh->o_vmapM,
                            bns->o_rhsq,
                            bns->o_fQM,
                            bns->o_q);
      occaTimerToc(mesh->device,"NonPmlTraceUpdateKernel");
    }

    if (mesh->MRABpmlNhaloElements[lev]){
      occaTimerTic(mesh->device,"PmlTraceUpdateKernel");
      bns->pmlTraceUpdateKernel(mesh->MRABpmlNhaloElements[lev],
                              mesh->o_MRABpmlHaloElementIds[lev],
                              mesh->o_MRABpmlHaloIds[lev],
                              bns->MRSAAB_C[lev-1], //
                              bns->MRAB_B[id+0], //
                              bns->MRAB_B[id+1],
                              bns->MRAB_B[id+2], //
                              bns->MRSAAB_B[id+0], //
                              bns->MRSAAB_B[id+1],
                              bns->MRSAAB_B[id+2], 
                              mesh->MRABshiftIndex[lev],
                              mesh->o_vmapM,
                              bns->o_rhsq,
                              bns->o_pmlrhsqx,
                              bns->o_pmlrhsqy,
                              bns->o_q,
                              bns->o_pmlqx,
                              bns->o_pmlqy,
                              bns->o_fQM);
      occaTimerToc(mesh->device,"PmlTraceUpdateKernel");
    }
   occaTimerToc(mesh->device,"TraceUpdateKernel");
  }



}

  
}

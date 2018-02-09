#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannMRSAABStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
                         dfloat * sendBuffer, dfloat *recvBuffer, char * options){

for (iint Ntick=0; Ntick < pow(2,mesh->MRABNlevels-1);Ntick++) {

    // intermediate stage time
    dfloat t = mesh->dt*(tstep*pow(2,mesh->MRABNlevels-1) + Ntick);

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

        // Halo exchange 

    if(mesh->totalHaloPairs>0){
        // extract halo on DEVICE
        #if ASYNC 
          mesh->device.setStream(dataStream);
        #endif
        
        iint Nentries = mesh->Nfp*mesh->Nfields*mesh->Nfaces;
        mesh->haloExtractKernel(mesh->totalHaloPairs,
                    Nentries,
                    mesh->o_haloElementList,
                    mesh->o_fQP,
                    mesh->o_haloBuffer);

        // copy extracted halo to HOST
        mesh->o_haloBuffer.asyncCopyTo(sendBuffer);

        #if ASYNC 
          mesh->device.setStream(defaultStream);
        #endif
      }



      for (iint l=0;l<lev;l++) {
        if (mesh->MRABNelements[l])
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

        if (mesh->MRABpmlNelements[l])
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
      }

     

    if(strstr(options, "CUBATURE")){ 


        for (iint l=0;l<lev;l++) {
        if (mesh->MRABNelements[l])
              mesh->relaxationKernel(mesh->MRABNelements[l],
                            mesh->o_MRABelementIds[l],
                            mesh->Nrhs,
                            mesh->MRABshiftIndex[l],
                            mesh->o_cubInterpT,
                            mesh->o_cubProjectT,
                            mesh->o_q,
                            mesh->o_rhsq);  

        if (mesh->MRABpmlNelements[l])
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
             // mesh->pmlRelaxationKernel(mesh->MRABpmlNelements[l],
             //                mesh->o_MRABpmlElementIds[l],
             //                mesh->o_MRABpmlIds[l],
             //                mesh->Nrhs,
             //                mesh->MRABshiftIndex[l],
             //                mesh->o_cubInterpT,
             //                mesh->o_cubProjectT,
             //                mesh->o_q,
             //                mesh->o_rhsq);
      }
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
        mesh->o_fQP.asyncCopyFrom(recvBuffer, haloBytes, offset);
        mesh->device.finish();        

        #if ASYNC 
          mesh->device.setStream(defaultStream);
        #endif
      }


      // SURFACE KERNELS for boltzmann Nodal DG
      for (iint l=0;l<lev;l++) {
        if (mesh->MRABNelements[l])
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
                              mesh->o_fQP,
                              mesh->o_rhsq);

        if (mesh->MRABpmlNelements[l])
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
                              mesh->o_fQP,
                              mesh->o_rhsq,
                              mesh->o_pmlrhsqx,
                              mesh->o_pmlrhsqy);
      }


      for (lev=0;lev<mesh->MRABNlevels;lev++)
        if ((Ntick+1) % (1<<lev) !=0) break; //find the max lev to update



     for (iint l=0; l<lev; l++) {
        
        const iint id = mrab_order*mesh->MRABNlevels*3 + l*3;


          if (mesh->MRABNelements[l])
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
                              mesh->o_fQP,
                              mesh->o_q);

          if (mesh->MRABpmlNelements[l])
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
                                  mesh->o_fQM,
                                  mesh->o_fQP);

          // //rotate index
          mesh->MRABshiftIndex[l] = (mesh->MRABshiftIndex[l]+1)%3;
        }


        if (lev<mesh->MRABNlevels) {

          // const iint id = mrab_order*mesh->MRABNlevels*3 + (lev-1)*3; // !!!!!
          const iint id = mrab_order*mesh->MRABNlevels*3 + (lev)*3;

          if (mesh->MRABNhaloElements[lev])
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
                                    mesh->o_fQP,
                                    mesh->o_q);

          if (mesh->MRABpmlNhaloElements[lev])
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
                                      mesh->o_fQM,
                                      mesh->o_fQP);
        }




}







  
}

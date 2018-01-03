#include "boltzmann2D.h"

void boltzmannRun2D(mesh2D *mesh, char *options){

   // For initial data
  #if 0
  boltzmannReport2D(mesh, 0 ,options);
  #endif


  // MPI send buffer
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  iint haloBytes = mesh->totalHaloPairs*mesh->Nfp*mesh->Nfields*mesh->Nfaces*sizeof(dfloat);
  if (haloBytes) {
    occa::memory o_sendBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
    occa::memory o_recvBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
    sendBuffer = (dfloat*) o_sendBufferPinned.getMappedPointer();
    recvBuffer = (dfloat*) o_recvBufferPinned.getMappedPointer();
  }

  

  // Populate Trace Buffer
   dfloat zero = 0.0;
  for (iint l=0; l<mesh->MRABNlevels; l++) {

    if(strstr(options,"MRAB")){
       if (mesh->MRABNelements[l])
      mesh->updateKernel(mesh->MRABNelements[l],
                        mesh->o_MRABelementIds[l],
                        zero,
                        zero, zero, zero,
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
                            zero,
                            zero, zero, zero,
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
    }

    else if(strstr(options,"MRSAAB")){
    
    if (mesh->MRABNelements[l])
      mesh->updateKernel(mesh->MRABNelements[l],
                        mesh->o_MRABelementIds[l],
                        zero,
                        zero, zero, zero,
                        zero, zero, zero,
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
                            zero,
                            zero, zero, zero,
                            zero, zero, zero,
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
    }

   
  }


occa::initTimer(mesh->device);


  // VOLUME KERNELS
    mesh->device.finish();
    occa::tic("Boltzmann Solver");

 for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){
   // for(iint tstep=0;tstep<100;++tstep){
 
      if(strstr(options, "MRAB")){
       boltzmannMRABStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer, options);
      }

      if(strstr(options, "MRSAAB")){
      boltzmannMRSAABStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer, options);
      }

      if(strstr(options, "SRAB")){
      boltzmannSRABStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer, options);
      }

      if(strstr(options, "SAAB")){
      boltzmannSAABStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer, options);
      }

       if(strstr(options, "LSERK")){
      boltzmannLSERKStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer, options);
      }

       if(strstr(options, "SARK")){
      boltzmannSARKStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer, options);
      }

       if(strstr(options, "LSIMEX")){
      boltzmannLSIMEXStep2D(mesh, tstep, haloBytes, sendBuffer, recvBuffer, options);
      }

      
    //  if(strstr(options, "REPORT")){
    //   if((tstep%mesh->errorStep)==0){
    //     boltzmannReport2D(mesh, tstep, options);
    //   }
    // }

  }



mesh->device.finish();
occa::toc("Boltzmann Solver");
  
//   // // For Final Time
boltzmannReport2D(mesh, mesh->NtimeSteps,options);

occa::printTimer();




  //













  
  
  // //Deallocate Halo MPI storage
  // free(recvBuffer);
  // free(sendBuffer);
}



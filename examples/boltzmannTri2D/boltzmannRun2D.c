#include "boltzmann2D.h"

void boltzmannRun2D(bns_t *bns, char *options){

  mesh2D *mesh = bns->mesh; 
    
  // MPI send buffer
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  int haloBytes;

  if(strstr(options,"MRAB") || strstr(options,"MRSAAB"))
    haloBytes = mesh->totalHaloPairs*mesh->Nfp*bns->Nfields*mesh->Nfaces*sizeof(dfloat);
  else
    haloBytes = mesh->totalHaloPairs*mesh->Np*bns->Nfields*sizeof(dfloat);
  
   


  if (haloBytes) {
    // // Allocate MPI send buffer for single rate integrators
    // sendBuffer = (dfloat*) malloc(haloBytes);
    // recvBuffer = (dfloat*) malloc(haloBytes);
    
    occa::memory o_sendBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
    occa::memory o_recvBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
    sendBuffer = (dfloat*) o_sendBufferPinned.getMappedPointer();
    recvBuffer = (dfloat*) o_recvBufferPinned.getMappedPointer();
  }

  

  // Populate Trace Buffer
   dfloat zero = 0.0;
  for (int l=0; l<mesh->MRABNlevels; l++) {

    if(strstr(options,"MRAB")){
       if (mesh->MRABNelements[l])
      bns->updateKernel(mesh->MRABNelements[l],
                        mesh->o_MRABelementIds[l],
                        zero,
                        zero, zero, zero,
                        mesh->MRABshiftIndex[l],
                        mesh->o_vmapM,
                        bns->o_rhsq,
                        bns->o_fQM,
                        bns->o_q);

    if (mesh->MRABpmlNelements[l])
      bns->pmlUpdateKernel(mesh->MRABpmlNelements[l],
                            mesh->o_MRABpmlElementIds[l],
                            mesh->o_MRABpmlIds[l],
                            zero,
                            zero, zero, zero,
                            mesh->MRABshiftIndex[l],
                            mesh->o_vmapM,
                            bns->o_rhsq,
                            bns->o_pmlrhsqx,
                            bns->o_pmlrhsqy,
                            bns->o_q,
                            bns->o_pmlqx,
                            bns->o_pmlqy,
                            bns->o_fQM);
    }

    else if(strstr(options,"MRSAAB")){
    
    if (mesh->MRABNelements[l])
      bns->updateKernel(mesh->MRABNelements[l],
                        mesh->o_MRABelementIds[l],
                        zero,
                        zero, zero, zero,
                        zero, zero, zero,
                        mesh->MRABshiftIndex[l],
                        mesh->o_vmapM,
                        bns->o_rhsq,
                        bns->o_fQM,
                        bns->o_q);

    if (mesh->MRABpmlNelements[l])
      bns->pmlUpdateKernel(mesh->MRABpmlNelements[l],
                            mesh->o_MRABpmlElementIds[l],
                            mesh->o_MRABpmlIds[l],
                            zero,
                            zero, zero, zero,
                            zero, zero, zero,
                            mesh->MRABshiftIndex[l],
                            mesh->o_vmapM,
                            bns->o_rhsq,
                            bns->o_pmlrhsqx,
                            bns->o_pmlrhsqy,
                            bns->o_q,
                            bns->o_pmlqx,
                            bns->o_pmlqy,
                            bns->o_fQM);
    }

   
  }

printf("N: %d Nsteps: %d dt: %.5e \n", mesh->N, bns->NtimeSteps, bns->dt);




double tic_tot = 0.f, elp_tot = 0.f; 
double tic_sol = 0.f, elp_sol = 0.f; 
double tic_out = 0.f, elp_out = 0.f;



occa::initTimer(mesh->device);
occaTimerTic(mesh->device,"BOLTZMANN");

tic_tot = MPI_Wtime();


if(!strstr(options,"DOPRI5")){

 for(int tstep=0;tstep<bns->NtimeSteps;++tstep){
      
   // for(int tstep=0;tstep<1;++tstep){
      tic_out = MPI_Wtime();

      if(strstr(options, "REPORT")){
        if((tstep%bns->errorStep)==0){
          boltzmannReport2D(bns, tstep, options);
        }
      }

      elp_out += (MPI_Wtime() - tic_out);

      
      tic_sol = MPI_Wtime();

      if(strstr(options, "MRAB")){
       occaTimerTic(mesh->device, "MRAB"); 
       boltzmannMRABStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
       occaTimerToc(mesh->device, "MRAB"); 
      }

      if(strstr(options, "MRSAAB")){
        occaTimerTic(mesh->device, "MRSAAB"); 
        boltzmannMRSAABStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "MRSAAB"); 
      }

      if(strstr(options, "SRAB")){
        occaTimerTic(mesh->device, "SRAB"); 
        boltzmannSRABStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "SRAB"); 
      }

      if(strstr(options, "SRSAAB")){
        occaTimerTic(mesh->device, "SRSAAB"); 
        boltzmannSAABStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "SRSAAB"); 
      }

      if(strstr(options, "LSERK")){
        occaTimerTic(mesh->device, "LSERK");  
        boltzmannLSERKStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "LSERK");  
      }

       if(strstr(options, "SARK")){
        occaTimerTic(mesh->device, "SARK");  
        boltzmannSARKStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "SARK");  
      }

      if(strstr(options, "LSIMEX")){
        occaTimerTic(mesh->device, "LSIMEX");  
        boltzmannLSIMEXStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "LSIMEX");  
      } 

      elp_sol += (MPI_Wtime() - tic_sol);

  }
}else{


  boltzmannRunDOPRI2D(bns, haloBytes, sendBuffer, recvBuffer, options);


}

  elp_tot += (MPI_Wtime() - tic_tot);    
  occaTimerToc(mesh->device, "BOLTZMANN");

  // compute maximum over all processes
  double gelp_tot  = 0.f, gelp_sol = 0.f, gelp_out = 0.f;

  MPI_Allreduce(&elp_tot, &gelp_tot, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&elp_out, &gelp_out, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&elp_sol, &gelp_sol, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  
  int rank, size; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(rank==0){
    printf("ORDER\tSIZE\tTOTAL_TIME\tSOLVER_TIME\tOUTPUT TIME\n");
    printf("%2d %2d %.5e %.5e %.5e\n", mesh->N, size, gelp_tot, gelp_sol, gelp_out); 

    char fname[BUFSIZ]; sprintf(fname, "boltzmannScaling2D.dat");
    FILE *fp; fp = fopen(fname, "a");
    fprintf(fp, "%2d %2d %.5e %.5e %.5e\n", mesh->N,size,gelp_tot, gelp_sol, gelp_out); 
    fclose(fp);
  }



 
printf("writing Final data\n");  
// For Final Time
//boltzmannReport2D(bns, bns->NtimeSteps,options);

occa::printTimer();
}



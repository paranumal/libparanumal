#include "bns.h"

void bnsRun(bns_t *bns, setupAide &options){

  mesh_t  *mesh = bns->mesh; 
    
  // MPI send buffer
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  int haloBytes;

  if(options.compareArgs("TIME INTEGRATOR","MRSAAB"))
    haloBytes = mesh->totalHaloPairs*mesh->Nfp*bns->Nfields*mesh->Nfaces*sizeof(dfloat);
  else
    haloBytes = mesh->totalHaloPairs*mesh->Np*bns->Nfields*sizeof(dfloat);
  
   


  if (haloBytes) {
    occa::memory o_sendBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
    occa::memory o_recvBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
    sendBuffer = (dfloat*) o_sendBufferPinned.getMappedPointer();
    recvBuffer = (dfloat*) o_recvBufferPinned.getMappedPointer();
  }

  
  if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){
  // // Populate Trace Buffer
  // dfloat zero = 0.0;
  // for (int l=0; l<mesh->MRABNlevels; l++) {    
  //   if (mesh->MRABNelements[l])
  //   bns->updateKernel(mesh->MRABNelements[l],
  //                     mesh->o_MRABelementIds[l],
  //                     zero,
  //                     zero, zero, zero,
  //                     zero, zero, zero,
  //                     mesh->MRABshiftIndex[l],
  //                     mesh->o_vmapM,
  //                     bns->o_rhsq,
  //                     bns->o_fQM,
  //                     bns->o_q);

  // if (mesh->MRABpmlNelements[l])
  //   bns->pmlUpdateKernel(mesh->MRABpmlNelements[l],
  //                       mesh->o_MRABpmlElementIds[l],
  //                       mesh->o_MRABpmlIds[l],
  //                       zero,
  //                       zero, zero, zero,
  //                       zero, zero, zero,
  //                       mesh->MRABshiftIndex[l],
  //                       mesh->o_vmapM,
  //                       bns->o_rhsq,
  //                       bns->o_pmlrhsqx,
  //                       bns->o_pmlrhsqy,
  //                       bns->o_q,
  //                       bns->o_pmlqx,
  //                       bns->o_pmlqy,
  //                       bns->o_fQM);
  // }

  }

printf("N: %d Nsteps: %d dt: %.5e \n", mesh->N, bns->NtimeSteps, bns->dt);




double tic_tot = 0.f, elp_tot = 0.f; 
double tic_sol = 0.f, elp_sol = 0.f; 
double tic_out = 0.f, elp_out = 0.f;



occa::initTimer(mesh->device);
occaTimerTic(mesh->device,"BOLTZMANN");

tic_tot = MPI_Wtime();


// if(!( strstr(options,"DOPRI5") || strstr(options,"XDOPRI") || strstr(options,"SAADRK") )){

// if( (fixed_dt==1) ){
 for(int tstep=0;tstep<bns->NtimeSteps;++tstep){
      
   // for(int tstep=0;tstep<1;++tstep){
      tic_out = MPI_Wtime();

      if(bns->reportFlag){
        if((tstep%bns->reportStep)==0){
          bnsReport(bns, tstep, options);
        }
      }

       if(bns->errorFlag){
        if((tstep%bns->errorStep)==0){
         bnsError(bns, tstep, options);
        }
      }
  

      elp_out += (MPI_Wtime() - tic_out);

      
      tic_sol = MPI_Wtime();

     
      if(options.compareArgs("TIME INTEGRATOR", "MRSAAB")){
        // occaTimerTic(mesh->device, "MRSAAB"); 
        // boltzmannMRSAABStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        // occaTimerToc(mesh->device, "MRSAAB"); 
      }

    
      if(options.compareArgs("TIME INTEGRATOR","LSERK")){
        occaTimerTic(mesh->device, "LSERK");  
        bnsLSERKStep(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "LSERK");  
      }

      elp_sol += (MPI_Wtime() - tic_sol);

    

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
  }



 
printf("writing Final data\n");  
// For Final Time
//bnsReport(bns, bns->NtimeSteps,options);

occa::printTimer();
}



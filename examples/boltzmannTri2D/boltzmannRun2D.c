#include "boltzmann2D.h"

void boltzmannRun2D(bns_t *bns, setupAide &options){

  mesh2D *mesh = bns->mesh; 
    
  // MPI send buffer
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  int haloBytes;

  if(options.compareArgs("TIME INTEGRATOR","MRAB") || options.compareArgs("TIME INTEGRATOR","MRSAAB"))
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

    if(options.compareArgs("TIME INTEGRATOR","MRAB")){
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

    else if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){
    
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

int fixed_dt = 0; options.getArgs("FIXED TIME STEP", fixed_dt);

printf(" Fixed dt: %d \n", fixed_dt);

// if(  (  options.compareArgs("TIME INTEGRATOR","DOPRI5") && (fixed_dt==1) ) && 
//     !(  options.compareArgs("TIME INTEGRATOR","XDOPRI") 
//      || options.compareArgs("TIME INTEGRATOR","SAADRK") 
//      || options.compareArgs("TIME INTEGRATOR","IMEXRK") )){
if(!( strstr(options,"DOPRI5") || strstr(options,"XDOPRI") || strstr(options,"SAADRK") )){

// if( (fixed_dt==1) ){
 for(int tstep=0;tstep<bns->NtimeSteps;++tstep){
      
   // for(int tstep=0;tstep<1;++tstep){
      tic_out = MPI_Wtime();

      if(bns->reportFlag){
        if((tstep%bns->reportStep)==0){
          boltzmannReport2D(bns, tstep, options);
        }
      }

       if(bns->errorFlag){
        if((tstep%bns->errorStep)==0){
         boltzmannError2D(bns, tstep, options);
        }
      }
  

      elp_out += (MPI_Wtime() - tic_out);

      
      tic_sol = MPI_Wtime();

      if(options.compareArgs("TIME INTEGRATOR","DOPRI5")){
        
        dfloat time = tstep*bns->dt; 

        // boltzmannIMEXRKStep2D(bns, time, haloBytes, sendBuffer, recvBuffer, options);
        boltzmannDOPRIStep2D(bns,bns->time, haloBytes,sendBuffer, recvBuffer, options);
        bns->o_q.copyFrom(bns->o_rkq);
        bns->o_pmlqx.copyFrom(bns->o_rkqx);
        bns->o_pmlqy.copyFrom(bns->o_rkqy);
       }

      if(options.compareArgs("TIME INTEGRATOR", "MRAB")){
       occaTimerTic(mesh->device, "MRAB"); 
       boltzmannMRABStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
       occaTimerToc(mesh->device, "MRAB"); 
      }

      if(options.compareArgs("TIME INTEGRATOR", "MRSAAB")){
        occaTimerTic(mesh->device, "MRSAAB"); 
        boltzmannMRSAABStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "MRSAAB"); 
      }

      if(options.compareArgs("TIME INTEGRATOR", "SRAB")){
        occaTimerTic(mesh->device, "SRAB"); 
        boltzmannSRABStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "SRAB"); 
      }

      if(options.compareArgs("TIME INTEGRATOR", "SRSAAB")){
        occaTimerTic(mesh->device, "SRSAAB"); 
        boltzmannSAABStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "SRSAAB"); 
      }

      if(options.compareArgs("TIME INTEGRATOR","LSERK")){
        occaTimerTic(mesh->device, "LSERK");  
        boltzmannLSERKStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "LSERK");  
      }

       if(options.compareArgs("TIME INTEGRATOR","SARK")){
        occaTimerTic(mesh->device, "SARK");  
        boltzmannSARKStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "SARK");  
      }

      if(options.compareArgs("TIME INTEGRATOR","LSIMEX")){
        occaTimerTic(mesh->device, "LSIMEX");  
        boltzmannLSIMEXStep2D(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "LSIMEX");  
      } 

      elp_sol += (MPI_Wtime() - tic_sol);

    

  }

// char fname[BUFSIZ]; sprintf(fname, "boltzmannAddaptiveDt2D_M_%1d_iTau_%.f_Re_%.0f_Ma_%.0f.dat", bns->tmethod, bns->tauInv, bns->Re, bns->Ma);
//   FILE *fp; fp = fopen(fname, "a");
//   fprintf(fp, "%.5e %.5e\n", bns->time, bns->dt); 
//   fclose(fp);

}

  else {

  printf("====================Running Embedded Scahemes==========================\n");
  boltzmannRunEmbedded2D(bns, haloBytes, sendBuffer, recvBuffer, options);


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

    // char fname[BUFSIZ]; sprintf(fname, "boltzmannTotalTime2D_M_%1d_iTau_%.f_Re_%.0f_Ma_%.0f.dat", bns->tmethod, bns->tauInv, bns->Re, bns->Ma);
    // FILE *fp; fp = fopen(fname, "a");
    // fprintf(fp, "%2d %2d %.5e %.5e %.5e %d %d %d %d\n", mesh->N,size,gelp_tot, gelp_sol, gelp_out, bns->NtimeSteps, bns->atstep, bns->tstep, bns->rtstep); 
    // fclose(fp);
  }



 
printf("writing Final data\n");  
// For Final Time
//boltzmannReport2D(bns, bns->NtimeSteps,options);

occa::printTimer();
}



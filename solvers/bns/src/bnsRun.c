/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

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
#if 0
    occa::memory o_sendBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
    occa::memory o_recvBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
    sendBuffer = (dfloat*) o_sendBufferPinned.getMappedPointer();
    recvBuffer = (dfloat*) o_recvBufferPinned.getMappedPointer();
#endif
    
    sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, haloBytes, NULL, bns->o_sendBufferPinned);
    recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, haloBytes, NULL, bns->o_recvBufferPinned);
  }

  if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){
    printf("Populating trace values\n");
    // Populate Trace Buffer
    dlong offset = mesh->Np*mesh->Nelements*bns->Nfields;
    for (int l=0; l<mesh->MRABNlevels; l++) {  
      const int id = 3*mesh->MRABNlevels*3 + 3*l;
      if (mesh->MRABNelements[l])
        bns->traceUpdateKernel(mesh->MRABNelements[l],
                               mesh->o_MRABelementIds[l],
                               offset,
                               mesh->MRABshiftIndex[l],
                               bns->MRSAAB_C[l-1], //
                               bns->MRAB_B[id+0], //
                               bns->MRAB_B[id+1],
                               bns->MRAB_B[id+2], //
                               bns->MRSAAB_B[id+0], //
                               bns->MRSAAB_B[id+1],
                               bns->MRSAAB_B[id+2],
                               mesh->o_vmapM,
                               bns->o_q,
                               bns->o_rhsq,
                               bns->o_fQM);
      // if(bns->pmlFlag){
        if (mesh->MRABpmlNelements[l])
          bns->traceUpdateKernel(mesh->MRABpmlNelements[l],
                                 mesh->o_MRABpmlElementIds[l],
                                 offset,
                                 mesh->MRABshiftIndex[l],
                                 bns->MRSAAB_C[l-1], //
                                 bns->MRAB_B[id+0], //
                                 bns->MRAB_B[id+1],
                                 bns->MRAB_B[id+2], //
                                 bns->MRSAAB_B[id+0], //
                                 bns->MRSAAB_B[id+1],
                                 bns->MRSAAB_B[id+2],
                                 mesh->o_vmapM,
                                 bns->o_q,
                                 bns->o_rhsq,
                                 bns->o_fQM);
      // }
    }
  }

  if(mesh->rank==0) printf("N: %d Nsteps: %d dt: %.5e \n", mesh->N, bns->NtimeSteps, bns->dt);


  double tic_tot = 0.f, elp_tot = 0.f; 
  double tic_sol = 0.f, elp_sol = 0.f; 
  double tic_out = 0.f, elp_out = 0.f;

  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"BOLTZMANN");

  tic_tot = MPI_Wtime();

  if( options.compareArgs("TIME INTEGRATOR", "MRSAAB")  || 
      options.compareArgs("TIME INTEGRATOR", "LSERK") ){

    // for(int tstep=0;tstep<1;++tstep){
    for(int tstep=0;tstep<bns->NtimeSteps;++tstep){
        tic_out = MPI_Wtime();

      if(bns->reportFlag){
        if((tstep%bns->reportStep)==0){ 
          dfloat time =0; 
          if(options.compareArgs("TIME INTEGRATOR", "MRSAAB"))
            time = bns->startTime + bns->dt*tstep*pow(2,(mesh->MRABNlevels-1));     
          else
            time = bns->startTime + tstep*bns->dt;

          bnsReport(bns, time, options);

          // Write a restart file
          if(bns->writeRestartFile){
            if(mesh->rank==0) printf("\nWriting Binary Restart File....");
              bnsRestartWrite(bns, options, time);
            if(mesh->rank==0) printf("done\n");
          }   
        }
      }

      if(bns->errorFlag){
        if((tstep%bns->errorStep)==0){
         dfloat time =0; 
          if(options.compareArgs("TIME INTEGRATOR", "MRSAAB"))
            time = bns->startTime + bns->dt*tstep*pow(2,(mesh->MRABNlevels-1));     
          else
            time = bns->startTime + tstep*bns->dt;
          
         bnsError(bns, tstep, options);
        }
      }
  

      elp_out += (MPI_Wtime() - tic_out);
      tic_sol = MPI_Wtime();

     
      if(options.compareArgs("TIME INTEGRATOR", "MRSAAB")){
        occaTimerTic(mesh->device, "MRSAAB"); 
        bnsMRSAABStep(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "MRSAAB"); 
      }

      if(options.compareArgs("TIME INTEGRATOR","LSERK")){
        occaTimerTic(mesh->device, "LSERK");  
        bnsLSERKStep(bns, tstep, haloBytes, sendBuffer, recvBuffer, options);
        occaTimerToc(mesh->device, "LSERK");

      if(!(tstep%100)){
        bns->o_q.copyTo(bns->q);
        for(int fld=0;fld<bns->Nfields;++fld){
          dfloat maxq = -1e9, minq = 1e9;
          for(int e=0;e<mesh->Nelements;++e){
            for(int n=0;n<mesh->Np;++n){
              hlong id = e*mesh->Np*mesh->Nfields + n;
              maxq = mymax(maxq, bns->q[id + fld*mesh->Np]);
              minq = mymin(minq, bns->q[id + fld*mesh->Np]);
            }
          }
        if(fld==0)
        printf("fld %d in [%g,%g]\n", fld, minq, maxq);
        }
      }
    }

      
      /*

      if(options.compareArgs("TIME INTEGRATOR","SARK")){
        occaTimerTic(mesh->device, "SARK");
        dfloat time = tstep*bns->dt;  
        bnsSARKStep(bns, time, haloBytes, sendBuffer, recvBuffer, options);
        bns->o_q.copyFrom(bns->o_rkq);
        if(mesh->pmlNelements){
          bns->o_pmlqx.copyFrom(bns->o_rkqx);
          bns->o_pmlqy.copyFrom(bns->o_rkqy);
          if(bns->dim==3)
            bns->o_pmlqz.copyFrom(bns->o_rkqz);
        }

        occaTimerToc(mesh->device, "SARK");  
      }
      */

      elp_sol += (MPI_Wtime() - tic_sol);
    }
  }else if( options.compareArgs("TIME INTEGRATOR", "SARK")){

    occaTimerTic(mesh->device, "SARK_TOTAL");
    bnsRunEmbedded(bns, haloBytes, sendBuffer, recvBuffer, options);
    occaTimerToc(mesh->device, "SARK_TOTAL");

  }else{
    printf("Wrong time stepper\n");
    exit(EXIT_FAILURE); 
  }

 


  elp_tot += (MPI_Wtime() - tic_tot);    
  occaTimerToc(mesh->device, "BOLTZMANN");

  // compute maximum over all processes
  double gelp_tot  = 0.f, gelp_sol = 0.f, gelp_out = 0.f;

  MPI_Allreduce(&elp_tot, &gelp_tot, 1, MPI_DOUBLE, MPI_MAX, mesh->comm);
  MPI_Allreduce(&elp_out, &gelp_out, 1, MPI_DOUBLE, MPI_MAX, mesh->comm);
  MPI_Allreduce(&elp_sol, &gelp_sol, 1, MPI_DOUBLE, MPI_MAX, mesh->comm);
  
  if(mesh->rank==0){
    printf("ORDER\tSIZE\tTOTAL_TIME\tSOLVER_TIME\tOUTPUT TIME\n");
    printf("%2d %2d %.5e %.5e %.5e\n", mesh->N, mesh->size, gelp_tot, gelp_sol, gelp_out); 
  }
 
  printf("writing Final data\n");  
  // For Final Time
  //bnsReport(bns, bns->NtimeSteps,options);

  occa::printTimer();
}



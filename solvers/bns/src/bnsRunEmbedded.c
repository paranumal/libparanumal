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

void bnsRunEmbedded(bns_t *bns, int haloBytes, dfloat * sendBuffer,
		    dfloat *recvBuffer, setupAide &options){

  mesh_t *mesh = bns->mesh;

  dfloat safe  = 0.95;   //safety factor
  //error control parameters
  dfloat beta       = 0.05; 
  dfloat factor1    = 0.25;
  dfloat factor2    = 10.0;
  dfloat exp1       = 1.0/bns->rkp -0.75*beta; 
  dfloat invfactor1 = 1.0/factor1;
  dfloat invfactor2 = 1.0/factor2;
  dfloat facold     = 1E-4;

  dfloat hmin = 1e9;
  for(dlong e=0;e<mesh->Nelements;++e){  
    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];
      dfloat hest = .5/(sJ*invJ);
      hmin = mymin(hmin, hest);
      }
  }


  // hard code this for the moment
  dfloat outputInterval;
  options.getArgs("OUTPUT INTERVAL", outputInterval);

  dfloat nextOutputTime = bns->startTime + outputInterval;
  dfloat outputNumber = 0;

  //initial time
  bns->time = bns->startTime;
  bns->tstep = 0;
  bns->atstep = 0; 
  bns->rtstep = 0;  

  // Compute Coefficients before starting loop
  bnsSAADRKCoefficients(bns, options);

  if(bns->reportFlag)
    bnsReport(bns, bns->time, options);

  int done =0;
  while (!done) {

    if (bns->dt<bns->dtMIN){
      printf("ERROR: Time step became too small at time step=%d\n", bns->tstep);
      exit (-1);
    }
    if (isnan(bns->dt)) {
      printf("ERROR: Solution became unstable at time step=%d\n", bns->tstep);
      exit (-1);
    }


    //check for final timestep
    if (bns->time+bns->dt > bns->finalTime){
      bns->dt = bns->finalTime-bns->time;
      done = 1;
    }

    occaTimerTic(mesh->device, "SARK_STEP"); 
    bnsSARKStep(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);
    occaTimerToc(mesh->device, "SARK_STEP"); 
    
    

    occaTimerTic(mesh->device, "SARK_ERROR"); 
    //Error estimation 
    dlong Ntotal = mesh->Nelements*mesh->Np*bns->Nfields;
    // printf("Ntotal: %d \t %d\n", Ntotal, bns->Nblock);
    bns->errorEstimateKernel(Ntotal, 
                            bns->ATOL,
                            bns->RTOL,
                            bns->o_q,
                            bns->o_rkq,
                            bns->o_rkerr,
                            bns->o_errtmp);

    bns->o_errtmp.copyTo(bns->errtmp);
    bns->o_rkerr.copyTo(bns->rkerr);


    dfloat maxerr = 0.f;
    for(dlong e=0; e<mesh->Nelements; e++){
      for(int n =0;n<mesh->Np;n++){
        const dlong id = e*mesh->Np + n; 
          maxerr= mymax(maxerr,bns->rkerr[id]);
       }
    }

    // printf("Max err:\t%.5e\n",maxerr);
    dfloat localerr = 0, err = 0;

    for(int n=0;n<bns->Nblock;++n){
      localerr += bns->errtmp[n];
    }

    MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
   
    // printf("Error\t:\t%.5e\t%.5e\t\n", err, bns->dt);
    err = sqrt(err/(bns->totalElements*mesh->Np*bns->Nfields));
    
    occaTimerToc(mesh->device, "SARK_ERROR"); 




    dfloat fac1 = pow(err,exp1);
    dfloat fac  = fac1/pow(facold,beta);

    fac = mymax(invfactor2, mymin(invfactor1,fac/safe));
    dfloat dtnew = bns->dt/fac;

    if(err<1.0){


      if(bns->reportFlag){
        occaTimerTic(mesh->device, "SARK_OUTPUT"); 
        // check for output during this step and do a mini-step
        if(bns->time<nextOutputTime && bns->time+bns->dt>nextOutputTime){

          // Write a restart file
          dfloat savedt = bns->dt;        
          // save rkq
          bns->o_saveq.copyFrom(bns->o_rkq);
          if(mesh->pmlNelements){
            bns->o_saveqx.copyFrom(bns->o_rkqx);
            bns->o_saveqy.copyFrom(bns->o_rkqy);
            if(bns->dim==3)
              bns->o_saveqz.copyFrom(bns->o_rkqz);
          }

          // change dt to match output
          bns->dt = nextOutputTime-bns->time;
          // print
          if(mesh->rank==0) printf("Taking output mini step: %g\n", bns->dt);
          
          // Compute new coefficients
          bnsSAADRKCoefficients(bns, options);

          // if(options.compareArgs("TIME INTEGRATOR","SARK"))  // SA Adaptive RK 
          bnsSARKStep(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);
          // shift for output
          bns->o_rkq.copyTo(bns->o_q);
          // output  (print from rkq)
          bnsReport(bns, nextOutputTime, options);
          
          // Write a restart file
          if(bns->writeRestartFile){
            if(mesh->rank==0) printf("\nWriting Binary Restart File....");
	    bnsRestartWrite(bns, options, nextOutputTime);
	    if(mesh->rank==0) printf("done\n");
          } 

          // restore time step
          bns->dt = savedt;

          // Go back to old coefficients
          bnsSAADRKCoefficients(bns, options);

          // increment next output time
          nextOutputTime += outputInterval;

          // accept saved rkq
          bns->o_q.copyFrom(bns->o_saveq);
          if(mesh->pmlNelements){
            bns->o_pmlqx.copyFrom(bns->o_saveqx);
            bns->o_pmlqy.copyFrom(bns->o_saveqy);
            if(bns->dim==3)
              bns->o_pmlqz.copyFrom(bns->o_saveqz);
          }
        } 
        else{
          // Accept rkq
          bns->o_q.copyFrom(bns->o_rkq);
          if(mesh->pmlNelements){
            bns->o_pmlqx.copyFrom(bns->o_rkqx);
            bns->o_pmlqy.copyFrom(bns->o_rkqy); 
            if(bns->dim==3)
              bns->o_pmlqz.copyFrom(bns->o_rkqz);
          } 
        }
        occaTimerToc(mesh->device, "SARK_OUTPUT"); 
      }
      else{
        // Accept rkq
        bns->o_q.copyFrom(bns->o_rkq);
        if(mesh->pmlNelements){
          bns->o_pmlqx.copyFrom(bns->o_rkqx);
          bns->o_pmlqy.copyFrom(bns->o_rkqy); 
          if(bns->dim==3)
            bns->o_pmlqz.copyFrom(bns->o_rkqz);
        }
      }

      if(bns->errorFlag){
        if((bns->tstep%bns->errorStep)==0){
            bnsError(bns, bns->time, options);
        } 
      }

      if(bns->outputForceStep){
        if(bns->tstep%bns->outputForceStep){
          bns->o_q.copyTo(bns->q);
          bnsForces(bns,bns->time,options);
        }
      }
      

      facold = mymax(err,1E-4);
      bns->time += bns->dt;

      if(mesh->rank==0) printf("\r time = %g (%d), dt = %g accepted (ratio dt/hmin = %g)               ", bns->time, bns->atstep, bns->dt, bns->dt/hmin);
      bns->tstep++;
    }
    else{
      bns->rtstep++; 
      dtnew = bns->dt/(mymax(invfactor1,fac1/safe));
      if(mesh->rank==0) printf("\r time = %g (%d), dt = %g rejected (ratio dt/min = %g), trying %g", bns->time,bns->atstep, bns->dt, bns->dt/hmin, dtnew);
      done =0;
    }
   
    bns->dt = dtnew;
    bns->atstep++;
    



    bnsSAADRKCoefficients(bns, options);
    #if 0
    char fname[BUFSIZ]; sprintf(fname, "boltzmannAddaptiveDt.dat");
    FILE *fp; fp = fopen(fname, "a");
    fprintf(fp, "%.5e %.5e\n", bns->time, bns->dt); 
    fclose(fp);
    #endif
  }


  

    printf("Total Step: %d, Rejected Step: %d, Accepted Step: %d \n", bns->atstep, bns->rtstep, bns->tstep);

}



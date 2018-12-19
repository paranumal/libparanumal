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

#include "cns.h"

void cnsRun(cns_t *cns, setupAide &options){

  mesh_t *mesh = cns->mesh;

  cnsReport(cns, 0, options);

  occa::timer timer;
  
  timer.initTimer(mesh->device);

  timer.tic("Run");
  
  if (options.compareArgs("TIME INTEGRATOR","DOPRI5")) {
    int Nregect = 0;

    // hard code this for the moment
    dfloat outputInterval;
    options.getArgs("TIME OUTPUT INTERVAL", outputInterval);
    bool timeIntervalFlag = (outputInterval > 0.);

    dfloat nextOutputTime = outputInterval;
    
    int outputTstepInterval;
    options.getArgs("TSTEP OUTPUT INTERVAL", outputTstepInterval);
    bool tstepIntervalFlag = (outputTstepInterval > 0);

    //initial time
    dfloat time = 0.0;
    int tstep=0, allStep = 0;

    int done =0;
    while (!done) {

      cns->advSwitch = 1;
      
      if (mesh->dt<cns->dtMIN){
        printf("ERROR: Time step became too small at time step=%d\n", tstep);
        exit (-1);
      }
      if (isnan(mesh->dt)) {
        printf("ERROR: Solution became unstable at time step=%d\n", tstep);
        exit (-1);
      }

      //check for final timestep
      if (time+mesh->dt > mesh->finalTime){
        mesh->dt = mesh->finalTime-time;
        done = 1;
      }

      // try a step with the current time step
      cnsDopriStep(cns, options, time);

      // compute Dopri estimator
      dfloat err = cnsDopriEstimate(cns);
                                         
      // build controller
      dfloat fac1 = pow(err,cns->exp1);
      dfloat fac = fac1/pow(cns->facold,cns->beta);

      fac = mymax(cns->invfactor2, mymin(cns->invfactor1,fac/cns->safe));
      dfloat dtnew = mymin(cns->dtMAX, mesh->dt/fac);

      if (err<1.0) { //dt is accepted

        // check for time interval output during this step
        if(timeIntervalFlag && time<nextOutputTime && time+mesh->dt>=nextOutputTime){
          cnsDopriOutputStep(cns, time,mesh->dt,nextOutputTime, cns->o_saveq);

          cns->o_saveq.copyTo(cns->o_q);
          
          // output  (print from rkq)
          cnsReport(cns, nextOutputTime, options);

          // increment next output time
          nextOutputTime += outputInterval;
        }
        
        // accept rkq
        cns->o_q.copyFrom(cns->o_rkq);

        time += mesh->dt;
        tstep++;

	if(0){
	  dfloat *maxStresses = (dfloat*) calloc(cns->Nstresses, sizeof(dfloat));
	  cns->o_viscousStresses.copyTo(cns->viscousStresses);
	  for(int e=0;e<mesh->Nelements;++e){
	    for(int fld=0;fld<cns->Nstresses;++fld){
	      for(int n=0;n<mesh->Np;++n){
		maxStresses[fld] =
		  mymax(fabs(cns->viscousStresses[e*mesh->Np*cns->Nstresses+fld*mesh->Np+n]), maxStresses[fld]);
	      }
	    }
	  }
	  for(int fld=0;fld<cns->Nstresses;++fld){
	    printf("maxStresses[%d] = %g\n", fld, maxStresses[fld]);
	  }
	  
	  dfloat *maxQ = (dfloat*) calloc(mesh->Nfields, sizeof(dfloat));
	  cns->o_q.copyTo(cns->q);
	  for(int e=0;e<mesh->Nelements;++e){
	    for(int fld=0;fld<mesh->Nfields;++fld){
	      for(int n=0;n<mesh->Np;++n){
		maxQ[fld] =
		  mymax(fabs(cns->q[e*mesh->Np*mesh->Nfields+fld*mesh->Np+n]), maxQ[fld]);
	      }
	    }
	  }
	  for(int fld=0;fld<mesh->Nfields;++fld){
	    printf("maxQ[%d] = %g\n", fld, maxQ[fld]);
	  }
	}      
	
	
	if(cns->outputForceStep){
	  if(tstep%cns->outputForceStep){
	    cns->o_q.copyTo(cns->q);
	    cns->mesh->device.finish();
	    cnsForces(cns,time);
	    
	  }
	}
	
	cns->facold = mymax(err,1E-4); // hard coded factor ?
	
	// check for time step interval output during this step
	if(tstepIntervalFlag && (tstep%outputTstepInterval==0)){
	  cns->o_q.copyTo(cns->q);//  ?????
          
	// output  (print from rkq)
	cnsReport(cns, nextOutputTime, options);

	// increment next output time
	nextOutputTime += outputInterval;
      }
    } else {
      dtnew = mesh->dt/(mymax(cns->invfactor1,fac1/cns->safe));
      Nregect++;

      done = 0;
    }

    mesh->dt = dtnew;
    allStep++;

    printf("\rTime = %.4e (%d). Average Dt = %.4e, Rejection rate = %.2g   ", time, tstep, time/(dfloat)tstep, Nregect/(dfloat) tstep); fflush(stdout);
  }
    
  mesh->device.finish();
    
  double elapsed  = timer.toc("Run");

  printf("\nRun took %lg seconds for %d accepted steps and %d total steps\n", elapsed, tstep, allStep);
    
} else if (options.compareArgs("TIME INTEGRATOR","LSERK4")) {

  for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

    dfloat time = tstep*mesh->dt;

    cnsLserkStep(cns, options, time);
      
    if(((tstep+1)%mesh->errorStep)==0){
      time += mesh->dt;
      cnsReport(cns, time, options);
    }
  }
 }
  
}

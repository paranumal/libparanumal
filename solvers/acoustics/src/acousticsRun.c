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

#include "acoustics.h"

void acousticsRun(acoustics_t *acoustics, setupAide &newOptions){

  mesh_t *mesh = acoustics->mesh;

  acousticsReport(acoustics, 0, newOptions);

  occa::timer timer;
  
  timer.initTimer(mesh->device);

  timer.tic("Run");
  
  if (newOptions.compareArgs("TIME INTEGRATOR","DOPRI5")) {
    
    // hard code this for the moment
    dfloat outputInterval;
    newOptions.getArgs("OUTPUT INTERVAL", outputInterval);
    
    dfloat nextOutputTime = outputInterval;
    dfloat outputNumber = 0;
    
    //initial time
    dfloat time = 0.0;
    int tstep=0, allStep = 0;

    int done =0;
    while (!done) {

      acoustics->advSwitch = 1;
      
      if (mesh->dt<acoustics->dtMIN){
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
      acousticsDopriStep(acoustics, newOptions, time);

      // compute Dopri estimator
      dfloat err = acousticsDopriEstimate(acoustics);
					 
      // build controller
      dfloat fac1 = pow(err,acoustics->exp1);
      dfloat fac = fac1/pow(acoustics->facold,acoustics->beta);

      fac = mymax(acoustics->invfactor2, mymin(acoustics->invfactor1,fac/acoustics->safe));
      dfloat dtnew = mesh->dt/fac;

      if (err<1.0) { //dt is accepted

	// check for output during this step and do a mini-step
	if(time<nextOutputTime && time+mesh->dt>nextOutputTime){
	  dfloat savedt = mesh->dt;
	  
	  // save rkq
	  acoustics->o_saveq.copyFrom(acoustics->o_rkq);

	  // change dt to match output
	  mesh->dt = nextOutputTime-time;

	  // print
	  printf("Taking output mini step: %g\n", mesh->dt);
	  
	  // time step to output
	  acousticsDopriStep(acoustics, newOptions, time);	  

	  // shift for output
	  acoustics->o_rkq.copyTo(acoustics->o_q);
	  
	  // output  (print from rkq)
	  acousticsReport(acoustics, nextOutputTime, newOptions);

	  // restore time step
	  mesh->dt = savedt;

	  // increment next output time
	  nextOutputTime += outputInterval;

	  // accept saved rkq
	  acoustics->o_q.copyFrom(acoustics->o_saveq);
	}
	else{
	  // accept rkq
	  acoustics->o_q.copyFrom(acoustics->o_rkq);
	}

        time += mesh->dt;

        acoustics->facold = mymax(err,1E-4); // hard coded factor ?

	printf("\r time = %g (%d), dt = %g accepted                      ", time, allStep,  mesh->dt);
        tstep++;
      } else {
        dtnew = mesh->dt/(mymax(acoustics->invfactor1,fac1/acoustics->safe));
	printf("\r time = %g (%d), dt = %g rejected, trying %g", time, allStep, mesh->dt, dtnew);

	done = 0;
      }
      mesh->dt = dtnew;
      allStep++;

    }
    
    mesh->device.finish();
    
    double elapsed  = timer.toc("Run");

    printf("run took %lg seconds for %d accepted steps and %d total steps\n", elapsed, tstep, allStep);
    
  } else if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")) {

    for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

      dfloat time = tstep*mesh->dt;

      acousticsLserkStep(acoustics, newOptions, time);
      
      if(((tstep+1)%mesh->errorStep)==0){
	time += mesh->dt;
        acousticsReport(acoustics, time, newOptions);
      }
    }
  }
  
}

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

#include "advection.h"

void advectionRun(advection_t *advection, setupAide &newOptions){

  mesh_t *mesh = advection->mesh;

  // advectionReport(advection, 0, newOptions);

  occa::timer timer;
  
  timer.initTimer(mesh->device);

  timer.tic("Run");

  int tstep=0, allStep = 0;

  MPI_Barrier(MPI_COMM_WORLD);

  mesh->device.finish();
  
  occa::streamTag start = mesh->device.tagStream();
  
  if (newOptions.compareArgs("TIME INTEGRATOR","DOPRI5")) {
    
    // hard code this for the moment
    dfloat outputInterval;
    newOptions.getArgs("OUTPUT INTERVAL", outputInterval);
    
    dfloat nextOutputTime = outputInterval;
    dfloat outputNumber = 0;
    
    //initial time
    dfloat time = 0.0;
  
    int done =0;
    while (!done) {

      advection->advSwitch = 1;
      
      if (mesh->dt<advection->dtMIN){
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
      advectionDopriStep(advection, newOptions, time);

      // compute Dopri estimator
      dfloat err = advectionDopriEstimate(advection);
					 
      // build controller
      dfloat fac1 = pow(err,advection->exp1);
      dfloat fac = fac1/pow(advection->facold,advection->beta);

      fac = mymax(advection->invfactor2, mymin(advection->invfactor1,fac/advection->safe));
      dfloat dtnew = mesh->dt/fac;

      if (err<1.0) { //dt is accepted

	// check for output during this step and do a mini-step
	if(time<nextOutputTime && time+mesh->dt>nextOutputTime){
	  dfloat savedt = mesh->dt;
	  
	  // save rkq
	  advection->o_saveq.copyFrom(advection->o_rkq);

	  // change dt to match output
	  mesh->dt = nextOutputTime-time;

	  // print
	  printf("Taking output mini step: %g\n", mesh->dt);
	  
	  // time step to output
	  advectionDopriStep(advection, newOptions, time);	  

	  // shift for output
	  advection->o_rkq.copyTo(advection->o_q);
	  
	  // output  (print from rkq)
	  advectionReport(advection, nextOutputTime, newOptions);

	  // restore time step
	  mesh->dt = savedt;

	  // increment next output time
	  nextOutputTime += outputInterval;

	  // accept saved rkq
	  advection->o_q.copyFrom(advection->o_saveq);
	}
	else{
	  // accept rkq
	  advection->o_q.copyFrom(advection->o_rkq);
	}

        time += mesh->dt;

        advection->facold = mymax(err,1E-4); // hard coded factor ?
	if(!(tstep%1000))
	  printf("\r time = %g (%d), dt = %g accepted                      ", time, allStep,  mesh->dt);
        tstep++;
      } else {
        dtnew = mesh->dt/(mymax(advection->invfactor1,fac1/advection->safe));

	if(!(tstep%1000))
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

      advectionLserkStep(advection, newOptions, time);
      
      if(((tstep+1)%mesh->errorStep)==0){
	time += mesh->dt;
        advectionReport(advection, time, newOptions);
      }
    }

    allStep = mesh->NtimeSteps;
  }
  
  occa::streamTag end = mesh->device.tagStream();

  mesh->device.finish();
  
  double elapsed = mesh->device.timeBetween(start, end);
  
  printf("%d %d %d %lg %lg %lg %lg %lg \%\%[ N, Nel, Nodes, elapsed, time/step, nodes/time, gnodes*Nsteps/time, gnodes*Nstages*Nsteps/time, %s ]\n",
	 mesh->N,
	 mesh->Nelements,
	 mesh->Np*mesh->Nelements,
	 elapsed,
	 elapsed/allStep,
	 (1.*mesh->Np*mesh->Nelements)/(elapsed),
	 (1.*mesh->Np*mesh->Nelements*allStep)/(1.e9*elapsed),
	 (1.*mesh->Np*mesh->Nelements*allStep*mesh->Nrk)/(1.e9*elapsed),
	 newOptions.getArgs("ADVECTION FORMULATION").c_str()
	 );
}

#include "cnsTri2D.h"

void cnsRunTri2D(cns_t *cns, setupAide &newOptions){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = cns->mesh;

  cnsReportTri2D(cns, 0, newOptions);

  occa::timer timer;
  
  timer.initTimer(mesh->device);

  timer.tic("Run");
  
  if (newOptions.compareArgs("TIME INTEGRATOR","DOPRI5")) {

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
    newOptions.getArgs("OUTPUT INTERVAL", outputInterval);
    
    dfloat nextOutputTime = outputInterval;
    dfloat outputNumber = 0;
    
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
      cnsDopriStepTri2D(cns, newOptions, time);
      
      // compute Dopri estimator
      dfloat err = cnsDopriEstimateTri2D(cns);
					 
      // build controller
      dfloat fac1 = pow(err,cns->exp1);
      dfloat fac = fac1/pow(cns->facold,cns->beta);

      fac = mymax(cns->invfactor2, mymin(cns->invfactor1,fac/cns->safe));
      dfloat dtnew = mesh->dt/fac;

      if (err<1.0) { //dt is accepted

	// check for output during this step and do a mini-step
	if(time<nextOutputTime && time+mesh->dt>nextOutputTime){
	  dfloat savedt = mesh->dt;
	  
	  // save rkq
	  cns->o_saveq.copyFrom(cns->o_rkq);

	  // change dt to match output
	  mesh->dt = nextOutputTime-time;

	  // print
	  printf("Taking output mini step: %g\n", mesh->dt);
	  
	  // time step to output
	  cnsDopriStepTri2D(cns, newOptions, time);	  

	  // shift for output
	  cns->o_rkq.copyTo(cns->o_q);
	  
	  // output  (print from rkq)
	  cnsReportTri2D(cns, nextOutputTime, newOptions);

	  // restore time step
	  mesh->dt = savedt;

	  // increment next output time
	  nextOutputTime += outputInterval;

	  // accept saved rkq
	  cns->o_q.copyFrom(cns->o_saveq);
	}
	else{
	  // accept rkq
	  cns->o_q.copyFrom(cns->o_rkq);
	}

        time += mesh->dt;

        cns->facold = mymax(err,1E-4); // hard coded factor ?

	printf("\r time = %g (%d), dt = %g accepted (ratio dt/hmin = %g)               ", time, allStep, mesh->dt, mesh->dt/hmin);
        tstep++;
      } else {
        dtnew = mesh->dt/(mymax(cns->invfactor1,fac1/cns->safe));
	printf("\r time = %g (%d), dt = %g rejected (ratio dt/min = %g), trying %g", time, allStep, mesh->dt, mesh->dt/hmin, dtnew);

	done = 0;
      }
      mesh->dt = dtnew;
      
      
      allStep++;

    }

#if 0
    cns->o_rkA.free();
    cns->o_rkE.free();
#endif
    
    mesh->device.finish();
    
    double elapsed  = timer.toc("Run");

    printf("run took %lg seconds for %d accepted steps and %d total steps\n", elapsed, tstep, allStep);
    
  } else if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")) {

    for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

      dfloat time = tstep*mesh->dt;

      cnsLserkStepTri2D(cns, newOptions, time);
      
      if(((tstep+1)%mesh->errorStep)==0){
	time += mesh->dt;
        cnsReportTri2D(cns, time, newOptions);
      }
    }
  }
  
}

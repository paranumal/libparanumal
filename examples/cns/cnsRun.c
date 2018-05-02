#include "cns.h"

void cnsRun(cns_t *cns, setupAide &options){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = cns->mesh;

  cnsReport(cns, 0, options);

  occa::timer timer;
  
  timer.initTimer(mesh->device);

  timer.tic("Run");
  
  if (options.compareArgs("TIME INTEGRATOR","DOPRI5")) {
    int Nregect = 0;

    // hard code this for the moment
    dfloat outputInterval;
    options.getArgs("OUTPUT INTERVAL", outputInterval);
    
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
      cnsDopriStep(cns, options, time);

      // compute Dopri estimator
      dfloat err = cnsDopriEstimate(cns);
                                         
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
          cnsDopriStep(cns, options, time);    

          // shift for output
          cns->o_rkq.copyTo(cns->o_q);
          
          // output  (print from rkq)
          cnsReport(cns, nextOutputTime, options);

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
        tstep++;

        cns->facold = mymax(err,1E-4); // hard coded factor ?
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

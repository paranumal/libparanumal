#include "advectionQuad3D.h"

void advectionRunDOPRIQuad3D(solver_t *solver){

  mesh_t *mesh = solver->mesh;
    
  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  dfloat * test_q = (dfloat *) calloc(mesh->Nfields*mesh->Nelements*mesh->Np,sizeof(dfloat));
    
  //kernel arguments
  dfloat alpha = 1./mesh->N;

  occa::timer timer;
  
  timer.initTimer(mesh->device);

  timer.tic("Run");
  
  mesh->filterKernelq0H(mesh->Nelements,
			mesh->o_dualProjMatrix,
			mesh->o_cubeFaceNumber,
			mesh->o_EToE,
			mesh->o_qpre,
			mesh->o_qPreFilter);
  
  mesh->filterKernelq0V(mesh->Nelements,
			alpha,
			mesh->o_dualProjMatrix,
			mesh->o_cubeFaceNumber,
			mesh->o_EToE,
			mesh->o_x,
			mesh->o_y,
			mesh->o_z,
			mesh->o_qPreFilter,
			mesh->o_qpre);
  
  int lastStep = 0;
  
  while (1) {
    if (mesh->dt<solver->dtmin){
      printf("ERROR: Time step became too small at time step=%d\n", solver->tstep);
      exit (-1);
    }
    if (isnan(mesh->dt)) {
      printf("ERROR: Solution became unstable at time step=%d\n", solver->tstep);
      exit (-1);
    }

    int outputStep = 0;

    // check for next output
    if((solver->time+mesh->dt > solver->nextOutputTime) && (solver->time<=solver->nextOutputTime)) {
      mesh->dt = solver->nextOutputTime-solver->time;
      outputStep = 1;
    }

    //check for final timestep
    if (solver->time+mesh->dt > solver->finalTime) {
      mesh->dt = solver->finalTime-solver->time;
      lastStep = 1;
    }
    
    for (iint rk = 0; rk < mesh->Nrk; ++rk) {

      dfloat t = solver->time + solver->rkc[rk]*mesh->dt;
      mesh->rkStageKernel(mesh->Nelements,
			  rk,
			  mesh->dt,
			  mesh->o_rka,
			  mesh->o_q,
			  mesh->o_resq,
			  mesh->o_rkq);

      // compute volume contribution to DG advection RHS
      mesh->volumeKernel(mesh->Nelements,
			 mesh->o_vgeo,
			 mesh->o_D,
			 mesh->o_x,
			 mesh->o_y,
			 mesh->o_z,
			 mesh->o_rkq,
			 mesh->o_prerhsq);

      mesh->surfaceKernel(mesh->Nelements,
			  mesh->o_sgeo,
			  mesh->o_LIFTT,
			  mesh->o_vmapM,
			  mesh->o_vmapP,
			  solver->tstep,
			  mesh->o_x,
			  mesh->o_y,
			  mesh->o_z,
			  mesh->o_rkq,
			  mesh->o_prerhsq);

      mesh->filterKernelq0H(mesh->Nelements,
			    mesh->o_dualProjMatrix,
			    mesh->o_cubeFaceNumber,
			    mesh->o_EToE,
			    mesh->o_rkq,
			    mesh->o_qPreFilter);
      
      mesh->filterKernelq0V(mesh->Nelements,
			    alpha,
			    mesh->o_dualProjMatrix,
			    mesh->o_cubeFaceNumber,
			    mesh->o_EToE,
			    mesh->o_x,
			    mesh->o_y,
			    mesh->o_z,
			    mesh->o_qPreFilter,
			    mesh->o_rkq);
      
      mesh->volumeCorrectionKernel(mesh->Nelements,
				   mesh->o_rkq,
				   mesh->o_qPreCorr);

      mesh->updateKernel(mesh->Nelements, 
                          rk,
                          mesh->dt,
                          mesh->o_rka,
                          mesh->o_rkb,
                          mesh->o_q,
                          mesh->o_rhsq,
			  mesh->o_qPreCorr,
                          mesh->o_resq, 
                          mesh->o_rkq,
                          mesh->o_rkerr);

    }
    mesh->rkErrorEstimateKernel(solver->Ntotal, 
			       solver->absTol,
			       solver->relTol,
			       mesh->o_q,
			       mesh->o_rkq,
			       mesh->o_rkerr,
			       mesh->o_errtmp);
    
    mesh->o_errtmp.copyTo(solver->errtmp);
    dfloat localerr = 0;
    dfloat err = 0;
    for(iint n=0;n<solver->Nblock;++n){
      localerr += solver->errtmp[n];
    }
    
    err = sqrt(localerr/mesh->Nelements);
    dfloat fac1 = pow(err,solver->exp1);
    dfloat fac = fac1/pow(solver->oldFactor,solver->beta);
    
    fac = mymax(solver->invfactor2, mymin(solver->invfactor1,fac/solver->safety));
    dfloat dtnew = mesh->dt/fac;

    if (err<1.0) { //dt is accepted
      solver->time += mesh->dt;

      solver->oldFactor = mymax(err,1E-4);

      mesh->o_q.copyFrom(mesh->o_rkq);
      //printf("dt = %g accepted\n", mesh->dt);
      if(outputStep){
	outputStep = 0;
	solver->nextOutputTime += solver->outputInterval;

	printf("tstep = %d, t = %g\n", solver->tstep, solver->time);
	fflush(stdout);

	mesh->o_q.copyTo(mesh->q);
	advectionPlotVTUQuad3DV2(mesh, "foo", solver->time/solver->outputInterval);
      }

      solver->tstep++;
      //if the last step was accepted, stop computing
      if (lastStep) break;
    }
    else {
      dtnew = mesh->dt/(mymax(solver->invfactor1,fac1/solver->safety));
      //printf("dtnew factors = %lf %lf %lf\n",solver->invfactor1,fac1,solver->safety);
      //printf("dt = %g rejected, trying %g\n", mesh->dt, dtnew);
    }
    mesh->dt = dtnew;
    solver->allStep++;
    mesh->filterKernelq0H(mesh->Nelements,
			  mesh->o_dualProjMatrix,
			  mesh->o_cubeFaceNumber,
			  mesh->o_EToE,
			  mesh->o_q,
			  mesh->o_qPreFilter);
      
    mesh->filterKernelq0V(mesh->Nelements,
			  alpha,
			  mesh->o_dualProjMatrix,
			  mesh->o_cubeFaceNumber,
			  mesh->o_EToE,
			  mesh->o_x,
			  mesh->o_y,
			  mesh->o_z,
			  mesh->o_qPreFilter,
			  mesh->o_q);
  }
  
  mesh->device.finish();
  
  double elapsed  = timer.toc("Run");
  
  printf("run took %lg seconds for %d accepted steps and %d total steps\n", elapsed, solver->tstep, solver->allStep);
}

#include "advectionQuad3D.h"

void advectionRunDOPRIQuad3D(solver_t *solver){

  mesh_t *mesh = solver->mesh;
    
  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*solver->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  dfloat * test_q = (dfloat *) calloc(solver->Nfields*mesh->Nelements*mesh->Np,sizeof(dfloat));
    
  //kernel arguments
  dfloat alpha = 1./mesh->N;

  occa::timer timer;
  
  timer.initTimer(solver->device);

  timer.tic("Run");
  
  solver->filterKernelH(mesh->Nelements,
		      solver->o_dualProjMatrix,
		      solver->o_cubeFaceNumber,
		      solver->o_EToE,
		      solver->o_q,
		      solver->o_qFilter);
  
  solver->filterKernelV(mesh->Nelements,
			alpha,
			solver->o_dualProjMatrix,
			solver->o_cubeFaceNumber,
			solver->o_EToE,
			solver->o_x,
			solver->o_y,
			solver->o_z,
			solver->o_q,
			solver->o_qFilter,
			solver->o_qFiltered);
  solver->o_qFiltered.copyTo(solver->o_q);
  
  int lastStep = 0;
  
  while (1) {
    if (solver->dt<solver->dtmin){
      printf("ERROR: Time step became too small at time step=%d\n", solver->tstep);
      exit (-1);
    }
    if (isnan(solver->dt)) {
      printf("ERROR: Solution became unstable at time step=%d\n", solver->tstep);
      exit (-1);
    }

    int outputStep = 0;

    // check for next output
    /*    if((solver->time+solver->dt > solver->nextOutputTime) && (solver->time<=solver->nextOutputTime)) {
      solver->dt = solver->nextOutputTime-solver->time;
      outputStep = 1;
      }*/

    //check for final timestep
    if (solver->time+solver->dt > solver->finalTime) {
      solver->dt = solver->finalTime-solver->time;
      lastStep = 1;
    }
    
    for (iint rk = 0; rk < solver->Nrk; ++rk) {

      dfloat t = solver->time + solver->rkC[rk]*solver->dt;
      solver->rkStageKernel(mesh->Nelements,
			  rk,
			  solver->dt,
			  solver->o_rkA,
			  solver->o_q,
			  solver->o_resq,
			  solver->o_rkq);

      
      // compute volume contribution to DG advection RHS
      solver->volumeKernel(mesh->Nelements,
			 solver->o_vgeo,
			 solver->o_D,
			 solver->o_x,
			 solver->o_y,
			 solver->o_z,
			 solver->o_rkq,
			 solver->o_rhsq);

      solver->surfaceKernel(mesh->Nelements,
			  solver->o_sgeo,
			  solver->o_LIFTT,
			  solver->o_vmapM,
			  solver->o_vmapP,
			  solver->tstep,
			  solver->o_x,
			  solver->o_y,
			  solver->o_z,
			  solver->o_rkq,
			  solver->o_rhsq);
      
      solver->filterKernelH(mesh->Nelements,
			    solver->o_dualProjMatrix,
			    solver->o_cubeFaceNumber,
			    solver->o_EToE,
			    solver->o_rhsq,
			    solver->o_qFilter);
      
      solver->filterKernelV(mesh->Nelements,
			    alpha,
			    solver->o_dualProjMatrix,
			    solver->o_cubeFaceNumber,
			    solver->o_EToE,
			    solver->o_x,
			    solver->o_y,
			    solver->o_z,
			    solver->o_rhsq,
			    solver->o_qFilter,
			    solver->o_qFiltered);
      solver->o_qFiltered.copyTo(solver->o_rhsq);
      
      solver->volumeCorrectionKernel(mesh->Nelements,
				     solver->o_rkq,
				     solver->o_qCorr);
	 
            
      solver->updateKernel(mesh->Nelements, 
                          rk,
                          solver->dt,
                          solver->o_rkA,
                          solver->o_rkE,
                          solver->o_q,
                          solver->o_rhsq,
			  solver->o_qCorr,
                          solver->o_resq, 
                          solver->o_rkq,
                          solver->o_rkerr);
	
    }

    /*    solver->o_rkerr.copyTo(test_q);
    for (iint e = 0; e < mesh->Nelements; ++e) {
      for (iint f = 0; f < solver->Nfields; ++f) {
	for (iint n = 0; n < mesh->Np; ++n) {
	  if (mesh->vgeo[e*mesh->Nvgeo*mesh->Np + 10*mesh->Np + n]*test_q[e*solver->Nfields*mesh->Np + f*mesh->Np + n]*test_q[e*solver->Nfields*mesh->Np + f*mesh->Np + n] > 0.00000000001) printf("outlier %d\n",mesh->cubeDistance[e]);
	}
      }
      }*/
    
    solver->rkErrorEstimateKernel(solver->Ntotal, 
				solver->absTol,
				solver->relTol,
				solver->o_vgeo,
				solver->o_q,
				solver->o_rkq,
				solver->o_rkerr,
				solver->o_errtmp);
    
    solver->o_errtmp.copyTo(solver->errtmp);
    dfloat localerr = 0;
    dfloat err = 0;

    for(iint n=0;n<solver->Nblock;++n){
      localerr += solver->errtmp[n];
    }
    
    err = sqrt(localerr/mesh->Nelements);
    dfloat fac1 = pow(err,solver->exp1);
    dfloat fac = fac1/pow(solver->oldFactor,solver->beta);
    
    fac = mymax(solver->invfactor2, mymin(solver->invfactor1,fac/solver->safety));
    dfloat dtnew = solver->dt/fac;

    if (err<1.0) { //dt is accepted
      solver->time += solver->dt;

      solver->oldFactor = mymax(err,1E-4);

      solver->o_q.copyFrom(solver->o_rkq);
      //printf("\r dt = %g accepted                                      ", solver->dt);
      /*if(outputStep){
	outputStep = 0;
	solver->nextOutputTime += solver->outputInterval;

	printf("tstep = %d, t = %g\n", solver->tstep, solver->time);
	fflush(stdout);

	solver->o_q.copyTo(mesh->q);
	advectionPlotVTUQuad3DV2(mesh, "foo", solver->time/solver->outputInterval);
	}*/

      solver->tstep++;
      //if the last step was accepted, stop computing
      if (lastStep) break;
    }
    else {
      dtnew = solver->dt/(mymax(solver->invfactor1,fac1/solver->safety));
      //printf("dtnew factors = %lf %lf %lf\n",solver->invfactor1,fac1,solver->safety);
      //printf("\r dt = %g rejected, trying %g                              ", solver->dt, dtnew);
    }
    solver->dt = dtnew;
    solver->allStep++;
    solver->filterKernelH(mesh->Nelements,
			  solver->o_dualProjMatrix,
			  solver->o_cubeFaceNumber,
			  solver->o_EToE,
			  solver->o_q,
			  solver->o_qFilter);
      
    solver->filterKernelV(mesh->Nelements,
			  alpha,
			  solver->o_dualProjMatrix,
			  solver->o_cubeFaceNumber,
			  solver->o_EToE,
			  solver->o_x,
			  solver->o_y,
			  solver->o_z,
			  solver->o_q,
			  solver->o_qFilter,
			  solver->o_qFiltered);
  }

  double elapsed  = timer.toc("Run");
  
  //mesh->device.finish();
  
  
  printf("\n run took %lg seconds for %d accepted steps and %d total steps\n", elapsed, solver->tstep, solver->allStep);
}

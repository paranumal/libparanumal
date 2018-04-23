#include "advectionQuad3D.h"

void advectionRunMRSAABQuad3D(solver_t *solver){

  mesh_t *mesh = solver->mesh;
  
  occa::initTimer(solver->device);
  
  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*solver->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  dfloat * test_q = (dfloat *) calloc(mesh->Nelements*mesh->Np*solver->Nfields*solver->Nrhs,sizeof(dfloat));
  
  /*for (int e = 0; e < mesh->Nelements; ++e) {
    for (int f = 0; f < mesh->Nfaces; ++f) {
    printf("%d ",mesh->EToF[e*mesh->Nfaces + f]);
    }
    printf("\n");
    }*/
  
  //kernel arguments
  dfloat alpha = 1./mesh->N;
  
  for(iint tstep=solver->Nrhs;tstep<solver->NtimeSteps;++tstep){
    
    for (iint Ntick=0; Ntick < pow(2,mesh->MRABNlevels-1);Ntick++) {
      
      iint mrab_order=solver->Nrhs-1;
      /*      if (tstep - mesh->Nrhs == 0) mrab_order = 0;
      else if (tstep - mesh->Nrhs == 1) mrab_order = 1;
      else mrab_order = 2;*/
      
      //synthesize actual stage time
      iint t = tstep*pow(2,mesh->MRABNlevels-1) + Ntick;

      iint lev;
      for (lev=0;lev<mesh->MRABNlevels;lev++)
	if (Ntick % (1<<lev) != 0) break;

      iint levS;
      for (levS=0;levS<mesh->MRABNlevels;levS++)
        if ((Ntick+1) % (1<<levS) !=0) break; //find the max lev to update
      
      for (iint l=0;l<lev;l++) {
	if (mesh->MRABNelements[l]) {
	  // compute volume contribution to DG boltzmann RHS
	  solver->volumeKernel(mesh->MRABNelements[l],
			     solver->o_MRABelementIds[l],
			     solver->MRABshiftIndex[l],
			     solver->o_vgeo,
			     solver->o_D,
			     solver->o_x,
			     solver->o_y,
			     solver->o_z,
			     solver->o_q,
			     solver->o_rhsq);
	}
      }
      
      occa::tic("surfaceKernel");
      
      for (iint l=0;l<lev;l++) {
	if (mesh->MRABNelements[l]) {

	  solver->surfaceKernel(mesh->MRABNelements[l],
			      solver->o_MRABelementIds[l],
			      solver->MRABshiftIndex[l],
			      solver->o_sgeo,
			      solver->o_LIFTT,
			      solver->o_vmapM,
			      solver->o_vmapP,
			      t,
			      solver->o_x,
			      solver->o_y,
			      solver->o_z,
			      solver->o_fQ,
			      solver->o_rhsq);
	  solver->lev_updates[l] = Ntick;
	}
      }
      occa::toc("surfaceKernel");
      
      solver->o_shift.copyFrom(solver->MRABshiftIndex);
      solver->o_lev_updates.copyFrom(solver->lev_updates);

      for (iint l = 0; l < lev; l++) {
		
      	solver->filterKernelH(mesh->MRABNelements[l],
			    solver->o_MRABelementIds[l],
			    solver->o_shift,
			    solver->o_dualProjMatrix,
			    solver->o_cubeFaceNumber,
			    solver->o_EToE,
			    solver->o_lev_updates,
			    solver->o_MRABlevels,
			    l,
			    solver->o_rhsq,
			    solver->o_qFilter);
      }
      for (iint l = 0; l < lev; l++) {
	
	solver->filterKernelV(mesh->MRABNelements[l],
			    solver->o_MRABelementIds[l],
			    solver->o_shift,
			    alpha,
			    solver->o_dualProjMatrix,
			    solver->o_cubeFaceNumber,
			    solver->o_EToE,
			    solver->o_x,
			    solver->o_y,
			    solver->o_z,
			    solver->o_lev_updates,
			    solver->o_MRABlevels,
			    l,
			    solver->o_rhsq,
			    solver->o_qFilter,
			    solver->o_qFiltered);	
			    }
      
      for (iint l=0;l<lev;l++) {
	if (mesh->MRABNelements[l]) {
	  solver->volumeCorrectionKernel(mesh->MRABNelements[l],
				       solver->o_MRABelementIds[l],
				       solver->MRABshiftIndex[l],
				       solver->o_q,
				       solver->o_qCorr);
	}
      }
      
      for (iint l = 0; l < levS; l++) {
	const iint id = mrab_order*mesh->MRABNlevels*solver->Nrhs + l*solver->Nrhs;
	occa::tic("updateKernel");
	
	if (mesh->MRABNelements[l]) {
	  solver->updateKernel(mesh->MRABNelements[l],
			     solver->o_MRABelementIds[l],
			     solver->MRSAAB_C[l],
			     solver->MRAB_A[id+0],
			     solver->MRAB_A[id+1],
			     solver->MRAB_A[id+2],
			     solver->MRAB_A[id+3],
			     solver->MRSAAB_A[id+0],
			     solver->MRSAAB_A[id+1],
			     solver->MRSAAB_A[id+2],
			     solver->MRSAAB_A[id+3],
			     solver->MRABshiftIndex[l],
			     solver->o_qFiltered,
			       //solver->o_rhsq,
			     solver->o_fQ,
			     solver->o_qCorr,
			     solver->o_q);
	
	  //we *must* use 2 here (n - 1), so rk coefficients point the right direction in time
	  solver->MRABshiftIndex[l] = (solver->MRABshiftIndex[l]+solver->Nrhs-1)%solver->Nrhs;
	}
      }
      
      occa::toc("updateKernel");
      
      if (levS<mesh->MRABNlevels) {
	const iint id = mrab_order*mesh->MRABNlevels*solver->Nrhs + levS*solver->Nrhs;
	
	if (mesh->MRABNhaloElements[levS]) {
	  solver->traceUpdateKernel(mesh->MRABNhaloElements[levS],
				  solver->o_MRABhaloIds[levS],
				  solver->MRSAAB_C[levS-1], //
				  solver->MRAB_B[id+0], //
				  solver->MRAB_B[id+1],
				  solver->MRAB_B[id+2],
				  solver->MRAB_B[id+3],//
				  solver->MRSAAB_B[id+0], //
				  solver->MRSAAB_B[id+1],
				  solver->MRSAAB_B[id+2],
				  solver->MRSAAB_B[id+3],
				  solver->MRABshiftIndex[levS],
				  solver->o_qFiltered,
				    //solver->o_rhsq,
				  solver->o_fQ,
				  solver->o_qPreFilter,
				  solver->o_qCorr,
				  solver->o_q);
      	}
      }
      
      for (iint l = 0; l < levS; ++l) {
	if (mesh->MRABNelements[l]) {
	  solver->filterKernelLevelsH(mesh->MRABNelements[l],
				solver->o_MRABelementIds[l],
				solver->o_dualProjMatrix,
				solver->o_cubeFaceNumber,
				solver->o_EToE,
				solver->o_fQ,
				solver->o_qPreFilter);
	}
      }
	
      for (iint l = 0; l < levS; ++l) {
	if (mesh->MRABNelements[l]) {
	  solver->filterKernelLevelsV(mesh->MRABNelements[l],
				    solver->o_MRABelementIds[l],
				    alpha,
				    solver->o_dualProjMatrix,
				    solver->o_cubeFaceNumber,
				    solver->o_EToE,
				    solver->o_x,
				    solver->o_y,
				    solver->o_z,
				    solver->o_qPreFilter,
				    solver->o_fQ,
				    solver->o_q);
	}
      }
    }

    /*    if (solver->NtimeSteps - (tstep + 1) < 20) {
      solver->o_q.copyTo(solver->q);
      dfloat t = solver->dt*(tstep+1)*pow(2,mesh->MRABNlevels-1);
      advectionErrorNormQuad3D(mesh,t,NULL,0);
      }*/
      
    
    // estimate maximum error
    if((((tstep+1)%solver->errorStep)==0)){
      //	dfloat t = (tstep+1)*solver->dt;
      dfloat t = solver->dt*((tstep+1)*pow(2,mesh->MRABNlevels-1));
      
      printf("tstep = %d, t = %g\n", tstep, t);
      fflush(stdout);
      // copy data back to host
      solver->o_q.copyTo(solver->q);

      // check for nans
      for(int n=0;n<solver->Nfields*mesh->Nelements*mesh->Np;++n){
	if(isnan(solver->q[n])){
	  printf("found nan\n");
	  exit(-1);
	}
      }
      /*
      advectionPlotNorms(mesh,"norms",tstep/solver->errorStep,solver->q);

      // output field files
      iint fld = 0;
      char fname[BUFSIZ];
      */
      //advectionPlotVTUQuad3DV2(mesh, "foo", tstep/solver->errorStep);
    }        
    occa::printTimer();
  }
  
  free(recvBuffer);
  free(sendBuffer);
}


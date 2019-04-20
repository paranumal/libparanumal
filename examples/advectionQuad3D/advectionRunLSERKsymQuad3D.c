#include "advectionQuad3D.h"

void advectionRunLSERKsymQuad3D(solver_t *solver,dfloat alpha_scale){

  mesh_t *mesh = solver->mesh;
    
  // MPI send buffer
  dfloat * test_q = (dfloat *) calloc(mesh->NgridElements*solver->Nfields*mesh->Np,sizeof(dfloat));
    
  //kernel arguments
  dfloat alpha = alpha_scale;

  iint Nboundary = mesh->NgridElements - mesh->Nelements;

  for(iint tstep=0;tstep < solver->NtimeSteps;++tstep){
	
      for (iint rk = 0; rk < solver->Nrk; ++rk) {
	
	//synthesize actual stage time
	dfloat t = tstep*solver->dt;

	
	// compute volume contribution to DG advection RHS
	/*	  	  solver->filterWeakKernelH(mesh->Nelements,
			      solver->o_dualTransMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_gridToE,
			      solver->o_cubeDistance,
			      solver->o_vgeo,
			      solver->o_qpre,
			      solver->o_qFilterw);
		
		  solver->filterWeakKernelV(mesh->Nelements,
			      alpha,
			      solver->o_dualTransMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_gridToE,
			      solver->o_x,
			      solver->o_y,
			      solver->o_z,
			      solver->o_cubeDistance,
			      solver->o_vgeo,
			      solver->o_qpre,
			      solver->o_qFilterw,
			      solver->o_qw);
	*/
		  
	// compute volume contribution to DG advection RHS
	solver->volumeKernel(mesh->Nelements,
			     solver->o_vgeo,
			     solver->o_weakD,
			     solver->o_x,
			     solver->o_y,
			     solver->o_z,
			       solver->o_mass,
			     solver->o_qpre,
			     //solver->o_qpre,
			     solver->o_qw,
			     solver->o_rhsqs,
			     solver->o_rhsqw
			     );
	  
	solver->surfaceKernel(mesh->Nelements,
			      solver->o_sgeo,
			      solver->o_vgeo,
			      solver->o_LIFTT,
			      solver->o_vmapM,
			      solver->o_vmapP,
			      t,
			      solver->o_x,
			      solver->o_y,
			      solver->o_z,
			      solver->o_qpre,
			      //solver->o_qpre,
			      solver->o_qw,
			      solver->o_rhsqs,
			      solver->o_rhsqw
			      );
	
	/*	solver->filterKernelH(mesh->Nelements,
			      solver->o_dualProjMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_gridToE,
			      solver->o_cubeDistance,
			      solver->o_vgeo,
			      solver->o_rhsqs,
			      solver->o_qFilters);
	
	solver->filterKernelV(mesh->Nelements,
			      alpha,
			      solver->o_dualProjMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_gridToE,
			      solver->o_x,
			      solver->o_y,
			      solver->o_z,
			      solver->o_cubeDistance,
			      solver->o_vgeo,
			      solver->o_rhsqs,
			      solver->o_qFilters,
			      solver->o_qs);
	*/	
	solver->volumeCorrectionKernel(mesh->Nelements,
				       solver->o_q,
				       solver->o_qCorr);

	/*solver->massMatrixKernel(mesh->Nelements,
				 solver->o_invmass,
				 solver->o_vgeo,
				 solver->o_rhsqs);
	*/
	solver->updateKernel(mesh->Nelements,
			     solver->dt,
			     solver->rka[rk],
			     solver->rkb[rk],
			     solver->o_rhsqs,
			     solver->o_rhsqw,
			     solver->o_qCorr,
			     solver->o_resq,
			     solver->o_qpre);

      }
	
      if (tstep == 0) {
	solver->o_qpre.copyTo(solver->q);
	advectionErrorNormQuad3D(solver,1*solver->dt,"start",0);
      }
  }
}

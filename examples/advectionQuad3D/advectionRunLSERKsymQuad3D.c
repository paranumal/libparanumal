#include "advectionQuad3D.h"

void advectionRunLSERKsymQuad3D(solver_t *solver,dfloat alpha_scale){

  mesh_t *mesh = solver->mesh;
    
  // MPI send buffer
  dfloat * test_q = (dfloat *) calloc(mesh->NgridElements*solver->Nfields*mesh->Np,sizeof(dfloat));
    
  //kernel arguments
  dfloat alpha = alpha_scale;

  iint Nboundary = mesh->NgridElements - mesh->Nelements;

  /*for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint n = 0; n < mesh->Np; ++n) {
      solver->q[e*mesh->Np*solver->Nfields + n] = mesh->x[e*mesh->Np + n];
    }
  }
  solver->o_qpre.copyFrom(solver->q);
  */

  /*if (Nboundary > 0) {
    solver->loadFilterGridKernel(Nboundary,
				 mesh->Nelements,
				 solver->o_rlocal,
				 solver->o_slocal,
				 solver->o_par_loc,
				 solver->o_perp_index,
				 solver->o_eInterp,
				 solver->o_overlapDirection,
				 solver->o_qpre);
  }
  
    solver->filterKernelH(mesh->Nelements,
			mesh->NgridElements,
			solver->o_dualProjMatrix,
			solver->o_cubeFaceNumber,
			solver->o_gridToE,
			solver->o_qpre,
			solver->o_qFilter);
  solver->filterKernelV(mesh->Nelements,
			mesh->NgridElements,
			alpha,
			solver->o_dualProjMatrix,
			solver->o_cubeFaceNumber,
			solver->o_gridToE,
			solver->o_x,
			solver->o_y,
			solver->o_z,
			solver->o_qpre,
			solver->o_qFilter,
			solver->o_q);
  
			solver->o_q.copyTo(solver->o_qpre);*/
  /*
  solver->o_q.copyTo(solver->q);
  advectionErrorNormQuad3D(solver,0,"basic",0);
  */

  /*  solver->filterKernelH(mesh->Nelements,
			solver->o_dualProjMatrix,
			solver->o_cubeFaceNumber,
			solver->o_EToE,
			solver->o_qpre,
			solver->o_qFilter);
      
  solver->filterKernelV(mesh->Nelements,
			alpha,
			solver->o_dualProjMatrix,
			solver->o_cubeFaceNumber,
			solver->o_EToE,
			solver->o_x,
			solver->o_y,
			solver->o_z,
			solver->o_qFilter,
			solver->o_qpre);
  */
  for(iint tstep=0;tstep < solver->NtimeSteps;++tstep){
	
      for (iint rk = 0; rk < solver->Nrk; ++rk) {
	
	//synthesize actual stage time
	dfloat t = tstep*solver->dt;
      
	// compute volume contribution to DG advection RHS
	solver->volumeKernel(mesh->Nelements,
			     solver->o_vgeo,
			     solver->o_D,
			     solver->o_weakD,
			     solver->o_x,
			     solver->o_y,
			     solver->o_z,
			     solver->o_qpre,
			     solver->o_qpre,
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
			      solver->o_qpre,
			      solver->o_rhsqs,
			      solver->o_rhsqw
			      );

	/*if (Nboundary > 0) {
	  solver->loadFilterGridKernel(Nboundary,
				       mesh->Nelements,
				       solver->o_rlocal,
				       solver->o_slocal,
				       solver->o_par_loc,
				       solver->o_perp_index,
				       solver->o_eInterp,
				       solver->o_overlapDirection,
				       solver->o_rhsq);
	}
	*/				     
	/*	solver->filterKernelH(mesh->Nelements,
			      solver->o_dualProjMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_gridToE,
			      solver->o_vgeo,
			      solver->o_cubeDistance,
			      solver->o_rhsqs,
			      //solver->o_rhsqw,
			      solver->o_qFilter);
		
	solver->filterKernelV(mesh->Nelements,
			      alpha,
			      solver->o_dualProjMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_gridToE,
			      solver->o_x,
			      solver->o_y,
			      solver->o_z,
			      solver->o_vgeo,
			      solver->o_cubeDistance,
			      solver->o_rhsqs,
			      //solver->o_rhsqw,
			      solver->o_qFilter,
			      solver->o_q);	
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

	/*	solver->filterKernelH(mesh->Nelements,
			      solver->o_dualProjMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_EToE,
			      solver->o_qpre,
			      solver->o_qFilter);
      
	solver->filterKernelV(mesh->Nelements,
			      alpha,
			      solver->o_dualProjMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_EToE,
			      solver->o_x,
			      solver->o_y,
			      solver->o_z,
			      solver->o_qpre,
			      solver->o_qFilter,
			      solver->o_q);
			      solver->o_q.copyTo(solver->o_qpre);*/
      }
	
      if (tstep == 0) {
	solver->o_qpre.copyTo(solver->q);
	advectionErrorNormQuad3D(solver,1*solver->dt,"start",0);
      }
  }
}

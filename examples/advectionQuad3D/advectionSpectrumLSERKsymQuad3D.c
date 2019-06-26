#include "advectionQuad3D.h"

void advectionSpectrumLSERKsymQuad3D(solver_t *solver,dfloat alpha_scale){

  mesh_t *mesh = solver->mesh;
    
  // MPI send buffer
  dfloat * test_qs = (dfloat *) calloc(mesh->NgridElements*solver->Nfields*mesh->Np,sizeof(dfloat));
  dfloat * test_qw = (dfloat *) calloc(mesh->NgridElements*solver->Nfields*mesh->Np,sizeof(dfloat));
    
  //kernel arguments
  dfloat alpha = alpha_scale;

  iint Nboundary = mesh->NgridElements - mesh->Nelements;

  FILE *matrix = fopen("./output_matrix","w");
  /*
    for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint n = 0; n < mesh->Np; ++n) {
      solver->q[e*mesh->Np*solver->Nfields + n] = 0;
    }
  }
  solver->o_qpre.copyFrom(solver->q);
  */
  
  //synthesize actual stage time
  dfloat t = 0;

  for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint f = 0; f < solver->Nfields; ++f) {
      for (iint n = 0; n < mesh->Np; ++n) {
	solver->q[e*mesh->Np*solver->Nfields + f*mesh->Np + n] = 0.;
      }
    }
  }
  solver->o_qpre.copyFrom(solver->q);
  
  for (iint e = 0; e < mesh->Nelements; ++e) {
      for (iint f = 0; f < 2; ++f) {
	  for (iint n = 0; n < mesh->Np; ++n) {
      
	  iint curr_pos = e*mesh->Np*solver->Nfields + f*mesh->Np + n;
	  solver->q[curr_pos] = 1.;
	  solver->o_qpre.copyFrom(solver->q);
	
	// compute volume contribution to DG advection RHS

	solver->filterWeakKernelH(mesh->Nelements,
			      solver->o_dualTransMatrix,
				  solver->o_invTransMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_gridToE,
			      solver->o_vgeo,
			      solver->o_cubeDistance,
			      solver->o_qpre,
			      solver->o_qFilter);
		
	solver->filterWeakKernelV(mesh->Nelements,
			      alpha,
			      solver->o_dualTransMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_gridToE,
			      solver->o_x,
			      solver->o_y,
			      solver->o_z,
			      solver->o_vgeo,
			      solver->o_cubeDistance,
			      solver->o_qpre,
			      solver->o_qFilter,
			      solver->o_qFiltered);
	  
	  solver->volumeKernel(mesh->Nelements,
			     solver->o_vgeo,
			     solver->o_D,
			     solver->o_weakD,
			     solver->o_x,
			     solver->o_y,
			     solver->o_z,
			     solver->o_qpre,
			     solver->o_qFiltered,
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
			      solver->o_qFiltered,
			      solver->o_rhsqs,
			      solver->o_rhsqw
			      );

	/*solver->loadFilterGridKernel(Nboundary,
				     mesh->Nelements,
				     solver->o_rlocal,
				     solver->o_slocal,
				     solver->o_par_loc,
				     solver->o_perp_index,
				     solver->o_eInterp,
				     solver->o_overlapDirection,
				     solver->o_rhsq);
	*/	
	solver->filterKernelH(mesh->Nelements,
			      solver->o_dualProjMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_gridToE,
			      solver->o_vgeo,
			      solver->o_cubeDistance,
			      solver->o_rhsqs,
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
			      solver->o_qFilter,
			      solver->o_q);
	
	
	/*	solver->massMatrixKernel(mesh->Nelements,
				 solver->o_invmass,
				 solver->o_vgeo,
				 solver->o_rhsqs);
	*/
	solver->q[curr_pos] = 0.;
	
	solver->o_q.copyTo(test_qs);
	solver->o_rhsqw.copyTo(test_qw);
	
	for (iint es = 0; es < mesh->Nelements; ++es) {
	  for (iint fs = 0; fs < 2; ++fs) {
	    for (iint ns = 0; ns < mesh->Np; ++ns) {
		fprintf(matrix,"%17.15g ",0.5*(test_qs[es*mesh->Np*solver->Nfields + fs*mesh->Np + ns]+test_qw[es*mesh->Np*solver->Nfields + fs*mesh->Np + ns]));
	    }
	  }
	}
	fprintf(matrix,"\n");
	  }
      }
  }
}
  

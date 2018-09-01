#include "advectionQuad3D.h"

void advectionSpectrumLSERKQuad3D(solver_t *solver,dfloat alpha_scale){

  mesh_t *mesh = solver->mesh;
    
  // MPI send buffer
  dfloat * test_q = (dfloat *) calloc(mesh->NgridElements*solver->Nfields*mesh->Np,sizeof(dfloat));
    
  //kernel arguments
  dfloat alpha = alpha_scale;

  iint Nboundary = mesh->NgridElements - mesh->Nelements;

  FILE *matrix = fopen("./output_matrix","w");
  double output_tol = 1e-18;
  
  /*for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint n = 0; n < mesh->Np; ++n) {
      solver->q[e*mesh->Np*solver->Nfields + n] = mesh->x[e*mesh->Np + n];
    }
  }
  solver->o_qpre.copyFrom(solver->q);
  */
  
  //synthesize actual stage time
  dfloat t = 0;

  for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint n = 0; n < mesh->Np; ++n) {
      solver->q[e*mesh->Np*solver->Nfields + n] = 0.;
    }
  }
  solver->o_qpre.copyFrom(solver->q);
  
  for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint n = 0; n < mesh->Np; ++n) {

      iint curr_pos = e*mesh->Nelements*solver->Nfields + n;
      solver->q[curr_pos] = 1.;
      solver->o_qpre.copyFrom(solver->q);
      
      // compute volume contribution to DG advection RHS
      solver->volumeKernel(mesh->Nelements,
			   solver->o_vgeo,
			   solver->o_D,
			   solver->o_x,
			   solver->o_y,
			   solver->o_z,
			   solver->o_qpre,
			   solver->o_rhsq);
      
      solver->surfaceKernel(mesh->Nelements,
			    solver->o_sgeo,
			    solver->o_LIFTT,
			    solver->o_vmapM,
			    solver->o_vmapP,
			    t,
			    solver->o_x,
			    solver->o_y,
			    solver->o_z,
			    solver->o_qpre,
			    solver->o_rhsq);
      
      /*      solver->loadFilterGridKernel(Nboundary,
                                   mesh->Nelements,
				   solver->o_rlocal,
				   solver->o_slocal,
				   solver->o_par_loc,
				   solver->o_perp_index,
				   solver->o_eInterp,
				   solver->o_overlapDirection,
				   solver->o_rhsq);
  
      solver->filterKernelH(mesh->Nelements,
                            mesh->NgridElements,
			    solver->o_dualProjMatrix,
			    solver->o_cubeFaceNumber,
			    solver->o_gridToE,
			    solver->o_rhsq,
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
			    solver->o_rhsq,
			    solver->o_qFilter,
			    solver->o_q);
      
      */
      solver->q[curr_pos] = 0.;

      solver->o_rhsq.copyTo(test_q);

      for (iint e = 0; e < mesh->Nelements; ++e) {
	for (iint n = 0; n < mesh->Np; ++n) {
	  if (test_q[e*mesh->Np*solver->Nfields + n] < output_tol) 
	    fprintf(matrix,"0");
	  else
	    fprintf(matrix,"%g",test_q[e*mesh->Np*solver->Nfields + n]);
	  fprintf(matrix," ");
	}
      }
      fprintf(matrix,"\n");
    }
  }
}

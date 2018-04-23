#include "advectionQuad3D.h"

void advectionRunLSERKbasicQuad3D(solver_t *solver){

  mesh_t *mesh = solver->mesh;
    
  // MPI send buffer
  dfloat * test_q = (dfloat *) calloc(mesh->Nelements*solver->Nfields*mesh->Np,sizeof(dfloat));
    
  //kernel arguments
  dfloat alpha = 1./mesh->N;

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
			solver->o_qFilter,
			solver->o_q);

  solver->o_q.copyTo(solver->o_qpre);
  
  for(iint tstep=0;tstep < solver->NtimeSteps;++tstep){
	
      for (iint rk = 0; rk < solver->Nrk; ++rk) {
	
	//synthesize actual stage time
	dfloat t = tstep*solver->dt;
      
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
			      solver->o_qFilter,
			      solver->o_rhsq);
	
	solver->volumeCorrectionKernel(mesh->Nelements,
				       solver->o_qpre,
				       solver->o_qCorr);
	
	solver->updateKernel(mesh->Nelements,
			     solver->dt,
			     solver->rka[rk],
			     solver->rkb[rk],
			     solver->o_rhsq,
			     solver->o_qCorr,
			     solver->o_resq,
			     solver->o_q,
			     solver->o_qpre);
      }
      solver->filterKernelH(mesh->Nelements,
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
      solver->o_q.copyFrom(solver->o_qpre);
  }
}

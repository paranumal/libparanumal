#include "advectionQuad3D.h"

void advectionSpectrumLSERKQuad3D(solver_t *solver,dfloat alpha_scale){

  mesh_t *mesh = solver->mesh;
    
  // MPI send buffer
  dfloat * test_qs = (dfloat *) calloc(mesh->NgridElements*solver->Nfields*mesh->Np,sizeof(dfloat));
  dfloat * test_qw = (dfloat *) calloc(mesh->NgridElements*solver->Nfields*mesh->Np,sizeof(dfloat));
    
  //kernel arguments
  dfloat alpha = alpha_scale;

  iint Nboundary = mesh->NgridElements - mesh->Nelements;

  FILE *matrix = fopen("./output_matrix","w");
  
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

  occa::kernel volumeKernel =
    solver->device.buildKernelFromSource("./advectionSpectrumKernelsQuad3D.okl",
					 "advectionVolumeQuad3D",
					 solver->kernelInfo);

  occa::kernel surfaceKernel =
    solver->device.buildKernelFromSource("./advectionSpectrumKernelsQuad3D.okl",
					 "advectionSurfaceQuad3D",
					 solver->kernelInfo);

  occa::kernel filterKernelH =
    solver->device.buildKernelFromSource("./advectionSpectrumKernelsQuad3D.okl",
					 "advectionFilterHsymQuad3D",
					 solver->kernelInfo);

  occa::kernel filterKernelV =
    solver->device.buildKernelFromSource("./advectionSpectrumKernelsQuad3D.okl",
					 "advectionFilterVsymQuad3D",
					 solver->kernelInfo);

  
  printf("Nfields = %d\n", solver->Nfields);

  int Nfields = 2;
  
  dfloat *h_q  = (dfloat*) calloc(mesh->Np*Nfields*mesh->Nelements,sizeof(dfloat));
  dfloat *h_Aq = (dfloat*) calloc(mesh->Np*Nfields*mesh->Nelements,sizeof(dfloat));
  dfloat *h_Fq = (dfloat*) calloc(mesh->Np*Nfields*mesh->Nelements,sizeof(dfloat));
  dfloat *h_Sq = (dfloat*) calloc(mesh->Np*Nfields*mesh->Nelements,sizeof(dfloat));
  dfloat *h_FHq = (dfloat*) calloc(mesh->Np*Nfields*mesh->Nelements,sizeof(dfloat));

  occa::memory o_q   = solver->device.malloc(mesh->Np*Nfields*mesh->Nelements*sizeof(dfloat), h_q);
  occa::memory o_Aq  = solver->device.malloc(mesh->Np*Nfields*mesh->Nelements*sizeof(dfloat), h_q);
  occa::memory o_Sq  = solver->device.malloc(mesh->Np*Nfields*mesh->Nelements*sizeof(dfloat), h_q);
  occa::memory o_Fq  = solver->device.malloc(mesh->Np*Nfields*mesh->Nelements*sizeof(dfloat), h_q);
  occa::memory o_FHq = solver->device.malloc(mesh->Np*Nfields*mesh->Nelements*sizeof(dfloat), h_q);

  FILE *fpA = fopen("transA.dat", "w");
  FILE *fpF = fopen("transF.dat", "w");
  FILE *fpS = fopen("transS.dat", "w");
  FILE *fpM = fopen("M.dat", "w");

  for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint f = 0; f < Nfields; ++f) {
      for (iint n = 0; n < mesh->Np; ++n) {
	
	iint curr_pos = e*mesh->Np*Nfields + f*mesh->Np + n;

	h_q[curr_pos] = 1.;
	o_q.copyFrom(h_q);
	h_q[curr_pos] = 0;
	o_Sq.copyFrom(h_q);
		
	// compute volume contribution to DG advection RHS
	volumeKernel(mesh->Nelements,
		     solver->o_vgeo,
		     solver->o_D,
		     solver->o_x,
		     solver->o_y,
		     solver->o_z,
		     o_q, 
		     o_Aq);
	
	surfaceKernel(mesh->Nelements,
		      solver->o_sgeo,
		      solver->o_LIFTT,
		      solver->o_vmapM,
		      solver->o_vmapP,
		      t,
		      solver->o_x,
		      solver->o_y,
		      solver->o_z,
		      o_q,
		      o_Aq,
		      o_Sq);

	// use filter action to compute filter matrix column
	filterKernelH(mesh->Nelements,
		      solver->o_dualProjMatrix,
		      solver->o_cubeFaceNumber,
		      solver->o_gridToE,
		      solver->o_vgeo,
		      solver->o_cubeDistance,
		      o_q,
		      o_FHq);
		
	filterKernelV(mesh->Nelements,
		      alpha,
		      solver->o_dualProjMatrix,
		      solver->o_cubeFaceNumber,
		      solver->o_gridToE,
			      solver->o_x,
		      solver->o_y,
		      solver->o_z,
		      solver->o_vgeo,
		      solver->o_cubeDistance,
		      o_q,
		      o_FHq,
		      o_Fq);
	
	o_Aq.copyTo(h_Aq);
	o_Sq.copyTo(h_Sq);
	o_Fq.copyTo(h_Fq);

	for(iint m=0;m<mesh->Nelements*mesh->Np*Nfields;++m){
	  fprintf(fpA, "%17.15lg ", h_Aq[m]);
	  fprintf(fpF, "%17.15lg ", h_Fq[m]);
	  fprintf(fpS, "%17.15lg ", h_Sq[m]);
	}
	fprintf(fpA, "\n");
	fprintf(fpF, "\n");
	fprintf(fpS, "\n");

	fprintf(fpM, "%17.15lg\n", mesh->vgeo[e*mesh->Np*mesh->Nvgeo+n+JWID*mesh->Np]);
      }
    }
  }
  
  fclose(fpA);
  fclose(fpF);
  fclose(fpM);

  exit(0);




  
  
  for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint f = 0; f < 2; ++f) {
      for (iint n = 0; n < mesh->Np; ++n) {
      
	iint curr_pos = e*mesh->Np*solver->Nfields + f*mesh->Np + n;
	solver->q[curr_pos] = 1.;
	solver->o_qpre.copyFrom(solver->q);
	
	// compute volume contribution to DG advection RHS
	solver->volumeKernel(mesh->Nelements,
			     solver->o_vgeo,
			     solver->o_D,
			     solver->o_weakD,
			     solver->o_x,
			     solver->o_y,
			     solver->o_z,
			     solver->o_mass,
			     solver->o_qpre,
			     solver->o_rhsqs,
			     solver->o_rhsqw
			     );
	
	solver->surfaceKernel(mesh->Nelements,
			      solver->o_vgeo,
			      solver->o_sgeo,
			      solver->o_LIFTT,
			      solver->o_vmapM,
			      solver->o_vmapP,
			      t,
			      solver->o_x,
			      solver->o_y,
			      solver->o_z,
			      solver->o_qpre,
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
	
	
	/*	solver->massMatrixKernel(mesh->Nelements,
		solver->o_invmass,
		solver->o_vgeo,
		solver->o_rhsqs);
	*/
	solver->q[curr_pos] = 0.;
	
	/*	solver->o_rhsqs.copyTo(test_qs);
		solver->o_rhsqw.copyTo(test_qw);*/

	solver->o_q.copyTo(test_qs);
	//solver->o_rhsqw.copyTo(test_qw);
	
	for (iint es = 0; es < mesh->Nelements; ++es) {
	  for (iint fs = 0; fs < 2; ++fs) {
	    for (iint ns = 0; ns < mesh->Np; ++ns) {
	      //fprintf(matrix,"%17.15g ",test_qs[es*mesh->Np*solver->Nfields + fs*mesh->Np + ns] + test_qw[es*mesh->Np*solver->Nfields + fs*mesh->Np + ns]);
	      fprintf(matrix,"%17.15g ",(test_qs[es*mesh->Np*solver->Nfields + fs*mesh->Np + ns]/*+test_qw[es*mesh->Np*solver->Nfields + fs*mesh->Np + ns]*/)/*/mesh->vgeo[es*mesh->Np*mesh->Nvgeo + 10*mesh->Np + ns]*//*mesh->vgeo[es*mesh->Np*mesh->Nvgeo + 9*mesh->Np + ns]*/);
	    }
	  }
	}
	fprintf(matrix,"\n");
	/*solver->o_rhsq.copyTo(solver->q);
	  for (iint es = 0; es < mesh->Nelements; ++es) {
	  for (iint ns = 0; ns < mesh->Np; ++ns) {
	  solver->q[es*mesh->Np*solver->Nfields + ns] -= 2*mesh->y[es*mesh->Np + ns];
	  }
	  }*/
	//advectionPlotVTUQuad3DV2(solver, "der", 0);
      }
    }
  }
}
  

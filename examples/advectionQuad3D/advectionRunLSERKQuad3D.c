#include "advectionQuad3D.h"

void advectionRunLSERKQuad3D(solver_t *solver){

  mesh_t *mesh = solver->mesh;
    
  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*solver->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  dfloat * test_q = (dfloat *) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
    
  //kernel arguments
  dfloat alpha = 1./mesh->N;

  solver->filterKernelq0H(mesh->Nelements,
			  solver->o_dualProjMatrix,
			  solver->o_cubeFaceNumber,
			  solver->o_EToE,
			  solver->o_qpre,
			  solver->o_qPreFilter);
  
    solver->filterKernelq0V(mesh->Nelements,
			    alpha,
			    solver->o_dualProjMatrix,
			    solver->o_cubeFaceNumber,
			    solver->o_EToE,
			    solver->o_x,
			    solver->o_y,
			    solver->o_z,
			    solver->o_qpre,
			    solver->o_qPreFilter,
			    solver->o_qPreFiltered);
    
    solver->o_qPreFiltered.copyTo(solver->o_qpre);
  
  for(iint tstep=0;tstep < solver->Nrhs;++tstep){
    for (iint Ntick=0; Ntick < pow(2,mesh->MRABNlevels-1);Ntick++) {      
      
      iint lev;
      for (lev=0;lev<mesh->MRABNlevels;lev++)
	if (Ntick % (1<<lev) !=0) break; //find the max lev to add to rhsq

      iint levS;
      for (levS=0;levS<mesh->MRABNlevels;levS++)
	if ((Ntick+1) % (1<<levS) !=0) break; //find the max lev to add to rhsq
	
      for (iint rk = 0; rk < solver->Nrk; ++rk) {
	
	//synthesize actual stage time
	dfloat t = tstep*pow(2,mesh->MRABNlevels-1) + Ntick;
      
	if(mesh->totalHaloPairs>0){
	  // extract halo on DEVICE
	  iint Nentries = mesh->Np*solver->Nfields;
	  
	  solver->haloExtractKernel(mesh->totalHaloPairs,
				  Nentries,
				  solver->o_haloElementList,
				  solver->o_qpre,
				  solver->o_haloBuffer);

	  // copy extracted halo to HOST 
	  solver->o_haloBuffer.copyTo(sendBuffer);      
	
	  // start halo exchange
	  meshHaloExchangeStart(mesh,
				mesh->Np*mesh->Nfields*sizeof(dfloat),
				sendBuffer,
				recvBuffer);
	}
	// compute volume contribution to DG advection RHS
	solver->volumePreKernel(mesh->Nelements,
			      solver->o_vgeo,
			      solver->o_D,
			      solver->o_x,
			      solver->o_y,
			      solver->o_z,
			      solver->o_qpre,
			      solver->o_prerhsq);

	if(mesh->totalHaloPairs>0){
	  // wait for halo data to arrive
	  meshHaloExchangeFinish(mesh);
	
	  // copy halo data to DEVICE
	  size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
	  solver->o_q.copyFrom(recvBuffer, haloBytes, offset);
	  solver->device.finish();
	}
	
	solver->surfacePreKernel(mesh->Nelements,
			       solver->o_sgeo,
			       solver->o_LIFTT,
			       solver->o_vmapM,
			       solver->o_vmapP,
			       t,
			       solver->o_x,
			       solver->o_y,
			       solver->o_z,
			       solver->o_qpre,
			       solver->o_prerhsq);
	
	/*	for (iint l = 0; l < mesh->MRABNlevels; ++l) {
	  iint saved = (l < lev)&&(rk == 0);
	  solver->filterKernelHLSERK(mesh->MRABNelements[l],
				   solver->o_MRABelementIds[l],
				   saved,
				   mesh->MRABshiftIndex[l],
				   solver->o_dualProjMatrix,
				   solver->o_cubeFaceNumber,
				   solver->o_EToE,
				   solver->o_prerhsq,
				   solver->o_qPreFilter,
				   solver->o_rhsq);
	}

	for (iint l = 0; l < mesh->MRABNlevels; ++l) {
	  iint saved = (l < lev)&&(rk == 0);
	  solver->filterKernelVLSERK(mesh->MRABNelements[l],
				   solver->o_MRABelementIds[l],
				   alpha,
				   saved,
				   solver->MRABshiftIndex[l],
				   solver->o_dualProjMatrix,
				   solver->o_cubeFaceNumber,
				   solver->o_EToE,
				   solver->o_x,
				   solver->o_y,
				   solver->o_z,
				   solver->o_qPreFilter,
				   solver->o_prerhsq,
				   solver->o_qFilter);
				   }
	*/
	for (iint l = 0; l < mesh->MRABNlevels; l++) {
	  iint saved = (l < lev)&&(rk == 0);
	  solver->volumeCorrPreKernel(mesh->MRABNelements[l],
				    solver->o_MRABelementIds[l],
				    saved,
				    solver->MRABshiftIndex[l],
				    solver->o_qpre,
				    solver->o_qPreCorr,
				    solver->o_qCorr);
	}
	
	for (iint l = 0; l < mesh->MRABNlevels; l++) {
	  iint saved = (l < lev)&&(rk == 0);
	  
	  if (mesh->MRABNelements[l]) {
	    solver->updatePreKernel(mesh->MRABNelements[l],
				  solver->o_MRABelementIds[l],
				  saved,
				  mesh->dt,
				  solver->rka[rk],
				  solver->rkb[rk],
				  solver->MRABshiftIndex[l],
				  solver->o_prerhsq,
				    //solver->o_qFiltered,
				  solver->o_rhsq,
				  solver->o_qPreCorr,
				  solver->o_resq,
				  solver->o_qpre);
	    
	    saved = (l < levS)&&(rk == solver->Nrk-1);
	    if (saved)
	      solver->MRABshiftIndex[l] = (solver->MRABshiftIndex[l]+solver->Nrhs-1)%solver->Nrhs;
	  }
	}
      }
      solver->filterKernelq0H(mesh->Nelements,
			      solver->o_dualProjMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_EToE,
			      solver->o_qpre,
			      solver->o_qPreFilter);
      
      solver->filterKernelq0V(mesh->Nelements,
			      alpha,
			      solver->o_dualProjMatrix,
			      solver->o_cubeFaceNumber,
			      solver->o_EToE,
			      solver->o_x,
			      solver->o_y,
			      solver->o_z,
			      solver->o_qpre,
			      solver->o_qPreFilter,
			      solver->o_qPreFiltered);

      solver->o_qPreFiltered.copyTo(solver->o_qpre);
      
    }
  }
  
  solver->traceUpdatePreKernel(mesh->Nelements,
			     solver->o_vmapM,
			     solver->o_fQ,
			     solver->o_qpre,
			     solver->o_q);

  
  free(recvBuffer);
  free(sendBuffer);
}

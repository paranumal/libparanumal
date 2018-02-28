#include "boltzmannQuad3D.h"

void boltzmannRunMRSAABQuad3D(solver_t *solver){

  mesh_t *mesh = solver->mesh;
  
  occa::initTimer(mesh->device);
  
  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  dfloat * test_q = (dfloat *) calloc(mesh->Nelements*mesh->Np*mesh->Nfields*mesh->Nrhs,sizeof(dfloat));
  
  /*for (int e = 0; e < mesh->Nelements; ++e) {
    for (int f = 0; f < mesh->Nfaces; ++f) {
      printf("%d ",mesh->EToF[e*mesh->Nfaces + f]);
    }
    printf("\n");
  }*/
  
  //kernel arguments
  dfloat alpha = 1./mesh->N;

  //filter the initial state
  mesh->filterKernelq0H(mesh->Nelements,
			alpha,
			mesh->o_dualProjMatrix,
			mesh->o_cubeFaceNumber,
			mesh->o_EToE,
			mesh->o_x,
			mesh->o_y,
			mesh->o_z,
			mesh->o_q,
			mesh->o_qFilter);
  mesh->filterKernelq0V(mesh->Nelements,
			alpha,
			mesh->o_dualProjMatrix,
			mesh->o_cubeFaceNumber,
			mesh->o_EToE,
			mesh->o_x,
			mesh->o_y,
			mesh->o_z,
			mesh->o_qFilter,
			mesh->o_q);
  
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){
     for (iint Ntick=0; Ntick < pow(2,mesh->MRABNlevels-1);Ntick++) {
       
      iint mrab_order = 0; 
      
      if(tstep==0)
	mrab_order = 0; // first order
      else if(tstep==1)
	mrab_order = 1; // second order
      else
	mrab_order = 2; // third order 

      //synthesize actual stage time
      iint t = tstep*pow(2,mesh->MRABNlevels-1) + Ntick;

      iint lev;
      for (lev=0;lev<mesh->MRABNlevels;lev++)
	if (Ntick % (1<<lev) != 0) break;

      if(mesh->totalHaloPairs>0){
	// extract halo on DEVICE
	iint Nentries = mesh->Np*mesh->Nfields;
	
	mesh->haloExtractKernel(mesh->totalHaloPairs,
				Nentries,
				mesh->o_haloElementList,
				mesh->o_q,
				mesh->o_haloBuffer);

	// copy extracted halo to HOST 
	mesh->o_haloBuffer.copyTo(sendBuffer);      
	
	// start halo exchange
	meshHaloExchangeStart(mesh,
			      mesh->Np*mesh->Nfields*sizeof(dfloat),
			      sendBuffer,
			      recvBuffer);
      }
      
      for (iint l=0;l<lev;l++) {
	if (mesh->MRABNelements[l]) {
	  // compute volume contribution to DG boltzmann RHS
	  mesh->volumeKernel(mesh->MRABNelements[l],
			     mesh->o_MRABelementIds[l],
			     mesh->MRABshiftIndex[l],
			     mesh->Nrhs,
			     mesh->o_vgeo,
			     mesh->o_D,
			     mesh->o_x,
			     mesh->o_y,
			     mesh->o_z,
			     mesh->o_q,
			     mesh->o_rhsq);
	}
      }

      if(mesh->totalHaloPairs>0){
	// wait for halo data to arrive
	meshHaloExchangeFinish(mesh);
	
	// copy halo data to DEVICE
	size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
	mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
      }
    
      mesh->device.finish();
      
      occa::tic("surfaceKernel");

      for (iint l=0;l<lev;l++) {
	if (mesh->MRABNelements[l]) {

	  mesh->surfaceKernel(mesh->MRABNelements[l],
			      mesh->o_MRABelementIds[l],
			      mesh->Nrhs,
			      mesh->MRABshiftIndex[l],
			      mesh->o_sgeo,
			      mesh->o_LIFTT,
			      mesh->o_vmapM,
			      mesh->o_mapP,
			      t,
			      mesh->o_x,
			      mesh->o_y,
			      mesh->o_z,
			      mesh->o_q,
			      mesh->o_fQM,
			      mesh->o_rhsq);
	  mesh->lev_updates[l] = Ntick;
	}
      }
      occa::toc("surfaceKernel");

            for (iint l = 0; l < lev; l++) {
	
       	mesh->o_shift.copyFrom(mesh->MRABshiftIndex);
	mesh->o_lev_updates.copyFrom(mesh->lev_updates);
	
      	mesh->filterKernelH(mesh->MRABNelements[l],
			    mesh->o_MRABelementIds[l],
			    mesh->Nrhs,
			    mesh->o_shift,
			    alpha,
			    mesh->o_dualProjMatrix,
			    mesh->o_cubeFaceNumber,
			    mesh->o_EToE,
			    mesh->o_x,
			    mesh->o_y,
			    mesh->o_z,
			    mesh->o_lev_updates,
			    mesh->o_MRABlevels,
			    l,
			    mesh->o_rhsq,
			    mesh->o_qFilter);
      }
      for (iint l = 0; l < lev; l++) {
	
	mesh->filterKernelV(mesh->MRABNelements[l],
			    mesh->o_MRABelementIds[l],
			    mesh->Nrhs,
			    mesh->o_shift,
			    alpha,
			    mesh->o_dualProjMatrix,
			    mesh->o_cubeFaceNumber,
			    mesh->o_EToE,
			    mesh->o_x,
			    mesh->o_y,
			    mesh->o_z,
			    mesh->o_lev_updates,
			    mesh->o_MRABlevels,
			    l,
			    mesh->o_qFilter,
			    mesh->o_rhsq);	
			    }
      
      for (lev=0;lev<mesh->MRABNlevels;lev++)
        if ((Ntick+1) % (1<<lev) !=0) break; //find the max lev to update
      
      for (iint l = 0; l < lev; l++) {
	const iint id = mrab_order*mesh->MRABNlevels*mesh->Nrhs + l*mesh->Nrhs;

	occa::tic("updateKernel");
	
	if (mesh->MRABNelements[l]) {
	  mesh->updateKernel(mesh->MRABNelements[l],
			     mesh->o_MRABelementIds[l],
			     mesh->MRSAAB_C[l],
			     mesh->MRAB_A[id+0],
			     mesh->MRAB_A[id+1],
			     mesh->MRAB_A[id+2],
			     mesh->MRSAAB_A[id+0],
			     mesh->MRSAAB_A[id+1],
			     mesh->MRSAAB_A[id+2],
			     mesh->MRABshiftIndex[l],
			     mesh->o_rhsq,
			     mesh->o_vmapM,
			     mesh->o_fQM,
			     mesh->o_q);

	  //we *must* use 2 here (n - 1), so rk coefficients point the right direction in time
	  mesh->MRABshiftIndex[l] = (mesh->MRABshiftIndex[l]+2)%mesh->Nrhs;
	}
	}
      
      occa::toc("updateKernel");
      
      if (lev<mesh->MRABNlevels) {
	
	const iint id = mrab_order*mesh->MRABNlevels*mesh->Nrhs + lev*mesh->Nrhs;
	
	if (mesh->MRABNhaloElements[lev]) {
	  //trace update using same kernel
	  mesh->traceUpdateKernel(mesh->MRABNhaloElements[lev],
				  mesh->o_MRABhaloIds[lev],
				  mesh->MRSAAB_C[lev-1], //
				  mesh->MRAB_B[id+0], //
				  mesh->MRAB_B[id+1],
				  mesh->MRAB_B[id+2], //
				  mesh->MRSAAB_B[id+0], //
				  mesh->MRSAAB_B[id+1],
				  mesh->MRSAAB_B[id+2], 
				  mesh->MRABshiftIndex[lev],
				  mesh->o_rhsq,
				  mesh->o_vmapM,
				  mesh->o_fQM,
				  mesh->o_q);
      	}
      }
     }

     // estimate maximum error
     if((tstep%mesh->errorStep)==0){
      //	dfloat t = (tstep+1)*mesh->dt;
      dfloat t = mesh->dt*((tstep+1)*pow(2,mesh->MRABNlevels-1));
	
      printf("tstep = %d, t = %g\n", tstep, t);
      fflush(stdout);
      // copy data back to host
      mesh->o_q.copyTo(mesh->q);
      boltzmannPlotNorms(mesh,"norms",tstep,mesh->q);

      // check for nans
      for(int n=0;n<mesh->Nfields*mesh->Nelements*mesh->Np;++n){
	if(isnan(mesh->q[n])){
	  printf("found nan\n");
	  exit(-1);
	}
      }

      // output field files
      iint fld = 0;
      char fname[BUFSIZ];

      boltzmannPlotVTUQuad3DV2(mesh, "foo", tstep/mesh->errorStep);
    }        
    occa::printTimer();
  }
  
  free(recvBuffer);
  free(sendBuffer);
}


#include "boltzmannQuad3D.h"

void boltzmannRunQuad3D(solver_t *solver){

  //
  mesh_t *mesh = solver->mesh;
  
  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  occa::initTimer(mesh->device);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){
    
    for(iint rk=0;rk<mesh->Nrk;++rk){
      
      // intermediate stage time
      dfloat t = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      printf("tstep = %d, rk = %d, t = %g\n", tstep, rk, t);
      
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


      // compute volume contribution to DG boltzmann RHS
      mesh->volumeKernel(mesh->Nelements,
			 mesh->o_vgeo,
			 mesh->o_D,
			 mesh->o_x,
			 mesh->o_y,
			 mesh->o_z,
			 mesh->o_q,
			 mesh->o_rhsq);
      
      if(mesh->totalHaloPairs>0){
	// wait for halo data to arrive
	meshHaloExchangeFinish(mesh);
	
	// copy halo data to DEVICE
	size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
	mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
      }
      
      mesh->device.finish();
      occa::tic("surfaceKernel");
      
      // compute surface contribution to DG boltzmann RHS
      mesh->surfaceKernel(mesh->Nelements,
			  mesh->o_sgeo,
			  mesh->o_LIFTT,
			  mesh->o_vmapM,
			  mesh->o_vmapP,
			  t,
			  mesh->o_x,
			  mesh->o_y,
			  mesh->o_z,
			  mesh->o_q,
			  mesh->o_rhsq);
      
      mesh->device.finish();
      occa::toc("surfaceKernel");
      
      dfloat tupdate = tstep*mesh->dt + mesh->dt*mesh->rkc[rk+1];

      mesh->device.finish();
      occa::tic("updateKernel");

      mesh->updateKernel(mesh->Nelements,
			 mesh->dt,
			 mesh->rka[rk],
			 mesh->rkb[rk],
			 mesh->o_rhsq,
			 mesh->o_resq,
			 mesh->o_q);

      mesh->device.finish();
      occa::toc("updateKernel");
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0){

      // copy data back to host
      mesh->o_q.copyTo(mesh->q);
      
      // do error stuff on host
      //      boltzmannErrorQuad2D(mesh, mesh->dt*(tstep+1));

      // compute vorticity
      //      boltzmannComputeVorticityQuad2D(mesh, mesh->q, 0, mesh->Nfields);
      
      // output field files
      iint fld = 0;
      char fname[BUFSIZ];
      sprintf(fname, "foo_T%04d.vtu", tstep/mesh->errorStep);
      boltzmannPlotVTUQuad3D(mesh, fname, fld);
    }
  }

  occa::printTimer();
  
  free(recvBuffer);
  free(sendBuffer);
}




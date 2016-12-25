#include <stdio.h>
#include <math.h>
#include "euler2D.h"


void eulerRun2D(mesh2D *mesh){

  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(iint rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat t = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      //      mesh->o_q.copyTo(mesh->q);
      //      meshBoltzmannError2D(mesh, mesh->dt*(tstep+1));
      
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
	meshHaloExchangeStart2D(mesh,
				mesh->Np*mesh->Nfields*sizeof(dfloat),
				sendBuffer,
				recvBuffer);
      }

#if 1
      // compute volume contribution to DG boltzmann RHS
      mesh->volumeKernel(mesh->Nelements,
			 mesh->o_vgeo,
			 mesh->o_cubInterpT,
			 mesh->o_cubDrWT,
			 mesh->o_cubDsWT,
			 mesh->o_q,
			 mesh->o_rhsq);
#endif
      
      if(mesh->totalHaloPairs>0){
	// wait for halo data to arrive
	meshHaloExchangeFinish2D(mesh);
	
	// copy halo data to DEVICE
	size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
	mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
      }

#if 1
      // compute surface contribution to DG boltzmann RHS
      mesh->surfaceKernel(mesh->Nelements,
			  mesh->o_sgeo,
			  mesh->o_intInterpT,
			  mesh->o_intLIFTT,
			  mesh->o_vmapM,
			  mesh->o_vmapP,
			  mesh->o_EToB,
			  t,
			  mesh->o_intx,
			  mesh->o_inty,
			  mesh->o_q,
			  mesh->o_rhsq);
#endif 
      // update solution using Runge-Kutta
      iint recombine = 0; (rk==mesh->Nrk-1); // recombine at end of RK step (q/2=>qx, q/2=>qy)
      mesh->updateKernel(mesh->Nelements*mesh->Np*mesh->Nfields,
			 mesh->dt,
			 mesh->rka[rk],
			 mesh->rkb[rk],
			 mesh->o_rhsq,
			 mesh->o_resq,
			 mesh->o_q);
      
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0){

      // copy data back to host
      mesh->o_q.copyTo(mesh->q);
      
      // do error stuff on host
      eulerError2D(mesh, mesh->dt*(tstep+1));

      // output field files
      iint fld = 0;
      char fname[BUFSIZ];
      sprintf(fname, "foo_T%04d", tstep/mesh->errorStep);
      meshPlotVTU2D(mesh, fname, fld);
    }
  }

  free(recvBuffer);
  free(sendBuffer);
}




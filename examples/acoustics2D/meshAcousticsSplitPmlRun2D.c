#include <stdio.h>
#include <math.h>
#include "mesh2D.h"


void meshAcousticsSplitPmlOccaRun2D(mesh2D *mesh){

  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(iint rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat t = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      
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

      // compute volume contribution to DG acoustics RHS
      mesh->volumeKernel(mesh->Nelements,
			 mesh->o_vgeo,
			 mesh->o_sigmax,
			 mesh->o_sigmay,
			 mesh->o_DrT,
			 mesh->o_DsT,
			 mesh->o_q,
			 mesh->o_pmlqx,
			 mesh->o_pmlqy,
			 mesh->o_rhspmlqx,
			 mesh->o_rhspmlqy);
      
      
      if(mesh->totalHaloPairs>0){
	// wait for halo data to arrive
	meshHaloExchangeFinish2D(mesh);
	
	// copy halo data to DEVICE
	size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
	mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
      }
      
      // compute surface contribution to DG acoustics RHS
      mesh->surfaceKernel(mesh->Nelements,
			  mesh->o_sgeo,
			  mesh->o_LIFTT,
			  mesh->o_vmapM,
			  mesh->o_vmapP,
			  mesh->o_EToB,
			  t,
			  mesh->o_x,
			  mesh->o_y,
			  mesh->o_q,
			  mesh->o_rhspmlqx,
			  mesh->o_rhspmlqy);
      
      
      if(!(tstep%5000)){

	mesh->o_rhspmlqx.copyTo(mesh->rhspmlqx);
	mesh->o_rhspmlqy.copyTo(mesh->rhspmlqy);
	
	
	dfloat maxrhsx = 0, maxrhsy = 0;
	for(int e=0;e<mesh->Nelements;++e){
	  for(int fld=0;fld<3;++fld){
	    for(int n=0;n<mesh->Np;++n){
	      maxrhsx = mymax(maxrhsx, fabs(mesh->rhspmlqx[n+fld*mesh->Np+mesh->Nfields*e*mesh->Np]));
	      maxrhsy = mymax(maxrhsy, fabs(mesh->rhspmlqy[n+fld*mesh->Np+mesh->Nfields*e*mesh->Np]));
	    }
	  }
	}
	printf("maxrhsx*dt=%lg, maxrhsqy*dt=%lg\n", mesh->dt*maxrhsx, mesh->dt*maxrhsy);
      }

      // update solution using Runge-Kutta
      iint recombine = 0; (rk==mesh->Nrk-1); // recombine at end of RK step (q/2=>qx, q/2=>qy)
      mesh->updateKernel(mesh->Nelements*mesh->Np*mesh->Nfields,
			 recombine,
			 mesh->dt,
			 mesh->rka[rk],
			 mesh->rkb[rk],
			 mesh->o_rhspmlqx,
			 mesh->o_rhspmlqy,
			 mesh->o_respmlqx,
			 mesh->o_respmlqy,
			 mesh->o_pmlqx,
			 mesh->o_pmlqy,
			 mesh->o_q);
      
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0){

      // copy data back to host
      mesh->o_q.copyTo(mesh->q);
      
      // do error stuff on host
      meshAcousticsError2D(mesh, mesh->dt*(tstep+1));

      // compute vorticity
      meshAcousticsComputeVorticity2D(mesh, mesh->q, 0, mesh->Nfields);
      
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




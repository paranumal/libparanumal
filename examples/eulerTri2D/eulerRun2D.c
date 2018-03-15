#include <stdio.h>
#include <math.h>
#include "euler2D.h"


void eulerRanges2D(mesh2D *mesh){

  mesh->o_q.copyTo(mesh->q);
  mesh->o_rhsq.copyTo(mesh->rhsq);

  dfloat minR = 1e9, maxR = -1e9;
  dfloat minRhsR = 1e9, maxRhsR = -1e9;

  dfloat minU = 1e9, maxU = -1e9;
  dfloat minRhsU = 1e9, maxRhsU = -1e9;

  dfloat minV = 1e9, maxV = -1e9;
  dfloat minRhsV = 1e9, maxRhsV = -1e9;

  
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      int id =	0+mesh->Nfields*n+mesh->Nfields*mesh->Np*e;
      minR = mymin(minR, mesh->q[id]);
      maxR = mymax(maxR, mesh->q[id]);
      minRhsR = mymin(minRhsR, mesh->rhsq[id]);
      maxRhsR = mymax(maxRhsR, mesh->rhsq[id]);

      ++id;
      minU = mymin(minU, mesh->q[id]);
      maxU = mymax(maxU, mesh->q[id]);
      minRhsU = mymin(minRhsU, mesh->rhsq[id]);
      maxRhsU = mymax(maxRhsU, mesh->rhsq[id]);

      ++id;
      minV = mymin(minV, mesh->q[id]);
      maxV = mymax(maxV, mesh->q[id]);
      minRhsV = mymin(minRhsV, mesh->rhsq[id]);
      maxRhsV = mymax(maxRhsV, mesh->rhsq[id]);
      
    }
  }
  printf("R in [%g,%g]\n", minR, maxR);
  printf("RHS R in [%g,%g]\n", minRhsR, maxRhsR);

  printf("U in [%g,%g]\n", minU, maxU);
  printf("RHS U in [%g,%g]\n", minRhsU, maxRhsU);

  printf("V in [%g,%g]\n", minV, maxV);
  printf("RHS V in [%g,%g]\n", minRhsV, maxRhsV);
  
}


void eulerRun2D(mesh2D *mesh){

  // MPI send buffer
  int haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(int rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat t = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      if(mesh->totalHaloPairs>0){
	// extract halo on DEVICE
	int Nentries = mesh->Np*mesh->Nfields;
	
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
	meshHaloExchangeFinish(mesh);
	
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
      int fld = 0;
      char fname[BUFSIZ];
      sprintf(fname, "foo_T%04d", tstep/mesh->errorStep);
      meshPlotVTU2D(mesh, fname, fld);
    }
  }

  free(recvBuffer);
  free(sendBuffer);
}




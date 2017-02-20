#include "acoustics3D.h"

void acousticsRun3D(mesh3D *mesh){

  // MPI send buffer
  dfloat *sendBuffer =
    (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields,
		     sizeof(dfloat));
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(iint rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat time = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      // halo exchange
      meshHaloExchange(mesh,
			 mesh->Np*mesh->Nfields,
			 mesh->q,
			 sendBuffer,
			 mesh->q+mesh->Np*mesh->Nfields*mesh->Nelements);
      
      // compute volume contribution to DG acoustics RHS
      acousticsVolume3D(mesh);

      // compute surface contribution to DG acoustics RHS
      acousticsSurface3D(mesh, time);

      // update solution using LSERK4
      acousticsUpdate3D(mesh, mesh->rka[rk], mesh->rkb[rk]);
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0)
      acousticsError3D(mesh, mesh->dt*(tstep+1));
  }

  // output field files
  iint fld = 3;
  meshPlotVTU3D(mesh, "foo", fld);

  free(sendBuffer);
}

void acousticsRun3Dbbdg(mesh3D *mesh){

  // MPI send buffer
  dfloat *sendBuffer =
    (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields,
         sizeof(dfloat));
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(iint rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat time = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      // halo exchange
      meshHaloExchange(mesh,
       mesh->Np*mesh->Nfields,
       mesh->q,
       sendBuffer,
       mesh->q+mesh->Np*mesh->Nfields*mesh->Nelements);
      
      // compute volume contribution to DG acoustics RHS
      acousticsVolume3Dbbdg(mesh);

      // compute surface contribution to DG acoustics RHS
      acousticsSurface3Dbbdg(mesh, time);

      // update solution using LSERK4
      acousticsUpdate3D(mesh, mesh->rka[rk], mesh->rkb[rk]);
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0) {
      
      //Transform to nodal basis
      dfloat qtmp[mesh->Nfields*mesh->Np];  
      for (iint e =0;e<mesh->Nelements;e++){
        iint id = e*mesh->Np*mesh->Nfields;
        
        for (iint n=0; n<mesh->Np; n++){
          qtmp[n*mesh->Nfields + 0] = mesh->q[id+n*mesh->Nfields+0];
          qtmp[n*mesh->Nfields + 1] = mesh->q[id+n*mesh->Nfields+1];
          qtmp[n*mesh->Nfields + 2] = mesh->q[id+n*mesh->Nfields+2];
          qtmp[n*mesh->Nfields + 3] = mesh->q[id+n*mesh->Nfields+3];
          mesh->q[id+n*mesh->Nfields+0] = 0.0;
          mesh->q[id+n*mesh->Nfields+1] = 0.0;
          mesh->q[id+n*mesh->Nfields+2] = 0.0;
          mesh->q[id+n*mesh->Nfields+3] = 0.0;
        }
        for (iint n=0;n<mesh->Np;n++){
          for (iint m=0; m<mesh->Np; m++){
            mesh->q[id+n*mesh->Nfields + 0] += mesh->VB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+0];
            mesh->q[id+n*mesh->Nfields + 1] += mesh->VB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+1];
            mesh->q[id+n*mesh->Nfields + 2] += mesh->VB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+2];
            mesh->q[id+n*mesh->Nfields + 3] += mesh->VB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+3];
          }
        }
      }
      
      // do error stuff on host
      acousticsError3D(mesh, mesh->dt*(tstep+1));
      
      // output field files
      iint fld = 3;
      meshPlotVTU3D(mesh, "foo", fld);
      
      //Transform to back to modal basis
      for (iint e =0;e<mesh->Nelements;e++){
        iint id = e*mesh->Np*mesh->Nfields;  

        for (iint n=0; n<mesh->Np; n++){
          qtmp[n*mesh->Nfields + 0] = mesh->q[id+n*mesh->Nfields+0];
          qtmp[n*mesh->Nfields + 1] = mesh->q[id+n*mesh->Nfields+1];
          qtmp[n*mesh->Nfields + 2] = mesh->q[id+n*mesh->Nfields+2];
          qtmp[n*mesh->Nfields + 3] = mesh->q[id+n*mesh->Nfields+3];
          mesh->q[id+n*mesh->Nfields+0] = 0.0;
          mesh->q[id+n*mesh->Nfields+1] = 0.0;
          mesh->q[id+n*mesh->Nfields+2] = 0.0;
          mesh->q[id+n*mesh->Nfields+3] = 0.0;
        }
        for (iint n=0;n<mesh->Np;n++){
          for (iint m=0; m<mesh->Np; m++){
            mesh->q[id+n*mesh->Nfields + 0] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+0];
            mesh->q[id+n*mesh->Nfields + 1] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+1];
            mesh->q[id+n*mesh->Nfields + 2] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+2];
            mesh->q[id+n*mesh->Nfields + 3] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+3];
          }
        }
      }
    }
  }

  free(sendBuffer);
}

void acousticsOccaRun3D(mesh3D *mesh){

  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(iint rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat time = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

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
      
      // compute volume contribution to DG acoustics RHS
      mesh->volumeKernel(mesh->Nelements,
			 mesh->o_vgeo,
			 mesh->o_DrT,
			 mesh->o_DsT,
			 mesh->o_DtT,
			 mesh->o_q,
			 mesh->o_rhsq);
      
      if(mesh->totalHaloPairs>0){
	// wait for halo data to arrive
	meshHaloExchangeFinish(mesh);
	
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
			  time,
			  mesh->o_x,
			  mesh->o_y,
			  mesh->o_z,
			  mesh->o_q,
			  mesh->o_rhsq);
      
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
      acousticsError3D(mesh, mesh->dt*(tstep+1));

      // output field files
      iint fld = 3;
      meshPlotVTU3D(mesh, "foo", fld);
    }
  }

  free(recvBuffer);
  free(sendBuffer);
}

void acousticsOccaRun3Dbbdg(mesh3D *mesh){

  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(iint rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat time = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

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
      
      // compute volume contribution to DG acoustics RHS
      mesh->volumeKernel(mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_D0ids,
       mesh->o_D1ids,
       mesh->o_D2ids,
       mesh->o_D3ids,
       mesh->o_Dvals,
       mesh->o_q,
       mesh->o_rhsq);
      
      if(mesh->totalHaloPairs>0){
  // wait for halo data to arrive
  meshHaloExchangeFinish(mesh);
  
  // copy halo data to DEVICE
  size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
  mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
      }
      
      // compute surface contribution to DG acoustics RHS
      mesh->surfaceKernel(mesh->Nelements,
        mesh->o_sgeo,
        mesh->o_L0ids,
        mesh->o_L0vals,
        mesh->o_ELids,
        mesh->o_ELvals,
        mesh->o_vmapM,
        mesh->o_vmapP,
        mesh->o_EToB,
        time,
        mesh->o_x,
        mesh->o_y,
        mesh->o_z,
        mesh->o_q,
        mesh->o_rhsq);
      
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

      //Transform to nodal basis
      dfloat qtmp[mesh->Nfields*mesh->Np];  
      for (iint e =0;e<mesh->Nelements;e++){
        iint id = e*mesh->Np*mesh->Nfields;
        
        for (iint n=0; n<mesh->Np; n++){
          qtmp[n*mesh->Nfields + 0] = mesh->q[id+n*mesh->Nfields+0];
          qtmp[n*mesh->Nfields + 1] = mesh->q[id+n*mesh->Nfields+1];
          qtmp[n*mesh->Nfields + 2] = mesh->q[id+n*mesh->Nfields+2];
          qtmp[n*mesh->Nfields + 3] = mesh->q[id+n*mesh->Nfields+3];
          mesh->q[id+n*mesh->Nfields+0] = 0.0;
          mesh->q[id+n*mesh->Nfields+1] = 0.0;
          mesh->q[id+n*mesh->Nfields+2] = 0.0;
          mesh->q[id+n*mesh->Nfields+3] = 0.0;
        }
        for (iint n=0;n<mesh->Np;n++){
          for (iint m=0; m<mesh->Np; m++){
            mesh->q[id+n*mesh->Nfields + 0] += mesh->VB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+0];
            mesh->q[id+n*mesh->Nfields + 1] += mesh->VB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+1];
            mesh->q[id+n*mesh->Nfields + 2] += mesh->VB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+2];
            mesh->q[id+n*mesh->Nfields + 3] += mesh->VB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+3];
          }
        }
      }

      // do error stuff on host
      acousticsError3D(mesh, mesh->dt*(tstep+1));
      
      // output field files
      iint fld = 3;
      meshPlotVTU3D(mesh, "foo", fld);
    }
  }

  free(recvBuffer);
  free(sendBuffer);
}

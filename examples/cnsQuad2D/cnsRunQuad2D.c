#include "cnsQuad2D.h"

void cnsRunQuad2D(cns_t *cns){

  mesh_t *mesh = cns->mesh;

  dfloat ramp = 1;
  
  // MPI send buffer
  int haloBytes = mesh->totalHaloPairs*mesh->Np*cns->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  int haloStressesBytes = mesh->totalHaloPairs*mesh->Np*cns->Nstresses*sizeof(dfloat);
  dfloat *sendStressesBuffer = (dfloat*) malloc(haloStressesBytes);
  dfloat *recvStressesBuffer = (dfloat*) malloc(haloStressesBytes);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(int rk=0;rk<mesh->Nrk;++rk){

      dfloat currentTime = tstep*mesh->dt + mesh->rkc[rk]*mesh->dt;
      
      // extract q halo on DEVICE
      if(mesh->totalHaloPairs>0){
	int Nentries = mesh->Np*cns->Nfields;
	
	mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_q, mesh->o_haloBuffer);
	
	// copy extracted halo to HOST 
	mesh->o_haloBuffer.copyTo(sendBuffer);      
	
	// start halo exchange
	meshHaloExchangeStart(mesh, mesh->Np*cns->Nfields*sizeof(dfloat), sendBuffer, recvBuffer);
      }

      // now compute viscous stresses
      cns->stressesVolumeKernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, cns->mu, cns->o_q, cns->o_viscousStresses);

      // wait for q halo data to arrive
      if(mesh->totalHaloPairs>0){
	meshHaloExchangeFinish(mesh);
	
	// copy halo data to DEVICE
	size_t offset = mesh->Np*cns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
	cns->o_q.copyFrom(recvBuffer, haloBytes, offset);
      }
      
      cns->stressesSurfaceKernel(mesh->Nelements, mesh->o_sgeo, cns->o_LIFTT,
				 mesh->o_vmapM, mesh->o_vmapP, mesh->o_EToB, currentTime,
				 mesh->o_x, mesh->o_y, ramp, cns->mu, cns->o_q, cns->o_viscousStresses);
      
      // extract stresses halo on DEVICE
      if(mesh->totalHaloPairs>0){
	int Nentries = mesh->Np*cns->Nstresses;
	
	mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_viscousStresses, cns->o_haloStressesBuffer);
	
	// copy extracted halo to HOST 
	mesh->o_haloBuffer.copyTo(sendBuffer);      
	
	// start halo exchange
	meshHaloExchangeStart(mesh, mesh->Np*cns->Nstresses*sizeof(dfloat), sendStressesBuffer, recvStressesBuffer);
      }
      
      // compute volume contribution to DG cns RHS
      cns->volumeKernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, cns->o_viscousStresses, cns->o_q, cns->o_rhsq);

      // wait for halo stresses data to arrive
      if(mesh->totalHaloPairs>0){
	meshHaloExchangeFinish(mesh);
	
	// copy halo data to DEVICE
	size_t offset = mesh->Np*cns->Nstresses*mesh->Nelements*sizeof(dfloat); // offset for halo data
	cns->o_viscousStresses.copyFrom(recvStressesBuffer, haloStressesBytes, offset);
      }
      
      // compute surface contribution to DG cns RHS (LIFTT ?)
      cns->surfaceKernel(mesh->Nelements, mesh->o_sgeo, cns->o_LIFTT, mesh->o_vmapM, mesh->o_vmapP, mesh->o_EToB,
			 currentTime, mesh->o_x, mesh->o_y, ramp, cns->mu, cns->o_q, cns->o_viscousStresses, cns->o_rhsq);
      
      // update solution using Runge-Kutta
      cns->updateKernel(mesh->Nelements, mesh->dt, mesh->rka[rk], mesh->rkb[rk], cns->o_rhsq, cns->o_resq, cns->o_q);
      
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0){
      
      // copy data back to host
      cns->o_q.copyTo(mesh->q);

      // do error stuff on host
      cnsError2D(mesh, mesh->dt*(tstep+1));

      // output field files
      int fld = 0;
      char fname[BUFSIZ];
      sprintf(fname, "foo_%05d.vtu", tstep/mesh->errorStep);
      cnsPlotVTUQuad2D(cns, fname);
    }
  }

  free(recvBuffer);
  free(sendBuffer);

  free(recvStressesBuffer);
  free(sendStressesBuffer);
}

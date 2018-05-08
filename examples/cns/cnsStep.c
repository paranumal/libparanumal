#include "cns.h"

void cnsDopriStep(cns_t *cns, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = cns->mesh;
  
  //RK step
  for(int rk=0;rk<cns->Nrk;++rk){
    
    // t_rk = t + C_rk*dt
    dfloat currentTime = time + cns->rkC[rk]*mesh->dt;
    
    //compute RK stage 
    // rkq = q + dt sum_{i=0}^{rk-1} a_{rk,i}*rhsq_i
    cns->rkStageKernel(mesh->Nelements,
		       rk,
		       mesh->dt,
		       cns->o_rkA,
		       cns->o_q,
		       cns->o_rkrhsq,
		       cns->o_rkq);
    
    //compute RHS
    // rhsq = F(currentTIme, rkq)

    // extract q halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*cns->Nfields;          
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_rkq, cns->o_haloBuffer);
      
      // copy extracted halo to HOST 
      cns->o_haloBuffer.copyTo(cns->sendBuffer);      
      
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*cns->Nfields*sizeof(dfloat), cns->sendBuffer, cns->recvBuffer);
    }

    // now compute viscous stresses
    cns->stressesVolumeKernel(mesh->Nelements, 
			      mesh->o_vgeo, 
			      mesh->o_Dmatrices,
			      cns->mu, 
			      cns->o_rkq, 
			      cns->o_viscousStresses);

    // wait for q halo data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
          
      // copy halo data to DEVICE
      size_t offset = mesh->Np*cns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      cns->o_rkq.copyFrom(cns->recvBuffer, cns->haloBytes, offset);
    }

    cns->stressesSurfaceKernel(mesh->Nelements, 
			       mesh->o_sgeo, 
			       mesh->o_LIFTT,
			       mesh->o_vmapM, 
			       mesh->o_vmapP, 
			       mesh->o_EToB, 
			       currentTime,
			       mesh->o_x, 
			       mesh->o_y,
			       mesh->o_z, 
			       cns->mu, 
			       cns->o_rkq, 
			       cns->o_viscousStresses);

    // extract stresses halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*cns->Nstresses;
      
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_viscousStresses, cns->o_haloStressesBuffer);
      
      // copy extracted halo to HOST 
      cns->o_haloStressesBuffer.copyTo(cns->sendStressesBuffer);      
      
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*cns->Nstresses*sizeof(dfloat), cns->sendStressesBuffer, cns->recvStressesBuffer);
    }

    // compute volume contribution to DG cns RHS
    if (newOptions.compareArgs("ADVECTION TYPE","CUBATURE")) {
      cns->cubatureVolumeKernel(mesh->Nelements, 
				cns->advSwitch, 
				mesh->o_vgeo,
				mesh->o_cubvgeo, 
				mesh->o_cubDWmatrices,
				mesh->o_cubInterpT,
				mesh->o_cubProjectT,
				cns->o_viscousStresses, 
				cns->o_rkq, 
				cns->o_rhsq);
    } else {
      cns->volumeKernel(mesh->Nelements, 
			cns->advSwitch, 
			mesh->o_vgeo, 
			mesh->o_Dmatrices,
			cns->o_viscousStresses, 
			cns->o_rkq, 
			cns->o_rhsq);
    }

    // wait for halo stresses data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
      
      // copy halo data to DEVICE
      size_t offset = mesh->Np*cns->Nstresses*mesh->Nelements*sizeof(dfloat); // offset for halo data
      cns->o_viscousStresses.copyFrom(cns->recvStressesBuffer, cns->haloStressesBytes, offset);
    }

    // compute surface contribution to DG cns RHS (LIFTT ?)
    if (newOptions.compareArgs("ADVECTION TYPE","CUBATURE")) {
      cns->cubatureSurfaceKernel(mesh->Nelements, 
				 cns->advSwitch,
				 mesh->o_vgeo, 
				 mesh->o_cubsgeo, 
				 mesh->o_vmapM, 
				 mesh->o_vmapP, 
				 mesh->o_EToB,
				 mesh->o_intInterpT,
				 mesh->o_intLIFTT, 
				 currentTime, 
				 mesh->o_intx, 
				 mesh->o_inty,
				 mesh->o_intz, 
				 cns->mu, 
				 cns->o_rkq, 
				 cns->o_viscousStresses, 
				 cns->o_rhsq);
    } else {
      cns->surfaceKernel(mesh->Nelements, 
			 cns->advSwitch, 
			 mesh->o_sgeo, 
			 mesh->o_LIFTT, 
			 mesh->o_vmapM, 
			 mesh->o_vmapP, 
			 mesh->o_EToB,
			 currentTime, 
			 mesh->o_x, 
			 mesh->o_y,
			 mesh->o_z, 
			 cns->mu, 
			 cns->o_rkq, 
			 cns->o_viscousStresses, 
			 cns->o_rhsq);
    }
    
    // update solution using Runge-Kutta
    // rkrhsq_rk = rhsq
    // if rk==6 
    //   q = q + dt*sum_{i=0}^{rk} rkA_{rk,i}*rkrhs_i
    //   rkerr = dt*sum_{i=0}^{rk} rkE_{rk,i}*rkrhs_i
    cns->rkUpdateKernel(mesh->Nelements, 
			rk,
			mesh->dt, 
			cns->o_rkA, 
			cns->o_rkE, 
			cns->o_q,
			cns->o_rhsq, 
			cns->o_rkrhsq, 
			cns->o_rkq,
			cns->o_rkerr);
  }
}


void cnsLserkStep(cns_t *cns, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = cns->mesh;
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  int advSwitch = 1;//(tstep>100);
    
  for(int rk=0;rk<mesh->Nrk;++rk){
      
    dfloat currentTime = time + mesh->rkc[rk]*mesh->dt;
      
    // extract q halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*cns->Nfields;
        
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_q, cns->o_haloBuffer);
        
      // copy extracted halo to HOST 
      cns->o_haloBuffer.copyTo(cns->sendBuffer);      
        
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*cns->Nfields*sizeof(dfloat), cns->sendBuffer, cns->recvBuffer);
    }
      
    // now compute viscous stresses
    cns->stressesVolumeKernel(mesh->Nelements, 
			      mesh->o_vgeo, 
			      mesh->o_Dmatrices, 
			      cns->mu, 
			      cns->o_q, 
			      cns->o_viscousStresses);
      
    // wait for q halo data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
        
      // copy halo data to DEVICE
      size_t offset = mesh->Np*cns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      cns->o_q.copyFrom(cns->recvBuffer, cns->haloBytes, offset);
    }
      
    cns->stressesSurfaceKernel(mesh->Nelements, 
			       mesh->o_sgeo, 
			       mesh->o_LIFTT,
			       mesh->o_vmapM, 
			       mesh->o_vmapP, 
			       mesh->o_EToB, 
			       currentTime,
			       mesh->o_x, 
			       mesh->o_y,
			       mesh->o_z, 
			       cns->mu, 
			       cns->o_q, 
			       cns->o_viscousStresses);
      
    // extract stresses halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*cns->Nstresses;
          
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_viscousStresses, cns->o_haloStressesBuffer);
        
      // copy extracted halo to HOST 
      cns->o_haloStressesBuffer.copyTo(cns->sendStressesBuffer);      
          
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*cns->Nstresses*sizeof(dfloat), cns->sendStressesBuffer, cns->recvStressesBuffer);
    }
      
    // compute volume contribution to DG cns RHS
    if (newOptions.compareArgs("ADVECTION TYPE","CUBATURE")) {

      cns->cubatureVolumeKernel(mesh->Nelements, 
				advSwitch, 
				mesh->o_vgeo,
				mesh->o_cubvgeo, 
				mesh->o_cubDWmatrices,
				mesh->o_cubInterpT,
				mesh->o_cubProjectT,
				cns->o_viscousStresses, 
				cns->o_q, 
				cns->o_rhsq);
    } else {
      cns->volumeKernel(mesh->Nelements, 
			advSwitch, 
			mesh->o_vgeo, 
			mesh->o_Dmatrices,
			cns->o_viscousStresses, 
			cns->o_q, 
			cns->o_rhsq);
    }

    // wait for halo stresses data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
        
      // copy halo data to DEVICE
      size_t offset = mesh->Np*cns->Nstresses*mesh->Nelements*sizeof(dfloat); // offset for halo data
      cns->o_viscousStresses.copyFrom(cns->recvStressesBuffer, cns->haloStressesBytes, offset);
    }
      
    // compute surface contribution to DG cns RHS (LIFTT ?)
    if (newOptions.compareArgs("ADVECTION TYPE","CUBATURE")) {
      cns->cubatureSurfaceKernel(mesh->Nelements, 
				 advSwitch,
				 mesh->o_vgeo, 
				 mesh->o_cubsgeo, 
				 mesh->o_intInterpT,
				 mesh->o_intLIFTT, 
				 mesh->o_vmapM, 
				 mesh->o_vmapP, 
				 mesh->o_EToB,
				 currentTime, 
				 mesh->o_intx, 
				 mesh->o_inty,
				 mesh->o_intz, 
				 cns->mu, 
				 cns->o_q, 
				 cns->o_viscousStresses, 
				 cns->o_rhsq);
    } else {
      cns->surfaceKernel(mesh->Nelements, 
			 advSwitch, 
			 mesh->o_sgeo, 
			 mesh->o_LIFTT, 
			 mesh->o_vmapM, 
			 mesh->o_vmapP, 
			 mesh->o_EToB,
			 currentTime, 
			 mesh->o_x, 
			 mesh->o_y,
			 mesh->o_z, 
			 cns->mu, 
			 cns->o_q, 
			 cns->o_viscousStresses, 
			 cns->o_rhsq);
    }
        
    // update solution using Runge-Kutta
    cns->updateKernel(mesh->Nelements, 
		      mesh->dt, 
		      mesh->rka[rk], 
		      mesh->rkb[rk], 
		      cns->o_rhsq, 
		      cns->o_resq, 
		      cns->o_q);
  }
}

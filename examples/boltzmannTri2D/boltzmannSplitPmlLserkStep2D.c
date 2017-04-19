#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannSplitPmlLserkStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char * options){

  // LSERK4 stages
  for(iint rk=0;rk<mesh->Nrk;++rk){

    // intermediate stage time
    dfloat t = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];
    
    if(mesh->totalHaloPairs>0){
      // EXCTRACT HALO  on DEVICE
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
    

    // COMPUTE RAMP FUNCTION 
    dfloat ramp, drampdt;
    boltzmannRampFunction2D(t, &ramp, &drampdt);
    //ramp = 1.0;  drampdt = 0.0;
    

    // VOLUME KERNELS
    mesh->device.finish();
    occa::tic("volumeKernel");
    
    // compute volume contribution to DG boltzmann RHS
    if(mesh->pmlNelements){	
      mesh->pmlVolumeKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    ramp, 
			    drampdt,
			    mesh->o_vgeo,
			    mesh->o_sigmax,
			    mesh->o_sigmay,
			    mesh->o_DrT,
			    mesh->o_DsT,
			    mesh->o_q,
			    mesh->o_pmlqx,
			    mesh->o_pmlqy,
			    mesh->o_rhsq,
			    mesh->o_rhspmlqx,
			    mesh->o_rhspmlqy);
    }

    // compute volume contribution to DG boltzmann RHS added d/dt (ramp(qbar)) to RHS
    if(mesh->nonPmlNelements){
      mesh->volumeKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 ramp, 
			 drampdt,
			 mesh->o_vgeo,
			 mesh->o_DrT,
			 mesh->o_DsT,
			 mesh->o_q,
			 mesh->o_rhsq);
	}
    
    
    mesh->device.finish();
    occa::toc("volumeKernel");

      

	if(strstr(options, "CUBATURE")){ 
	// VOLUME KERNELS
    mesh->device.finish();
    occa::tic("relaxationKernel");
		// compute relaxation terms using cubature integration
		if(mesh->pmlNelements){
		  mesh->pmlRelaxationKernel(mesh->pmlNelements,
				     mesh->o_pmlElementIds,
				     mesh->o_cubInterpT,
				     mesh->o_cubProjectT,
				     mesh->o_q,
				     mesh->o_rhsq,
				     mesh->o_rhspmlqx,
				     mesh->o_rhspmlqy);
		}

		// compute relaxation terms using cubature
		if(mesh->nonPmlNelements){
		  mesh->relaxationKernel(mesh->nonPmlNelements,
				     mesh->o_nonPmlElementIds,
				     mesh->o_cubInterpT,
				     mesh->o_cubProjectT,
				     mesh->o_q,
				     mesh->o_rhsq);   
		}
		 // VOLUME KERNELS
    mesh->device.finish();
    occa::toc("relaxationKernel");
	}

	 


    // COMPLETE HALO EXCHANGE
    if(mesh->totalHaloPairs>0){
      // wait for halo data to arrive
      meshHaloExchangeFinish(mesh);
      
      // copy halo data to DEVICE
      size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
    }


    // SURFACE KERNELS
    mesh->device.finish();
    occa::tic("surfaceKernel");

     if(mesh->pmlNelements){
    // compute surface contribution to DG boltzmann RHS
    mesh->pmlSurfaceKernel(mesh->pmlNelements,
			   mesh->o_pmlElementIds,
			   mesh->o_sgeo,
			   mesh->o_LIFTT,
			   mesh->o_vmapM,
			   mesh->o_vmapP,
			   mesh->o_EToB,
			   t,
			   mesh->o_x,
			   mesh->o_y,
			   ramp,
			   mesh->o_q,
			   mesh->o_rhsq,
			   mesh->o_rhspmlqx,
			   mesh->o_rhspmlqy);
    }
    
    if(mesh->nonPmlNelements){
      mesh->surfaceKernel(mesh->nonPmlNelements,
			  mesh->o_nonPmlElementIds,
			  mesh->o_sgeo,
			  mesh->o_LIFTT,
			  mesh->o_vmapM,
			  mesh->o_vmapP,
			  mesh->o_EToB,
			  t,
			  mesh->o_x,
			  mesh->o_y,
			  ramp,
			  mesh->o_q,
			  mesh->o_rhsq);
    }
    
    mesh->device.finish();
    occa::toc("surfaceKernel");
    
    // ramp function for flow at next RK stage
    dfloat tupdate = tstep*mesh->dt + mesh->dt*mesh->rkc[rk+1];
    dfloat rampUpdate, drampdtUpdate;
    boltzmannRampFunction2D(tupdate, &rampUpdate, &drampdtUpdate);
    //rampUpdate = 1.0f; drampdtUpdate = 0.0f;


    //UPDATE
    mesh->device.finish();
    occa::tic("updateKernel");
    
    //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    if (mesh->pmlNelements){   
      mesh->pmlUpdateKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    mesh->dt,
			    mesh->rka[rk],
			    mesh->rkb[rk],
			    rampUpdate,
			    mesh->o_rhsq,
			    mesh->o_rhspmlqx,
			    mesh->o_rhspmlqy,
			    mesh->o_resq,
			    mesh->o_respmlqx,
			    mesh->o_respmlqy,
			    mesh->o_pmlqx,
			    mesh->o_pmlqy,
			    mesh->o_q);
   }
    
    if(mesh->nonPmlNelements){
      mesh->updateKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->rka[rk],
			 mesh->rkb[rk],
			 mesh->o_rhsq,
			 mesh->o_resq,
			 mesh->o_q);
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");      
    
  }
}

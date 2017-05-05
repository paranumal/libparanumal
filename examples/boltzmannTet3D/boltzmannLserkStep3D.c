#include "boltzmann3D.h"

// complete a time step using LSERK4
void boltzmannLserkStep3D(mesh3D *mesh, iint tstep, iint haloBytes,
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
    boltzmannRampFunction3D(t, &ramp, &drampdt);
    
    // VOLUME KERNELS
    mesh->device.finish();
    occa::tic("volumeKernel");
    
    // compute volume contribution to DG boltzmann RHS
    if(mesh->pmlNelements){	
    	mesh->device.finish();
       occa::tic("PML_volumeKernel");

      mesh->pmlVolumeKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    ramp, 
			    drampdt,
			    mesh->o_vgeo,
			    mesh->o_sigmax,
			    mesh->o_sigmay,
			    mesh->o_sigmaz,
			    mesh->o_DrT,
			    mesh->o_DsT,
			    mesh->o_DtT,
			    mesh->o_q,
                mesh->o_pmlq,
			    mesh->o_rhsq,
			    mesh->o_rhspmlq);
      //
       mesh->device.finish();
      occa::toc("PML_volumeKernel");	

    }

    // compute volume contribution to DG boltzmann RHS added d/dt (ramp(qbar)) to RHS
    if(mesh->nonPmlNelements){

      mesh->device.finish();
      occa::tic("NONPML_volumeKernel");	

      mesh->volumeKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 ramp, 
			 drampdt,
			 mesh->o_vgeo,
			 mesh->o_DrT,
			 mesh->o_DsT,
			 mesh->o_DtT,
			 mesh->o_q,
			 mesh->o_rhsq);

      mesh->device.finish();
      occa::toc("NONPML_volumeKernel");

	}
    mesh->device.finish();
    occa::toc("volumeKernel");

      

	if(strstr(options, "CUBATURE")){ 
	// VOLUME KERNELS
    mesh->device.finish();
    occa::tic("relaxationKernel");
		// // compute relaxation terms using cubature integration
		// if(mesh->pmlNelements){
		// 	mesh->device.finish();
  //          occa::tic("PML_relaxationKernel");

		//   mesh->pmlRelaxationKernel(mesh->pmlNelements,
		// 		     mesh->o_pmlElementIds,
		// 		     mesh->o_cubInterpT,
		// 		     mesh->o_cubProjectT,
		// 		     mesh->o_q,
		// 		     mesh->o_rhsq,
		// 		     mesh->o_rhspmlqx,
		// 		     mesh->o_rhspmlqy);

		//    mesh->device.finish();
  //          occa::toc("PML_relaxationKernel");
		// }

		// compute relaxation terms using cubature
		if(mesh->nonPmlNelements){
		  mesh->device.finish();
           occa::tic("NONPML_relaxationKernel");
           	
		  mesh->relaxationKernel(mesh->nonPmlNelements,
				     mesh->o_nonPmlElementIds,
				     mesh->o_cubInterpT,
				     mesh->o_cubProjectT,
				     mesh->o_q,
				     mesh->o_rhsq);  

	      mesh->device.finish();
           occa::toc("NONPML_relaxationKernel");			      
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
     
    mesh->device.finish();
    occa::tic("PML_surfaceKernel"); 

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
			   mesh->o_z,
			   ramp,
			   mesh->o_q,
			   mesh->o_rhsq,
			   mesh->o_rhspmlq);

    mesh->device.finish();
    occa::toc("PML_surfaceKernel"); 
    }
    
    if(mesh->nonPmlNelements){
    	  mesh->device.finish();
    occa::tic("NONPML_surfaceKernel"); 

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
			  mesh->o_z,
			  ramp,
			  mesh->o_q,
			  mesh->o_rhsq);

    mesh->device.finish();
    occa::toc("NONPML_surfaceKernel"); 
    }
    
    mesh->device.finish();
    occa::toc("surfaceKernel");
    
    // ramp function for flow at next RK stage
    dfloat tupdate = tstep*mesh->dt + mesh->dt*mesh->rkc[rk+1];
    dfloat rampUpdate, drampdtUpdate;
    boltzmannRampFunction3D(tupdate, &rampUpdate, &drampdtUpdate);


    //UPDATE
    mesh->device.finish();
    occa::tic("updateKernel");
    
    //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    if (mesh->pmlNelements){   
    mesh->device.finish();
      occa::tic("PML_updateKernel"); 
      mesh->pmlUpdateKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    mesh->dt,
			    mesh->rka[rk],
			    mesh->rkb[rk],
			    rampUpdate,
			    mesh->o_rhsq,
			    mesh->o_rhspmlq,
			    mesh->o_resq,
			    mesh->o_respmlq,
			    mesh->o_q,
			    mesh->o_pmlq);
       mesh->device.finish();
      occa::toc("PML_updateKernel"); 
   }
    
    if(mesh->nonPmlNelements){
    	 mesh->device.finish();
      occa::tic("NONPML_updateKernel");   
      mesh->updateKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->rka[rk],
			 mesh->rkb[rk],
			 mesh->o_rhsq,
			 mesh->o_resq,
			 mesh->o_q);
       mesh->device.finish();
      occa::toc("NONPML_updateKernel");   
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");      
    
  }
}

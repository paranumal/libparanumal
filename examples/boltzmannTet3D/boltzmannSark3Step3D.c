#include "boltzmann3D.h"

// complete a time step using LSERK4
void boltzmannSark3Step3D(mesh3D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer,char * options){


    // STAGE1
    dfloat t = tstep*mesh->dt + mesh->dt*mesh->rk3c[0];

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

    dfloat ramp, drampdt;
    boltzmannRampFunction3D(t, &ramp, &drampdt);

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

    // compute volume contribution to DG boltzmann RHS
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
	    
	    if(mesh->pmlNelements){

	       mesh->device.finish();
           occa::tic("PML_relaxationKernel");
		   mesh->pmlRelaxationKernel(mesh->pmlNelements,
				     mesh->o_pmlElementIds,
				     mesh->o_cubInterpT,
				     mesh->o_cubProjectT,
				     mesh->o_q,
				     mesh->o_rhsq); 
		   mesh->device.finish();
           occa::toc("PML_relaxationKernel");
		}
	  
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




    // complete halo exchange
    if(mesh->totalHaloPairs>0){
      // wait for halo data to arrive
      meshHaloExchangeFinish(mesh);
      
      // copy halo data to DEVICE
      size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
    }
    
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


    mesh->device.finish();
    occa::tic("updateKernel");
    
   //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    if (mesh->pmlNelements){   
      mesh->device.finish();
      occa::tic("PML_updateKernel");  
      mesh->pmlUpdateStageKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    mesh->dt,
			    mesh->sarke[0],
			    mesh->rk3a[1][0], //a21
			    mesh->sarka[1][0], //a21
			    mesh->rk3a[0][0], // 0.
			    mesh->sarka[0][0], // 0. 
			    ramp,
			    mesh->o_rhsq,
			    mesh->o_rhspmlq,
			    mesh->o_rhsq2,
			    mesh->o_rhspmlq2,
			    mesh->o_qold,
			    mesh->o_pmlqold,
			    mesh->o_q,
			    mesh->o_pmlq);
        mesh->device.finish();
      occa::toc("PML_updateKernel"); 
    }
    
    if(mesh->nonPmlNelements){

    	 mesh->device.finish();
    	  occa::tic("NONPML_updateKernel");  
      mesh->updateStageKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->sarke[0],
			 mesh->rk3a[1][0], //a21
			 mesh->sarka[1][0], //a21
			 mesh->rk3a[0][0], // 0.
			 mesh->sarka[0][0], // 0. 
			 mesh->o_rhsq,
			 mesh->o_rhsq2,
			 mesh->o_qold,
			 mesh->o_q);
       mesh->device.finish();
      occa::toc("NONPML_updateKernel");  
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");      
   


   // Stage 2

      t = tstep*mesh->dt + mesh->dt*mesh->rk3c[1];

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




    boltzmannRampFunction3D(t, &ramp, &drampdt);

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
			    mesh->o_rhsq2,
			    mesh->o_rhspmlq2);
       	mesh->device.finish();
       occa::toc("PML_volumeKernel");
    } 

    // compute volume contribution to DG boltzmann RHS
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
			 mesh->o_rhsq2);
      mesh->device.finish();
      occa::toc("NONPML_volumeKernel");
	}
    
    
    mesh->device.finish();
    occa::toc("volumeKernel");


   if(strstr(options, "CUBATURE")){ 
	// VOLUME KERNELS
    mesh->device.finish();
    occa::tic("relaxationKernel");
	    
	    if(mesh->pmlNelements){

	       mesh->device.finish();
           occa::tic("PML_relaxationKernel");
		  mesh->pmlRelaxationKernel(mesh->pmlNelements,
				     mesh->o_pmlElementIds,
				     mesh->o_cubInterpT,
				     mesh->o_cubProjectT,
				     mesh->o_q,
				     mesh->o_rhsq2); 
		  	mesh->device.finish();
           occa::toc("PML_relaxationKernel");
		}
	  
	    // compute relaxation terms using cubature
	    if(mesh->nonPmlNelements){
	    	mesh->device.finish();
           occa::tic("NONPML_relaxationKernel");
	      mesh->relaxationKernel(mesh->nonPmlNelements,
				     mesh->o_nonPmlElementIds,
				     mesh->o_cubInterpT,
				     mesh->o_cubProjectT,
				     mesh->o_q,
				     mesh->o_rhsq2); 
		    mesh->device.finish();
           occa::toc("NONPML_relaxationKernel");  
	    }

	  // VOLUME KERNELS
    mesh->device.finish();
    occa::toc("relaxationKernel");
	}



    // complete halo exchange
    if(mesh->totalHaloPairs>0){
      // wait for halo data to arrive
      meshHaloExchangeFinish(mesh);
      
      // copy halo data to DEVICE
      size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
    }
    
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
			   mesh->o_rhsq2,
			   mesh->o_rhspmlq2);
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
			  mesh->o_rhsq2);

    	 mesh->device.finish();
    occa::toc("NONPML_surfaceKernel"); 
    }

    mesh->device.finish();
    occa::toc("surfaceKernel");

     mesh->device.finish();
    occa::tic("updateKernel");
    
   //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    if (mesh->pmlNelements){ 
    mesh->device.finish();
      occa::tic("PML_updateKernel");    
      mesh->pmlUpdateStageKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    mesh->dt,
			    mesh->sarke[1],
				mesh->rk3a[2][0], // a31
			    mesh->sarka[2][0], // a31
			    mesh->rk3a[2][1], // a32
			    mesh->sarka[2][1], // a32
			    ramp,
			    mesh->o_rhsq,
			    mesh->o_rhspmlq,
			    mesh->o_rhsq2,
			    mesh->o_rhspmlq2,
			    mesh->o_qold,
			    mesh->o_pmlqold,
			    mesh->o_q,
			    mesh->o_pmlq);
      mesh->device.finish();
      occa::toc("PML_updateKernel");  
    }
    
    if(mesh->nonPmlNelements){
    	 mesh->device.finish();
      occa::tic("NONPML_updateKernel");   
      mesh->updateStageKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->sarke[1],
			 mesh->rk3a[2][0], // a31
			 mesh->sarka[2][0], // a31
			 mesh->rk3a[2][1], // a32
			 mesh->sarka[2][1], // a32 
			 mesh->o_rhsq,
			 mesh->o_rhsq2,
			 mesh->o_qold,
			 mesh->o_q);
      mesh->device.finish();
      occa::toc("NONPML_updateKernel");   
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");      


    


    // Stage 3

     t = tstep*mesh->dt + mesh->dt*mesh->rk3c[2];

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

    boltzmannRampFunction3D(t, &ramp, &drampdt);

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
			    mesh->o_rhsq3,
			    mesh->o_rhspmlq3);
      //
       mesh->device.finish();
      occa::toc("PML_volumeKernel");	
    } 

    // compute volume contribution to DG boltzmann RHS
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
			 mesh->o_rhsq3);
      mesh->device.finish();
      occa::toc("NONPML_volumeKernel");
	}
    
    mesh->device.finish();
    occa::toc("volumeKernel");




     if(strstr(options, "CUBATURE")){ 
	// VOLUME KERNELS
    mesh->device.finish();
    occa::tic("relaxationKernel");
	    
	    if(mesh->pmlNelements){
	 mesh->device.finish();
           occa::tic("PML_relaxationKernel");

		  mesh->pmlRelaxationKernel(mesh->pmlNelements,
				     mesh->o_pmlElementIds,
				     mesh->o_cubInterpT,
				     mesh->o_cubProjectT,
				     mesh->o_q,
				     mesh->o_rhsq3); 

		   mesh->device.finish();
           occa::toc("PML_relaxationKernel");
		}
	  
	    // compute relaxation terms using cubature
	    if(mesh->nonPmlNelements){
	    	mesh->device.finish();
           occa::tic("NONPML_relaxationKernel");
	      mesh->relaxationKernel(mesh->nonPmlNelements,
				     mesh->o_nonPmlElementIds,
				     mesh->o_cubInterpT,
				     mesh->o_cubProjectT,
				     mesh->o_q,
				     mesh->o_rhsq3);  

		   mesh->device.finish();
           occa::toc("NONPML_relaxationKernel"); 
	    }

	  // VOLUME KERNELS
    mesh->device.finish();
    occa::toc("relaxationKernel");
	}

    // complete halo exchange
    if(mesh->totalHaloPairs>0){
      // wait for halo data to arrive
      meshHaloExchangeFinish(mesh);
      
      // copy halo data to DEVICE
      size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
    }
    
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
			   mesh->o_rhsq3,
			   mesh->o_rhspmlq3);

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
			  mesh->o_rhsq3);
        mesh->device.finish();
    occa::toc("NONPML_surfaceKernel"); 
    }
    
    mesh->device.finish();
    occa::toc("surfaceKernel");


    mesh->device.finish();
    occa::tic("updateKernel");
    
    if (mesh->pmlNelements){   
    	mesh->device.finish();
      occa::tic("PML_updateKernel");  
      mesh->pmlUpdateKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    mesh->dt,
			    mesh->sarke[2],
			 	mesh->rk3b[0], // b1
			 	mesh->sarkb[0], // b1
			 	mesh->rk3b[1], // b2
				mesh->sarkb[1], //b2
				mesh->rk3b[2], // b3
				mesh->sarkb[2], // b3
			    ramp,
			    mesh->o_rhsq,
			    mesh->o_rhspmlq,
			    mesh->o_rhsq2,
			    mesh->o_rhspmlq2,
			    mesh->o_rhsq3,
			    mesh->o_rhspmlq3,
			    mesh->o_qold,
			    mesh->o_pmlqold,
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
			 mesh->sarke[2],
			 mesh->rk3b[0], // b1
			 mesh->sarkb[0], // b1
			 mesh->rk3b[1], // b2
			 mesh->sarkb[1], //b2
			 mesh->rk3b[2], // b3
			 mesh->sarkb[2], // b3
			 mesh->o_rhsq,
			 mesh->o_rhsq2,
			 mesh->o_rhsq3,
			 mesh->o_qold,
			 mesh->o_q);
       mesh->device.finish();
      occa::toc("NONPML_updateKernel");   
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");      

  
   
}

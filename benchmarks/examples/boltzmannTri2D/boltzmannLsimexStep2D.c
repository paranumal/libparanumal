#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannLsimexStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				   dfloat * sendBuffer, dfloat *recvBuffer, char * options){


  for(iint k=0;k<mesh->Nimex;++k)
    {

      // intermediate stage time
      dfloat t = tstep*mesh->dt + mesh->dt*mesh->LsimexC[k];
    
      dfloat ramp, drampdt;
      boltzmannRampFunction2D(t, &ramp, &drampdt);

      // ramp = 1.0; 
      // drampdt = 0.0; 

      // RESIDUAL UPDATE, i.e. Y = Q+ (a(k,k-1)-b(k-1))_ex *Y + (a(k,k-1)-b(k-1))_im *Z
      mesh->device.finish();
      occa::tic("residualUpdateKernel");
      // compute volume contribution to DG boltzmann RHS
      if(mesh->pmlNelements){

      	mesh->device.finish();
        occa::tic("PML_residualUpdateKernel");
		//	printf("pmlNel = %d\n", mesh->pmlNelements);
		mesh->pmlResidualUpdateKernel(mesh->pmlNelements,
					      mesh->o_pmlElementIds,
					      mesh->dt,
					      ramp,
						  mesh->LsimexABi[k],
						  mesh->LsimexABe[k],
						  mesh->o_q,
						  mesh->o_pmlqx,
						  mesh->o_pmlqy,
						  mesh->o_qZ,
						  mesh->o_qY,
                          mesh->o_qYx,
                          mesh->o_qYy);
		mesh->device.finish();
        occa::toc("PML_residualUpdateKernel");
      }
      if(mesh->nonPmlNelements){
      	mesh->device.finish();
        occa::tic("NONPML_residualUpdateKernel");
		//	printf("nonPmlNel = %d\n", mesh->nonPmlNelements);
		mesh->residualUpdateKernel(mesh->nonPmlNelements,
					   mesh->o_nonPmlElementIds,
					   mesh->dt,
					   mesh->LsimexABi[k],
					   mesh->LsimexABe[k],
					   mesh->o_q,
					   mesh->o_qZ,
					   mesh->o_qY);
		mesh->device.finish();
        occa::toc("NONPML_residualUpdateKernel");
      }
      mesh->device.finish();
      occa::toc("residualUpdateKernel");
    

      

      mesh->device.finish();
      occa::tic("PML_Implicit");
      // Compute Implicit Part of Boltzmann, node based no communication
      if(mesh->pmlNelements){
      	//Implicit Solve Satge
      mesh->device.finish();
      occa::tic("PML_ImplicitSolve");
	  mesh->pmlImplicitSolveKernel(mesh->pmlNelements,
				     mesh->o_pmlElementIds,
				     mesh->dt,
					 mesh->LsimexAd[k],
					 mesh->o_cubInterpT,
					 mesh->o_cubProjectT,
					 mesh->o_qY,
					 mesh->o_qZ); 
	  mesh->device.finish();
      occa::toc("PML_ImplicitSolve");
     
      //
      mesh->device.finish();
      occa::tic("PML_ImplicitUpdate");
	// No surface term for implicit part
	mesh->pmlImplicitUpdateKernel(mesh->pmlNelements,
				      mesh->o_pmlElementIds,
				      mesh->dt,
				      ramp,
				      mesh->LsimexAd[k],
				      mesh->o_qY,
				      mesh->o_qYx,
				      mesh->o_qYy,
				      mesh->o_qZ,
				      mesh->o_pmlqx,
				      mesh->o_pmlqy,
				      mesh->o_qSx,
				      mesh->o_qSy,
				      mesh->o_qS,
				      mesh->o_q);
	  mesh->device.finish();
      occa::toc("PML_ImplicitUpdate");
      }
    
      mesh->device.finish();
      occa::toc("PML_Implicit");

      // compute volume contribution to DG boltzmann RHS
      
      mesh->device.finish();
      occa::tic("NONPML_Implicit");
      if(mesh->nonPmlNelements){
      mesh->device.finish();
      occa::tic("NONPML_ImplicitSolve");

         mesh->implicitSolveKernel(mesh->nonPmlNelements,
				  mesh->o_nonPmlElementIds,
				  mesh->dt,
				  mesh->LsimexAd[k],
				  mesh->o_cubInterpT,
				  mesh->o_cubProjectT,
				  mesh->o_qY,
				  mesh->o_qZ); 
       mesh->device.finish();
      occa::toc("NONPML_ImplicitSolve");

      mesh->device.finish();
      occa::tic("NONPML_ImplicitUpdate");
       
	//No surface term for implicit part
	mesh->implicitUpdateKernel(mesh->nonPmlNelements,
				   mesh->o_nonPmlElementIds,
				   mesh->dt,
				   mesh->LsimexAd[k],
				   mesh->o_qZ,
				   mesh->o_qY,
				   mesh->o_q,
				   mesh->o_qS);	
   	  mesh->device.finish();
      occa::toc("NONPML_ImplicitUpdate");
      }
      
      mesh->device.finish();
      occa::toc("NONPML_Implicit");


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
			    mesh->o_DrT,
			    mesh->o_DsT,
			    mesh->o_q,
			    mesh->o_pmlqx,
			    mesh->o_pmlqy,
			    mesh->o_qY,
			    mesh->o_qYx,
			    mesh->o_qYy);
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
			 mesh->o_q,
			 mesh->o_qY);

      mesh->device.finish();
      occa::toc("NONPML_volumeKernel");
	}
    
    
    mesh->device.finish();
    occa::toc("volumeKernel");



      // complete halo exchange
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

    // SURFACE KERNELS
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
			   ramp,
			   mesh->o_q,
			   mesh->o_qY,
			   mesh->o_qYx,
			   mesh->o_qYy);
    
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
			  ramp,
			  mesh->o_q,
			  mesh->o_qY);
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

	mesh->pmlUpdateKernel(mesh->pmlNelements,
			      mesh->o_pmlElementIds,
			      mesh->dt,
			      mesh->LsimexB[k],
			      ramp,
			      mesh->o_qZ,
			      mesh->o_qY,
			      mesh->o_qYx,
			      mesh->o_qYy,
			      mesh->o_qS,
			      mesh->o_qSx,
			      mesh->o_qSy,
			      mesh->o_pmlqx,
			      mesh->o_pmlqy,
			      mesh->o_q);

	  mesh->device.finish();
      occa::toc("PML_updateKernel");  
      }
    
      if(mesh->nonPmlNelements){

      mesh->device.finish();
      occa::tic("NONPML_updateKernel");   	
	mesh->updateKernel(mesh->nonPmlNelements,
			   mesh->o_nonPmlElementIds,
			   mesh->dt,
			   mesh->LsimexB[k],
			   mesh->o_qZ,
			   mesh->o_qY,
			   mesh->o_qS,
			   mesh->o_q);

	  mesh->device.finish();
      occa::toc("NONPML_updateKernel");  
      }
    
      mesh->device.finish();
      occa::toc("updateKernel");      
    
    }
}

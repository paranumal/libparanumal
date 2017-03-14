#include "boltzmann2D.h"

// complete a time step using Semi-Analytic Adams-Bahforth 
void boltzmannSplitPmlSaabStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer){


 // time
 dfloat t = tstep*mesh->dt;

if(mesh->totalHaloPairs>0){
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
    boltzmannRampFunction2D(t, &ramp, &drampdt);

    mesh->device.finish();
    occa::tic("volumeKernel");
    
    // compute volume contribution to DG boltzmann RHS
    if(mesh->pmlNelements){	
      mesh->pmlVolumeKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    mesh->o_vgeo,
			    mesh->o_sigmax,
			    mesh->o_sigmay,
			    mesh->o_DrT,
			    mesh->o_DsT,
			    mesh->o_q,
			    mesh->o_pmlqx,
			    mesh->o_pmlqy,
			    mesh->o_pmlNT,
			    mesh->o_rhspmlqx,
			    mesh->o_rhspmlqy,
			    mesh->o_rhspmlNT);
    }


   
    
    // compute volume contribution to DG boltzmann RHS
    // added d/dt (ramp(qbar)) to RHS
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

#if CUBATURE_ENABLED
    // compute relaxation terms using cubature
    if(mesh->pmlNelements)
      mesh->pmlRelaxationKernel(mesh->pmlNelements,
			     mesh->o_pmlElementIds,
			     mesh->o_cubInterpT,
			     mesh->o_cubProjectT,
			     mesh->o_q,
			     mesh->o_rhspmlqx,
			     mesh->o_rhspmlqy,
			     mesh->o_rhspmlNT);
  
    // compute relaxation terms using cubature
    if(mesh->nonPmlNelements)
      mesh->relaxationKernel(mesh->nonPmlNelements,
			     mesh->o_nonPmlElementIds,
			     mesh->o_cubInterpT,
			     mesh->o_cubProjectT,
			     mesh->o_q,
			     mesh->o_rhsq);
    
#endif

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
    
  mesh->device.finish();
    occa::tic("updateKernel"); 

 if(tstep==0){

	// if (mesh->pmlNelements)   
	// 	      mesh->pmlUpdateFirstOrderKernel(mesh->pmlNelements,
	// 				mesh->o_pmlElementIds,
	// 			    mesh->dt,
	// 			    ramp,
	// 			    mesh->o_rhspmlqx,
	// 			    mesh->o_rhspmlqy,
	// 			    mesh->o_rhspmlNT,
	// 			    mesh->o_rhspmlqx2,
	// 			    mesh->o_rhspmlqy2,
	// 			    mesh->o_rhspmlNT2,
	// 			    mesh->o_rhspmlqx3,
	// 			    mesh->o_rhspmlqy3,
	// 			    mesh->o_rhspmlNT3,
	// 			    mesh->o_pmlqx,
	// 			    mesh->o_pmlqy,
	// 			    mesh->o_pmlNT,
	// 			    mesh->o_q);

	if(mesh->nonPmlNelements)
		    if(mesh->nonPmlNelements){
				dfloat abc1 = mesh->dt;
		    	dfloat abc2 = 0.0;
		    	dfloat abc3 = 0.0;
		    	
		    	dfloat cc = -mesh->tauInv; 
		    	dfloat dtc = mesh->dt; 
		    	dfloat saabc1 = (pow(cc,3)*pow(dtc,4))/24 + (pow(cc,2)*pow(dtc,3))/6 
		    	               + (cc*pow(dtc,2))/2 + dtc;
		    	dfloat saabc2 = 0.0;
		    	dfloat saabc3 = 0.0; 
		     mesh->updateKernel(mesh->nonPmlNelements,
				 mesh->o_nonPmlElementIds,
				 mesh->dt,
				 mesh->saabexp,
				 abc1,
				 abc2,
				 abc3,
				 saabc1,
				 saabc2,
				 saabc3,
				 mesh->o_rhsq3,
				 mesh->o_rhsq2,
				 mesh->o_rhsq,
				 mesh->o_q);
		}


}
else if(tstep==1){


	// if (mesh->pmlNelements)   
 //  			mesh->pmlUpdateSecondOrderKernel(mesh->pmlNelements,
	// 				mesh->o_pmlElementIds,
	// 			    mesh->dt,
	// 			    ramp,
	// 			    mesh->o_rhspmlqx,
	// 			    mesh->o_rhspmlqy,
	// 			    mesh->o_rhspmlNT,
	// 			    mesh->o_rhspmlqx2,
	// 			    mesh->o_rhspmlqy2,
	// 			    mesh->o_rhspmlNT2,
	// 			    mesh->o_rhspmlqx3,
	// 			    mesh->o_rhspmlqy3,
	// 			    mesh->o_rhspmlNT3,
	// 			    mesh->o_pmlqx,
	// 			    mesh->o_pmlqy,
	// 			    mesh->o_pmlNT,
	// 			    mesh->o_q);


   	if(mesh->nonPmlNelements)
		    if(mesh->nonPmlNelements){
				dfloat abc1 = 3.0*mesh->dt/2.0;
    	        dfloat abc2 = -1.0*mesh->dt/2.0;
    	        dfloat abc3 = 0.0;
		    	
		    	dfloat cc = -mesh->tauInv; 
		    	dfloat dtc = mesh->dt; 
		    	dfloat saabc1 = (pow(cc,3)*pow(dtc,4))/20. + (5.*pow(cc,2)*pow(dtc,3))/24. 
		    	               + (2.*cc*pow(dtc,2))/3. + 3.*dtc/2.;

		    	dfloat saabc2 = -(pow(cc,3)*pow(dtc,4))/120. - (pow(cc,2)*pow(dtc,3))/24. 
		    	               - (cc*pow(dtc,2))/6. - dtc/2.;
		    	dfloat saabc3 = 0.0; 
		     mesh->updateKernel(mesh->nonPmlNelements,
				 mesh->o_nonPmlElementIds,
				 mesh->dt,
				 mesh->saabexp,
				 abc1,
				 abc2,
				 abc3,
				 saabc1,
				 saabc2,
				 saabc3,
				 mesh->o_rhsq3,
				 mesh->o_rhsq2,
				 mesh->o_rhsq,
				 mesh->o_q);
		}


 }

 else{

 // if (mesh->pmlNelements)   
 //  			mesh->pmlUpdateKernel(mesh->pmlNelements,
	// 	    mesh->o_pmlElementIds,
	// 	    mesh->dt,
	// 	    ramp,
	// 	    mesh->o_rhspmlqx,
	// 	    mesh->o_rhspmlqy,
	// 	    mesh->o_rhspmlNT,
	// 	    mesh->o_rhspmlqx2,
	// 	    mesh->o_rhspmlqy2,
	// 	    mesh->o_rhspmlNT2,
	// 	    mesh->o_rhspmlqx3,
	// 	    mesh->o_rhspmlqy3,
	// 	    mesh->o_rhspmlNT3,
	// 	    mesh->o_pmlqx,
	// 	    mesh->o_pmlqy,
	// 	    mesh->o_pmlNT,
	// 	    mesh->o_q);

	if(mesh->nonPmlNelements){
				
		     mesh->updateKernel(mesh->nonPmlNelements,
				 mesh->o_nonPmlElementIds,
				 mesh->dt,
				 mesh->saabexp,
				 mesh->mrab[0],
				 mesh->mrab[1],
				 mesh->mrab[2],
				 mesh->saab[0],
				 mesh->saab[1],
				 mesh->saab[2],
				 mesh->o_rhsq3,
				 mesh->o_rhsq2,
				 mesh->o_rhsq,
				 mesh->o_q);
		}






 }






  
  
}

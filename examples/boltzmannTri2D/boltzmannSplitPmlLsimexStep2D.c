#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannSplitPmlLsimexStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer){

  // LSERK4 stages
  for(iint k=0;k<mesh->Nimex;++k)
  {

    // intermediate stage time
     dfloat t = tstep*mesh->dt + mesh->dt*mesh->LsimexC[k];
   
    
    
    dfloat ramp, drampdt;
    boltzmannRampFunction2D(t, &ramp, &drampdt);
    
   
    mesh->device.finish();
    occa::tic("residualUpdateKernel");
    
    // compute volume contribution to DG boltzmann RHS
    if(mesh->pmlNelements){
      mesh->pmlResidualUpdateKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    mesh->dt,
			    ramp,
			    mesh->LsimexABi[k],
			    mesh->LsimexABe[k],
			    mesh->o_pmlqx,
			    mesh->o_pmlqy,
			    mesh->o_pmlNT,
			    mesh->o_qZx,
			    mesh->o_qZy,
			    mesh->o_qZnt,
			    mesh->o_qYx,
			    mesh->o_qYy,
			    mesh->o_qYnt,
			    mesh->o_qY);
    }
    if(mesh->nonPmlNelements){
      mesh->residualUpdateKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->LsimexABi[k],
			 mesh->LsimexABe[k],
			 mesh->o_q,
			 mesh->o_qZ,
			 mesh->o_qY);
    }
    
     
    
    mesh->device.finish();
    occa::toc("residualUpdateKernel");
    


    //Implicit Solve Satge

     mesh->device.finish();
    occa::tic("implicitSolve");
    
    // compute volume contribution to DG boltzmann RHS
    if(mesh->pmlNelements){
      mesh->pmlNRIterationKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    mesh->dt,
			    ramp,
			    mesh->LsimexAd[k],
			    mesh->o_cubInterpT,
		        mesh->o_cubProjectT,
			    mesh->o_qYx,
			    mesh->o_qYy,
			    mesh->o_qYnt,
			    mesh->o_qZx,
			    mesh->o_qZy,
			    mesh->o_qZnt,
			    mesh->o_qZ);

      mesh->pmlImplicitVolumeKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    ramp,
			    mesh->o_cubInterpT,
			    mesh->o_cubProjectT,
			    mesh->o_qZx,
			    mesh->o_qZy,
			    mesh->o_qZnt,
			    mesh->o_qZ);
     // No surface term for implicit part
      mesh->pmlImplicitUpdateKernel(mesh->pmlNelements,
			  mesh->o_pmlElementIds,
             mesh->dt,
             ramp,
             mesh->LsimexAd[k],
             mesh->o_qZx,
			 mesh->o_qZy,
			 mesh->o_qZnt,
			 mesh->o_qYx,
			 mesh->o_qYy,
			 mesh->o_qYnt,
			 mesh->o_pmlqx,
			 mesh->o_pmlqy,
			 mesh->o_pmlNT,
			 mesh->o_qSx,
			 mesh->o_qSy,
			 mesh->o_qSnt,
			 mesh->o_q);
    }
    
    // compute volume contribution to DG boltzmann RHS
    // added d/dt (ramp(qbar)) to RHS
    if(mesh->nonPmlNelements)
    {
      mesh->NRIterationKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->LsimexAd[k],
			 mesh->o_cubInterpT,
		     mesh->o_cubProjectT,
		     mesh->o_qY,
		     mesh->o_qZ);

      mesh->implicitVolumeKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->o_cubInterpT,
			 mesh->o_cubProjectT,
			 mesh->o_qZ);	
      //No surface term for implicit part
      mesh->implicitUpdateKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
             mesh->dt,
             mesh->LsimexAd[k],
             mesh->o_qZ,
             mesh->o_qY,
             mesh->o_q,
             mesh->o_qS);		
    }
    
    
    mesh->device.finish();
    occa::toc("implicitSolve");



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

   

    mesh->device.finish();
    occa::tic("volumeKernel");
    
    //compute volume contribution to DG boltzmann RHS
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
			    mesh->o_qYx,
			    mesh->o_qYy,
			    mesh->o_qYnt);
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
			 mesh->o_qY);
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
    
    mesh->device.finish();
    occa::tic("surfaceKernel");
    if (mesh->pmlNelements){  
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
			   mesh->o_qYx,
			   mesh->o_qYy);
}
    
    if(mesh->nonPmlNelements)
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
    occa::toc("surfaceKernel");
   
    
    mesh->device.finish();
    occa::tic("updateKernel");
    
    //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    if (mesh->pmlNelements){   
      mesh->pmlUpdateKernel(mesh->pmlNelements,
			    mesh->o_pmlElementIds,
			    mesh->dt,
			    mesh->LsimexB[k],
			    ramp,
			    mesh->o_qZx,
			    mesh->o_qZy,
			    mesh->o_qZnt,
			    mesh->o_qYx,
			    mesh->o_qYy,
			    mesh->o_qYnt,
			    mesh->o_qSx,
			    mesh->o_qSy,
			    mesh->o_qSnt,
			    mesh->o_pmlqx,
			    mesh->o_pmlqy,
			    mesh->o_pmlNT,
			    mesh->o_q);
  }
    
    if(mesh->nonPmlNelements){
      mesh->updateKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->LsimexB[k],
			 mesh->o_qZ,
			 mesh->o_qY,
			 mesh->o_qS,
			 mesh->o_q);
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");      
    
  }
}

#include "boltzmann2D.h"

// complete a time step using LSERK4
void boltzmannSplitPmlSark54Step2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer){


    // STAGE1
    dfloat t = tstep*mesh->dt + mesh->dt*mesh->rkc[0];

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
			 mesh->o_rhsq1);
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
			     mesh->o_rhsq1);
    
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
     
     if(mesh->pmlNelements)
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
			  mesh->o_rhsq1);
    
    mesh->device.finish();
    occa::toc("surfaceKernel");


    mesh->device.finish();
    occa::tic("updateKernel");
    
    //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    // if (mesh->pmlNelements)   
    //   mesh->pmlUpdateKernel(mesh->pmlNelements,
			 //    mesh->o_pmlElementIds,
			 //    mesh->dt,
			 //    mesh->rka[rk],
			 //    mesh->rkb[rk],
			 //    rampUpdate,
			 //    mesh->o_rhspmlqx,
			 //    mesh->o_rhspmlqy,
			 //    mesh->o_rhspmlNT,
			 //    mesh->o_respmlqx,
			 //    mesh->o_respmlqy,
			 //    mesh->o_respmlNT,
			 //    mesh->o_pmlqx,
			 //    mesh->o_pmlqy,
			 //    mesh->o_pmlNT,
			 //    mesh->o_q);
    
    if(mesh->nonPmlNelements){
      mesh->updateStageKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->sarke[1],
			 mesh->rk4a[1][0], //a21
			 mesh->sarka[1][0], //a21
			 mesh->rk4a[0][0], 
			 mesh->sarka[0][0],
			 mesh->rk4a[0][0], 
			 mesh->sarka[0][0],
			 mesh->rk4a[0][0], 
			 mesh->sarka[0][0],
			 mesh->o_rhsq1,
			 mesh->o_rhsq2,
			 mesh->o_rhsq3,
			 mesh->o_rhsq4,
			 mesh->o_qold,
			 mesh->o_q);
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");      
   


   


   // Stage 2

      t = tstep*mesh->dt + mesh->dt*mesh->rkc[1];

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
			 mesh->o_rhsq2);
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
			     mesh->o_rhsq2);
    
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
     
     if(mesh->pmlNelements)
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
			  mesh->o_rhsq2);
    
    mesh->device.finish();
    occa::toc("surfaceKernel");


    mesh->device.finish();
    occa::tic("updateKernel");
    
    //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    // if (mesh->pmlNelements)   
    //   mesh->pmlUpdateKernel(mesh->pmlNelements,
			 //    mesh->o_pmlElementIds,
			 //    mesh->dt,
			 //    mesh->rka[rk],
			 //    mesh->rkb[rk],
			 //    rampUpdate,
			 //    mesh->o_rhspmlqx,
			 //    mesh->o_rhspmlqy,
			 //    mesh->o_rhspmlNT,
			 //    mesh->o_respmlqx,
			 //    mesh->o_respmlqy,
			 //    mesh->o_respmlNT,
			 //    mesh->o_pmlqx,
			 //    mesh->o_pmlqy,
			 //    mesh->o_pmlNT,
			 //    mesh->o_q);
    
     if(mesh->nonPmlNelements){
      mesh->updateStageKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->sarke[2],
			 mesh->rk4a[2][0], //a21
			 mesh->sarka[2][0], //a21
			 mesh->rk4a[2][1], 
			 mesh->sarka[2][1],
			 mesh->rk4a[0][0], 
			 mesh->sarka[0][0],
			 mesh->rk4a[0][0], 
			 mesh->sarka[0][0],
			 mesh->o_rhsq1,
			 mesh->o_rhsq2,
			 mesh->o_rhsq3,
			 mesh->o_rhsq4,
			 mesh->o_qold,
			 mesh->o_q);
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");      
   


    
    // Stage 3
      t = tstep*mesh->dt + mesh->dt*mesh->rkc[2];

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
			 mesh->o_rhsq3);
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
			     mesh->o_rhsq3);
    
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
     
     if(mesh->pmlNelements)
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
			  mesh->o_rhsq3);
    
    mesh->device.finish();
    occa::toc("surfaceKernel");


    mesh->device.finish();
    occa::tic("updateKernel");
    
    //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    // if (mesh->pmlNelements)   
    //   mesh->pmlUpdateKernel(mesh->pmlNelements,
			 //    mesh->o_pmlElementIds,
			 //    mesh->dt,
			 //    mesh->rka[rk],
			 //    mesh->rkb[rk],
			 //    rampUpdate,
			 //    mesh->o_rhspmlqx,
			 //    mesh->o_rhspmlqy,
			 //    mesh->o_rhspmlNT,
			 //    mesh->o_respmlqx,
			 //    mesh->o_respmlqy,
			 //    mesh->o_respmlNT,
			 //    mesh->o_pmlqx,
			 //    mesh->o_pmlqy,
			 //    mesh->o_pmlNT,
			 //    mesh->o_q);
    
     if(mesh->nonPmlNelements){
      mesh->updateStageKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->sarke[3],
			 mesh->rk4a[3][0], //a21
			 mesh->sarka[3][0], //a21
			 mesh->rk4a[3][1], 
			 mesh->sarka[3][1],
			 mesh->rk4a[3][2], 
			 mesh->sarka[3][2],
			 mesh->rk4a[0][0], 
			 mesh->sarka[0][0],
			 mesh->o_rhsq1,
			 mesh->o_rhsq2,
			 mesh->o_rhsq3,
			 mesh->o_rhsq4,
			 mesh->o_qold,
			 mesh->o_q);
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");   


    
    // Stage 4
      t = tstep*mesh->dt + mesh->dt*mesh->rkc[3];

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
			 mesh->o_rhsq4);
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
			     mesh->o_rhsq4);
    
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
     
     if(mesh->pmlNelements)
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
			  mesh->o_rhsq4);
    
    mesh->device.finish();
    occa::toc("surfaceKernel");


    mesh->device.finish();
    occa::tic("updateKernel");
    
    //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    // if (mesh->pmlNelements)   
    //   mesh->pmlUpdateKernel(mesh->pmlNelements,
			 //    mesh->o_pmlElementIds,
			 //    mesh->dt,
			 //    mesh->rka[rk],
			 //    mesh->rkb[rk],
			 //    rampUpdate,
			 //    mesh->o_rhspmlqx,
			 //    mesh->o_rhspmlqy,
			 //    mesh->o_rhspmlNT,
			 //    mesh->o_respmlqx,
			 //    mesh->o_respmlqy,
			 //    mesh->o_respmlNT,
			 //    mesh->o_pmlqx,
			 //    mesh->o_pmlqy,
			 //    mesh->o_pmlNT,
			 //    mesh->o_q);
    
     if(mesh->nonPmlNelements){
      mesh->updateStageKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->sarke[4],
			 mesh->rk4a[4][0], //a21
			 mesh->sarka[4][0], //a21
			 mesh->rk4a[4][1], 
			 mesh->sarka[4][1],
			 mesh->rk4a[4][2], 
			 mesh->sarka[4][2],
			 mesh->rk4a[4][3], 
			 mesh->sarka[4][3],
			 mesh->o_rhsq1,
			 mesh->o_rhsq2,
			 mesh->o_rhsq3,
			 mesh->o_rhsq4,
			 mesh->o_qold,
			 mesh->o_q);
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");   




    // Stage 4
      t = tstep*mesh->dt + mesh->dt*mesh->rkc[4];

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
			 mesh->o_rhsq5);
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
			     mesh->o_rhsq5);
    
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
     
     if(mesh->pmlNelements)
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
			  mesh->o_rhsq5);
    
    mesh->device.finish();
    occa::toc("surfaceKernel");


    mesh->device.finish();
    occa::tic("updateKernel");
    
    //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    // if (mesh->pmlNelements)   
    //   mesh->pmlUpdateKernel(mesh->pmlNelements,
			 //    mesh->o_pmlElementIds,
			 //    mesh->dt,
			 //    mesh->rka[rk],
			 //    mesh->rkb[rk],
			 //    rampUpdate,
			 //    mesh->o_rhspmlqx,
			 //    mesh->o_rhspmlqy,
			 //    mesh->o_rhspmlNT,
			 //    mesh->o_respmlqx,
			 //    mesh->o_respmlqy,
			 //    mesh->o_respmlNT,
			 //    mesh->o_pmlqx,
			 //    mesh->o_pmlqy,
			 //    mesh->o_pmlNT,
			 //    mesh->o_q);
    
     if(mesh->nonPmlNelements){
      mesh->updateKernel(mesh->nonPmlNelements,
			 mesh->o_nonPmlElementIds,
			 mesh->dt,
			 mesh->sarke[5],
			 mesh->rk4b[0], //a21
			 mesh->sarkb[0], //a21
			 mesh->rk4b[1], //a21
			 mesh->sarkb[1], //a21
			 mesh->rk4b[2], //a21
			 mesh->sarkb[2], //a21
			 mesh->rk4b[3], //a21
			 mesh->sarkb[3], //a21
			 mesh->rk4b[4], //a21
			 mesh->sarkb[4], //a21
			 mesh->o_rhsq1,
			 mesh->o_rhsq2,
			 mesh->o_rhsq3,
			 mesh->o_rhsq4,
			 mesh->o_rhsq5,
			 mesh->o_qold,
			 mesh->o_q);
    }
    
    mesh->device.finish();
    occa::toc("updateKernel");   
















  
   
}

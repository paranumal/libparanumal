#include "acoustics2D.h"

void acousticsRun2D(mesh2D *mesh){

  // MPI send buffer
  int haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(int rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat t = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];
      
      if(mesh->totalHaloPairs>0){
	// extract halo node data
	meshHaloExtract(mesh, mesh->Np*mesh->Nfields*sizeof(dfloat),
			  mesh->q, sendBuffer);
	
	// start halo exchange
	meshHaloExchangeStart(mesh,
			      mesh->Np*mesh->Nfields*sizeof(dfloat),
			      sendBuffer,
			      recvBuffer);
      }
      
      // compute volume contribution to DG acoustics RHS
      acousticsVolume2D(mesh);

      // compute volume contributions to from PML
      acousticsPml2D(mesh);
      
      if(mesh->totalHaloPairs>0){
	// wait for halo data to arrive
	meshHaloExchangeFinish(mesh);
	
	// copy data to halo zone of q
	memcpy(mesh->q+mesh->Np*mesh->Nfields*mesh->Nelements, recvBuffer, haloBytes);
      }
      
      // compute surface contribution to DG acoustics RHS
      acousticsSurface2D(mesh, t);
      
      // update solution using Euler Forward
      acousticsUpdate2D(mesh, mesh->rka[rk], mesh->rkb[rk]);

      // update pml variables
      acousticsPmlUpdate2D(mesh, mesh->rka[rk], mesh->rkb[rk]);
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0){
      
      // do error stuff on host
      acousticsError2D(mesh, mesh->dt*(tstep+1));
      
      // output field files
      int fld = 2;
      meshPlotVTU2D(mesh, "foo", fld);
    }
    
  }

  free(recvBuffer);
  free(sendBuffer);
}


void acousticsOccaRun2D(mesh2D *mesh){

  // MPI send buffer
  int haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(int rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat t = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      //printf("mesh->totalHaloPairs = %d\n",mesh->totalHaloPairs);
      if(mesh->totalHaloPairs>0){
	// extract halo on DEVICE
	int Nentries = mesh->Np*mesh->Nfields;
	
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
			  t,
			  mesh->o_x,
			  mesh->o_y,
			  mesh->o_q,
			  mesh->o_rhsq);
      
      // update solution using Runge-Kutta
#if 0
      
      mesh->updateKernel(mesh->Nelements*mesh->Np*mesh->Nfields,
			 mesh->dt,
			 mesh->rka[rk],
			 mesh->rkb[rk],
			 mesh->o_rhsq,
			 mesh->o_resq,
			 mesh->o_q);
      
#else // WADG
      //printf("running update on step %d\n",tstep);
      
      mesh->updateKernel(mesh->Nelements,
			 mesh->dt,
			 mesh->rka[rk],
			 mesh->rkb[rk],
			 mesh->cubNp,
			 mesh->o_cubInterpT,
			 mesh->o_cubProjectT,
			 mesh->o_c2,
			 mesh->o_rhsq,
			 mesh->o_resq,
			 mesh->o_q);
#endif
      
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0){

      // copy data back to host
      mesh->o_q.copyTo(mesh->q);

      // do error stuff on host
      acousticsError2D(mesh, mesh->dt*(tstep+1));

      // output field files
      int fld = 2;
      meshPlotVTU2D(mesh, "foo", fld);
    }
  }

  free(recvBuffer);
  free(sendBuffer);
}


// async HOST<> DEVICE copies contributed by Nigel Nunn
void acousticsOccaAsyncRun2D(mesh2D *mesh){

  // MPI send buffer
  int haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);

  // offset for halo data
  size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements * sizeof(dfloat);
  int Nentries = mesh->Np*mesh->Nfields;
  
  // using mapped memory for host-device transfers
  occa::memory o_pinned_send_Q, o_pinned_recv_Q;
  dfloat *H_send_Q=NULL, *H_recv_Q=NULL;

  if (haloBytes>0) {

    // OpenCL ties mapped memory to "current" stream,
    // so here we make the secondary stream current:
    mesh->device.setStream(mesh->stream1);

    // allocate mapped memory
    o_pinned_send_Q = mesh->device.mappedAlloc(haloBytes, NULL);
    o_pinned_recv_Q = mesh->device.mappedAlloc(haloBytes, NULL);

    // get pointers
    H_send_Q = (dfloat*) o_pinned_send_Q.getMappedPointer();
    H_recv_Q = (dfloat*) o_pinned_recv_Q.getMappedPointer();

    // flush current stream
    mesh->device.finish();

    // revert to main stream
    mesh->device.setStream(mesh->stream0);
    mesh->device.finish();
  }

  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(int rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat t = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      if(mesh->totalHaloPairs>0){

	// make sure that the main stream is synced before halo exchange 
	mesh->device.finish();
	
        // switch to secondary stream
        mesh->device.setStream(mesh->stream1);
	
	// extract halo on DEVICE
	mesh->haloExtractKernel(mesh->totalHaloPairs,
				Nentries,
				mesh->o_haloElementList,
				mesh->o_q,
				mesh->o_haloBuffer);

        // (async) copy extracted halo to HOST (2nd stream)
        mesh->o_haloBuffer.asyncCopyTo(H_send_Q);

	// switch back to main stream
        mesh->device.setStream(mesh->stream0);
      }

      // compute volume contribution to DG acoustics RHS
      mesh->volumeKernel(mesh->Nelements,
			 mesh->o_vgeo,
			 mesh->o_DrT,
			 mesh->o_DsT,
			 mesh->o_q,
			 mesh->o_rhsq);
      
      if(mesh->totalHaloPairs>0){

        // NBN: switch to secondary stream
        mesh->device.setStream(mesh->stream1);

        // NBN: ensure "send" buffer is loaded
        mesh->device.finish();
        
        // start halo exchange (NBN: moved from above)
        meshHaloExchangeStart(mesh,
			      mesh->Np*mesh->Nfields * sizeof(dfloat),
			      H_send_Q,
			      H_recv_Q);

        // wait for halo data to arrive
        meshHaloExchangeFinish(mesh);

        // copy halo data to DEVICE
        // NBN: making this copy "async" means the copy will be 
        // queued on stream0, then surf and update kernels will 
        // immediately be queued on the same stream.
        mesh->o_q.asyncCopyFrom(H_recv_Q,haloBytes,offset);
	
	// NBN: ensure async copy has finished
        mesh->device.finish();

        // NBN: run remaining work on stream0
        mesh->device.setStream(mesh->stream0);  
      }   
      
      // compute surface contribution to DG acoustics RHS
      mesh->surfaceKernel(mesh->Nelements,
			  mesh->o_sgeo,
			  mesh->o_LIFTT,
			  mesh->o_vmapM,
			  mesh->o_vmapP,
			  mesh->o_EToB,
			  t,
			  mesh->o_x,
			  mesh->o_y,
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
      acousticsError2D(mesh, mesh->dt*(tstep+1));

      // output field files
      int fld = 2;
      meshPlotVTU2D(mesh, "foo", fld);
    }
  }

  if (haloBytes>0){  // NBN: clean-up

    // OpenCL: mapped memory is tied to 2nd stream
    mesh->device.setStream(mesh->stream1);
    mesh->device.finish();

    o_pinned_send_Q.free();
    o_pinned_recv_Q.free();      // release mapped memory

    H_recv_Q = NULL;              // invalidate pointers
    H_send_Q = NULL;

    mesh->device.setStream(mesh->stream0);
  }
}


// async HOST<> DEVICE copies contributed by Nigel Nunn
// modified to split surface into interior and not-interior elements
void acousticsSplitSurfaceOccaAsyncRun2D(mesh2D *mesh){

  // MPI send buffer
  int haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);

  // offset for halo data
  size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements * sizeof(dfloat);
  int Nentries = mesh->Np*mesh->Nfields;
  
  // using mapped memory for host-device transfers
  occa::memory o_pinned_send_Q, o_pinned_recv_Q;
  dfloat *H_send_Q=NULL, *H_recv_Q=NULL;

  if (haloBytes>0) {

    // OpenCL ties mapped memory to "current" stream,
    // so here we make the secondary stream current:
    mesh->device.setStream(mesh->stream1);

    // allocate mapped memory
    o_pinned_send_Q = mesh->device.mappedAlloc(haloBytes, NULL);
    o_pinned_recv_Q = mesh->device.mappedAlloc(haloBytes, NULL);

    // get pointers
    H_send_Q = (dfloat*) o_pinned_send_Q.getMappedPointer();
    H_recv_Q = (dfloat*) o_pinned_recv_Q.getMappedPointer();

    // flush current stream
    mesh->device.finish();

    // revert to main stream
    mesh->device.setStream(mesh->stream0);
    mesh->device.finish();
  }

  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(int rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat t = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      if(mesh->totalHaloPairs>0){

	// make sure that the main stream is synced before halo exchange 
	mesh->device.finish();
	
        // switch to secondary stream
        mesh->device.setStream(mesh->stream1);
	
	// extract halo on DEVICE
	mesh->haloExtractKernel(mesh->totalHaloPairs,
				Nentries,
				mesh->o_haloElementList,
				mesh->o_q,
				mesh->o_haloBuffer);

        // (async) copy extracted halo to HOST (2nd stream)
        mesh->o_haloBuffer.asyncCopyTo(H_send_Q);

	// switch back to main stream
        mesh->device.setStream(mesh->stream0);
      }

      // compute volume contribution to DG acoustics RHS
      mesh->volumeKernel(mesh->Nelements,
			 mesh->o_vgeo,
			 mesh->o_DrT,
			 mesh->o_DsT,
			 mesh->o_q,
			 mesh->o_rhsq);
      
      if(mesh->pmlNelements){
	// add PML contributions
	mesh->pmlKernel(mesh->pmlNelements,
			mesh->o_pmlElementList,
			mesh->o_pmlSigmaX,
			mesh->o_pmlSigmaY,
			mesh->o_q,
				 mesh->o_pmlq,
			mesh->o_rhsq,
			mesh->o_pmlrhsq);
	
	// update PML variables
	mesh->pmlUpdateKernel(mesh->pmlNelements*mesh->Np*mesh->pmlNfields,
			      mesh->dt,
			      mesh->rka[rk],
			      mesh->rkb[rk],
			      mesh->o_pmlrhsq,
			      mesh->o_pmlresq,
			      mesh->o_pmlq);
      }
      
      // compute surface contribution to DG acoustics RHS
      mesh->partialSurfaceKernel(mesh->NinternalElements,
				 mesh->o_internalElementIds,
				 mesh->o_sgeo,
				 mesh->o_LIFTT,
				 mesh->o_vmapM,
				 mesh->o_vmapP,
				 mesh->o_EToB,
				 t,
				 mesh->o_x,
				 mesh->o_y,
				 mesh->o_q,
				 mesh->o_rhsq);

      if(mesh->totalHaloPairs>0){

        // NBN: switch to secondary stream
        mesh->device.setStream(mesh->stream1);

        // NBN: ensure "send" buffer is loaded
        mesh->device.finish();
        
        // start halo exchange (NBN: moved from above)
        meshHaloExchangeStart(mesh,
			      mesh->Np*mesh->Nfields * sizeof(dfloat),
			      H_send_Q,
			      H_recv_Q);

        // wait for halo data to arrive
        meshHaloExchangeFinish(mesh);

        // copy halo data to DEVICE on stream 1
        mesh->o_q.asyncCopyFrom(H_recv_Q,haloBytes,offset);
	
	// NBN: ensure async copy has finished
        mesh->device.finish();
	
	// NBN: run remaining work on stream0
        mesh->device.setStream(mesh->stream0);  
      }   
      
      if(mesh->NnotInternalElements){
	// compute surface contribution to DG acoustics RHS on stream 0
	mesh->partialSurfaceKernel(mesh->NnotInternalElements,
				   mesh->o_notInternalElementIds,
				   mesh->o_sgeo,
				   mesh->o_LIFTT,
				   mesh->o_vmapM,
				   mesh->o_vmapP,
				   mesh->o_EToB,
				   t,
				   mesh->o_x,
				   mesh->o_y,
				   mesh->o_q,
				   mesh->o_rhsq);
      }
      
      
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
      acousticsError2D(mesh, mesh->dt*(tstep+1));

      // output field files
      int fld = 2;
      meshPlotVTU2D(mesh, "foo", fld);
    }
  }

  if (haloBytes>0){  // NBN: clean-up

    // OpenCL: mapped memory is tied to 2nd stream
    mesh->device.setStream(mesh->stream1);
    mesh->device.finish();

    o_pinned_send_Q.free();
    o_pinned_recv_Q.free();      // release mapped memory

    H_recv_Q = NULL;              // invalidate pointers
    H_send_Q = NULL;

    mesh->device.setStream(mesh->stream0);
  }
}

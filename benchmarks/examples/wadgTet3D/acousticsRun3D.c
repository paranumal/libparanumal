#include "acoustics3D.h"

void acousticsRun3D(mesh3D *mesh){

  // MPI send buffer
  dfloat *sendBuffer =
    (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields,
		     sizeof(dfloat));
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(iint rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat time = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      // halo exchange
      meshHaloExchange(mesh,
		       mesh->Np*mesh->Nfields,
		       mesh->q,
		       sendBuffer,
		       mesh->q+mesh->Np*mesh->Nfields*mesh->Nelements);
      
      // compute volume contribution to DG acoustics RHS
      acousticsVolume3D(mesh);

      // compute surface contribution to DG acoustics RHS
      acousticsSurface3D(mesh, time);

      // update solution using LSERK4
      acousticsUpdate3D(mesh, mesh->rka[rk], mesh->rkb[rk]);
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0)
      acousticsError3D(mesh, mesh->dt*(tstep+1));
  }

  free(sendBuffer);
}

void acousticsOccaRun3D(mesh3D *mesh){

  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(iint rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat time = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      //printf("mesh->totalHaloPairs = %d\n",mesh->totalHaloPairs);
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
      
      // compute volume contribution to DG acoustics RHS
      mesh->volumeKernel(mesh->Nelements,
			 mesh->o_vgeo,
			 mesh->o_DrT,
			 mesh->o_DsT,
			 mesh->o_DtT,
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
			  time,
			  mesh->o_x,
			  mesh->o_y,
			  mesh->o_z,
			  mesh->o_q,
			  mesh->o_rhsq);
      
#if 0
      // update solution using Runge-Kutta
      mesh->updateKernel(mesh->Nelements*mesh->Np*mesh->Nfields,
			 mesh->dt,
			 mesh->rka[rk],
			 mesh->rkb[rk],
			 mesh->o_rhsq,
			 mesh->o_resq,
			 mesh->o_q);

#else // wadg
      // update solution using Runge-Kutta
      mesh->updateKernel(mesh->Nelements,
			 mesh->dt,
			 mesh->rka[rk],
			 mesh->rkb[rk],
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
      acousticsError3D(mesh, mesh->dt*(tstep+1));

      // output field files
      iint fld = 2;
      meshPlotVTU3D(mesh, "foo", fld);
    }
  }

  free(recvBuffer);
  free(sendBuffer);
}


void acousticsTimeWadg3D(mesh3D *mesh){

  occa::kernelInfo kernelInfo;
  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_N", mesh->N);
  kernelInfo.addDefine("p_Np", mesh->Np);
  kernelInfo.addDefine("p_cubNp", mesh->cubNp);  
  kernelInfo.addDefine("p_Nfp", mesh->Nfp);
  kernelInfo.addDefine("p_Nfaces", mesh->Nfaces);
  kernelInfo.addDefine("p_NfacesNfp", mesh->Nfp*mesh->Nfaces);
  kernelInfo.addDefine("p_Nvgeo", mesh->Nvgeo);
  kernelInfo.addDefine("p_Nsgeo", mesh->Nsgeo);
  
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);
  
  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  if(sizeof(dfloat)==4){
    kernelInfo.addDefine("dfloat","float");
    kernelInfo.addDefine("dfloat4","float4");
  }
  if(sizeof(dfloat)==8){
    kernelInfo.addDefine("dfloat","double");
    kernelInfo.addDefine("dfloat4","double4");
  }

  if(sizeof(iint)==4){
    kernelInfo.addDefine("iint","int");
  }
  if(sizeof(iint)==8){
    kernelInfo.addDefine("iint","long long int");
  }

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo.addCompilerFlag("--ftz=true");
    kernelInfo.addCompilerFlag("--prec-div=false");
    kernelInfo.addCompilerFlag("--prec-sqrt=false");
    kernelInfo.addCompilerFlag("--use_fast_math");
    kernelInfo.addCompilerFlag("--fmad=true"); // compiler option for cuda
  }


  
  dfloat * invMc = (dfloat*) calloc(mesh->Np*mesh->Np*mesh->Nelements,sizeof(dfloat));  
  for(iint e=0;e<mesh->Nelements;++e){
    for (int j = 0; j < mesh->Np; ++j){
      for (int i = 0; i < mesh->Np; ++i){
	// store arbitrary values for testing
	invMc[i + j*mesh->Np + e*mesh->Np*mesh->Np] =
	  (dfloat) i+j+e;
      }
    }
  }  
  occa::memory o_invMc =
    mesh->device.malloc(mesh->Np*mesh->Np*mesh->Nelements*sizeof(dfloat),
			invMc);  
  
  occa::kernel applyWADG =
    mesh->device.buildKernelFromSource("./examples/wadgTet3D/wadgKernels.okl",
				       "applyWADG",
				       kernelInfo);  

  occa::kernel applyMassInv =
    mesh->device.buildKernelFromSource("./examples/wadgTet3D/wadgKernels.okl",
				       "applyMassInv",
				       kernelInfo);  
  
  // load/store comparison
  double wadgTime = 0.0, massInvTime = 0.0;
  
  int Nsteps = 10;
  for(iint i = 0; i < Nsteps; ++i){

    clock_t begin = clock();
    
    applyWADG(mesh->Nelements,
	      mesh->cubNp,
	      mesh->o_cubInterpT,
	      mesh->o_cubProjectT,
	      mesh->o_c2, // weight
	      mesh->o_q);
    mesh->device.finish();
    clock_t end = clock();
    wadgTime += (double)(end - begin) / CLOCKS_PER_SEC;

    begin = clock();
    applyMassInv(mesh->Nelements,		 
		 o_invMc,
		 mesh->o_q);      
    mesh->device.finish();
    end = clock();
    massInvTime += (double)(end - begin) / CLOCKS_PER_SEC;
  }

  printf("wadg time = %g, mass inv time = %g\n",wadgTime,massInvTime);
}


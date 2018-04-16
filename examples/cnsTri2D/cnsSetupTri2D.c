#include "cnsTri2D.h"

cns_t *cnsSetupTri2D(mesh2D *mesh, char *options, char* boundaryHeaderFileName){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,MPI_COMM_WORLD);

  int deviceID = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) deviceID++;
  }

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", deviceID);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 0");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");  

  cns_t *cns = (cns_t*) calloc(1, sizeof(cns_t));

  mesh->Nfields = 3;
  cns->Nfields = mesh->Nfields;
  
  cns->Nstresses = 3;
  cns->mesh = mesh;

  dlong Ntotal = mesh->Nelements*mesh->Np*mesh->Nfields;
  cns->Nblock = (Ntotal+blockSize-1)/blockSize;
  
  hlong localElements = (hlong) mesh->Nelements;
  MPI_Allreduce(&localElements, &(cns->totalElements), 1, MPI_HLONG, MPI_SUM, MPI_COMM_WORLD);

  // mean flow
  cns->rbar = 1;
  cns->ubar = 0.2;
  cns->vbar = 0;

  // viscosity
  dfloat Re = 5000;
  cns->mu = cns->ubar/Re; // assumes reference length 1

  // speed of sound (assuming isothermal unit bulk flow) = sqrt(RT)
  dfloat mach = 0.17;
  cns->RT = cns->ubar*cns->ubar/(mach*mach);
  
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  cns->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  
  if (strstr(options,"LSERK")) {
    cns->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
		  		sizeof(dfloat));
  }

  if (strstr(options,"DOPRI5")) {
    int NrkStages = 7;
    cns->rkq  = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
          sizeof(dfloat));
    cns->rkrhsq = (dfloat*) calloc(NrkStages*mesh->Nelements*mesh->Np*mesh->Nfields,
          sizeof(dfloat));
    cns->rkerr  = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
          sizeof(dfloat));

    cns->errtmp = (dfloat*) calloc(cns->Nblock, sizeof(dfloat));
  }

  cns->viscousStresses = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*cns->Nstresses,
					   sizeof(dfloat));

  cns->Vort = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  // fix this later (initial conditions)
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];

      dlong qbase = e*mesh->Np*mesh->Nfields + n;

#if 0
      cnsGaussianPulse2D(x, y, t,
			 mesh->q+qbase,
			 mesh->q+qbase+mesh->Np,
			 mesh->q+qbase+2*mesh->Np);
#else
      mesh->q[qbase+0*mesh->Np] = cns->rbar;
      //      mesh->q[qbase+1*mesh->Np] = cns->rbar*cns->ubar*y*(6-y)/9.;
      const dfloat alpha = 20;
      mesh->q[qbase+1*mesh->Np] = cns->rbar*cns->ubar; //*tanh(alpha*y)*tanh(alpha*(6-y));
      mesh->q[qbase+2*mesh->Np] = cns->rbar*cns->vbar;
#endif
    }
  }

  // set penalty parameter
  mesh->Lambda2 = 0.5;
  
  // set time step
  dfloat hmin = 1e9;
  for(dlong e=0;e<mesh->Nelements;++e){  

    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      if(invJ<0) printf("invJ = %g\n", invJ);
      
      // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
      // h = 0.5/(sJ/J)
      
      dfloat hest = .5/(sJ*invJ);

      hmin = mymin(hmin, hest);
    }
  }

  // need to change cfl and defn of dt
  dfloat cfl = 0.5; // depends on the stability region size

  dfloat dtAdv  = hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(cns->RT));
  dfloat dtVisc = pow(hmin, 2)/(pow(mesh->N+1,4)*cns->mu);

  dfloat dt = cfl*mymin(dtAdv, dtVisc);
  dt = cfl*dtAdv;
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  
  //
  mesh->finalTime = 2;
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  if(strstr(options,"LSERK")) {
    mesh->dt = mesh->finalTime/mesh->NtimeSteps;
  }

  if (rank ==0) printf("dtAdv = %lg (before cfl), dtVisc = %lg (before cfl), dt = %lg\n",
   dtAdv, dtVisc, dt);

  cns->frame = 0;
  // errorStep
  mesh->errorStep = 1000;

  if (rank ==0) printf("dt = %g\n", mesh->dt);

  // OCCA build stuff
  
  occa::kernelInfo kernelInfo;
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  //add boundary data to kernel info
  kernelInfo.addInclude(boundaryHeaderFileName);
 
  cns->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);
  
  cns->o_viscousStresses =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*cns->Nstresses*sizeof(dfloat),
			cns->viscousStresses);
  
  cns->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), cns->rhsq);

  if (strstr(options,"LSERK")) {
    cns->o_resq =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), cns->resq);
  }

  if (strstr(options,"DOPRI5")) {
    int NrkStages = 7;
    cns->o_rkq =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), cns->rkq);
    cns->o_rkrhsq =
      mesh->device.malloc(NrkStages*mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), cns->rkrhsq);
    cns->o_rkerr =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), cns->rkerr);
  
    cns->o_errtmp = mesh->device.malloc(cns->Nblock*sizeof(dfloat), cns->errtmp);
  }

  
  cns->o_Vort = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), cns->Vort);
  

  if(mesh->totalHaloPairs>0){
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));

    cns->o_haloStressesBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*cns->Nstresses*sizeof(dfloat));
  
    // MPI send buffer
    cns->haloBytes = mesh->totalHaloPairs*mesh->Np*cns->Nfields*sizeof(dfloat);
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(cns->haloBytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(cns->haloBytes, NULL);
    cns->o_haloBuffer = mesh->device.malloc(cns->haloBytes);
    cns->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    cns->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();

    cns->haloStressesBytes = mesh->totalHaloPairs*mesh->Np*cns->Nstresses*sizeof(dfloat);
    occa::memory o_sendStressesBuffer = mesh->device.mappedAlloc(cns->haloStressesBytes, NULL);
    occa::memory o_recvStressesBuffer = mesh->device.mappedAlloc(cns->haloStressesBytes, NULL);
    cns->o_haloStressesBuffer = mesh->device.malloc(cns->haloStressesBytes);
    cns->sendStressesBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    cns->recvStressesBuffer = (dfloat*) o_recvBuffer.getMappedPointer();
  }

  //  p_RT, p_rbar, p_ubar, p_vbar
  // p_half, p_two, p_third, p_Nstresses
  
  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_Nstresses", cns->Nstresses);

  kernelInfo.addDefine("p_RT", cns->RT);

  dfloat sqrtRT = sqrt(cns->RT);
  kernelInfo.addDefine("p_sqrtRT", sqrtRT);
  
  kernelInfo.addDefine("p_rbar", cns->rbar);
  kernelInfo.addDefine("p_ubar", cns->ubar);
  kernelInfo.addDefine("p_vbar", cns->vbar);
  

  const dfloat p_one = 1.0, p_two = 2.0, p_half = 1./2., p_third = 1./3., p_zero = 0;

  kernelInfo.addDefine("p_two", p_two);
  kernelInfo.addDefine("p_one", p_one);
  kernelInfo.addDefine("p_half", p_half);
  kernelInfo.addDefine("p_third", p_third);
  kernelInfo.addDefine("p_zero", p_zero);
  
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  int cubMaxNodes = mymax(mesh->Np, (mesh->intNfp*mesh->Nfaces));
  kernelInfo.addDefine("p_cubMaxNodes", cubMaxNodes);

  kernelInfo.addDefine("p_Lambda2", 0.5f);

  kernelInfo.addDefine("p_blockSize", blockSize);


  kernelInfo.addParserFlag("automate-add-barriers", "disabled");
  
  cns->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsVolumeTri2D.okl",
				       "cnsVolumeTri2D",
				       kernelInfo);

  cns->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsSurfaceTri2D.okl",
				       "cnsSurfaceTri2D",
				       kernelInfo);

  cns->cubatureVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsCubatureVolumeTri2D.okl",
               "cnsCubatureVolumeTri2D",
               kernelInfo);

  cns->cubatureSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsCubatureSurfaceTri2D.okl",
               "cnsCubatureSurfaceTri2D",
               kernelInfo);

  cns->stressesVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsVolumeTri2D.okl",
				       "cnsStressesVolumeTri2D",
				       kernelInfo);

  cns->stressesSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsSurfaceTri2D.okl",
				       "cnsStressesSurfaceTri2D",
				       kernelInfo);
  
  cns->vorticityKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsVorticityTri2D.okl",
               "cnsVorticityTri2D",
               kernelInfo);


  cns->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsUpdate2D.okl",
				       "cnsUpdate2D",
				       kernelInfo);

  cns->rkUpdateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsUpdate2D.okl",
               "cnsRkUpdate2D",
               kernelInfo);
  cns->rkStageKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsUpdate2D.okl",
               "cnsRkStage2D",
               kernelInfo);
  cns->rkErrorEstimateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsUpdate2D.okl",
               "cnsErrorEstimate2D",
               kernelInfo);

  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
				       "meshHaloExtract2D",
				       kernelInfo);

  return cns;
}

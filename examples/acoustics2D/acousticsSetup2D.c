#include "acoustics2D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS
// #undef USE_2_STREAMS

void acousticsSetup2D(mesh2D *mesh){

  mesh->Nfields = 4;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  
  // fix this later (initial conditions)
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];

#if 0
      acousticsCavitySolution2D(x, y, t,
				mesh->q+cnt, mesh->q+cnt+1, mesh->q+cnt+2);

      acousticsGaussianPulse2D(x, y, t,
			       mesh->q+cnt, mesh->q+cnt+1, mesh->q+cnt+2);
#endif
      mesh->q[cnt+2] = 1.;
      
      cnt += mesh->Nfields;
    }
  }

  // set penalty parameter
  mesh->Lambda2 = 0.5;
  
  // set time step
  dfloat hmin = 1e9;
  for(iint e=0;e<mesh->Nelements;++e){  

    for(iint f=0;f<mesh->Nfaces;++f){
      iint sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
      // h = 0.5/(sJ/J)
      
      dfloat hest = .5/(sJ*invJ);

      hmin = mymin(hmin, hest);
    }
  }
  dfloat cfl = .1; // depends on the stability region size

  dfloat dt = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*mesh->Lambda2);

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //
  mesh->finalTime = 4;
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  mesh->dt = mesh->finalTime/mesh->NtimeSteps;

  // errorStep
  mesh->errorStep = 100;

  printf("dt = %g\n", mesh->dt);

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", rank);
  //  sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //  sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);

  occa::kernelInfo kernelInfo;

  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);
  
#if 0  
  mesh->device.setup(deviceConfig);


  // build Dr, Ds, LIFT transposes
  dfloat *DrT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *DsT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      DrT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
      DsT[n+m*mesh->Np] = mesh->Ds[n*mesh->Np+m];
    }
  }

  dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      LIFTT[n+m*mesh->Np] = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
    }
  }

  // find elements that have all neighbors on this process
  iint *internalElementIds = (iint*) calloc(mesh->Nelements, sizeof(iint));
  iint *notInternalElementIds = (iint*) calloc(mesh->Nelements, sizeof(iint));

  iint Ninterior = 0, NnotInterior = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    iint flag = 0;
    for(iint f=0;f<mesh->Nfaces;++f)
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1)
	flag = 1;
    if(!flag)
      internalElementIds[Ninterior++] = e;
    else
      notInternalElementIds[NnotInterior++] = e;
  }

  printf("NinteriorElements = %d, NnotInternalElements = %d\n", Ninterior, NnotInterior);
  
  mesh->NinternalElements = Ninterior;
  mesh->NnotInternalElements = NnotInterior;
  mesh->o_internalElementIds    = mesh->device.malloc(Ninterior*sizeof(iint), internalElementIds);
  if(NnotInterior>0)
    mesh->o_notInternalElementIds = mesh->device.malloc(NnotInterior*sizeof(iint), notInternalElementIds);
  
  // OCCA allocate device memory (remember to go back for halo)
  mesh->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);
  mesh->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
  mesh->o_resq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->resq);

  mesh->o_Dr = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				   mesh->Dr);

  mesh->o_Ds = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				   mesh->Ds);

  mesh->o_DrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				    DrT);

  mesh->o_DsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				    DsT);

  mesh->o_LIFT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			mesh->LIFT);

  mesh->o_LIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			LIFTT);

  
  mesh->o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*sizeof(dfloat),
			mesh->vgeo);
  
  mesh->o_sgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nsgeo*sizeof(dfloat),
			mesh->sgeo);

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
			mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
			mesh->vmapP);

  mesh->o_EToB =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),
			mesh->EToB);

  mesh->o_x =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
			mesh->x);

  mesh->o_y =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
			mesh->y);
  
  if(mesh->totalHaloPairs>0){
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      mesh->device.malloc(mesh->totalHaloPairs*sizeof(iint), mesh->haloElementList);
    
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));
  }

  //-------------------------------------
  // NBN: 2 streams for async MPI updates
  // {Vol, Surf, update}  run on q[0]
  // {halo-get, copy} run on q[1]
  //-------------------------------------
  mesh->stream0 = mesh->device.getStream();
#ifdef USE_2_STREAMS
  mesh->stream1 = mesh->device.createStream();  // NBN: second stream
#else
  mesh->stream1 = mesh->stream0;                // NBN: stream1 == stream0
#endif
  mesh->device.setStream(mesh->stream0);
  //-------------------------------------
  
  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_N", mesh->N);
  kernelInfo.addDefine("p_Np", mesh->Np);
  kernelInfo.addDefine("p_Nfp", mesh->Nfp);
  kernelInfo.addDefine("p_Nfaces", mesh->Nfaces);
  kernelInfo.addDefine("p_NfacesNfp", mesh->Nfp*mesh->Nfaces);
  kernelInfo.addDefine("p_Nvgeo", mesh->Nvgeo);
  kernelInfo.addDefine("p_Nsgeo", mesh->Nsgeo);

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
#endif
  
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);
  
  kernelInfo.addDefine("p_Lambda2", 0.5f);


  
  mesh->volumeKernel =
    mesh->device.buildKernelFromSource("okl/acousticsVolume2D.okl",
				       "acousticsVolume2D_o0",
				       kernelInfo);

  mesh->surfaceKernel =
    mesh->device.buildKernelFromSource("okl/acousticsSurface2D.okl",
				       "acousticsSurface2D_s0",
				       kernelInfo);

  mesh->partialSurfaceKernel =
    mesh->device.buildKernelFromSource("okl/acousticsPartialSurface2D.okl",
				       "acousticsPartialSurface2D_s0",
				       kernelInfo);

  mesh->updateKernel =
    mesh->device.buildKernelFromSource("okl/acousticsUpdate2D.okl",
				       "acousticsUpdate2D",
				       kernelInfo);
  
  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource("okl/meshHaloExtract2D.okl",
				       "meshHaloExtract2D",
				       kernelInfo);

  mesh->pmlKernel =
    mesh->device.buildKernelFromSource("okl/acousticsPml2D.okl",
				       "acousticsPml2D",
				       kernelInfo);

  mesh->pmlUpdateKernel =
    mesh->device.buildKernelFromSource("okl/acousticsPmlUpdate2D.okl",
				       "acousticsPmlUpdate2D",
				       kernelInfo);

  
  int Ntests = 10;

#define maxNkernels 100

  int NvolumeKernels = 7;
  occa::kernel *acousticsVolumeKernels = new occa::kernel[maxNkernels];
  char kernelNames[maxNkernels][BUFSIZ];
  double bestElapsed = 1e9;
  
  for(int ker=0;ker<NvolumeKernels;++ker){
    sprintf(kernelNames[ker], "acousticsVolume2D_o%d", ker);
    
    acousticsVolumeKernels[ker] =
      mesh->device.buildKernelFromSource("okl/acousticsVolume2D.okl", kernelNames[ker], kernelInfo);
    
    mesh->device.finish();
    occa::tic(kernelNames[ker]);
    for(int test=0;test<Ntests;++test)
      acousticsVolumeKernels[ker](mesh->Nelements,
				      mesh->o_vgeo,
				      mesh->o_DrT,
				      mesh->o_DsT,
				      mesh->o_q,
				      mesh->o_rhsq);
    mesh->device.finish();
    double elapsed = occa::toc(kernelNames[ker]);
    if(elapsed<bestElapsed){
      mesh->volumeKernel = acousticsVolumeKernels[ker];
      printf("promoting kernel: %d (time %g)\n", ker, elapsed);
      bestElapsed = elapsed;
    }
    else{
      printf("not promoting kernel: %d (time %g)\n", ker, elapsed);
    }
  }

  int NsurfaceKernels = 6;
  char surfaceKernelNames[maxNkernels][BUFSIZ];
  occa::kernel *acousticsSurfaceKernels = new occa::kernel[maxNkernels];
  bestElapsed = 1e9;
  
  for(int ker=0;ker<NsurfaceKernels;++ker){
    sprintf(surfaceKernelNames[ker], "acousticsSurface2D_s%d", ker);
    
    acousticsSurfaceKernels[ker] =
      mesh->device.buildKernelFromSource("okl/acousticsSurface2D.okl",
					 surfaceKernelNames[ker], kernelInfo);
    
    mesh->device.finish();
    occa::tic(surfaceKernelNames[ker]);
    for(int test=0;test<Ntests;++test){
      dfloat t  = 0;
      acousticsSurfaceKernels[ker](mesh->Nelements,
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
    mesh->device.finish();
    double elapsed = occa::toc(surfaceKernelNames[ker]);
    if(elapsed<bestElapsed){
      mesh->surfaceKernel = acousticsSurfaceKernels[ker];
      printf("promoting kernel: %d (time %g)\n", ker, elapsed);
      bestElapsed = elapsed;
    }
    else{
      printf("not promoting kernel: %d (time %g)\n", ker, elapsed);
    }
  }

  int NpartialSurfaceKernels = 6;
  char partialSurfaceKernelNames[maxNkernels][BUFSIZ];
  occa::kernel *acousticsPartialSurfaceKernels = new occa::kernel[maxNkernels];
  bestElapsed = 1e9;
  
  for(int ker=0;ker<NpartialSurfaceKernels;++ker){
    sprintf(partialSurfaceKernelNames[ker], "acousticsPartialSurface2D_s%d", ker);
    
    acousticsPartialSurfaceKernels[ker] =
      mesh->device.buildKernelFromSource("okl/acousticsPartialSurface2D.okl",
					 partialSurfaceKernelNames[ker], kernelInfo);
    
    mesh->device.finish();
    occa::tic(partialSurfaceKernelNames[ker]);
    for(int test=0;test<Ntests;++test){
      dfloat t  = 0;
      acousticsPartialSurfaceKernels[ker](mesh->NinternalElements,
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
    }
    mesh->device.finish();
    double elapsed = occa::toc(partialSurfaceKernelNames[ker]);
    if(elapsed<bestElapsed){
      mesh->partialSurfaceKernel = acousticsPartialSurfaceKernels[ker];
      printf("promoting kernel: %d (time %g)\n", ker, elapsed);
      bestElapsed = elapsed;
    }
    else{
      printf("not promoting kernel: %d (time %g)\n", ker, elapsed);
    }
  }

  // set up perfectly matched layer
  dfloat xmin = -.8, xmax = .8, ymin = -.8, ymax = .8;
  dfloat xsigma = 40, ysigma = 40;
  acousticsPmlSetup2D(mesh, xmin, xmax, ymin, ymax, xsigma, ysigma);
 

}

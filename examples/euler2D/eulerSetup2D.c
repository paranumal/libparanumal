#include "mpi.h"
#include <math.h>
#include "euler2D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS
// #undef USE_2_STREAMS

void eulerSetup2D(mesh2D *mesh){
  
  mesh->Nfields = 4;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  // set temperature, gas constant, wave speeds
  mesh->RT = 1.;
  mesh->sqrtRT = 1.; // sqrt(mesh->RT);  
  
  // initial conditions
  // uniform flow
  dfloat rho = 1, u = 0, v = 0; //u = 1.f/sqrt(2.f), v = 1.f/sqrt(2.f); 
  dfloat Rbar = rho;
  dfloat Ubar = rho*u;
  dfloat Vbar = rho*v;

  printf("%17.15lf %17.15lf %17.15lf\n", Rbar, Ubar, Vbar);
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];

      mesh->q[cnt+0] = Rbar*(1.+.1*exp(-20*(x*x+y*y))); // uniform density, zero flow
      mesh->q[cnt+1] = Ubar;
      mesh->q[cnt+2] = Vbar;
    
      cnt += mesh->Nfields;

    }
  }

  // set time step
  dfloat hmin = 1e9, hmax = 0;
  for(iint e=0;e<mesh->Nelements;++e){  

    for(iint f=0;f<mesh->Nfaces;++f){
      iint sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
      // h = 0.5/(sJ/J)
      
      dfloat hest = .5/(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }
  }

  dfloat cfl = .4; // depends on the stability region size (was .4)

  // dt ~ cfl (h/(N+1)^2)/(Lambda^2*fastest wave speed)
  dfloat dt = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*(sqrt(Ubar*Ubar+Vbar*Vbar)/Rbar+mesh->sqrtRT));

  printf("hmin = %g\n", hmin);
  printf("hmax = %g\n", hmax);
  printf("cfl = %g\n", cfl);
  printf("dt = %g\n", dt);
  printf("max wave speed = %g\n", mesh->sqrtRT);
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //
  mesh->finalTime = 100;
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
  //  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank+1)%3);
  //  sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //  sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  sprintf(deviceConfig, "mode = Serial");
  mesh->device.setup(deviceConfig);

  // build Dr, Ds, LIFT transposes
  dfloat *DrT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *DsT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Np;++m){
      DrT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
      DsT[n+m*mesh->Np] = mesh->Ds[n*mesh->Np+m];
    }
  }

  dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      LIFTT[n+m*mesh->Np] = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
    }
  }

  // build volume cubature matrix transposes
  dfloat *cubDrWT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  dfloat *cubDsWT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  dfloat *cubInterpT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->cubNp;++m){
      cubDrWT[n+m*mesh->Np] = mesh->cubDrW[n*mesh->cubNp+m];
      cubDsWT[n+m*mesh->Np] = mesh->cubDsW[n*mesh->cubNp+m];
      cubInterpT[m+n*mesh->cubNp] = mesh->cubInterp[m*mesh->Np+n];
      printf("%g @ ", cubInterpT[m+n*mesh->cubNp]);
    }
  }

  // build surface integration matrix transposes
  dfloat *intLIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  dfloat *intInterpT = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Nfaces*mesh->intNfp;++m){
      intLIFTT[n+m*mesh->Np] = mesh->intLIFT[n*mesh->intNfp*mesh->Nfaces+m];
    }
  }
  for(int n=0;n<mesh->intNfp*mesh->Nfaces;++n){
    for(int m=0;m<mesh->Nfp;++m){
      intInterpT[n+m*mesh->Nfaces*mesh->intNfp] = mesh->intInterp[n*mesh->Nfp + m];
    }
  }

  mesh->intx = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  mesh->inty = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      for(iint n=0;n<mesh->intNfp;++n){
	dfloat ix = 0, iy = 0;
	for(iint m=0;m<mesh->Nfp;++m){
	  iint vid = mesh->vmapM[m+f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces];
	  dfloat xm = mesh->x[vid];
	  dfloat ym = mesh->y[vid];
	  dfloat Inm = mesh->intInterp[n+f*mesh->intNfp+m*mesh->intNfp*mesh->Nfaces];
	  ix += Inm*xm;
	  iy += Inm*ym;
	}
	iint id = n + f*mesh->intNfp + e*mesh->Nfaces*mesh->intNfp;
	mesh->intx[id] = ix;
	mesh->inty[id] = iy;
      }
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


  mesh->o_intx =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			mesh->intx);

  mesh->o_inty =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			mesh->inty);

  
  if(mesh->totalHaloPairs>0){
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      mesh->device.malloc(mesh->totalHaloPairs*sizeof(iint), mesh->haloElementList);
    
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));
  }

  mesh->o_cubInterpT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubInterpT);

  mesh->o_cubDrWT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubDrWT);
  
  mesh->o_cubDsWT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubDsWT);

  mesh->o_intInterpT =
    mesh->device.malloc(mesh->Nfp*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			intInterpT);

  mesh->o_intLIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			intLIFTT);
  
  
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
  
  occa::kernelInfo kernelInfo;

  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_N", mesh->N);
  kernelInfo.addDefine("p_Np", mesh->Np);
  kernelInfo.addDefine("p_Nfp", mesh->Nfp);
  kernelInfo.addDefine("p_Nfaces", mesh->Nfaces);
  kernelInfo.addDefine("p_NfacesNfp", mesh->Nfp*mesh->Nfaces);
  kernelInfo.addDefine("p_Nvgeo", mesh->Nvgeo);
  kernelInfo.addDefine("p_Nsgeo", mesh->Nsgeo);

  printf("cubNp = %d, intNfp = %d\n", mesh->cubNp, mesh->intNfp);
  
  kernelInfo.addDefine("p_cubNp", mesh->cubNp);
  kernelInfo.addDefine("p_intNfp", mesh->intNfp);
  kernelInfo.addDefine("p_intNfpNfaces", mesh->intNfp*mesh->Nfaces);

  int maxVolumeNodes = mymax(mesh->Np, mesh->cubNp);
  kernelInfo.addDefine("p_maxVolumeNodes", maxVolumeNodes);

  int maxSurfaceNodes = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo.addDefine("p_maxSurfaceNodes", maxSurfaceNodes);
  printf("maxSurfaceNodes=%d\n", maxSurfaceNodes);
  
  int NblockV = 256/maxVolumeNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 256/maxSurfaceNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  // physics 
  kernelInfo.addDefine("p_RT", mesh->RT);
  kernelInfo.addDefine("p_sqrtRT", mesh->sqrtRT);
  kernelInfo.addDefine("p_isqrtRT", 1.f/mesh->sqrtRT);

  kernelInfo.addDefine("p_Rbar", Rbar);
  kernelInfo.addDefine("p_Ubar", Ubar);
  kernelInfo.addDefine("p_Vbar", Vbar);

  if(sizeof(dfloat)==4){
    kernelInfo.addDefine("dfloat","float");
    kernelInfo.addDefine("dfloat4","float4");
    kernelInfo.addDefine("dfloat8","float8");
  }
  if(sizeof(dfloat)==8){
    kernelInfo.addDefine("dfloat","double");
    kernelInfo.addDefine("dfloat4","double4");
    kernelInfo.addDefine("dfloat8","double8");
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

  mesh->volumeKernel =
    mesh->device.buildKernelFromSource("okl/eulerVolume2D.okl",
				       "eulerVolume2D",
				       kernelInfo);
  mesh->surfaceKernel =
    mesh->device.buildKernelFromSource("okl/eulerSurface2D.okl",
				       "eulerSurface2D",
				       kernelInfo);
  mesh->updateKernel =
    mesh->device.buildKernelFromSource("okl/eulerUpdate2D.okl",
				       "eulerUpdate2D",
				       kernelInfo);

  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource("okl/meshHaloExtract2D.okl",
				       "meshHaloExtract2D",
				       kernelInfo);
  
}

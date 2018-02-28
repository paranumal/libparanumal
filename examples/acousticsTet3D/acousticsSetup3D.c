#include "acoustics3D.h"

void acousticsSetup3D(mesh3D *mesh){

  mesh->Nfields = 4;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  // fix this later (initial conditions)
  int cnt = 0;
  dfloat time = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];
      
      acousticsCavitySolution3D(x, y, z, time,
				mesh->q+cnt,
				mesh->q+cnt+1,
				mesh->q+cnt+2,
				mesh->q+cnt+3);

      cnt += mesh->Nfields;
    }
    printf("\n");
  }

  // set penalty parameter
  mesh->Lambda2 = 0.5;
  
  // set time step
  dfloat hmin = 1e9;
  for(int e=0;e<mesh->Nelements;++e){  

    for(int f=0;f<mesh->Nfaces;++f){
      int sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
      // h = 0.5/(sJ/J)
      
      dfloat hest = .5/(sJ*invJ);

      hmin = mymin(hmin, hest);
    }
  }
  dfloat cfl = 1; // depends on the stability region size

  dfloat dt = cfl*hmin/((mesh->N+1)*(mesh->N+1)*mesh->Lambda2);

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //
  mesh->finalTime = .4;
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  mesh->dt = mesh->finalTime/mesh->NtimeSteps;

  // errorStep
  mesh->errorStep = 10;

  printf("dt = %g\n", mesh->dt);

  // output mesh
  meshVTU3D(mesh, "foo.vtu");

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", 0);
  //  sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //  sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);

 
  occa::kernelInfo kernelInfo;

  // generic occa device set up
  meshOccaSetup3D(mesh, deviceConfig, kernelInfo);

  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  
  mesh->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsVolume3D.okl",
				       "acousticsVolume3D_o0",
				       kernelInfo);

  mesh->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsSurface3D.okl",
				       "acousticsSurface3D_s0",
				       kernelInfo);

  mesh->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsUpdate3D.okl",
				       "acousticsUpdate3D",
				       kernelInfo);

  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract3D.okl",
				       "meshHaloExtract3D",
				       kernelInfo);

  int Ntests = 10;

#define maxNkernels 100

  int NvolumeKernels = 7;
  occa::kernel *acousticsVolumeKernels = new occa::kernel[maxNkernels];
  char kernelNames[maxNkernels][BUFSIZ];
  double bestElapsed = 1e9;
  
  for(int ker=0;ker<NvolumeKernels;++ker){
    sprintf(kernelNames[ker], "acousticsVolume3D_o%d", ker);
    
    acousticsVolumeKernels[ker] =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsVolume3D.okl", kernelNames[ker], kernelInfo);
    
    mesh->device.finish();
    occa::tic(kernelNames[ker]);
    for(int test=0;test<Ntests;++test)
      acousticsVolumeKernels[ker](mesh->Nelements,
				      mesh->o_vgeo,
				      mesh->o_DrT,
				      mesh->o_DsT,
				      mesh->o_DtT,
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

  int NsurfaceKernels = 5;
  char surfaceKernelNames[maxNkernels][BUFSIZ];
  occa::kernel *acousticsSurfaceKernels = new occa::kernel[maxNkernels];
  bestElapsed = 1e9;
  
  for(int ker=0;ker<NsurfaceKernels;++ker){
    sprintf(surfaceKernelNames[ker], "acousticsSurface3D_s%d", ker);
    
    acousticsSurfaceKernels[ker] =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsSurface3D.okl",
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
				       mesh->o_z,
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


}

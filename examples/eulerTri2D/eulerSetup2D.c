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
  mesh->finalTime = 10;
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
  //sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank+1)%3);
  sprintf(deviceConfig, "mode = OpenCL, deviceID = 1, platformID = 0");
  //  sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //  sprintf(deviceConfig, "mode = Serial");

  occa::kernelInfo kernelInfo;

  // generic set up
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);
  
  // physics 
  kernelInfo.addDefine("p_RT", mesh->RT);
  kernelInfo.addDefine("p_sqrtRT", mesh->sqrtRT);
  kernelInfo.addDefine("p_isqrtRT", 1.f/mesh->sqrtRT);

  kernelInfo.addDefine("p_Rbar", Rbar);
  kernelInfo.addDefine("p_Ubar", Ubar);
  kernelInfo.addDefine("p_Vbar", Vbar);

  // euler kernel specific quantities
  int maxVolumeNodes = mymax(mesh->Np, mesh->cubNp);
  kernelInfo.addDefine("p_maxVolumeNodes", maxVolumeNodes);

  int maxSurfaceNodes = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo.addDefine("p_maxSurfaceNodes", maxSurfaceNodes);
  printf("maxSurfaceNodes=%d\n", maxSurfaceNodes);
  
  int NblockV = 128/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 128/maxSurfaceNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  // build euler specific kernels
  mesh->volumeKernel =
    mesh->device.buildKernelFromSource("okl/eulerVolume2D.okl",
				       "eulerVolume2D_c1",
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

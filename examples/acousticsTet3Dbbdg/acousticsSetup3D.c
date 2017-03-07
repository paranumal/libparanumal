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
  iint cnt = 0;
  dfloat time = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];
      
      cnt = e*mesh->Np*mesh->Nfields + n*mesh->Nfields;
      acousticsCavitySolution3D(x, y, z, time,
				mesh->q+cnt,
				mesh->q+cnt+1,
				mesh->q+cnt+2,
				mesh->q+cnt+3);
    }
  }

  //Transform to BB modal space
  #if USE_BERN
    dfloat qtmp[mesh->Nfields*mesh->Np];
    for (iint e =0;e<mesh->Nelements;e++){
      cnt = e*mesh->Np*mesh->Nfields;
      
      for (iint n=0; n<mesh->Np; n++){
        qtmp[n*mesh->Nfields + 0] = mesh->q[cnt+n*mesh->Nfields+0];
        qtmp[n*mesh->Nfields + 1] = mesh->q[cnt+n*mesh->Nfields+1];
        qtmp[n*mesh->Nfields + 2] = mesh->q[cnt+n*mesh->Nfields+2];
        qtmp[n*mesh->Nfields + 3] = mesh->q[cnt+n*mesh->Nfields+3];
        mesh->q[cnt+n*mesh->Nfields+0] = 0.0;
        mesh->q[cnt+n*mesh->Nfields+1] = 0.0;
        mesh->q[cnt+n*mesh->Nfields+2] = 0.0;
        mesh->q[cnt+n*mesh->Nfields+3] = 0.0;
      }
      for (iint n=0;n<mesh->Np;n++){
        for (iint m=0; m<mesh->Np; m++){
          mesh->q[cnt+n*mesh->Nfields + 0] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+0];
          mesh->q[cnt+n*mesh->Nfields + 1] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+1];
          mesh->q[cnt+n*mesh->Nfields + 2] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+2];
          mesh->q[cnt+n*mesh->Nfields + 3] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+3];
        }
      }
    }
  #endif

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
  dfloat cfl = 0.5; // depends on the stability region size

  dfloat dt = cfl*hmin/((mesh->N+1)*(mesh->N+1)*mesh->Lambda2);

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //
  mesh->finalTime = .1;
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  mesh->dt = mesh->finalTime/mesh->NtimeSteps;

  // errorStep
  mesh->errorStep = 100;

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
  //  sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 0);

  mesh->device.setup(deviceConfig);

  occa::kernelInfo kernelInfo;

  // generic occa device set up
  meshOccaSetup3D(mesh, deviceConfig, kernelInfo);


  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);
  
  kernelInfo.addDefine("p_Lambda2", 0.5f);

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

  mesh->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgVolume3D.okl",
				       "acousticsVolume3Dbbdg",
				       kernelInfo);

  mesh->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgSurface3D.okl",
				       "acousticsSurface3Dbbdg",
				       kernelInfo);
  
  mesh->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsUpdate3D.okl",
				       "acousticsUpdate3D",
				       kernelInfo);

  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract3D.okl",
				       "meshHaloExtract3D",
				       kernelInfo);

}

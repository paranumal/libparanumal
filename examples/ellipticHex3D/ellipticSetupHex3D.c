#include "ellipticHex3D.h"

void ellipticSetupHex3D(mesh3D *mesh){

  mesh->Nfields = 1;
  
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
      
      cnt += mesh->Nfields;
    }
  }

  // set penalty parameter
  mesh->Lambda2 = 0.5;

  printf("Nelements = %d, Nfaces = %d\n", mesh->Nelements, mesh->Nfaces);
  
  // set time step
  dfloat hmin = 1e9;
  for(iint e=0;e<mesh->Nelements;++e){  

    for(iint n=0;n<mesh->Nfaces*mesh->Nfp;++n){
      iint sid = mesh->Nsgeo*(mesh->Nfp*mesh->Nfaces*e +  n);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
      // h = 0.5/(sJ/J)
      
      dfloat hest = .5/(sJ*invJ);
      
      hmin = mymin(hmin, hest);
    }
  }

  printf("hmin = %g\n", hmin);
  
  dfloat cfl = 1; // depends on the stability region size

  dfloat dt = cfl*hmin/((mesh->N+1)*(mesh->N+1)*mesh->Lambda2);

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //
  mesh->finalTime = .4;
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
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);

  occa::kernelInfo kernelInfo;

  void meshOccaSetup3D(mesh3D *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo);
  meshOccaSetup3D(mesh, deviceConfig, kernelInfo);
  
  kernelInfo.addDefine("p_Lambda2", 0.5f);

  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource("okl/meshHaloExtract3D.okl",
				       "meshHaloExtract3D",
				       kernelInfo);

  mesh->gatherKernel =
    mesh->device.buildKernelFromSource("okl/gather.okl",
				       "gather",
				       kernelInfo);

  mesh->scatterKernel =
    mesh->device.buildKernelFromSource("okl/scatter.okl",
				       "scatter",
				       kernelInfo);

  mesh->getKernel =
    mesh->device.buildKernelFromSource("okl/get.okl",
				       "get",
				       kernelInfo);

  mesh->putKernel =
    mesh->device.buildKernelFromSource("okl/put.okl",
				       "put",
				       kernelInfo);


  mesh->AxKernel =
    mesh->device.buildKernelFromSource("okl/ellipticAxHex3D.okl",
				       "ellipticAxHex3D",
				       kernelInfo);

  mesh->weightedInnerProduct1Kernel =
    mesh->device.buildKernelFromSource("okl/weightedInnerProduct1.okl",
				       "weightedInnerProduct1",
				       kernelInfo);

  mesh->weightedInnerProduct2Kernel =
    mesh->device.buildKernelFromSource("okl/weightedInnerProduct2.okl",
				       "weightedInnerProduct2",
				       kernelInfo);

  mesh->scaledAddKernel =
      mesh->device.buildKernelFromSource("okl/scaledAdd.okl",
					 "scaledAdd",
					 kernelInfo);

  mesh->dotMultiplyKernel =
      mesh->device.buildKernelFromSource("okl/dotMultiply.okl",
					 "dotMultiply",
					 kernelInfo);

  mesh->dotDivideKernel =
      mesh->device.buildKernelFromSource("okl/dotDivide.okl",
					 "dotDivide",
					 kernelInfo);

  
  // set up gslib MPI gather-scatter and OCCA gather/scatter arrays
  meshParallelGatherScatterSetup3D(mesh,
				   mesh->Np*mesh->Nelements,
				   sizeof(dfloat),
				   mesh->gatherLocalIds,
				   mesh->gatherBaseIds, 
				   mesh->gatherBaseRanks,
				   mesh->gatherMaxRanks,
				   mesh->NuniqueBases,
				   mesh->o_gatherNodeOffsets,
				   mesh->o_gatherLocalNodes,
				   mesh->o_gatherTmp,
				   mesh->NnodeHalo,
				   mesh->o_nodeHaloIds,
				   mesh->o_subGatherTmp,
				   (void**)&(mesh->subGatherTmp),
				   (void**)&(mesh->gsh));


  // find maximum degree
  {
    for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
      mesh->rhsq[n] = 1;
    }
    mesh->o_rhsq.copyFrom(mesh->rhsq);
    
    void meshParallelGatherScatter3D(mesh3D *mesh, occa::memory &o_v, occa::memory &o_gsv, const char *type);
    meshParallelGatherScatter3D(mesh, mesh->o_rhsq, mesh->o_rhsq, dfloatString);

    mesh->o_rhsq.copyTo(mesh->rhsq);
    
    dfloat maxDegree = 0, minDegree = 1e9;
    dfloat gatherMaxDegree = 0, gatherMinDegree = 1e9;
    for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
      maxDegree = mymax(maxDegree, mesh->rhsq[n]);
      minDegree = mymin(minDegree, mesh->rhsq[n]);
    }

    MPI_Allreduce(&maxDegree, &gatherMaxDegree, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minDegree, &gatherMinDegree, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

    if(rank==0){
      printf("max degree = " dfloatFormat "\n", gatherMaxDegree);
      printf("min degree = " dfloatFormat "\n", gatherMinDegree);
    }
  }

  // build weights for continuous SEM L2 project --->
  iint Ntotal = mesh->Nelements*mesh->Np;
  dfloat *localMM = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat wJ = mesh->ggeo[e*mesh->Np*mesh->Nggeo + n + GWJID*mesh->Np];
      localMM[n+e*mesh->Np] = wJ;
    }
  }

  occa::memory o_localMM = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
  occa::memory o_MM      = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);

  // sum up all contributions at base nodes and scatter back
  meshParallelGatherScatter3D(mesh, o_localMM, o_MM, dfloatString);

  mesh->o_projectL2 = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
  mesh->dotDivideKernel(Ntotal, o_localMM, o_MM, mesh->o_projectL2);

  free(localMM); o_MM.free(); o_localMM.free();
  // <------
  

  
}

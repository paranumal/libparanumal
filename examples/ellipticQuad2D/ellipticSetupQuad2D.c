#include "ellipticQuad2D.h"

void ellipticSetupQuad2D(mesh2D *mesh, ogs_t **ogs, precon_t **precon, dfloat lambda, const char *options){

  mesh->Nfields = 1;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", 0);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //  sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);

  occa::kernelInfo kernelInfo;

  void meshOccaSetup2D(mesh2D *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo);
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);
  
  kernelInfo.addDefine("p_NqP", (mesh->Nq+2));
  kernelInfo.addDefine("p_NpP", (mesh->NqP*mesh->NqP));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);

  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax); 

  int NblockV = 512/mesh->Np; // get close to 256 threads
  kernelInfo.addDefine("p_NblockV", NblockV);
  
  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
				       "meshHaloExtract2D",
				       kernelInfo);

  mesh->gatherKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/gather.okl",
				       "gather",
				       kernelInfo);

  mesh->scatterKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/scatter.okl",
				       "scatter",
				       kernelInfo);

  mesh->getKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/get.okl",
				       "get",
				       kernelInfo);

  mesh->putKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/put.okl",
				       "put",
				       kernelInfo);


  mesh->AxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxQuad2D.okl",
				       "ellipticAxQuad2D_e0",
				       kernelInfo);

  mesh->weightedInnerProduct1Kernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/weightedInnerProduct1.okl",
				       "weightedInnerProduct1",
				       kernelInfo);

  mesh->weightedInnerProduct2Kernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/weightedInnerProduct2.okl",
				       "weightedInnerProduct2",
				       kernelInfo);

  mesh->innerProductKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/innerProduct.okl",
				       "innerProduct",
				       kernelInfo);
  
  mesh->scaledAddKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/scaledAdd.okl",
					 "scaledAdd",
					 kernelInfo);

  mesh->dotMultiplyKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/dotMultiply.okl",
					 "dotMultiply",
					 kernelInfo);

  mesh->dotDivideKernel = 
      mesh->device.buildKernelFromSource(DHOLMES "/okl/dotDivide.okl",
					 "dotDivide",
					 kernelInfo);


  mesh->gradientKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientQuad2D.okl",
				       "ellipticGradientQuad2D",
					 kernelInfo);


  mesh->ipdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgQuad2D.okl",
				       "ellipticAxIpdgQuad2D",
				       kernelInfo);  


  // set up gslib MPI gather-scatter and OCCA gather/scatter arrays
  *ogs = meshParallelGatherScatterSetup(mesh,
					mesh->Np*mesh->Nelements,
					sizeof(dfloat),
					mesh->gatherLocalIds,
					mesh->gatherBaseIds, 
					mesh->gatherHaloFlags);

  *precon = ellipticPreconditionerSetupQuad2D(mesh, *ogs, lambda, options);

  (*precon)->preconKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticOasPreconQuad2D.okl",
				       "ellipticOasPreconQuad2D",
				       kernelInfo);
  
  (*precon)->restrictKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconRestrictQuad2D.okl",
				       "ellipticFooQuad2D",
				       kernelInfo);

  (*precon)->coarsenKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconCoarsen.okl",
				       "ellipticPreconCoarsen",
				       kernelInfo);

  (*precon)->prolongateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconProlongate.okl",
				       "ellipticPreconProlongate",
				       kernelInfo);


  // find maximum degree
  {
    for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
      mesh->rhsq[n] = 1;
    }
    mesh->o_rhsq.copyFrom(mesh->rhsq);
    
    ellipticParallelGatherScatter2D(mesh, *ogs, mesh->o_rhsq, mesh->o_rhsq, dfloatString, "add");

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
  ellipticParallelGatherScatter2D(mesh, *ogs, o_localMM, o_MM, dfloatString, "add");

  mesh->o_projectL2 = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
  mesh->dotDivideKernel(Ntotal, o_localMM, o_MM, mesh->o_projectL2);

  free(localMM); o_MM.free(); o_localMM.free();
  // <------

  // coarse grid preconditioner (only continous elements)
  // do this here since we need A*x
  ellipticCoarsePreconditionerSetupQuad2D(mesh, *precon, *ogs, lambda, options);


  
}

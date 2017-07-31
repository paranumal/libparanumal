#include "ellipticTet3D.h"

void ellipticComputeDegreeVector(mesh3D *mesh, iint Ntotal, ogs_t *ogs, dfloat *deg){

  // build degree vector
  for(iint n=0;n<Ntotal;++n)
    deg[n] = 1;

  occa::memory o_deg = mesh->device.malloc(Ntotal*sizeof(dfloat), deg);
  
  o_deg.copyFrom(deg);
  
  ellipticParallelGatherScatterTet3D(mesh, ogs, o_deg, o_deg, dfloatString, "add");
  
  o_deg.copyTo(deg);

  mesh->device.finish();
  o_deg.free();
  
}

solver_t *ellipticSolveSetupTet3D(mesh_t *mesh, dfloat tau, dfloat lambda, iint *BCType, 
                          occa::kernelInfo &kernelInfo, const char *options){

  iint rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint Ntotal = mesh->Np*mesh->Nelements;
  iint NtotalP = mesh->NpP*mesh->Nelements;
  iint Nblock = (Ntotal+blockSize-1)/blockSize;
  iint Nhalo = mesh->Np*mesh->totalHaloPairs;
  iint Nall   = Ntotal + Nhalo;
  iint NallP  = NtotalP;

  int NblockV = mymax(1,1024/mesh->Np); // works for CUDA
  int NblockS = mymax(1,1024/(mesh->Nfp*mesh->Nfaces)); // works for CUDA
  
  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));

  solver->tau = tau;

  solver->mesh = mesh;

  solver->p   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->z   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->zP  = (dfloat*) calloc(NallP,  sizeof(dfloat));
  solver->Ax  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->Ap  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));

  solver->grad = (dfloat*) calloc(Nall*4, sizeof(dfloat));

  solver->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
  solver->o_zP  = mesh->device.malloc(NallP*sizeof(dfloat),solver->zP); // CAUTION
  solver->o_Ax  = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), solver->Ap);
  solver->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), solver->tmp);

  solver->o_grad  = mesh->device.malloc(Nall*4*sizeof(dfloat), solver->grad);
  
  //setup async halo stream
  solver->defaultStream = mesh->device.getStream();
  solver->dataStream = mesh->device.createStream();
  mesh->device.setStream(solver->defaultStream);
  
  iint Nbytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  if(Nbytes>0){
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
    
    solver->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    solver->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();
  }else{
    solver->sendBuffer = NULL;
    solver->recvBuffer = NULL;
  }
  mesh->device.setStream(solver->defaultStream);

  solver->type = strdup(dfloatString);

  solver->Nblock = Nblock;

  //fill geometric factors in halo
  if(mesh->totalHaloPairs){
    iint Nlocal = mesh->Np*mesh->Nelements;
    iint Nhalo  = mesh->Np*mesh->totalHaloPairs;
    dfloat *vgeoSendBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Nvgeo, sizeof(dfloat));
    
    // import geometric factors from halo elements
    mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (mesh->Nelements+mesh->totalHaloPairs)*mesh->Nvgeo*sizeof(dfloat));
    
    meshHaloExchange(mesh,
         mesh->Nvgeo*sizeof(dfloat),
         mesh->vgeo,
         vgeoSendBuffer,
         mesh->vgeo + mesh->Nelements*mesh->Nvgeo);
    
    mesh->o_vgeo =
      mesh->device.malloc((mesh->Nelements + mesh->totalHaloPairs)*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);
  }

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");
  //kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");

  kernelInfo.addDefine("p_blockSize", blockSize);
  kernelInfo.addDefine("p_NblockV", NblockV);
  kernelInfo.addDefine("p_NblockS", NblockS);
  
  // add custom defines
  kernelInfo.addDefine("p_NpP", (mesh->Np+mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);

  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax); 

  //int NblockV = 256/mesh->Np; // get close to 256 threads
  //kernelInfo.addDefine("p_NblockV", NblockV);
  
  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract3D.okl",
				       "meshHaloExtract3D",
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
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTet3D.okl",
               "ellipticAxTet3D",
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
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientTet3D.okl", 
				       "ellipticGradientTet3D",
					 kernelInfo);


  mesh->ipdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgTet3D.okl", 
				       "ellipticAxIpdgTet3D",
				       kernelInfo);  

  solver->rhsBCIpdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticRhsBCIpdgTet3D.okl",
               "ellipticRhsBCIpdgTet3D",
               kernelInfo);
  
  // set up gslib MPI gather-scatter and OCCA gather/scatter arrays
  occaTimerTic(mesh->device,"GatherScatterSetup");
  solver->ogs = meshParallelGatherScatterSetup(mesh,
					       mesh->Np*mesh->Nelements,
					       sizeof(dfloat),
					       mesh->gatherLocalIds,
					       mesh->gatherBaseIds, 
					       mesh->gatherHaloFlags);
  occaTimerToc(mesh->device,"GatherScatterSetup"); 
  
  occaTimerTic(mesh->device,"PreconditionerSetup");
  solver->precon = ellipticPreconditionerSetupTet3D(mesh, solver->ogs, tau, lambda, BCType, options);
  
  solver->precon->preconKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticOasPreconTet3D.okl", 
				       "ellipticOasPreconTet3D",
				       kernelInfo);
  
  solver->precon->restrictKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconRestrictTet3D.okl", 
				       "ellipticFooTet3D",
				       kernelInfo);

  solver->precon->coarsenKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconCoarsen.okl",
				       "ellipticPreconCoarsen",
				       kernelInfo);

  solver->precon->prolongateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconProlongate.okl",
				       "ellipticPreconProlongate",
				       kernelInfo);

  solver->precon->blockJacobiKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBlockJacobiPreconTet3D.okl",
               "ellipticBlockJacobiPreconTet3D",
               kernelInfo);
  occaTimerToc(mesh->device,"PreconditionerSetup");

  occaTimerTic(mesh->device,"DegreeVectorSetup");
  dfloat *invDegree = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *degree = (dfloat*) calloc(Ntotal, sizeof(dfloat));

  solver->o_invDegree = mesh->device.malloc(Ntotal*sizeof(dfloat), invDegree);
  
  ellipticComputeDegreeVector(mesh, Ntotal, solver->ogs, degree);

  for(iint n=0;n<Ntotal;++n){ // need to weight inner products{
    if(degree[n] == 0) printf("WARNING!!!!\n");
    invDegree[n] = 1./degree[n];
  }
  
  solver->o_invDegree.copyFrom(invDegree);
  occaTimerToc(mesh->device,"DegreeVectorSetup");


  if (strstr(options,"MATRIXFREE")) { 
    //set matrix free A in parAlmond
    void **args = (void **) calloc(2,sizeof(void *));
    dfloat *vlambda = (dfloat *) calloc(1,sizeof(dfloat));
    
    *vlambda = lambda;

    args[0] = (void *) solver;
    args[1] = (void *) vlambda;

    parAlmondSetMatFreeAX(solver->precon->parAlmond,ellipticMatrixFreeAx,args);
  }

  return solver;
}

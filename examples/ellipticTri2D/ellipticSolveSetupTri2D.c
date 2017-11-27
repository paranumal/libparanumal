#include "ellipticTri2D.h"

void ellipticComputeDegreeVector(mesh2D *mesh, iint Ntotal, ogs_t *ogs, dfloat *deg){

  // build degree vector
  for(iint n=0;n<Ntotal;++n)
    deg[n] = 1;

  occa::memory o_deg = mesh->device.malloc(Ntotal*sizeof(dfloat), deg);

  o_deg.copyFrom(deg);

  ellipticParallelGatherScatterTri2D(mesh, ogs, o_deg, o_deg, dfloatString, "add");

  o_deg.copyTo(deg);

  mesh->device.finish();
  o_deg.free();

}

solver_t *ellipticSolveSetupTri2D(mesh_t *mesh, dfloat tau, dfloat lambda, iint*BCType,
                      occa::kernelInfo &kernelInfo, const char *options, const char *parAlmondOptions){

  iint Ntotal = mesh->Np*mesh->Nelements;
  iint Nblock = (Ntotal+blockSize-1)/blockSize;
  iint Nhalo = mesh->Np*mesh->totalHaloPairs;
  iint Nall   = Ntotal + Nhalo;

  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));

  solver->tau = tau;

  solver->mesh = mesh;

  solver->p   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->z   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->Ax  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->Ap  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));

  solver->grad = (dfloat*) calloc(Nall*4, sizeof(dfloat));

  solver->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);

  solver->o_res = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
  solver->o_Sres = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
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
    iint Nlocal = mesh->Nelements*mesh->Np;
    iint Nhalo  = mesh->totalHaloPairs*mesh->Np;
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

  //check all the bounaries for a Dirichlet
  bool allNeumann = (lambda==0) ? true :false;
  solver->allNeumannPenalty = 1;
  iint totalElements = 0;
  MPI_Allreduce(&(mesh->Nelements), &totalElements, 1, MPI_IINT, MPI_SUM, MPI_COMM_WORLD);
  solver->allNeumannScale = 1.0/sqrt(mesh->Np*totalElements);

  solver->EToB = (int *) calloc(mesh->Nelements*mesh->Nfaces,sizeof(int));
  for (iint e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if (bc>0) {
        int BC = BCType[bc]; //get the type of boundary
        solver->EToB[e*mesh->Nfaces+f] = BC; //record it
        if (BC!=2) allNeumann = false; //check if its a Dirchlet
      }
    }
  }
  MPI_Allreduce(&allNeumann, &(solver->allNeumann), 1, MPI::BOOL, MPI_LAND, MPI_COMM_WORLD);
  printf("allNeumann = %d \n", solver->allNeumann);

  solver->o_EToB = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int), solver->EToB);

  //add standard boundary functions
  char *boundaryHeaderFileName;
  boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTri2D/ellipticBoundary2D.h");
  kernelInfo.addInclude(boundaryHeaderFileName);

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");
  }

  kernelInfo.addDefine("p_blockSize", blockSize);

  // add custom defines
  kernelInfo.addDefine("p_NpP", (mesh->Np+mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);

  //sizes for the coarsen and prolongation kernels. degree N to degree 1
  kernelInfo.addDefine("p_NpFine", mesh->Np);
  kernelInfo.addDefine("p_NpCoarse", mesh->Nverts);

  kernelInfo.addDefine("p_NpFEM", mesh->NpFEM);

  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax);

  int NblockV = 256/mesh->Np; // get close to 256 threads
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockP = 512/(4*mesh->Np); // get close to 256 threads
  kernelInfo.addDefine("p_NblockP", NblockP);

  int NblockG;
  if(mesh->Np<=32) NblockG = ( 32/mesh->Np );
  else NblockG = 256/mesh->Np;
  kernelInfo.addDefine("p_NblockG", NblockG);

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

  mesh->gatherScatterKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/gatherScatter.okl",
               "gatherScatter",
               kernelInfo);

  mesh->getKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/get.okl",
				       "get",
				       kernelInfo);

  mesh->putKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/put.okl",
				       "put",
				       kernelInfo);

  mesh->sumKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/sum.okl",
               "sum",
               kernelInfo);

  mesh->addScalarKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/addScalar.okl",
               "addScalar",
               kernelInfo);

  solver->AxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTri2D.okl",
               "ellipticAxTri2D",
               kernelInfo);

  solver->partialAxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTri2D.okl",
               "ellipticPartialAxTri2D",
               kernelInfo);

  solver->weightedInnerProduct1Kernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/weightedInnerProduct1.okl",
				       "weightedInnerProduct1",
				       kernelInfo);

  solver->weightedInnerProduct2Kernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/weightedInnerProduct2.okl",
				       "weightedInnerProduct2",
				       kernelInfo);

  solver->innerProductKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/innerProduct.okl",
				       "innerProduct",
				       kernelInfo);

  solver->scaledAddKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/scaledAdd.okl",
					 "scaledAdd",
					 kernelInfo);

  solver->dotMultiplyKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/dotMultiply.okl",
					 "dotMultiply",
					 kernelInfo);

  solver->dotDivideKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/dotDivide.okl",
					 "dotDivide",
					 kernelInfo);

  solver->gradientKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientTri2D.okl",
				       "ellipticGradientTri2D",
					 kernelInfo);

  solver->partialGradientKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientTri2D.okl",
               "ellipticPartialGradientTri2D",
                kernelInfo);

  solver->ipdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgTri2D.okl",
				       "ellipticAxIpdgTri2D",
				       kernelInfo);

  solver->partialIpdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgTri2D.okl",
               "ellipticPartialAxIpdgTri2D",
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

  solver->precon = (precon_t*) calloc(1, sizeof(precon_t));

  #if 0

  solver->precon->overlappingPatchKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticOasPreconTri2D.okl",
               "ellipticOasPreconTri2D",
               kernelInfo);


  #endif

  solver->precon->restrictKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconRestrictTri2D.okl",
               "ellipticFooTri2D",
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
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBlockJacobiPreconTri2D.okl",
				       "ellipticBlockJacobiPreconTri2D",
				       kernelInfo);

  solver->precon->approxPatchSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
               "ellipticApproxPatchSolver2D",
               kernelInfo);

  solver->precon->exactPatchSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
               "ellipticExactPatchSolver2D",
               kernelInfo);

  solver->precon->patchGatherKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchGather.okl",
               "ellipticPatchGather",
               kernelInfo);

  solver->precon->approxFacePatchSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
               "ellipticApproxFacePatchSolver2D",
               kernelInfo);

  solver->precon->exactFacePatchSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
               "ellipticExactFacePatchSolver2D",
               kernelInfo);

  solver->precon->facePatchGatherKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchGather.okl",
               "ellipticFacePatchGather",
               kernelInfo);

  solver->precon->approxBlockJacobiSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
               "ellipticApproxBlockJacobiSolver2D",
               kernelInfo);

  solver->precon->exactBlockJacobiSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
               "ellipticExactBlockJacobiSolver2D",
               kernelInfo);

  solver->precon->SEMFEMInterpKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticSEMFEMInterpTri2D.okl",
               "ellipticSEMFEMInterpTri2D",
               kernelInfo);

  solver->precon->SEMFEMAnterpKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticSEMFEMAnterpTri2D.okl",
               "ellipticSEMFEMAnterpTri2D",
               kernelInfo);

  long long int pre = mesh->device.memoryAllocated();

  occaTimerTic(mesh->device,"PreconditionerSetup");
  ellipticPreconditionerSetupTri2D(solver, solver->ogs, tau, lambda, BCType,  options, parAlmondOptions);
  occaTimerToc(mesh->device,"PreconditionerSetup");

  long long int usedBytes = mesh->device.memoryAllocated()-pre;

  solver->precon->preconBytes = usedBytes;

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



  // set up separate gather scatter infrastructure for halo and non halo nodes
  ellipticParallelGatherScatterSetup(mesh,
                                     mesh->Np*mesh->Nelements,
                                     sizeof(dfloat),
                                     mesh->gatherLocalIds,
                                     mesh->gatherBaseIds,
                                     mesh->gatherHaloFlags,
                                     &(solver->halo),
                                     &(solver->nonHalo));

  // count elements that contribute to global C0 gather-scatter
  iint globalCount = 0;
  iint localCount = 0;
  iint *localHaloFlags = (iint*) calloc(mesh->Np*mesh->Nelements, sizeof(int));

  for(iint n=0;n<mesh->Np*mesh->Nelements;++n)
    localHaloFlags[mesh->gatherLocalIds[n]] += mesh->gatherHaloFlags[n];

  for(iint e=0;e<mesh->Nelements;++e){
    iint isHalo = 0;
    for(iint n=0;n<mesh->Np;++n){
      if(localHaloFlags[e*mesh->Np+n]>0){
        isHalo = 1;
      }
      if(localHaloFlags[e*mesh->Np+n]<0){
        printf("found halo flag %d\n", localHaloFlags[e*mesh->Np+n]);
      }
    }
    globalCount += isHalo;
    localCount += 1-isHalo;
  }

  iint *globalGatherElementList    = (iint*) calloc(globalCount, sizeof(iint));
  iint *localGatherElementList = (iint*) calloc(localCount, sizeof(iint));

  globalCount = 0;
  localCount = 0;

  for(iint e=0;e<mesh->Nelements;++e){
    iint isHalo = 0;
    for(iint n=0;n<mesh->Np;++n){
      if(localHaloFlags[e*mesh->Np+n]>0){
        isHalo = 1;
      }
    }
    if(isHalo){
      globalGatherElementList[globalCount++] = e;
    }
    else{
      localGatherElementList[localCount++] = e;
    }
  }
  printf("local = %d, global = %d\n", localCount, globalCount);

  solver->NglobalGatherElements = globalCount;
  solver->NlocalGatherElements = localCount;

  if(globalCount)
    solver->o_globalGatherElementList =
      mesh->device.malloc(globalCount*sizeof(iint), globalGatherElementList);

  if(localCount)
    solver->o_localGatherElementList =
      mesh->device.malloc(localCount*sizeof(iint), localGatherElementList);

  free(localHaloFlags);

  return solver;
}

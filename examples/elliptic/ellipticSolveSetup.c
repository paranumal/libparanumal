#include "elliptic.h"

void ellipticSolveSetup(elliptic_t *elliptic, dfloat lambda, occa::kernelInfo &kernelInfo, setupAide options){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = elliptic->mesh;

  dlong Ntotal = mesh->Np*mesh->Nelements;
  dlong Nblock = (Ntotal+blockSize-1)/blockSize;
  dlong Nhalo = mesh->Np*mesh->totalHaloPairs;
  dlong Nall   = Ntotal + Nhalo;

  //tau
  elliptic->tau = 2.0*(mesh->N+1)*(mesh->N+2)/2.0;

  elliptic->p   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  elliptic->z   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  elliptic->Ax  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  elliptic->Ap  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  elliptic->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));

  elliptic->grad = (dfloat*) calloc(Nall*4, sizeof(dfloat));

  elliptic->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->p);
  elliptic->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), elliptic->p);
  elliptic->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->z);

  elliptic->o_res = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->z);
  elliptic->o_Sres = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->z);
  elliptic->o_Ax  = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->p);
  elliptic->o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->Ap);
  elliptic->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), elliptic->tmp);

  elliptic->o_grad  = mesh->device.malloc(Nall*4*sizeof(dfloat), elliptic->grad);

  //setup async halo stream
  elliptic->defaultStream = mesh->defaultStream;
  elliptic->dataStream = mesh->dataStream;

  dlong Nbytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  if(Nbytes>0){
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(Nbytes, NULL);

    elliptic->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    elliptic->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();

    occa::memory o_gradSendBuffer = mesh->device.mappedAlloc(2*Nbytes, NULL);
    occa::memory o_gradRecvBuffer = mesh->device.mappedAlloc(2*Nbytes, NULL);

    elliptic->gradSendBuffer = (dfloat*) o_gradSendBuffer.getMappedPointer();
    elliptic->gradRecvBuffer = (dfloat*) o_gradRecvBuffer.getMappedPointer();
  }else{
    elliptic->sendBuffer = NULL;
    elliptic->recvBuffer = NULL;
  }
  mesh->device.setStream(elliptic->defaultStream);

  elliptic->type = strdup(dfloatString);

  elliptic->Nblock = Nblock;

  //fill geometric factors in halo
  if(mesh->totalHaloPairs){
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
    free(vgeoSendBuffer);
  }


  //build inverse of mass matrix
  mesh->invMM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  for (int n=0;n<mesh->Np*mesh->Np;n++)
    mesh->invMM[n] = mesh->MM[n];
  matrixInverse(mesh->Np,mesh->invMM);


  //check all the bounaries for a Dirichlet
  bool allNeumann = (lambda==0) ? true :false;
  elliptic->allNeumannPenalty = 1;
  hlong localElements = (hlong) mesh->Nelements;
  hlong totalElements = 0;
  MPI_Allreduce(&localElements, &totalElements, 1, MPI_HLONG, MPI_SUM, MPI_COMM_WORLD);
  elliptic->allNeumannScale = 1.0/sqrt(mesh->Np*totalElements);

  elliptic->EToB = (int *) calloc(mesh->Nelements*mesh->Nfaces,sizeof(int));
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if (bc>0) {
        int BC = elliptic->BCType[bc]; //get the type of boundary
        elliptic->EToB[e*mesh->Nfaces+f] = BC; //record it
        if (BC!=2) allNeumann = false; //check if its a Dirchlet
      }
    }
  }
  MPI_Allreduce(&allNeumann, &(elliptic->allNeumann), 1, MPI::BOOL, MPI_LAND, MPI_COMM_WORLD);
  if (rank==0&& options.compareArgs("VERBOSE","TRUE")) printf("allNeumann = %d \n", elliptic->allNeumann);

  //set surface mass matrix for continuous boundary conditions
  mesh->sMT = (dfloat *) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp,sizeof(dfloat));
  for (int n=0;n<mesh->Np;n++) {
    for (int m=0;m<mesh->Nfp*mesh->Nfaces;m++) {
      dfloat MSnm = 0;
      for (int i=0;i<mesh->Np;i++){
        MSnm += mesh->MM[n+i*mesh->Np]*mesh->LIFT[m+i*mesh->Nfp*mesh->Nfaces];
      }
      mesh->sMT[n+m*mesh->Np]  = MSnm;
    }
  }
  mesh->o_sMT = mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat), mesh->sMT);

  /* sparse basis setup */
  //build inverse vandermonde matrix
  mesh->invSparseV = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  for (int n=0;n<mesh->Np*mesh->Np;n++)
    mesh->invSparseV[n] = mesh->sparseV[n];
  matrixInverse(mesh->Np,mesh->invSparseV);

  int paddedRowSize = 4*((mesh->SparseNnzPerRow+3)/4); //make the nnz per row a multiple of 4

  char* IndTchar = (char*) calloc(paddedRowSize*mesh->Np,sizeof(char));
  for (int m=0;m<paddedRowSize/4;m++) {
    for (int n=0;n<mesh->Np;n++) {
      for (int k=0;k<4;k++) {
        if (k+4*m < mesh->SparseNnzPerRow) {
          IndTchar[k+4*n+m*4*mesh->Np] = mesh->sparseStackedNZ[n+(k+4*m)*mesh->Np];        
        } else {
          IndTchar[k+4*n+m*4*mesh->Np] = 0;
        }
      }
    }
  }
  
  mesh->o_sparseSrrT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSrrT);
  mesh->o_sparseSrsT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSrsT);
  mesh->o_sparseSssT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSssT);
  
  mesh->o_sparseStackedNZ = mesh->device.malloc(mesh->Np*paddedRowSize*sizeof(char), IndTchar);
  free(IndTchar);

  //copy boundary flags
  elliptic->o_EToB = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int), elliptic->EToB);

  if (rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    occa::setVerboseCompilation(true);
  else 
    occa::setVerboseCompilation(false);

  //add standard boundary functions
  char *boundaryHeaderFileName;
  boundaryHeaderFileName = strdup(DHOLMES "/examples/elliptic/data/ellipticBoundary2D.h");
  kernelInfo.addInclude(boundaryHeaderFileName);

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");
  }

  if(mesh->device.mode()=="Serial")
    kernelInfo.addCompilerFlag("-g");

  kernelInfo.addDefine("p_blockSize", blockSize);

  // add custom defines
  kernelInfo.addDefine("p_NpP", (mesh->Np+mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);
  kernelInfo.addDefine("p_SparseNnzPerRow", paddedRowSize);


  //sizes for the coarsen and prolongation kernels. degree N to degree 1
  kernelInfo.addDefine("p_NpFine", mesh->Np);
  kernelInfo.addDefine("p_NpCoarse", mesh->Nverts);

  kernelInfo.addDefine("p_NpFEM", mesh->NpFEM);

  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax);

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 256/mesh->Np; // works for CUDA
  int NnodesV = 1; //hard coded for now
  kernelInfo.addDefine("p_NblockV", NblockV);
  kernelInfo.addDefine("p_NnodesV", NnodesV);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  int NblockP = 256/(4*mesh->Np); // get close to 256 threads
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

  mesh->maskKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/mask.okl",
               "mask",
               kernelInfo);

  
  elliptic->AxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTri2D.okl",
               "ellipticAxTri2D",
               kernelInfo);

  elliptic->partialAxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTri2D.okl",
               "ellipticPartialAxTri2D",
               kernelInfo);

  elliptic->weightedInnerProduct1Kernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/weightedInnerProduct1.okl",
				       "weightedInnerProduct1",
				       kernelInfo);

  elliptic->weightedInnerProduct2Kernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/weightedInnerProduct2.okl",
				       "weightedInnerProduct2",
				       kernelInfo);

  elliptic->innerProductKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/innerProduct.okl",
				       "innerProduct",
				       kernelInfo);

  elliptic->scaledAddKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/scaledAdd.okl",
					 "scaledAdd",
					 kernelInfo);

  elliptic->dotMultiplyKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/dotMultiply.okl",
					 "dotMultiply",
					 kernelInfo);

  elliptic->dotDivideKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/dotDivide.okl",
					 "dotDivide",
					 kernelInfo);

  if (options.compareArgs("BASIS","BERN")) {

    elliptic->gradientKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientBBTri2D.okl",
               "ellipticGradientBBTri2D",
           kernelInfo);

    elliptic->partialGradientKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientBBTri2D.okl",
                 "ellipticPartialGradientBBTri2D",
                  kernelInfo);
  
    elliptic->ipdgKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgBBTri2D.okl",
                 "ellipticAxIpdgBBTri2D",
                 kernelInfo);

    elliptic->partialIpdgKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgBBTri2D.okl",
                 "ellipticPartialAxIpdgBBTri2D",
                 kernelInfo);
      
  } else if (options.compareArgs("BASIS","NODAL")) {

    elliptic->gradientKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientTri2D.okl",
               "ellipticGradientTri2D_v0",
           kernelInfo);

    elliptic->partialGradientKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientTri2D.okl",
             "ellipticPartialGradientTri2D_v0",
           kernelInfo);
 
    elliptic->ipdgKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgTri2D.okl",
                 "ellipticAxIpdgTri2D",
                 kernelInfo);

    elliptic->partialIpdgKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgTri2D.okl",
                 "ellipticPartialAxIpdgTri2D",
                 kernelInfo);
  }


  //on-host version of gather-scatter
  int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;
  mesh->hostGsh = gsParallelGatherScatterSetup(mesh->Nelements*mesh->Np, mesh->globalIds,verbose);

  // set up separate gather scatter infrastructure for halo and non halo nodes
  ellipticParallelGatherScatterSetup(elliptic);

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  elliptic->mapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) elliptic->mapB[n+e*mesh->Np] = 1E9;
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh->Nfp;n++) {
          int BCFlag = elliptic->BCType[bc];
          int fid = mesh->faceNodes[n+f*mesh->Nfp];
          elliptic->mapB[fid+e*mesh->Np] = mymin(BCFlag,elliptic->mapB[fid+e*mesh->Np]);
        }
      }
    }
  }
  gsParallelGatherScatter(mesh->hostGsh, elliptic->mapB, "int", "min"); 

  //use the bc flags to find masked ids
  elliptic->Nmasked = 0;
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (elliptic->mapB[n] == 1E9) {
      elliptic->mapB[n] = 0.;
    } else if (elliptic->mapB[n] == 1) { //Dirichlet boundary
      elliptic->Nmasked++;
    }
  }
  elliptic->o_mapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), elliptic->mapB);
  
  elliptic->maskIds = (dlong *) calloc(elliptic->Nmasked, sizeof(dlong));
  elliptic->Nmasked =0; //reset
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (elliptic->mapB[n] == 1) elliptic->maskIds[elliptic->Nmasked++] = n;
  }
  if (elliptic->Nmasked) elliptic->o_maskIds = mesh->device.malloc(elliptic->Nmasked*sizeof(dlong), elliptic->maskIds);

  /*preconditioner setup */
  elliptic->precon = (precon_t*) calloc(1, sizeof(precon_t));

  #if 0
  elliptic->precon->overlappingPatchKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticOasPreconTri2D.okl",
               "ellipticOasPreconTri2D",
               kernelInfo);
  #endif

  elliptic->precon->coarsenKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconCoarsen.okl",
				       "ellipticPreconCoarsen",
				       kernelInfo);

  elliptic->precon->prolongateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconProlongate.okl",
				       "ellipticPreconProlongate",
				       kernelInfo);


  elliptic->precon->blockJacobiKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBlockJacobiPreconTri2D.okl",
				       "ellipticBlockJacobiPreconTri2D",
				       kernelInfo);

  elliptic->precon->partialblockJacobiKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBlockJacobiPreconTri2D.okl",
               "ellipticPartialBlockJacobiPreconTri2D",
               kernelInfo);

  elliptic->precon->approxBlockJacobiSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
               "ellipticApproxBlockJacobiSolver2D",
               kernelInfo);

  elliptic->precon->exactBlockJacobiSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
               "ellipticExactBlockJacobiSolver2D",
               kernelInfo);

  elliptic->precon->SEMFEMInterpKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticSEMFEMInterp.okl",
               "ellipticSEMFEMInterp",
               kernelInfo);

  elliptic->precon->SEMFEMAnterpKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticSEMFEMAnterp.okl",
               "ellipticSEMFEMAnterp",
               kernelInfo);

  long long int pre = mesh->device.memoryAllocated();

  occaTimerTic(mesh->device,"PreconditionerSetup");
  ellipticPreconditionerSetup(elliptic, elliptic->ogs, lambda);
  occaTimerToc(mesh->device,"PreconditionerSetup");

  long long int usedBytes = mesh->device.memoryAllocated()-pre;

  elliptic->precon->preconBytes = usedBytes;
}

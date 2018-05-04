#include "elliptic.h"

void ellipticSolveSetup(elliptic_t *elliptic, dfloat lambda, occa::kernelInfo &kernelInfo){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  //sanity checking
  if (options.compareArgs("BASIS","BERN") && elliptic->elementType!=TRIANGLES) {
    printf("ERROR: BERN basis is only available for triangular elements\n");
    MPI_Finalize();
    exit(-1);
  }
  if (options.compareArgs("PRECONDITIONER","MASSMATRIX") && elliptic->elementType!=TRIANGLES 
                                                         && elliptic->elementType!=TETRAHEDRA ) {
    printf("ERROR: MASSMATRIX preconditioner is only available for triangle and tetrhedra elements. Use JACOBI instead.\n");
    MPI_Finalize();
    exit(-1);
  }
  if (options.compareArgs("PRECONDITIONER","MASSMATRIX") && lambda==0) {
    printf("ERROR: MASSMATRIX preconditioner is unavailble when lambda=0. \n");
    exit(-1);
  }

  dlong Ntotal = mesh->Np*mesh->Nelements;
  dlong Nblock = (Ntotal+blockSize-1)/blockSize;
  dlong Nhalo = mesh->Np*mesh->totalHaloPairs;
  dlong Nall   = Ntotal + Nhalo;

  //tau
  if (elliptic->elementType==TRIANGLES || elliptic->elementType==QUADRILATERALS)
    elliptic->tau = 2.0*(mesh->N+1)*(mesh->N+2)/2.0;
  else 
    elliptic->tau = 2.0*(mesh->N+1)*(mesh->N+3);

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
    dlong Nlocal = mesh->Nelements;
    dlong Nhalo = mesh->totalHaloPairs;
    size_t Nbytes = mesh->Nvgeo*sizeof(dfloat);

    if (elliptic->elementType==QUADRILATERALS || elliptic->elementType==HEXAHEDRA) {
      Nlocal *= mesh->Np;
      Nhalo *= mesh->Np;
      Nbytes *= mesh->Np;
    }

    dfloat *vgeoSendBuffer = (dfloat*) calloc(Nhalo*mesh->Nvgeo, sizeof(dfloat));

    // import geometric factors from halo elements
    mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat));

    meshHaloExchange(mesh,
         Nbytes,
         mesh->vgeo,
         vgeoSendBuffer,
         mesh->vgeo + Nlocal*mesh->Nvgeo);

    mesh->o_vgeo =
      mesh->device.malloc((Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);
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

  //copy boundary flags
  elliptic->o_EToB = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int), elliptic->EToB);

  if (rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    occa::setVerboseCompilation(true);
  else 
    occa::setVerboseCompilation(false);

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");
  }

  if(mesh->device.mode()=="Serial")
    kernelInfo.addCompilerFlag("-g");

  // set kernel name suffix
  char *suffix;
  
  if(elliptic->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(elliptic->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(elliptic->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(elliptic->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];


  for (int r=0;r<size;r++) {
    if (r==rank) {

      //mesh kernels 
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


      mesh->addScalarKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/addScalar.okl",
                   "addScalar",
                   kernelInfo);

      mesh->maskKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/mask.okl",
                   "mask",
                   kernelInfo);


      kernelInfo.addDefine("p_blockSize", blockSize);
      

      mesh->sumKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/sum.okl",
                   "sum",
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
      
      // add custom defines
      kernelInfo.addDefine("p_NpP", (mesh->Np+mesh->Nfp*mesh->Nfaces));
      kernelInfo.addDefine("p_Nverts", mesh->Nverts);

      //sizes for the coarsen and prolongation kernels. degree N to degree 1
      kernelInfo.addDefine("p_NpFine", mesh->Np);
      kernelInfo.addDefine("p_NpCoarse", mesh->Nverts);

      if (elliptic->elementType==QUADRILATERALS || elliptic->elementType==HEXAHEDRA) {
        kernelInfo.addDefine("p_NqFine", mesh->N+1);
        kernelInfo.addDefine("p_NqCoarse", 2);
      }

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

      //add standard boundary functions
      char *boundaryHeaderFileName;
      if (elliptic->dim==2)
        boundaryHeaderFileName = strdup(DHOLMES "/examples/elliptic/data/ellipticBoundary2D.h");
      else if (elliptic->dim==3)
        boundaryHeaderFileName = strdup(DHOLMES "/examples/elliptic/data/ellipticBoundary3D.h");
      kernelInfo.addInclude(boundaryHeaderFileName);

      sprintf(fileName, "okl/ellipticAx%s.okl", suffix);
      sprintf(kernelName, "ellipticAx%s", suffix);
      elliptic->AxKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

      sprintf(kernelName, "ellipticPartialAx%s", suffix);
      elliptic->partialAxKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);


      if (options.compareArgs("BASIS","BERN")) {

        sprintf(fileName, "okl/ellipticGradientBB%s.okl", suffix);
        sprintf(kernelName, "ellipticGradientBB%s", suffix);

        elliptic->gradientKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialGradientBB%s", suffix);
        elliptic->partialGradientKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
      
        sprintf(fileName, "okl/ellipticAxIpdgBB%s.okl", suffix);
        sprintf(kernelName, "ellipticAxIpdgBB%s", suffix);
        elliptic->ipdgKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialAxIpdgBB%s", suffix);
        elliptic->partialIpdgKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
          
      } else if (options.compareArgs("BASIS","NODAL")) {

        sprintf(fileName, "okl/ellipticGradient%s.okl", suffix);
        sprintf(kernelName, "ellipticGradient%s", suffix);

        elliptic->gradientKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialGradient%s", suffix);
        elliptic->partialGradientKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

        sprintf(fileName, "okl/ellipticAxIpdg%s.okl", suffix);
        sprintf(kernelName, "ellipticAxIpdg%s", suffix);
        elliptic->ipdgKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialAxIpdg%s", suffix);
        elliptic->partialIpdgKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
      }
    }
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


  for (int r=0;r<size;r++) {
    if (r==rank) {

      sprintf(fileName, "okl/ellipticPreconCoarsen%s.okl", suffix);
      sprintf(kernelName, "ellipticPreconCoarsen%s", suffix);
      elliptic->precon->coarsenKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

      sprintf(fileName, "okl/ellipticPreconProlongate%s.okl", suffix);
      sprintf(kernelName, "ellipticPreconProlongate%s", suffix);
      elliptic->precon->prolongateKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

      sprintf(fileName, "okl/ellipticBlockJacobiPrecon.okl");
      sprintf(kernelName, "ellipticBlockJacobiPrecon");
      elliptic->precon->blockJacobiKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

      sprintf(kernelName, "ellipticPartialBlockJacobiPrecon");
      elliptic->precon->partialblockJacobiKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

      sprintf(fileName, "okl/ellipticPatchSolver.okl");
      sprintf(kernelName, "ellipticApproxBlockJacobiSolver");
      elliptic->precon->approxBlockJacobiSolverKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

      if (   elliptic->elementType == TRIANGLES 
          || elliptic->elementType == TETRAHEDRA) {
        elliptic->precon->SEMFEMInterpKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticSEMFEMInterp.okl",
                     "ellipticSEMFEMInterp",
                     kernelInfo);

        elliptic->precon->SEMFEMAnterpKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticSEMFEMAnterp.okl",
                     "ellipticSEMFEMAnterp",
                     kernelInfo);
      }
    }
  }

  long long int pre = mesh->device.memoryAllocated();

  occaTimerTic(mesh->device,"PreconditionerSetup");
  ellipticPreconditionerSetup(elliptic, elliptic->ogs, lambda);
  occaTimerToc(mesh->device,"PreconditionerSetup");

  long long int usedBytes = mesh->device.memoryAllocated()-pre;

  elliptic->precon->preconBytes = usedBytes;
}

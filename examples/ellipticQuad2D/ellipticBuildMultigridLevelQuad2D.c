#include "ellipticQuad2D.h"

void matrixInverse(int N, dfloat *A);

// create solver and mesh structs for multigrid levels
solver_t *ellipticBuildMultigridLevelQuad2D(solver_t *baseSolver, int Nc, int Nf, int *BCType, const char *options){

  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));
  memcpy(solver,baseSolver,sizeof(solver_t));

  //populate the mini-mesh using the mesh struct
  mesh2D *mesh = (mesh2D*) calloc(1,sizeof(mesh2D));
  memcpy(mesh,baseSolver->mesh,sizeof(mesh2D));

  solver->mesh = mesh;

  // load reference (r,s) element nodes
  meshLoadReferenceNodesQuad2D(mesh, Nc);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesQuad2D(mesh);

  // compute geometric factors
  meshGeometricFactorsQuad2D(mesh);

  // create halo extension for x,y arrays
  dlong totalHaloNodes = mesh->totalHaloPairs*mesh->Np;
  dlong localNodes     = mesh->Nelements*mesh->Np;
  // temporary send buffer
  dfloat *sendBuffer = (dfloat*) calloc(totalHaloNodes, sizeof(dfloat));

  // extend x,y arrays to hold coordinates of node coordinates of elements in halo
  mesh->x = (dfloat*) realloc(mesh->x, (localNodes+totalHaloNodes)*sizeof(dfloat));
  mesh->y = (dfloat*) realloc(mesh->y, (localNodes+totalHaloNodes)*sizeof(dfloat));
  meshHaloExchange(mesh, mesh->Np*sizeof(dfloat), mesh->x, sendBuffer, mesh->x + localNodes);
  meshHaloExchange(mesh, mesh->Np*sizeof(dfloat), mesh->y, sendBuffer, mesh->y + localNodes);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes2D(mesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsQuad2D(mesh);

  // global nodes
  meshParallelConnectNodes(mesh);

  //dont need these once vmap is made
  free(mesh->x);
  free(mesh->y);
  free(sendBuffer);

  dlong Ntotal = mesh->Np*mesh->Nelements;
  dlong Nblock = (Ntotal+blockSize-1)/blockSize;

  solver->Nblock = Nblock;

  mesh->o_D = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);
  
  mesh->o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*mesh->Np*sizeof(dfloat),
      mesh->vgeo);
  mesh->o_sgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nfp*mesh->Nsgeo*sizeof(dfloat),
      mesh->sgeo);
  mesh->o_ggeo =
    mesh->device.malloc(mesh->Nelements*mesh->Np*mesh->Nggeo*sizeof(dfloat),
      mesh->ggeo);

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(int),
      mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(int),
      mesh->vmapP);


  //fill geometric factors in halo
  if(mesh->totalHaloPairs){
    dlong Nlocal = mesh->Np*mesh->Nelements;
    dlong Nhalo = mesh->totalHaloPairs*mesh->Np;
    dfloat *vgeoSendBuffer = (dfloat*) calloc(Nhalo*mesh->Nvgeo, sizeof(dfloat));

    // import geometric factors from halo elements
    mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat));

    meshHaloExchange(mesh,
         mesh->Nvgeo*mesh->Np*sizeof(dfloat),
         mesh->vgeo,
         vgeoSendBuffer,
         mesh->vgeo + Nlocal*mesh->Nvgeo);

    mesh->o_vgeo =
      mesh->device.malloc((Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);
    free(vgeoSendBuffer);
  }


  //set the normalization constant for the allNeumann Poisson problem on this coarse mesh
  hlong localElements = (hlong) mesh->Nelements;
  hlong totalElements = 0;
  MPI_Allreduce(&localElements, &totalElements, 1, MPI_HLONG, MPI_SUM, MPI_COMM_WORLD);
  solver->allNeumannScale = 1.0/sqrt(mesh->Np*totalElements);

  // info for kernel construction
  occa::kernelInfo kernelInfo;

  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_N", mesh->N);
  kernelInfo.addDefine("p_Nq", mesh->N+1);
  kernelInfo.addDefine("p_Np", mesh->Np);
  kernelInfo.addDefine("p_Nfp", mesh->Nfp);
  kernelInfo.addDefine("p_Nfaces", mesh->Nfaces);
  kernelInfo.addDefine("p_NfacesNfp", mesh->Nfp*mesh->Nfaces);
  kernelInfo.addDefine("p_Nvgeo", mesh->Nvgeo);
  kernelInfo.addDefine("p_Nsgeo", mesh->Nsgeo);
  kernelInfo.addDefine("p_Nggeo", mesh->Nggeo);

  kernelInfo.addDefine("p_NXID", NXID);
  kernelInfo.addDefine("p_NYID", NYID);
  kernelInfo.addDefine("p_SJID", SJID);
  kernelInfo.addDefine("p_IJID", IJID);
  kernelInfo.addDefine("p_WSJID", WSJID);
  kernelInfo.addDefine("p_IHID", IHID);

  kernelInfo.addDefine("p_max_EL_nnz", mesh->max_EL_nnz); // for Bernstein Bezier lift

  kernelInfo.addDefine("p_cubNp", mesh->cubNp);
  kernelInfo.addDefine("p_intNfp", mesh->intNfp);
  kernelInfo.addDefine("p_intNfpNfaces", mesh->intNfp*mesh->Nfaces);

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

  if(sizeof(int)==4){
    kernelInfo.addDefine("dlong","int");
  }
  if(sizeof(int)==8){
    kernelInfo.addDefine("dlong","long long int");
  }

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo.addCompilerFlag("--ftz=true");
    kernelInfo.addCompilerFlag("--prec-div=false");
    kernelInfo.addCompilerFlag("--prec-sqrt=false");
    kernelInfo.addCompilerFlag("--use_fast_math");
    kernelInfo.addCompilerFlag("--fmad=true"); // compiler option for cuda
    kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");
  }

  kernelInfo.addDefine("p_G00ID", G00ID);
  kernelInfo.addDefine("p_G01ID", G01ID);
  kernelInfo.addDefine("p_G11ID", G11ID);
  kernelInfo.addDefine("p_GWJID", GWJID);


  kernelInfo.addDefine("p_RXID", RXID);
  kernelInfo.addDefine("p_SXID", SXID);

  kernelInfo.addDefine("p_RYID", RYID);
  kernelInfo.addDefine("p_SYID", SYID);

  kernelInfo.addDefine("p_JID", JID);
  kernelInfo.addDefine("p_JWID", JWID);


  //add standard boundary functions
  char *boundaryHeaderFileName;
  boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticQuad2D/ellipticBoundary2D.h");
  kernelInfo.addInclude(boundaryHeaderFileName);

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  kernelInfo.addDefine("p_blockSize", blockSize);

  // add custom defines
  kernelInfo.addDefine("p_NqP", (mesh->Nq+2));
  kernelInfo.addDefine("p_NpP", (mesh->NqP*mesh->NqP));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);


  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax);

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 256/mesh->Np; // get close to 256 threads
  kernelInfo.addDefine("p_NblockV", NblockV);

  int one = 1; //set to one for now. TODO: try optimizing over these
  kernelInfo.addDefine("p_NnodesV", one);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  int NblockP = mymax(256/(5*mesh->Np),1); // get close to 256 threads
  kernelInfo.addDefine("p_NblockP", NblockP);

  int NblockG;
  if(mesh->Np<=32) NblockG = ( 32/mesh->Np );
  else NblockG = 256/mesh->Np;
  kernelInfo.addDefine("p_NblockG", NblockG);

  solver->scaledAddKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/scaledAdd.okl",
           "scaledAdd",
           kernelInfo);

  solver->dotMultiplyKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/dotMultiply.okl",
           "dotMultiply",
           kernelInfo);

  solver->AxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxQuad2D.okl",
               "ellipticAxQuad2D_e0",
               kernelInfo);

  solver->partialAxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxQuad2D.okl",
               "ellipticPartialAxQuad2D",
               kernelInfo);

  solver->gradientKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientQuad2D.okl",
               "ellipticGradientQuad2D",
           kernelInfo);

  solver->partialGradientKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientQuad2D.okl",
               "ellipticPartialGradientQuad2D",
                kernelInfo);

  solver->ipdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgQuad2D.okl",
               "ellipticAxIpdgQuad2D",
               kernelInfo);

  solver->partialIpdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgQuad2D.okl",
               "ellipticPartialAxIpdgQuad2D",
               kernelInfo);

  //new precon struct
  solver->precon = (precon_t *) calloc(1,sizeof(precon_t));

  solver->precon->overlappingPatchKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticOasPreconQuad2D.okl",
               "ellipticOasPreconQuad2D",
               kernelInfo);

  solver->precon->restrictKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconRestrictQuad2D.okl",
               "ellipticFooQuad2D",
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

  int NpFine   = (Nf+1)*(Nf+1);
  int NpCoarse = (Nc+1)*(Nc+1);
  int NqFine   = (Nf+1)*(Nf+1);
  int NqCoarse = (Nc+1)*(Nc+1);
  kernelInfo.addDefine("p_NpFine", NpFine);
  kernelInfo.addDefine("p_NpCoarse", NpCoarse);
  kernelInfo.addDefine("p_NqFine", Nf+1);
  kernelInfo.addDefine("p_NqCoarse", Nc+1);
  
  solver->precon->coarsenKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconCoarsenQuad2D.okl",
               "ellipticPreconCoarsenQuad2D",
               kernelInfo);

  solver->precon->prolongateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconProlongateQuad2D.okl",
               "ellipticPreconProlongateQuad2D",
               kernelInfo);

  //on host gather-scatter
  int verbose = strstr(options,"VERBOSE") ? 1:0;
  mesh->hostGsh = gsParallelGatherScatterSetup(mesh->Nelements*mesh->Np, mesh->globalIds, verbose);

  // set up separate gather scatter infrastructure for halo and non halo nodes
  ellipticParallelGatherScatterSetup(solver,options);

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  solver->mapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) solver->mapB[n+e*mesh->Np] = 1E9;
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh->Nfp;n++) {
          int BCFlag = BCType[bc];
          int fid = mesh->faceNodes[n+f*mesh->Nfp];
          solver->mapB[fid+e*mesh->Np] = mymin(BCFlag,solver->mapB[fid+e*mesh->Np]);
        }
      }
    }
  }
  gsParallelGatherScatter(mesh->hostGsh, solver->mapB, "int", "min"); 

  //use the bc flags to find masked ids
  solver->Nmasked = 0;
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (solver->mapB[n] == 1E9) {
      solver->mapB[n] = 0.;
    } else if (solver->mapB[n] == 1) { //Dirichlet boundary
      solver->Nmasked++;
    }
  }
  solver->o_mapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), solver->mapB);
  
  solver->maskIds = (dlong *) calloc(solver->Nmasked, sizeof(dlong));
  solver->Nmasked =0; //reset
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (solver->mapB[n] == 1) solver->maskIds[solver->Nmasked++] = n;
  }
  if (solver->Nmasked) solver->o_maskIds = mesh->device.malloc(solver->Nmasked*sizeof(dlong), solver->maskIds);

  return solver;
}

#include "elliptic.h"

// create elliptic and mesh structs for multigrid levels
elliptic_t *ellipticBuildMultigridLevel(elliptic_t *baseElliptic, int Nc, int Nf){

  elliptic_t *elliptic = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  memcpy(elliptic,baseElliptic,sizeof(elliptic_t));

  //populate the mini-mesh using the mesh struct
  mesh2D *mesh = (mesh2D*) calloc(1,sizeof(mesh2D));
  memcpy(mesh,baseElliptic->mesh,sizeof(mesh2D));

  elliptic->mesh = mesh;

  setupAide options = elliptic->options;

  // load reference (r,s) element nodes
  meshLoadReferenceNodesTri2D(mesh, Nc);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesTri2D(mesh);

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

  // global nodes
  meshParallelConnectNodes(mesh);

  //dont need these once vmap is made
  free(mesh->x);
  free(mesh->y);
  free(sendBuffer);

  dlong Ntotal = mesh->Np*mesh->Nelements;
  dlong Nblock = (Ntotal+blockSize-1)/blockSize;

  elliptic->Nblock = Nblock;

  // build Dr, Ds, LIFT transposes
  dfloat *DrT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *DsT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      DrT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
      DsT[n+m*mesh->Np] = mesh->Ds[n*mesh->Np+m];
    }
  }

  dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      LIFTT[n+m*mesh->Np] = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
    }
  }

  //build element stiffness matrices
  dfloat *SrrT, *SrsT, *SsrT, *SssT;
  mesh->Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  mesh->Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  mesh->Ssr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  mesh->Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  for (int n=0;n<mesh->Np;n++) {
    for (int m=0;m<mesh->Np;m++) {
      for (int k=0;k<mesh->Np;k++) {
        for (int l=0;l<mesh->Np;l++) {
          mesh->Srr[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
          mesh->Srs[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
          mesh->Ssr[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
          mesh->Sss[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
        }
      } 
    }
  }
  SrrT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  SrsT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  SsrT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  SssT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  for (int n=0;n<mesh->Np;n++) {
    for (int m=0;m<mesh->Np;m++) {  
      SrrT[m+n*mesh->Np] = mesh->Srr[n+m*mesh->Np];
      SrsT[m+n*mesh->Np] = mesh->Srs[n+m*mesh->Np];
      SsrT[m+n*mesh->Np] = mesh->Ssr[n+m*mesh->Np];
      SssT[m+n*mesh->Np] = mesh->Sss[n+m*mesh->Np];
    }
  }

  // deriv operators: transpose from row major to column major
  int *D1ids = (int*) calloc(mesh->Np*3,sizeof(int));
  int *D2ids = (int*) calloc(mesh->Np*3,sizeof(int));
  int *D3ids = (int*) calloc(mesh->Np*3,sizeof(int));
  dfloat *Dvals = (dfloat*) calloc(mesh->Np*3,sizeof(dfloat));

  dfloat *VBq = (dfloat*) calloc(mesh->Np*mesh->cubNp,sizeof(dfloat));
  dfloat *PBq = (dfloat*) calloc(mesh->Np*mesh->cubNp,sizeof(dfloat));

  dfloat *L0vals = (dfloat*) calloc(mesh->Nfp*3,sizeof(dfloat)); // tridiag
  int *ELids = (int*) calloc(1+mesh->Np*mesh->max_EL_nnz,sizeof(int));
  dfloat *ELvals = (dfloat*) calloc(1+mesh->Np*mesh->max_EL_nnz,sizeof(dfloat));


  for (int i = 0; i < mesh->Np; ++i){
    for (int j = 0; j < 3; ++j){
      D1ids[i+j*mesh->Np] = mesh->D1ids[j+i*3];
      D2ids[i+j*mesh->Np] = mesh->D2ids[j+i*3];
      D3ids[i+j*mesh->Np] = mesh->D3ids[j+i*3];
      Dvals[i+j*mesh->Np] = mesh->Dvals[j+i*3];
    }
  }

  for (int i = 0; i < mesh->cubNp; ++i){
    for (int j = 0; j < mesh->Np; ++j){
      VBq[i+j*mesh->cubNp] = mesh->VBq[j+i*mesh->Np];
      PBq[j+i*mesh->Np] = mesh->PBq[i+j*mesh->cubNp];
    }
  }


  for (int i = 0; i < mesh->Nfp; ++i){
    for (int j = 0; j < 3; ++j){
      L0vals[i+j*mesh->Nfp] = mesh->L0vals[j+i*3];
    }
  }

  for (int i = 0; i < mesh->Np; ++i){
    for (int j = 0; j < mesh->max_EL_nnz; ++j){
      ELids[i + j*mesh->Np] = mesh->ELids[j+i*mesh->max_EL_nnz];
      ELvals[i + j*mesh->Np] = mesh->ELvals[j+i*mesh->max_EL_nnz];
    }
  }

  //BB mass matrix
  mesh->BBMM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  for (int n = 0; n < mesh->Np; ++n){
    for (int m = 0; m < mesh->Np; ++m){
      for (int i = 0; i < mesh->Np; ++i){
        for (int j = 0; j < mesh->Np; ++j){
          mesh->BBMM[n+m*mesh->Np] += mesh->VB[m+j*mesh->Np]*mesh->MM[i+j*mesh->Np]*mesh->VB[n+i*mesh->Np];
        }
      }
    }
  }

  mesh->o_Dr = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
           mesh->Dr);

  mesh->o_Ds = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
           mesh->Ds);

  mesh->o_DrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
            DrT);

  mesh->o_DsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
            DsT);

  mesh->o_LIFT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
      mesh->LIFT);

  mesh->o_LIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
      LIFTT);

  mesh->o_MM =
    mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
      mesh->MM);

  mesh->o_SrrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SrrT);
  mesh->o_SrsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SrsT);
  mesh->o_SsrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SsrT);
  mesh->o_SssT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SssT);

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(int),
      mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(int),
      mesh->vmapP);

  mesh->o_D1ids = mesh->device.malloc(mesh->Np*3*sizeof(int),D1ids);
  mesh->o_D2ids = mesh->device.malloc(mesh->Np*3*sizeof(int),D2ids);
  mesh->o_D3ids = mesh->device.malloc(mesh->Np*3*sizeof(int),D3ids);
  mesh->o_Dvals = mesh->device.malloc(mesh->Np*3*sizeof(dfloat),Dvals);

  mesh->o_BBMM = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),mesh->BBMM);

  mesh->o_VBq = mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),VBq);
  mesh->o_PBq = mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),PBq);

  mesh->o_L0vals = mesh->device.malloc(mesh->Nfp*3*sizeof(dfloat),L0vals);
  mesh->o_ELids =
    mesh->device.malloc(mesh->Np*mesh->max_EL_nnz*sizeof(int),ELids);
  mesh->o_ELvals =
    mesh->device.malloc(mesh->Np*mesh->max_EL_nnz*sizeof(dfloat),ELvals);

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

  
  //set the normalization constant for the allNeumann Poisson problem on this coarse mesh
  hlong localElements = (hlong) mesh->Nelements;
  hlong totalElements = 0;
  MPI_Allreduce(&localElements, &totalElements, 1, MPI_HLONG, MPI_SUM, MPI_COMM_WORLD);
  elliptic->allNeumannScale = 1.0/sqrt(mesh->Np*totalElements);

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

  if(sizeof(dlong)==4){
    kernelInfo.addDefine("dlong","int");
  }
  if(sizeof(dlong)==8){
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
  boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTri2D/ellipticBoundary2D.h");
  kernelInfo.addInclude(boundaryHeaderFileName);

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  kernelInfo.addDefine("p_blockSize", blockSize);

  // add custom defines
  kernelInfo.addDefine("p_NpP", (mesh->Np+mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);
  kernelInfo.addDefine("p_SparseNnzPerRow", paddedRowSize);

  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax);

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int one = 1; //set to one for now. TODO: try optimizing over these
  kernelInfo.addDefine("p_NnodesV", one);

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

  elliptic->scaledAddKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/scaledAdd.okl",
           "scaledAdd",
           kernelInfo);

  elliptic->dotMultiplyKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/dotMultiply.okl",
           "dotMultiply",
           kernelInfo);


  elliptic->AxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTri2D.okl",
               "ellipticAxTri2D",
               kernelInfo);

  elliptic->partialAxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTri2D.okl",
               "ellipticPartialAxTri2D",
               kernelInfo);

  if (options.compareArgs("BASIS", "BERN")) {

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

  } else if (options.compareArgs("BASIS", "NODAL")) {

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

  //new precon struct
  elliptic->precon = (precon_t *) calloc(1,sizeof(precon_t));

#if 0
  elliptic->precon->overlappingPatchKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticOasPreconTri2D.okl",
               "ellipticOasPreconTri2D",
               kernelInfo);
#endif

  elliptic->precon->blockJacobiKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBlockJacobiPreconTri2D.okl",
               "ellipticBlockJacobiPreconTri2D",
               kernelInfo);

  elliptic->precon->approxBlockJacobiSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
               "ellipticApproxBlockJacobiSolver2D",
               kernelInfo);

  elliptic->precon->exactBlockJacobiSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
               "ellipticExactBlockJacobiSolver2D",
               kernelInfo);

  //sizes for the coarsen and prolongation kernels. degree NFine to degree N
  int NpFine   = (Nf+1)*(Nf+2)/2;
  int NpCoarse = (Nc+1)*(Nc+2)/2;
  kernelInfo.addDefine("p_NpFine", NpFine);
  kernelInfo.addDefine("p_NpCoarse", NpCoarse);

  elliptic->precon->coarsenKernel =
    elliptic->mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconCoarsen.okl",
             "ellipticPreconCoarsen",
             kernelInfo);

  elliptic->precon->prolongateKernel =
    elliptic->mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconProlongate.okl",
             "ellipticPreconProlongate",
             kernelInfo);

  //on host gather-scatter
  int verbose = options.compareArgs("VERBOSE", "TRUE") ? 1:0;
  mesh->hostGsh = gsParallelGatherScatterSetup(mesh->Nelements*mesh->Np, mesh->globalIds, verbose);

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

  free(DrT); free(DsT); free(LIFTT);
  free(SrrT); free(SrsT); free(SsrT); free(SssT);
  free(D1ids); free(D2ids); free(D3ids); free(Dvals);

  free(VBq); free(PBq);
  free(L0vals); free(ELids ); free(ELvals);

  return elliptic;
}

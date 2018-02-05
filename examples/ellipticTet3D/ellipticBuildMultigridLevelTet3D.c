#include "ellipticTet3D.h"

void matrixInverse(int N, dfloat *A);

// create solver and mesh structs for multigrid levels
solver_t *ellipticBuildMultigridLevelTet3D(solver_t *baseSolver, int Nc, int Nf, int *BCType, const char *options){

  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));
  memcpy(solver,baseSolver,sizeof(solver_t));

  //populate the mini-mesh using the mesh struct
  mesh3D *mesh = (mesh3D*) calloc(1,sizeof(mesh3D));
  memcpy(mesh,baseSolver->mesh,sizeof(mesh3D));

  solver->mesh = mesh;

  // load reference (r,s) element nodes
  meshLoadReferenceNodesTet3D(mesh, Nc);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesTet3D(mesh);

  // create halo extension for x,y arrays
  iint totalHaloNodes = mesh->totalHaloPairs*mesh->Np;
  iint localNodes     = mesh->Nelements*mesh->Np;
  // temporary send buffer
  dfloat *sendBuffer = (dfloat*) calloc(totalHaloNodes, sizeof(dfloat));

  // extend x,y arrays to hold coordinates of node coordinates of elements in halo
  mesh->x = (dfloat*) realloc(mesh->x, (localNodes+totalHaloNodes)*sizeof(dfloat));
  mesh->y = (dfloat*) realloc(mesh->y, (localNodes+totalHaloNodes)*sizeof(dfloat));
  mesh->z = (dfloat*) realloc(mesh->z, (localNodes+totalHaloNodes)*sizeof(dfloat));
  meshHaloExchange(mesh, mesh->Np*sizeof(dfloat), mesh->x, sendBuffer, mesh->x + localNodes);
  meshHaloExchange(mesh, mesh->Np*sizeof(dfloat), mesh->y, sendBuffer, mesh->y + localNodes);
  meshHaloExchange(mesh, mesh->Np*sizeof(dfloat), mesh->z, sendBuffer, mesh->z + localNodes);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes3D(mesh);

  // global nodes
  meshParallelConnectNodes(mesh); 

  //dont need these once vmap is made
  free(mesh->x);
  free(mesh->y);
  free(mesh->z);
  free(sendBuffer);

  iint Ntotal = mesh->Np*mesh->Nelements;
  iint Nblock = (Ntotal+blockSize-1)/blockSize;
  iint Nhalo = mesh->Np*mesh->totalHaloPairs;
  iint Nall   = Ntotal + Nhalo;

  solver->Nblock = Nblock;

  // build Dr, Ds, LIFT transposes
  dfloat *DrT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *DsT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *DtT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Np;++m){
      DrT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
      DsT[n+m*mesh->Np] = mesh->Ds[n*mesh->Np+m];
      DtT[n+m*mesh->Np] = mesh->Dt[n*mesh->Np+m];
    }
  }

  dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      LIFTT[n+m*mesh->Np] = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
    }
  }

  //build element stiffness matrices
  if (mesh->Nverts==4) {
    mesh->Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Srt = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Ssr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Sst = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Str = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Sts = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Stt = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {
        for (iint k=0;k<mesh->Np;k++) {
          for (iint l=0;l<mesh->Np;l++) {
            mesh->Srr[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
            mesh->Srs[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
            mesh->Srt[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dt[m+k*mesh->Np];
            mesh->Ssr[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
            mesh->Sss[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
            mesh->Sst[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dt[m+k*mesh->Np];
            mesh->Str[m+n*mesh->Np] += mesh->Dt[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
            mesh->Sts[m+n*mesh->Np] += mesh->Dt[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
            mesh->Stt[m+n*mesh->Np] += mesh->Dt[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dt[m+k*mesh->Np];
          }
        }
      }
    }
  }

  mesh->o_Dr = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
           mesh->Dr);

  mesh->o_Ds = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
           mesh->Ds);

  mesh->o_Dt = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
           mesh->Dt);

  mesh->o_DrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
            DrT);

  mesh->o_DsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
            DsT);

  mesh->o_DtT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
            DtT);

  mesh->o_LIFT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
      mesh->LIFT);

  mesh->o_LIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
      LIFTT);

  mesh->o_MM =
    mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
      mesh->MM);

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
      mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
      mesh->vmapP);


  //set the normalization constant for the allNeumann Poisson problem on this coarse mesh
  iint totalElements = 0;
  MPI_Allreduce(&(mesh->Nelements), &totalElements, 1, MPI_IINT, MPI_SUM, MPI_COMM_WORLD);
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

  kernelInfo.addDefine("p_max_EL_nnz", mesh->max_EL_nnz); // for Bernstein Bezier lift

  kernelInfo.addDefine("p_NXID", NXID);
  kernelInfo.addDefine("p_NYID", NYID);
  kernelInfo.addDefine("p_NZID", NZID);
  kernelInfo.addDefine("p_SJID", SJID);
  kernelInfo.addDefine("p_IJID", IJID);
  kernelInfo.addDefine("p_WSJID", WSJID);
  kernelInfo.addDefine("p_IHID", IHID);


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
    kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");
  }

  kernelInfo.addDefine("p_G00ID", G00ID);
  kernelInfo.addDefine("p_G01ID", G01ID);
  kernelInfo.addDefine("p_G02ID", G02ID);
  kernelInfo.addDefine("p_G11ID", G11ID);
  kernelInfo.addDefine("p_G12ID", G12ID);
  kernelInfo.addDefine("p_G22ID", G22ID);
  kernelInfo.addDefine("p_GWJID", GWJID);


  kernelInfo.addDefine("p_RXID", RXID);
  kernelInfo.addDefine("p_SXID", SXID);
  kernelInfo.addDefine("p_TXID", TXID);

  kernelInfo.addDefine("p_RYID", RYID);
  kernelInfo.addDefine("p_SYID", SYID);
  kernelInfo.addDefine("p_TYID", TYID);

  kernelInfo.addDefine("p_RZID", RZID);
  kernelInfo.addDefine("p_SZID", SZID);
  kernelInfo.addDefine("p_TZID", TZID);

  kernelInfo.addDefine("p_JID", JID);
  kernelInfo.addDefine("p_JWID", JWID);


  //add standard boundary functions
  char *boundaryHeaderFileName;
  boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTet3D/ellipticBoundary3D.h");
  kernelInfo.addInclude(boundaryHeaderFileName);

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  kernelInfo.addDefine("p_blockSize", blockSize);

  // add custom defines
  kernelInfo.addDefine("p_NpP", (mesh->Np+mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);

  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax);

  int NblockV = mymax(1,256/mesh->Np); // get close to 256 threads
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockP = mymax(1,512/(5*mesh->Np)); // get close to 256 threads
  kernelInfo.addDefine("p_NblockP", NblockP);

  int NblockG;
  if(mesh->Np<=32) NblockG = ( 32/mesh->Np );
  else NblockG = mymax(1,256/mesh->Np);
  kernelInfo.addDefine("p_NblockG", NblockG);

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

  solver->scaledAddKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/scaledAdd.okl",
           "scaledAdd",
           kernelInfo);

  solver->dotMultiplyKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/dotMultiply.okl",
           "dotMultiply",
           kernelInfo);


  solver->AxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTet3D.okl",
               "ellipticAxTet3D",
               kernelInfo);

  solver->partialAxKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTet3D.okl",
                 "ellipticPartialAxTet3D",
                 kernelInfo);


  solver->gradientKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientTet3D.okl",
               "ellipticGradientTet3D",
           kernelInfo);

  solver->partialGradientKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientTet3D.okl",
                 "ellipticPartialGradientTet3D",
                  kernelInfo);

  solver->ipdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgTet3D.okl",
               "ellipticAxIpdgTet3D",
               kernelInfo);

  solver->partialIpdgKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgTet3D.okl",
                 "ellipticPartialAxIpdgTet3D",
                 kernelInfo);

  //new precon struct
  solver->precon = (precon_t *) calloc(1,sizeof(precon_t));

#if 0
  solver->precon->overlappingPatchKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticOasPreconTet3D.okl",
               "ellipticOasPreconTet3D",
               kernelInfo);
#endif

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

  solver->precon->approxPatchSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver3D.okl",
               "ellipticApproxPatchSolver3D",
               kernelInfo);

  solver->precon->exactPatchSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver3D.okl",
               "ellipticExactPatchSolver3D",
               kernelInfo);

  solver->precon->patchGatherKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchGather.okl",
               "ellipticPatchGather",
               kernelInfo);

  solver->precon->approxFacePatchSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver3D.okl",
               "ellipticApproxFacePatchSolver3D",
               kernelInfo);

  solver->precon->exactFacePatchSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver3D.okl",
               "ellipticExactFacePatchSolver3D",
               kernelInfo);

  solver->precon->facePatchGatherKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchGather.okl",
               "ellipticFacePatchGather",
               kernelInfo);

  solver->precon->approxBlockJacobiSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver3D.okl",
               "ellipticApproxBlockJacobiSolver3D",
               kernelInfo);

  solver->precon->exactBlockJacobiSolverKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver3D.okl",
               "ellipticExactBlockJacobiSolver3D",
               kernelInfo);

  //sizes for the coarsen and prolongation kernels. degree NFine to degree N
  int NpFine   = (Nf+1)*(Nf+2)*(Nf+3)/6;
  int NpCoarse = (Nc+1)*(Nc+2)*(Nf+3)/6;
  kernelInfo.addDefine("p_NpFine", NpFine);
  kernelInfo.addDefine("p_NpCoarse", NpCoarse);

  solver->precon->coarsenKernel =
    solver->mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconCoarsen.okl",
             "ellipticPreconCoarsen",
             kernelInfo);

  solver->precon->prolongateKernel =
    solver->mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconProlongate.okl",
             "ellipticPreconProlongate",
             kernelInfo);

  //on host gather-scatter
  mesh->hostGsh = gsParallelGatherScatterSetup(mesh->Nelements*mesh->Np, mesh->globalIds);

  // setup occa gather scatter
  mesh->ogs = meshParallelGatherScatterSetup(mesh,Ntotal,
                                             mesh->gatherLocalIds,
                                             mesh->gatherBaseIds,
                                             mesh->gatherBaseRanks,
                                             mesh->gatherHaloFlags);
  solver->o_invDegree = mesh->ogs->o_invDegree;

  // set up separate gather scatter infrastructure for halo and non halo nodes
  ellipticParallelGatherScatterSetup(solver);

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)  
  mesh->mapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (iint e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) mesh->mapB[n+e*mesh->Np] = 1E9;
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh->Nfp;n++) {
          int BCFlag = BCType[bc];
          int fid = mesh->faceNodes[n+f*mesh->Nfp];
          mesh->mapB[fid+e*mesh->Np] = BCFlag;
        }
      }
    }
  }
  gsParallelGatherScatter(mesh->hostGsh, mesh->mapB, "int", "min"); 

  mesh->Nmasked = 0;
  for (iint n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (mesh->mapB[n] == 1E9) {
      mesh->mapB[n] = 0.;
    } else if (mesh->mapB[n] == 1) { //Dirichlet boundary
      mesh->mapB[n] = 1.;
      mesh->Nmasked++;
    }
  }
  mesh->o_mapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), mesh->mapB);
  
  mesh->maskIds = (iint *) calloc(mesh->Nmasked, sizeof(iint));
  mesh->Nmasked =0; //reset
  for (iint n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (mesh->mapB[n] == 1) mesh->maskIds[mesh->Nmasked++] = n;
  }
  if (mesh->Nmasked) mesh->o_maskIds = mesh->device.malloc(mesh->Nmasked*sizeof(iint), mesh->maskIds);

  free(DrT); free(DsT); free(DtT); free(LIFTT);

  return solver;
}

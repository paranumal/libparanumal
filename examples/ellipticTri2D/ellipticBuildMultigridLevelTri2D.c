#include "ellipticTri2D.h"

// create solver and mesh structs for multigrid levels
solver_t *ellipticBuildMultigridLevelTri2D(solver_t *baseSolver, int* levelDegree, int n, const char *options){

  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));
  memcpy(solver,baseSolver,sizeof(solver_t));

  //populate the mini-mesh using the mesh struct
  mesh2D *mesh = (mesh2D*) calloc(1,sizeof(mesh2D));
  memcpy(mesh,baseSolver->mesh,sizeof(mesh2D));

  solver->mesh = mesh;

  int N = levelDegree[n];
  int NFine = levelDegree[n-1];
  int NpFine = (NFine+1)*(NFine+2)/2;
  int NpCoarse = (N+1)*(N+2)/2;

  // load reference (r,s) element nodes
  meshLoadReferenceNodesTri2D(mesh, N);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesTri2D(mesh);

  // create halo extension for x,y arrays
  iint totalHaloNodes = mesh->totalHaloPairs*mesh->Np;
  iint localNodes     = mesh->Nelements*mesh->Np;
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
  //meshParallelConnectNodes(mesh); //not needed for Ipdg, but possibly re-enable later

  //dont need these once vmap is made
  free(mesh->x);
  free(mesh->y);
  free(sendBuffer);

  // build Dr, Ds, LIFT transposes
  dfloat *DrT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *DsT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Np;++m){
      DrT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
      DsT[n+m*mesh->Np] = mesh->Ds[n*mesh->Np+m];
    }
  }

  dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      LIFTT[n+m*mesh->Np] = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
    }
  }

  //build element stiffness matrices
  if (mesh->Nverts == 3) {
    mesh->Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Ssr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {
        for (iint k=0;k<mesh->Np;k++) {
          for (iint l=0;l<mesh->Np;l++) {
            mesh->Srr[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
            mesh->Srs[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
            mesh->Ssr[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
            mesh->Sss[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
          }
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

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
      mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
      mesh->vmapP);

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
  kernelInfo.addDefine("p_G11ID", G11ID);
  kernelInfo.addDefine("p_GWJID", GWJID);


  kernelInfo.addDefine("p_RXID", RXID);
  kernelInfo.addDefine("p_SXID", SXID);

  kernelInfo.addDefine("p_RYID", RYID);
  kernelInfo.addDefine("p_SYID", SYID);

  kernelInfo.addDefine("p_JID", JID);
  kernelInfo.addDefine("p_JWID", JWID);


  iint Ntotal = mesh->Np*mesh->Nelements;
  iint Nblock = (Ntotal+blockSize-1)/blockSize;
  iint Nhalo = mesh->Np*mesh->totalHaloPairs;
  iint Nall   = Ntotal + Nhalo;

  solver->Nblock = Nblock;

  //add standard boundary functions
  char *boundaryHeaderFileName;
  boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTri2D/ellipticBoundary2D.h");
  kernelInfo.addInclude(boundaryHeaderFileName);

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  kernelInfo.addDefine("p_blockSize", blockSize);

  // add custom defines
  kernelInfo.addDefine("p_NpP", (mesh->Np+mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);

  //sizes for the coarsen and prolongation kernels. degree NFine to degree N
  kernelInfo.addDefine("p_NpFine", NpFine);
  kernelInfo.addDefine("p_NpCoarse", NpCoarse);

  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax);

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  int NblockP = 256/(4*mesh->Np); // get close to 256 threads
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
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTri2D.okl",
               "ellipticAxTri2D",
               kernelInfo);

  solver->partialAxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTri2D.okl",
               "ellipticPartialAxTri2D",
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

  #if USE_BERN
    solver->BRGradientVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBassiRebayBBTri2D.okl",
                 "ellipticBBBRGradientVolume2D",
                 kernelInfo);

    solver->BRGradientSurfaceKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBassiRebayBBTri2D.okl",
                 "ellipticBBBRGradientSurface2D",
                 kernelInfo);

    solver->BRDivergenceVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBassiRebayBBTri2D.okl",
                 "ellipticBBBRDivergenceVolume2D",
                 kernelInfo);

    solver->BRDivergenceSurfaceKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBassiRebayBBTri2D.okl",
                 "ellipticBBBRDivergenceSurface2D",
                 kernelInfo);
  #else
    solver->BRGradientVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBassiRebayTri2D.okl",
                 "ellipticBRGradientVolume2D",
                 kernelInfo);

    solver->BRGradientSurfaceKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBassiRebayTri2D.okl",
                 "ellipticBRGradientSurface2D",
                 kernelInfo);

    solver->BRDivergenceVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBassiRebayTri2D.okl",
                 "ellipticBRDivergenceVolume2D",
                 kernelInfo);

    solver->BRDivergenceSurfaceKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticBassiRebayTri2D.okl",
                 "ellipticBRDivergenceSurface2D",
                 kernelInfo);
  #endif

  //new precon struct
  solver->precon = (precon_t *) calloc(1,sizeof(precon_t));

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

  free(DrT); free(DsT); free(LIFTT);

  return solver;
}

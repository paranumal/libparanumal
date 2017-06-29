#include "ellipticTri2D.h"

// create mini-meshes for each polynomial degree
void ellipticBuildAmgLevelsTri2D(solver_t *solver){

  // build some levels
  mesh2D *mesh = solver->mesh;

  // start with all polynomial levels { just for draft }
  int Nlevels = mesh->N;
  mesh2D **meshLevels = (mesh2D**) calloc(Nlevels, sizeof(mesh2D*));
  solver2D **solverLevels = (mesh2D**) calloc(Nlevels, sizeof(mesh2D*));

  // info for kernel construction
  occa::kernelInfo kernelInfo;
  
  kernelInfo.addDefine("p_Nfaces", mesh->Nfaces);
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);

  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_Nvgeo", mesh->Nvgeo);
  kernelInfo.addDefine("p_Nsgeo", mesh->Nsgeo);
  kernelInfo.addDefine("p_Nggeo", mesh->Nggeo);

  kernelInfo.addDefine("p_NXID", NXID);
  kernelInfo.addDefine("p_NYID", NYID);
  kernelInfo.addDefine("p_SJID", SJID);
  kernelInfo.addDefine("p_IJID", IJID);
  kernelInfo.addDefine("p_WSJID", WSJID);
  kernelInfo.addDefine("p_IHID", IHID);

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
  
  for(int level=Nlevels-1;level>=1;--level){ // hard coded for all degrees at the moment
    
    // hard code degree for this level
    iint levelN = level; 
    
    mesh2D *meshL = meshLevels[level];
    solver2D *solverL = solverLevels[level];
    
    // copy from original mesh (important to capture geofacs, occa device, ... WILL BREAK FOR QUADS)
    memcpy(meshL, mesh, sizeof(mesh2D*));


    // reload for new degree
    meshLoadReferenceNodesTri2D(meshL, levelN);

    // set up halo exchange info for MPI (do before connect face nodes)
    meshHaloSetup(meshL);
    
    // connect face nodes (find trace indices)
    meshConnectFaceNodes2D(meshL);

    // global nodes ( do we need this ? )
    meshParallelConnectNodes(meshL);

    // ----------------------------------------------------------------------
    // specialized for matrix-free IPDG: DrT, DsT, LIFTT, MM
    // build Dr, Ds, LIFT transposes
    dfloat *DrT = (dfloat*) calloc(meshL->Np*meshL->Np, sizeof(dfloat));
    dfloat *DsT = (dfloat*) calloc(meshL->Np*meshL->Np, sizeof(dfloat));
    dfloat *LIFTT = (dfloat*) calloc(meshL->Np*meshL->Nfaces*meshL->Nfp, sizeof(dfloat));
    
    for(iint n=0;n<meshL->Np;++n){
      for(iint m=0;m<meshL->Np;++m){
	DrT[n+m*meshL->Np] = meshL->Dr[n*meshL->Np+m];
	DsT[n+m*meshL->Np] = meshL->Ds[n*meshL->Np+m];
      }
    }
    
    for(iint n=0;n<meshL->Np;++n){
      for(iint m=0;m<meshL->Nfaces*meshL->Nfp;++m){
	LIFTT[n+m*meshL->Np] = meshL->LIFT[n*meshL->Nfp*meshL->Nfaces+m];
      }
    }

    // only set up essentials on DEVICE
    meshL->o_Dr = meshL->device.malloc(meshL->Np*meshL->Np*sizeof(dfloat),  meshL->Dr);
    meshL->o_Ds = meshL->device.malloc(meshL->Np*meshL->Np*sizeof(dfloat),  meshL->Ds);
    
    meshL->o_DrT = meshL->device.malloc(meshL->Np*meshL->Np*sizeof(dfloat), DrT);
    meshL->o_DsT = meshL->device.malloc(meshL->Np*meshL->Np*sizeof(dfloat), DsT);
    
    meshL->o_LIFT  = meshL->device.malloc(meshL->Np*meshL->Nfaces*meshL->Nfp*sizeof(dfloat), meshL->LIFT);
    meshL->o_LIFTT = meshL->device.malloc(meshL->Np*meshL->Nfaces*meshL->Nfp*sizeof(dfloat), LIFTT);

    meshL->o_vmapM = meshL->device.malloc(meshL->Nelements*meshL->Nfp*meshL->Nfaces*sizeof(iint), meshL->vmapM);
    meshL->o_vmapP = meshL->device.malloc(meshL->Nelements*meshL->Nfp*meshL->Nfaces*sizeof(iint), meshL->vmapP);
    
    meshL->o_x = meshL->device.malloc(meshL->Nelements*meshL->Np*sizeof(dfloat), meshL->x);
    meshL->o_y = meshL->device.malloc(meshL->Nelements*meshL->Np*sizeof(dfloat), meshL->y);
    
    // ----------------------------------------------------------------------

    // build level specific kernels
    occa::kernelInfo kernelInfoL = kernelInfo;
    
    kernelInfoL.addDefine("p_N", meshL->N);
    kernelInfoL.addDefine("p_Nq", meshL->N+1);
    kernelInfoL.addDefine("p_Np", meshL->Np);
    kernelInfoL.addDefine("p_Nfp", meshL->Nfp);
    kernelInfoL.addDefine("p_NfacesNfp", meshL->Nfp*meshL->Nfaces);
    
    // add custom defines
    kernelInfoL.addDefine("p_NpP", (meshL->Np+meshL->Nfp*meshL->Nfaces));

    int Nmax = mymax(meshL->Np, meshL->Nfaces*meshL->Nfp);
    kernelInfoL.addDefine("p_Nmax", Nmax);
    
    int NblockV = 256/meshL->Np; // get close to 256 threads
    kernelInfoL.addDefine("p_NblockV", NblockV);
    
    int NblockP = 512/(4*meshL->Np); // get close to 256 threads
    kernelInfoL.addDefine("p_NblockP", NblockP);
    
    meshL->haloExtractKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
					  "meshHaloExtract2D",
					  kernelInfo);
    
    meshL->gatherKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/gather.okl",
					  "gather",
					  kernelInfo);
    
    meshL->scatterKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/scatter.okl",
					  "scatter",
					  kernelInfo);
    
    meshL->getKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/get.okl",
					  "get",
					  kernelInfo);
    
    meshL->putKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/put.okl",
					  "put",
				       kernelInfo);
    
    solverL->AxKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxTri2D.okl",
					  "ellipticAxTri2D",
					  kernelInfo);
    
    solverL->gradientKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientTri2D.okl",
					  "ellipticGradientTri2D",
					  kernelInfo);
    
    solverL->ipdgKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgTri2D.okl",
					  "ellipticAxIpdgTri2D",
					  kernelInfo);
    
    occaTimerTic(meshL->device,"PreconditionerSetup");
    solverL->precon = ellipticPreconditionerSetupTri2D(solver, solverL->ogs, tau, lambda, BCType,  options);
    occaTimerToc(meshL->device,"PreconditionerSetup");
    
    solverL->precon->preconKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/ellipticOasPreconTri2D.okl",
					  "ellipticOasPreconTri2D",
					  kernelInfo);
    
    solverL->precon->coarsenKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconCoarsen.okl",
					  "ellipticPreconCoarsen",
					  kernelInfo);
    
    solverL->precon->prolongateKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconProlongate.okl",
					  "ellipticPreconProlongate",
					  kernelInfo);
    
    
    solverL->precon->blockJacobiKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/ellipticBlockJacobiPreconTri2D.okl",
					  "ellipticBlockJacobiPreconTri2D",
					  kernelInfo);
    
    solverL->precon->patchSolverKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
					  "ellipticPatchSolver2D",
					  kernelInfo);
    
    solverL->precon->approxPatchSolverKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
					  "ellipticApproxPatchSolver2D",
					  kernelInfo);
    
    solverL->precon->localPatchSolverKernel =
      meshL->device.buildKernelFromSource(DHOLMES "/okl/ellipticPatchSolver2D.okl",
					  "ellipticLocalPatchSolver2D",
					  kernelInfo);

    free(DrT); free(DsT); free(LIFTT);
  }
  
}

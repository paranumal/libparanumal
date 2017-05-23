#include "ellipticQuad2D.h"

void ellipticComputeDegreeVector(mesh2D *mesh, iint Ntotal, ogs_t *ogs, dfloat *deg){

  // build degree vector
  for(iint n=0;n<Ntotal;++n)
    deg[n] = 1;

  occa::memory o_deg = mesh->device.malloc(Ntotal*sizeof(dfloat), deg);
  
  o_deg.copyFrom(deg);
  
  ellipticParallelGatherScatterQuad2D(mesh, ogs, o_deg, o_deg, dfloatString, "add");
  
  o_deg.copyTo(deg);

  mesh->device.finish();
  o_deg.free();
  
}

solver_t *ellipticSolveSetupQuad2D(mesh_t *mesh, dfloat lambda, occa::kernelInfo &kernelInfo, const char *options) {

	iint rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint Ntotal = mesh->Np*mesh->Nelements;
  iint NtotalP = mesh->NqP*mesh->NqP*mesh->Nelements;
  iint Nblock = (Ntotal+blockSize-1)/blockSize;
  iint Nhalo = mesh->Np*mesh->totalHaloPairs;
  iint Nall   = Ntotal + Nhalo;
  iint NallP  = NtotalP;
  
  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));

  solver->mesh = mesh;

  solver->p   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->z   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->zP  = (dfloat*) calloc(NallP,  sizeof(dfloat));
  solver->Ax   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->Ap  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));
  solver->grad = (dfloat*) calloc(Nall*4, sizeof(dfloat));
	
  // need to rename o_r, o_x to avoid confusion
  solver->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), solver->r);
  solver->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
  solver->o_zP  = mesh->device.malloc(NallP*sizeof(dfloat),solver->zP); // CAUTION
  solver->o_Ax  = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), solver->Ap);
  solver->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), solver->tmp);
  solver->o_grad= mesh->device.malloc(Nall*4*sizeof(dfloat), solver->grad);
  
  solver->sendBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(dfloat));
  solver->recvBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(dfloat));
 
  solver->Nblock = Nblock;

  kernelInfo.addDefine("p_NqP", (mesh->Nq+2));
  kernelInfo.addDefine("p_NpP", (mesh->NqP*mesh->NqP));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);

  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax); 

  int NblockV = 256/mesh->Np; // get close to 256 threads
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

  occaTimerTic(mesh->device,"GatherScatterSetup");
  // set up gslib MPI gather-scatter and OCCA gather/scatter arrays
  solver->ogs = meshParallelGatherScatterSetup(mesh,
					mesh->Np*mesh->Nelements,
					sizeof(dfloat),
					mesh->gatherLocalIds,
					mesh->gatherBaseIds, 
					mesh->gatherHaloFlags);
  occaTimerToc(mesh->device,"GatherScatterSetup");

  occaTimerTic(mesh->device,"PreconditionerSetup");
  solver->precon = ellipticPreconditionerSetupQuad2D(mesh, solver->ogs, lambda, options);

  
  solver->precon->preconKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticOasPreconQuad2D.okl",
				       "ellipticOasPreconQuad2D",
				       kernelInfo);
  
  solver->precon->restrictKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconRestrictQuad2D.okl",
				       "ellipticFooQuad2D",
				       kernelInfo);

  solver->precon->coarsenKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconCoarsen.okl",
				       "ellipticPreconCoarsen",
				       kernelInfo);

  solver->precon->prolongateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconProlongate.okl",
				       "ellipticPreconProlongate",
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

  //fill geometric factors in halo
  if(mesh->totalHaloPairs){
    iint Nlocal = mesh->Np*mesh->Nelements;
    iint Nhalo = mesh->totalHaloPairs*mesh->Np;
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
  }

  // build weights for continuous SEM L2 project --->
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
  ellipticParallelGatherScatterQuad2D(mesh, solver->ogs, o_localMM, o_MM, dfloatString, "add");

  mesh->o_projectL2 = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
  mesh->dotDivideKernel(Ntotal, o_localMM, o_MM, mesh->o_projectL2);

  free(localMM); o_MM.free(); o_localMM.free();

  //set matrix free function pointers
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
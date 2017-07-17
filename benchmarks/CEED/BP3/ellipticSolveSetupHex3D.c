#include "ellipticHex3D.h"

void ellipticComputeDegreeVector(mesh3D *mesh, iint Ntotal, ogs_t *ogs, dfloat *deg){

  // build degree vector
  for(iint n=0;n<Ntotal;++n)
    deg[n] = 1;

  occa::memory o_deg = mesh->device.malloc(Ntotal*sizeof(dfloat), deg);
  
  o_deg.copyFrom(deg);
  
  ellipticParallelGatherScatter(mesh, ogs, o_deg, o_deg, dfloatString, "add");
  
  o_deg.copyTo(deg);

  mesh->device.finish();
  o_deg.free();
  
}

solver_t *ellipticSolveSetupHex3D(mesh_t *mesh, dfloat lambda, occa::kernelInfo &kernelInfo, const char *options) {

	iint rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  int NblockV = 1024/mesh->Np; // works for CUDA
  int NblockS = 1024/maxNodes; // works for CUDA
  int NblockG;

  iint gNq = mesh->Nq+1;
  iint gNp = gNq*gNq*gNq;
  iint gNq2 = gNq*gNq;
  if(gNq2<=32) NblockG = ( 512/gNq2 );
  else NblockG = 1; 
  NblockG = 512/gNq2;

  iint Ntotal = mesh->Np*mesh->Nelements;
  iint NtotalP = mesh->NqP*mesh->NqP*mesh->NqP*mesh->Nelements;

  iint Nblock = (Ntotal+blockSize-1)/blockSize;
  printf("Nblock = %d\n", Nblock);


  iint Nhalo = mesh->Np*mesh->totalHaloPairs;
  iint Nall   = Ntotal + Nhalo;
  iint NallP  = NtotalP;

  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));

  solver->mesh = mesh;

  solver->p   = (dfloat*) calloc(Nall, sizeof(dfloat));
  solver->r   = (dfloat*) calloc(Nall, sizeof(dfloat));
  solver->z   = (dfloat*) calloc(Nall, sizeof(dfloat));
  solver->zP  = (dfloat*) calloc(NallP, sizeof(dfloat));
  solver->Ax  = (dfloat*) calloc(Nall, sizeof(dfloat));
  solver->Ap  = (dfloat*) calloc(Nall, sizeof(dfloat));
  printf("Nblkck = %d\n", Nblock);
  solver->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));
  solver->grad = (dfloat*) calloc(4*(Ntotal+Nhalo), sizeof(dfloat));
  
  solver->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), solver->r);
  solver->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
  solver->o_zP  = mesh->device.malloc(NallP*sizeof(dfloat), solver->zP); // CAUTION
  solver->o_Ax  = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), solver->Ap);
  solver->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), solver->tmp);
  solver->o_grad  = mesh->device.malloc(Nall*4*sizeof(dfloat), solver->grad);
  solver->o_pAp  = mesh->device.malloc(sizeof(dfloat));

  iint Nbytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  
#if 0
  solver->sendBuffer = (dfloat*) calloc(Nbytes/sizeof(dfloat), sizeof(dfloat));
  solver->recvBuffer = (dfloat*) calloc(Nbytes/sizeof(dfloat), sizeof(dfloat));
#else
  solver->defaultStream = mesh->device.getStream();
  solver->dataStream = mesh->device.createStream();
  mesh->device.setStream(solver->defaultStream);
  
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
#endif
  solver->Nblock = Nblock;

  // BP3 specific stuff starts here
  dfloat *gD = (dfloat*) calloc(gNq*mesh->Nq, sizeof(dfloat));
  dfloat *gI = (dfloat*) calloc(gNq*mesh->Nq, sizeof(dfloat));
  dfloat *gggeo = (dfloat*) calloc(gNp*mesh->Nelements*mesh->Nggeo, sizeof(dfloat));

  srand48(32);
  for(iint n=0;n<gNq*mesh->Nq;++n){
    gD[n] = drand48();
    gI[n] = drand48();
  }

  for(iint n=0;n<gNp*mesh->Nelements*mesh->Nggeo;++n){
    // gggeo[n] = drand48();
    gggeo[n] = .1;
  }
  
  solver->o_gD = mesh->device.malloc(gNq*mesh->Nq*sizeof(dfloat), gD);
  solver->o_gI = mesh->device.malloc(gNq*mesh->Nq*sizeof(dfloat), gI);
  solver->o_gggeo = mesh->device.malloc(mesh->Nggeo*gNp*mesh->Nelements*sizeof(dfloat), gggeo);
  // BP3 specific stuff ends here 

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");
  //  kernelInfo.addCompilerFlag("-G");


  kernelInfo.addDefine("p_maxNodes", maxNodes);
  kernelInfo.addDefine("p_Nmax", maxNodes);

  kernelInfo.addDefine("p_NblockV", NblockV);
  kernelInfo.addDefine("p_NblockS", NblockS);

  kernelInfo.addDefine("p_NblockG", NblockG);
  printf("NblockG = %d\n", NblockG);

  kernelInfo.addDefine("p_Lambda2", 0.5f);
  kernelInfo.addDefine("p_NqP", (mesh->Nq+2));
  kernelInfo.addDefine("p_NpP", (mesh->NqP*mesh->NqP*mesh->NqP));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);

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

  solver->AxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxHex3D.okl",
               "ellipticAxHex3D_e3",
               kernelInfo);

  solver->partialAxKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxHex3D.okl",
				       "ellipticPartialAxHex3D_e3",
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

  solver->gradientKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientHex3D.okl",
				       "ellipticGradientHex3D",
				       kernelInfo);
  
  solver->partialGradientKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticGradientHex3D.okl",
				       "ellipticPartialGradientHex3D",
				       kernelInfo);


  solver->ipdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgHex3D.okl",
               "ellipticAxIpdgHex3D",
               kernelInfo);

  solver->partialIpdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAxIpdgHex3D.okl",
				       "ellipticPartialAxIpdgHex3D",
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
    iint Nlocal = mesh->Nelements*mesh->Np;
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

  ellipticParallelGatherScatter(mesh, solver->ogs, o_localMM, o_MM, dfloatString, "add");

  mesh->o_projectL2 = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
  mesh->dotDivideKernel(Ntotal, o_localMM, o_MM, mesh->o_projectL2);

  free(localMM); o_MM.free(); o_localMM.free();

  if(rank==0)
    printf("starting elliptic parallel gather scatter setup\n");

  // set up separate gather scatter infrastructure for halo and non halo nodes
  //  mesh->device.setStream(solver->dataStream);
  ellipticParallelGatherScatterSetup(mesh,
				     mesh->Np*mesh->Nelements,
				     sizeof(dfloat),
				     mesh->gatherLocalIds,
				     mesh->gatherBaseIds, 
				     mesh->gatherHaloFlags,
				     &(solver->halo),
				     &(solver->nonHalo));
  //  mesh->device.setStream(solver->defaultStream);

  
  // count elements that contribute to global C0 gather-scatter
  iint globalCount = 0;
  iint localCount = 0;
  iint *localHaloFlags = (iint*) calloc(mesh->Np*mesh->Nelements, sizeof(int));

  for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
    localHaloFlags[mesh->gatherLocalIds[n]] += mesh->gatherHaloFlags[n];
  }
  
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
  
  printf("local = %d, global = %d\n", localCount, globalCount);
  
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

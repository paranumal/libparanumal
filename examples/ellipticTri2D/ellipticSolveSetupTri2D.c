#include "ellipticTri2D.h"

solver_t *ellipticSolveSetupTri2D(mesh_t *mesh, dfloat lambda, occa::kernelInfo &kernelInfo){

  iint Ntotal = mesh->Np*mesh->Nelements;
  iint NtotalP = mesh->NpP*mesh->Nelements;
  iint Nblock = (Ntotal+blockSize-1)/blockSize;
  iint Nhalo = mesh->Np*mesh->totalHaloPairs;
  iint Nall   = Ntotal + Nhalo;
  iint NallP  = NtotalP;
  
  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));

  solver->mesh = mesh;

  solver->p   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->z   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->zP  = (dfloat*) calloc(NallP,  sizeof(dfloat));
  solver->Ax  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->Ap  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));

  solver->grad = (dfloat*) calloc(Nall*4, sizeof(dfloat));

  solver->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
  solver->o_zP  = mesh->device.malloc(NallP*sizeof(dfloat),solver->zP); // CAUTION
  solver->o_Ax  = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), solver->Ap);
  solver->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), solver->tmp);

  solver->o_grad  = mesh->device.malloc(Nall*4*sizeof(dfloat), solver->grad);
  
  // use this for OAS precon pairwise halo exchange
  solver->sendBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(dfloat));
  solver->recvBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(dfloat));

  solver->type = strdup(dfloatString);

  solver->Nblock = Nblock;

  // add custom defines
  kernelInfo.addDefine("p_NpP", (mesh->Np+mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);

  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax); 

  int NblockV = 256/mesh->Np; // get close to 256 threads
  kernelInfo.addDefine("p_NblockV", NblockV);
  
  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource("okl/meshHaloExtract2D.okl",
				       "meshHaloExtract2D",
				       kernelInfo);

  mesh->gatherKernel =
    mesh->device.buildKernelFromSource("okl/gather.okl",
				       "gather",
				       kernelInfo);

  mesh->scatterKernel =
    mesh->device.buildKernelFromSource("okl/scatter.okl",
				       "scatter",
				       kernelInfo);

  mesh->getKernel =
    mesh->device.buildKernelFromSource("okl/get.okl",
				       "get",
				       kernelInfo);

  mesh->putKernel =
    mesh->device.buildKernelFromSource("okl/put.okl",
				       "put",
				       kernelInfo);

  mesh->weightedInnerProduct1Kernel =
    mesh->device.buildKernelFromSource("okl/weightedInnerProduct1.okl",
				       "weightedInnerProduct1",
				       kernelInfo);

  mesh->weightedInnerProduct2Kernel =
    mesh->device.buildKernelFromSource("okl/weightedInnerProduct2.okl",
				       "weightedInnerProduct2",
				       kernelInfo);

  mesh->innerProductKernel =
    mesh->device.buildKernelFromSource("okl/innerProduct.okl",
				       "innerProduct",
				       kernelInfo);
  
  mesh->scaledAddKernel =
      mesh->device.buildKernelFromSource("okl/scaledAdd.okl",
					 "scaledAdd",
					 kernelInfo);

  mesh->dotMultiplyKernel =
      mesh->device.buildKernelFromSource("okl/dotMultiply.okl",
					 "dotMultiply",
					 kernelInfo);

  mesh->dotDivideKernel = 
      mesh->device.buildKernelFromSource("okl/dotDivide.okl",
					 "dotDivide",
					 kernelInfo);


  mesh->gradientKernel = 
    mesh->device.buildKernelFromSource("okl/ellipticGradientTri2D.okl",
				       "ellipticGradientTri2D",
					 kernelInfo);


  mesh->ipdgKernel =
    mesh->device.buildKernelFromSource("okl/ellipticAxIpdgTri2D.okl",
				       "ellipticAxIpdgTri2D",
				       kernelInfo);  

  
  // set up gslib MPI gather-scatter and OCCA gather/scatter arrays
  solver->ogs = meshParallelGatherScatterSetup(mesh,
					       mesh->Np*mesh->Nelements,
					       sizeof(dfloat),
					       mesh->gatherLocalIds,
					       mesh->gatherBaseIds, 
					       mesh->gatherHaloFlags);
  
  
  solver->precon = ellipticPreconditionerSetupTri2D(mesh, solver->ogs, lambda);
  
  solver->precon->preconKernel = 
    mesh->device.buildKernelFromSource("okl/ellipticOasPreconTri2D.okl",
				       "ellipticOasPreconTri2D",
				       kernelInfo);
  
  solver->precon->restrictKernel =
    mesh->device.buildKernelFromSource("okl/ellipticPreconRestrictTri2D.okl",
				       "ellipticFooTri2D",
				       kernelInfo);

  solver->precon->coarsenKernel =
    mesh->device.buildKernelFromSource("okl/ellipticPreconCoarsen.okl",
				       "ellipticPreconCoarsen",
				       kernelInfo);

  solver->precon->prolongateKernel =
    mesh->device.buildKernelFromSource("okl/ellipticPreconProlongate.okl",
				       "ellipticPreconProlongate",
				       kernelInfo);

  // probably should relocate this
  // build weights for continuous SEM L2 project --->                                                        
  dfloat *localMM = (dfloat*) calloc(Ntotal, sizeof(dfloat));

  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat J = mesh->vgeo[e*mesh->Nvgeo + JID];
      localMM[n+e*mesh->Np] = J;
    }
  }

  occa::memory o_localMM = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
  occa::memory o_MM      = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);

  // sum up all contributions at base nodes and scatter back                                                 
  ellipticParallelGatherScatterTri2D(mesh, solver->ogs, o_localMM, o_MM, dfloatString, "add");

  mesh->o_projectL2 = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
  mesh->dotDivideKernel(Ntotal, o_localMM, o_MM, mesh->o_projectL2);

  free(localMM); o_MM.free(); o_localMM.free();
  // <------                                            

  
  return solver;
}

/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "adaptive.h"


void adaptiveSolveSetup(adaptive_t *adaptive, dfloat lambda, occa::properties &kernelInfo){

  mesh_t *mesh = adaptive->mesh;
  setupAide options = adaptive->options;

  //sanity checking

  if (options.compareArgs("PRECONDITIONER","MASSMATRIX")){
    printf("ERROR: MASSMATRIX preconditioner is only available for triangle and tetrhedra elements. Use JACOBI instead.\n");
    MPI_Finalize();
    exit(-1);
  }
  if (options.compareArgs("PRECONDITIONER","MASSMATRIX") && lambda==0) {
    printf("ERROR: MASSMATRIX preconditioner is unavailble when lambda=0. \n");
    MPI_Finalize();
    exit(-1);
  }

  dlong Ntotal = mesh->Np*mesh->Nelements;

  dlong Nhalo = mesh->Np*mesh->totalHaloPairs;
  dlong Nall   = Ntotal + Nhalo;

  dlong Nblock  = mymax(1,(Ntotal+blockSize-1)/blockSize);
  dlong Nblock2 = mymax(1,(Nblock+blockSize-1)/blockSize);

  dlong NthreadsUpdatePCG = 256;
  dlong NblocksUpdatePCG = mymin((Ntotal+NthreadsUpdatePCG-1)/NthreadsUpdatePCG, 160);
 
  adaptive->NthreadsUpdatePCG = NthreadsUpdatePCG;
  adaptive->NblocksUpdatePCG = NblocksUpdatePCG;
  
  //tau
  adaptive->tau = 2.0*(mesh->N+1)*(mesh->N+3);

  adaptive->p   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  adaptive->z   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  adaptive->Ax  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  adaptive->Ap  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  adaptive->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));

  adaptive->grad = (dfloat*) calloc(Nall*4, sizeof(dfloat));

  adaptive->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), adaptive->p);
  adaptive->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), adaptive->p);
  adaptive->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), adaptive->z);

  adaptive->o_res = mesh->device.malloc(Nall*sizeof(dfloat), adaptive->z);
  adaptive->o_Sres = mesh->device.malloc(Nall*sizeof(dfloat), adaptive->z);
  adaptive->o_Ax  = mesh->device.malloc(Nall*sizeof(dfloat), adaptive->p);
  adaptive->o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), adaptive->Ap);
  adaptive->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), adaptive->tmp);
  adaptive->o_tmp2 = mesh->device.malloc(Nblock2*sizeof(dfloat), adaptive->tmp);

  adaptive->tmpNormr = (dfloat*) calloc(adaptive->NblocksUpdatePCG,sizeof(dfloat));
  adaptive->o_tmpNormr = mesh->device.malloc(adaptive->NblocksUpdatePCG*sizeof(dfloat), adaptive->tmpNormr);

  
  adaptive->o_grad  = mesh->device.malloc(Nall*4*sizeof(dfloat), adaptive->grad);

  int useFlexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");
  
  if(options.compareArgs("KRYLOV SOLVER", "NONBLOCKING")){
    
    int Nwork = (useFlexible) ? 9: 6; 

    adaptive->o_pcgWork = new occa::memory[Nwork];
    for(int n=0;n<Nwork;++n){
      adaptive->o_pcgWork[n]  = mesh->device.malloc(Nall*sizeof(dfloat), adaptive->z);
    }

    adaptive->tmppdots = (dfloat*) calloc(adaptive->NblocksUpdatePCG,sizeof(dfloat));
    adaptive->o_tmppdots = mesh->device.malloc(adaptive->NblocksUpdatePCG*sizeof(dfloat), adaptive->tmppdots);
    
    adaptive->tmprdotz = (dfloat*) calloc(adaptive->NblocksUpdatePCG,sizeof(dfloat));
    adaptive->o_tmprdotz = mesh->device.malloc(adaptive->NblocksUpdatePCG*sizeof(dfloat), adaptive->tmprdotz);
    
    adaptive->tmpzdotz = (dfloat*) calloc(adaptive->NblocksUpdatePCG,sizeof(dfloat));
    adaptive->o_tmpzdotz = mesh->device.malloc(adaptive->NblocksUpdatePCG*sizeof(dfloat), adaptive->tmpzdotz);
    
    adaptive->tmprdotr = (dfloat*) calloc(adaptive->NblocksUpdatePCG,sizeof(dfloat));
    adaptive->o_tmprdotr = mesh->device.malloc(adaptive->NblocksUpdatePCG*sizeof(dfloat), adaptive->tmprdotr);

    adaptive->tmpudotr = (dfloat*) calloc(adaptive->NblocksUpdatePCG,sizeof(dfloat));
    adaptive->o_tmpudotr = mesh->device.malloc(adaptive->NblocksUpdatePCG*sizeof(dfloat), adaptive->tmpudotr);
    
    adaptive->tmpudots = (dfloat*) calloc(adaptive->NblocksUpdatePCG,sizeof(dfloat));
    adaptive->o_tmpudots = mesh->device.malloc(adaptive->NblocksUpdatePCG*sizeof(dfloat), adaptive->tmpudots);
    
    adaptive->tmpudotw = (dfloat*) calloc(adaptive->NblocksUpdatePCG,sizeof(dfloat));
    adaptive->o_tmpudotw = mesh->device.malloc(adaptive->NblocksUpdatePCG*sizeof(dfloat), adaptive->tmpudotw);
  }
  
  //setup async halo stream
  adaptive->defaultStream = mesh->defaultStream;
  adaptive->dataStream = mesh->dataStream;

  dlong Nbytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  if(Nbytes>0){
    adaptive->sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, Nbytes, NULL, adaptive->o_sendBuffer, adaptive->h_sendBuffer);
    adaptive->recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, Nbytes, NULL, adaptive->o_recvBuffer, adaptive->h_recvBuffer);
    adaptive->gradSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, 2*Nbytes, NULL, adaptive->o_gradSendBuffer, adaptive->h_gradSendBuffer);
    adaptive->gradRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, 2*Nbytes, NULL, adaptive->o_gradRecvBuffer, adaptive->h_gradRecvBuffer);

  }else{
    adaptive->sendBuffer = NULL;
    adaptive->recvBuffer = NULL;
  }
  mesh->device.setStream(adaptive->defaultStream);

  adaptive->type = strdup(dfloatString);

  adaptive->Nblock = Nblock;
  adaptive->Nblock2 = Nblock2;

  //fill geometric factors in halo
  if(mesh->totalHaloPairs){
    dlong Nlocal = mesh->Nelements;
    dlong Nhalo = mesh->totalHaloPairs;
    size_t Nbytes = mesh->Nvgeo*sizeof(dfloat);

    if (adaptive->elementType==QUADRILATERALS || adaptive->elementType==HEXAHEDRA) {
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

  // count total number of elements
  hlong NelementsLocal = mesh->Nelements;
  hlong NelementsGlobal = 0;

  MPI_Allreduce(&NelementsLocal, &NelementsGlobal, 1, MPI_HLONG, MPI_SUM, mesh->comm);

  adaptive->NelementsGlobal = NelementsGlobal;

  //check all the bounaries for a Dirichlet
  bool allNeumann = (lambda==0) ? true :false;
  adaptive->allNeumannPenalty = 1.;
  hlong localElements = (hlong) mesh->Nelements;
  hlong totalElements = 0;
  MPI_Allreduce(&localElements, &totalElements, 1, MPI_HLONG, MPI_SUM, mesh->comm);
  adaptive->allNeumannScale = 1./sqrt((dfloat)mesh->Np*totalElements);

  adaptive->EToB = (int *) calloc(mesh->Nelements*mesh->Nfaces,sizeof(int));
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if (bc>0) {
        int BC = adaptive->BCType[bc]; //get the type of boundary
        adaptive->EToB[e*mesh->Nfaces+f] = BC; //record it
        if (BC!=2) allNeumann = false; //check if its a Dirchlet
      }
    }
  }


  // !!!!!! Removed MPI::BOOL since some mpi versions complains about it !!!!!
  int lallNeumann, gallNeumann;
  lallNeumann = allNeumann ? 0:1;
  MPI_Allreduce(&lallNeumann, &gallNeumann, 1, MPI_INT, MPI_SUM, mesh->comm);
  adaptive->allNeumann = (gallNeumann>0) ? false: true;

  // MPI_Allreduce(&allNeumann, &(adaptive->allNeumann), 1, MPI::BOOL, MPI_LAND, mesh->comm);
  //  if (mesh->rank==0&& options.compareArgs("VERBOSE","TRUE"))
    printf("allNeumann = %d \n", adaptive->allNeumann);

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
  adaptive->o_EToB = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int), adaptive->EToB);

#if 0
  if (mesh->rank==0 && options.compareArgs("VERBOSE","TRUE"))
    occa::setVerboseCompilation(true);
  else
    occa::setVerboseCompilation(false);
#endif

  //setup an unmasked gs handle
  int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;
  meshParallelGatherScatterSetup(mesh, Ntotal, mesh->globalIds, mesh->comm, verbose);

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  adaptive->mapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  int largeNumber = 1<<20;
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) adaptive->mapB[n+e*mesh->Np] = largeNumber;
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh->Nfp;n++) {
          int BCFlag = adaptive->BCType[bc];
          int fid = mesh->faceNodes[n+f*mesh->Nfp];
          adaptive->mapB[fid+e*mesh->Np] = mymin(BCFlag,adaptive->mapB[fid+e*mesh->Np]);
        }
      }
    }
  }
  ogsGatherScatter(adaptive->mapB, ogsInt, ogsMin, mesh->ogs);

  //use the bc flags to find masked ids
  adaptive->Nmasked = 0;
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (adaptive->mapB[n] == largeNumber) {
      adaptive->mapB[n] = 0.;
    } else if (adaptive->mapB[n] == 1) { //Dirichlet boundary
      adaptive->Nmasked++;
    }
  }

  
  adaptive->o_mapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), adaptive->mapB);

  adaptive->maskIds = (dlong *) calloc(adaptive->Nmasked, sizeof(dlong));
  adaptive->Nmasked =0; //reset
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (adaptive->mapB[n] == 1) adaptive->maskIds[adaptive->Nmasked++] = n;
  }
  if (adaptive->Nmasked) adaptive->o_maskIds = mesh->device.malloc(adaptive->Nmasked*sizeof(dlong), adaptive->maskIds);

  //make a masked version of the global id numbering
  mesh->maskedGlobalIds = (hlong *) calloc(Ntotal,sizeof(hlong));
  memcpy(mesh->maskedGlobalIds, mesh->globalIds, Ntotal*sizeof(hlong));
  for (dlong n=0;n<adaptive->Nmasked;n++)
    mesh->maskedGlobalIds[adaptive->maskIds[n]] = 0;

  //use the masked ids to make another gs handle
  adaptive->ogs = ogsSetup(Ntotal, mesh->maskedGlobalIds, mesh->comm, verbose, mesh->device);
  adaptive->o_invDegree = adaptive->ogs->o_invDegree;

  // count number of actual DOFS
  // TW: HERE

  /*preconditioner setup */
  //  adaptive->precon = new precon_t[1];
  adaptive->precon = (precon_t*) calloc(1, sizeof(precon_t));

  //  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo["compiler_flags"] += "-Xptxas -dlcm=ca";
  }

  if(mesh->device.mode()=="Serial")
    kernelInfo["compiler_flags"] += "-g";

  // set kernel name suffix
  char *suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  for (int r=0;r<2;r++){
    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {      

      //mesh kernels
      mesh->haloExtractKernel =
        mesh->device.buildKernel(DHOLMES "/okl/meshHaloExtract2D.okl",
                                       "meshHaloExtract2D",
                                       kernelInfo);

      mesh->addScalarKernel =
        mesh->device.buildKernel(DHOLMES "/okl/addScalar.okl",
                   "addScalar",
                   kernelInfo);

      mesh->maskKernel =
        mesh->device.buildKernel(DHOLMES "/okl/mask.okl",
                   "mask",
                   kernelInfo);


      kernelInfo["defines/" "p_blockSize"]= blockSize;


      mesh->sumKernel =
        mesh->device.buildKernel(DHOLMES "/okl/sum.okl",
                   "sum",
                   kernelInfo);

      adaptive->weightedInnerProduct1Kernel =
        mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct1.okl",
                                       "weightedInnerProduct1",
                                       kernelInfo);

      adaptive->weightedInnerProduct2Kernel =
        mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct2.okl",
                                       "weightedInnerProduct2",
                                       kernelInfo);

      adaptive->innerProductKernel =
        mesh->device.buildKernel(DHOLMES "/okl/innerProduct.okl",
                                       "innerProduct",
                                       kernelInfo);

      adaptive->weightedNorm2Kernel =
        mesh->device.buildKernel(DHOLMES "/okl/weightedNorm2.okl",
                                           "weightedNorm2",
                                           kernelInfo);

      adaptive->norm2Kernel =
        mesh->device.buildKernel(DHOLMES "/okl/norm2.okl",
                                           "norm2",
                                           kernelInfo);

      adaptive->scaledAddKernel =
          mesh->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
                                         "scaledAdd",
                                         kernelInfo);

      adaptive->dotMultiplyKernel =
          mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
                                         "dotMultiply",
                                         kernelInfo);

      adaptive->dotMultiplyAddKernel =
          mesh->device.buildKernel(DHOLMES "/okl/dotMultiplyAdd.okl",
                                         "dotMultiplyAdd",
                                         kernelInfo);

      adaptive->dotDivideKernel =
          mesh->device.buildKernel(DHOLMES "/okl/dotDivide.okl",
                                         "dotDivide",
                                         kernelInfo);

      // add custom defines
      kernelInfo["defines/" "p_NpP"]= (mesh->Np+mesh->Nfp*mesh->Nfaces);
      kernelInfo["defines/" "p_Nverts"]= mesh->Nverts;

      //sizes for the coarsen and prolongation kernels. degree N to degree 1
      kernelInfo["defines/" "p_NpFine"]= mesh->Np;
      kernelInfo["defines/" "p_NpCoarse"]= mesh->Nverts;


      if (adaptive->elementType==QUADRILATERALS || adaptive->elementType==HEXAHEDRA) {
        kernelInfo["defines/" "p_NqFine"]= mesh->N+1;
        kernelInfo["defines/" "p_NqCoarse"]= 2;
      }

      kernelInfo["defines/" "p_NpFEM"]= mesh->NpFEM;

      int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
      kernelInfo["defines/" "p_Nmax"]= Nmax;

      int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
      kernelInfo["defines/" "p_maxNodes"]= maxNodes;

      int NblockV = mymax(1,maxNthreads/mesh->Np); // works for CUDA
      int NnodesV = 1; //hard coded for now
      kernelInfo["defines/" "p_NblockV"]= NblockV;
      kernelInfo["defines/" "p_NnodesV"]= NnodesV;
      kernelInfo["defines/" "p_NblockVFine"]= NblockV;
      kernelInfo["defines/" "p_NblockVCoarse"]= NblockV;

      int NblockS = mymax(1,maxNthreads/maxNodes); // works for CUDA
      kernelInfo["defines/" "p_NblockS"]= NblockS;

      int NblockP = mymax(1,maxNthreads/(4*mesh->Np)); // get close to maxNthreads threads
      kernelInfo["defines/" "p_NblockP"]= NblockP;

      int NblockG;
      if(mesh->Np<=32) NblockG = ( 32/mesh->Np );
      else NblockG = maxNthreads/mesh->Np;
      kernelInfo["defines/" "p_NblockG"]= NblockG;

      kernelInfo["defines/" "p_halfC"]= (int)((mesh->cubNq+1)/2);
      kernelInfo["defines/" "p_halfN"]= (int)((mesh->Nq+1)/2);

      kernelInfo["defines/" "p_NthreadsUpdatePCG"] = (int) NthreadsUpdatePCG; // WARNING SHOULD BE MULTIPLE OF 32
      kernelInfo["defines/" "p_NwarpsUpdatePCG"] = (int) (NthreadsUpdatePCG/32); // WARNING: CUDA SPECIFIC
      
      //      cout << kernelInfo ;
      
      //add standard boundary functions
      char *boundaryHeaderFileName;
      if (adaptive->dim==2)
        boundaryHeaderFileName = strdup(DADAPTIVE "/data/adaptiveBoundary2D.h");
      else if (adaptive->dim==3)
        boundaryHeaderFileName = strdup(DADAPTIVE "/data/adaptiveBoundary3D.h");
      kernelInfo["includes"] += boundaryHeaderFileName;


      sprintf(fileName,  DADAPTIVE "/okl/adaptiveAx%s.okl", suffix);
      sprintf(kernelName, "adaptiveAx%s", suffix);

      occa::properties dfloatKernelInfo = kernelInfo;
      occa::properties floatKernelInfo = kernelInfo;
      floatKernelInfo["defines/" "pfloat"]= "float";
      dfloatKernelInfo["defines/" "pfloat"]= dfloatString;
      
      adaptive->AxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);

      if(adaptive->elementType!=HEXAHEDRA){
        sprintf(kernelName, "adaptivePartialAx%s", suffix);
      }
      else{
        if(adaptive->options.compareArgs("ELEMENT MAP", "TRILINEAR")){
          sprintf(kernelName, "adaptivePartialAxTrilinear%s", suffix);
        }else{
          sprintf(kernelName, "adaptivePartialAx%s", suffix);
        }
      }

      adaptive->partialAxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);
      adaptive->partialAxKernel2 = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);
      adaptive->partialFloatAxKernel = mesh->device.buildKernel(fileName,kernelName,floatKernelInfo);
      
      // only for Hex3D - cubature Ax
      if(adaptive->elementType==HEXAHEDRA){
	//	printf("BUILDING partialCubatureAxKernel\n");
	sprintf(fileName,  DADAPTIVE "/okl/adaptiveCubatureAx%s.okl", suffix);

	sprintf(kernelName, "adaptiveCubaturePartialAx%s", suffix);
	adaptive->partialCubatureAxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);
      }

      // combined PCG update and r.r kernel
      
      adaptive->updatePCGKernel =
	mesh->device.buildKernel(DADAPTIVE "/okl/adaptiveUpdatePCG.okl",
				 "adaptiveUpdatePCG", dfloatKernelInfo);


      // combined update for Non-blocking PCG
      adaptive->update1NBPCGKernel =
	mesh->device.buildKernel(DADAPTIVE "/okl/adaptiveUpdateNBPCG.okl",
				 "adaptiveUpdate1NBPCG", dfloatKernelInfo);

      adaptive->update2NBPCGKernel =
	mesh->device.buildKernel(DADAPTIVE "/okl/adaptiveUpdateNBPCG.okl",
				 "adaptiveUpdate2NBPCG", dfloatKernelInfo);


      // combined update for Non-blocking flexible PCG
      adaptive->update0NBFPCGKernel =
	mesh->device.buildKernel(DADAPTIVE "/okl/adaptiveUpdateNBFPCG.okl",
				 "adaptiveUpdate0NBFPCG", dfloatKernelInfo);

      adaptive->update1NBFPCGKernel =
	mesh->device.buildKernel(DADAPTIVE "/okl/adaptiveUpdateNBFPCG.okl",
				 "adaptiveUpdate1NBFPCG", dfloatKernelInfo);
      
      
      // Not implemented for Quad3D !!!!!
      if (options.compareArgs("BASIS","BERN")) {

        sprintf(fileName, DADAPTIVE "/okl/adaptiveGradientBB%s.okl", suffix);
        sprintf(kernelName, "adaptiveGradientBB%s", suffix);

        adaptive->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "adaptivePartialGradientBB%s", suffix);
        adaptive->partialGradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(fileName, DADAPTIVE "/okl/adaptiveAxIpdgBB%s.okl", suffix);
        sprintf(kernelName, "adaptiveAxIpdgBB%s", suffix);
        adaptive->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "adaptivePartialAxIpdgBB%s", suffix);
        adaptive->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      } else if (options.compareArgs("BASIS","NODAL")) {

        sprintf(fileName, DADAPTIVE "/okl/adaptiveGradient%s.okl", suffix);
        sprintf(kernelName, "adaptiveGradient%s", suffix);

        adaptive->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "adaptivePartialGradient%s", suffix);
        adaptive->partialGradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(fileName, DADAPTIVE "/okl/adaptiveAxIpdg%s.okl", suffix);
        sprintf(kernelName, "adaptiveAxIpdg%s", suffix);
        adaptive->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "adaptivePartialAxIpdg%s", suffix);
        adaptive->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
      }

      // Use the same kernel with quads for the following kenels
      sprintf(fileName, DADAPTIVE "/okl/adaptivePreconCoarsen%s.okl", suffix);
      sprintf(kernelName, "adaptivePreconCoarsen%s", suffix);
      adaptive->precon->coarsenKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(fileName, DADAPTIVE "/okl/adaptivePreconProlongate%s.okl", suffix);
      sprintf(kernelName, "adaptivePreconProlongate%s", suffix);
      adaptive->precon->prolongateKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      

      sprintf(fileName, DADAPTIVE "/okl/adaptiveBlockJacobiPrecon.okl");
      sprintf(kernelName, "adaptiveBlockJacobiPrecon");
      adaptive->precon->blockJacobiKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(kernelName, "adaptivePartialBlockJacobiPrecon");
      adaptive->precon->partialblockJacobiKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

    }

    MPI_Barrier(mesh->comm);
  }

  // TW: WARNING C0 appropriate only
  mesh->sumKernel(mesh->Nelements*mesh->Np, adaptive->o_invDegree, adaptive->o_tmp);
  adaptive->o_tmp.copyTo(adaptive->tmp);

  dfloat nullProjectWeightLocal = 0;
  dfloat nullProjectWeightGlobal = 0;
  for(dlong n=0;n<adaptive->Nblock;++n)
    nullProjectWeightLocal += adaptive->tmp[n];

  MPI_Allreduce(&nullProjectWeightLocal, &nullProjectWeightGlobal, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
  
  adaptive->nullProjectWeightGlobal = 1./nullProjectWeightGlobal;
  

  
  long long int pre = mesh->device.memoryAllocated();

  adaptivePreconditionerSetup(adaptive, adaptive->ogs, lambda, kernelInfo);

  long long int usedBytes = mesh->device.memoryAllocated()-pre;

  adaptive->precon->preconBytes = usedBytes;

}

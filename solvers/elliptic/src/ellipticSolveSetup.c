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

#include "elliptic.h"


void ellipticSolveSetup(elliptic_t *elliptic, dfloat lambda, occa::properties &kernelInfo){

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
    MPI_Finalize();
    exit(-1);
  }

  dlong Ntotal = mesh->Np*mesh->Nelements;
  dlong Nblock = mymax(1,(Ntotal+blockSize-1)/blockSize);
  dlong Nhalo = mesh->Np*mesh->totalHaloPairs;
  dlong Nall   = Ntotal + Nhalo;

  dlong Nblock2 = mymax(1,(Nblock+blockSize-1)/blockSize);

  //tau
  if (elliptic->elementType==TRIANGLES || 
      elliptic->elementType==QUADRILATERALS || 
      (elliptic->dim==3 && elliptic->elementType==QUADRILATERALS))
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
  elliptic->o_tmp2 = mesh->device.malloc(Nblock2*sizeof(dfloat), elliptic->tmp);

  elliptic->o_grad  = mesh->device.malloc(Nall*4*sizeof(dfloat), elliptic->grad);

  //setup async halo stream
  elliptic->defaultStream = mesh->defaultStream;
  elliptic->dataStream = mesh->dataStream;

  dlong Nbytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  if(Nbytes>0){
#if 0
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(Nbytes, NULL);

    elliptic->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    elliptic->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();


    occa::memory o_gradSendBuffer = mesh->device.mappedAlloc(2*Nbytes, NULL);
    occa::memory o_gradRecvBuffer = mesh->device.mappedAlloc(2*Nbytes, NULL);

    elliptic->gradSendBuffer = (dfloat*) o_gradSendBuffer.getMappedPointer();
    elliptic->gradRecvBuffer = (dfloat*) o_gradRecvBuffer.getMappedPointer();
#endif

    elliptic->sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, Nbytes, NULL, elliptic->o_sendBuffer);
    elliptic->recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, Nbytes, NULL, elliptic->o_recvBuffer);
    elliptic->gradSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, 2*Nbytes, NULL, elliptic->o_gradSendBuffer);
    elliptic->gradRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, 2*Nbytes, NULL, elliptic->o_gradRecvBuffer);

  }else{
    elliptic->sendBuffer = NULL;
    elliptic->recvBuffer = NULL;
  }
  mesh->device.setStream(elliptic->defaultStream);

  elliptic->type = strdup(dfloatString);

  elliptic->Nblock = Nblock;
  elliptic->Nblock2 = Nblock2;

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
  elliptic->allNeumannPenalty = 1.;
  hlong localElements = (hlong) mesh->Nelements;
  hlong totalElements = 0;
  MPI_Allreduce(&localElements, &totalElements, 1, MPI_HLONG, MPI_SUM, mesh->comm);
  elliptic->allNeumannScale = 1./sqrt((dfloat)mesh->Np*totalElements);

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


  // !!!!!! Removed MPI::BOOL since some mpi versions complains about it !!!!!
  int lallNeumann, gallNeumann;
  lallNeumann = allNeumann ? 0:1;
  MPI_Allreduce(&lallNeumann, &gallNeumann, 1, MPI_INT, MPI_SUM, mesh->comm);
  elliptic->allNeumann = (gallNeumann>0) ? false: true;

  // MPI_Allreduce(&allNeumann, &(elliptic->allNeumann), 1, MPI::BOOL, MPI_LAND, mesh->comm);
  if (mesh->rank==0&& options.compareArgs("VERBOSE","TRUE")) printf("allNeumann = %d \n", elliptic->allNeumann);

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
  ogsGatherScatter(elliptic->mapB, ogsInt, ogsMin, mesh->ogs);

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

  //make a masked version of the global id numbering
  mesh->maskedGlobalIds = (hlong *) calloc(Ntotal,sizeof(hlong));
  memcpy(mesh->maskedGlobalIds, mesh->globalIds, Ntotal*sizeof(hlong));
  for (dlong n=0;n<elliptic->Nmasked;n++)
    mesh->maskedGlobalIds[elliptic->maskIds[n]] = 0;

  //use the masked ids to make another gs handle
  elliptic->ogs = ogsSetup(Ntotal, mesh->maskedGlobalIds, mesh->comm, verbose, mesh->device);
  elliptic->o_invDegree = elliptic->ogs->o_invDegree;

  /*preconditioner setup */
  elliptic->precon = (precon_t*) calloc(1, sizeof(precon_t));

  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo["compiler_flags"] += "-Xptxas -dlcm=ca";
  }

  if(mesh->device.mode()=="Serial")
    kernelInfo["compiler_flags"] += "-g";

  // set kernel name suffix
  char *suffix;

  if(elliptic->elementType==TRIANGLES){
    if(elliptic->dim==2)
      suffix = strdup("Tri2D");
    else
      suffix = strdup("Tri3D");
  }
  if(elliptic->elementType==QUADRILATERALS){
    if(elliptic->dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D");
  }
  if(elliptic->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(elliptic->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];


  for (int r=0;r<mesh->size;r++) {
    if (r==mesh->rank) {

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

      elliptic->weightedInnerProduct1Kernel =
        mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct1.okl",
                                       "weightedInnerProduct1",
                                       kernelInfo);

      elliptic->weightedInnerProduct2Kernel =
        mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct2.okl",
                                       "weightedInnerProduct2",
                                       kernelInfo);

      elliptic->innerProductKernel =
        mesh->device.buildKernel(DHOLMES "/okl/innerProduct.okl",
                                       "innerProduct",
                                       kernelInfo);

      elliptic->weightedNorm2Kernel =
        mesh->device.buildKernel(DHOLMES "/okl/weightedNorm2.okl",
                                           "weightedNorm2",
                                           kernelInfo);

      elliptic->norm2Kernel =
        mesh->device.buildKernel(DHOLMES "/okl/norm2.okl",
                                           "norm2",
                                           kernelInfo);


      elliptic->scaledAddKernel =
          mesh->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
                                         "scaledAdd",
                                         kernelInfo);

      elliptic->dotMultiplyKernel =
          mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
                                         "dotMultiply",
                                         kernelInfo);

      elliptic->dotDivideKernel =
          mesh->device.buildKernel(DHOLMES "/okl/dotDivide.okl",
                                         "dotDivide",
                                         kernelInfo);

      // add custom defines
      kernelInfo["defines/" "p_NpP"]= (mesh->Np+mesh->Nfp*mesh->Nfaces);
      kernelInfo["defines/" "p_Nverts"]= mesh->Nverts;

      //sizes for the coarsen and prolongation kernels. degree N to degree 1
      kernelInfo["defines/" "p_NpFine"]= mesh->Np;
      kernelInfo["defines/" "p_NpCoarse"]= mesh->Nverts;


      if (elliptic->elementType==QUADRILATERALS || elliptic->elementType==HEXAHEDRA) {
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

      //add standard boundary functions
      char *boundaryHeaderFileName;
      if (elliptic->dim==2)
        boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary2D.h");
      else if (elliptic->dim==3)
        boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary3D.h");
      kernelInfo["includes"] += boundaryHeaderFileName;


      sprintf(fileName,  DELLIPTIC "/okl/ellipticAx%s.okl", suffix);
      sprintf(kernelName, "ellipticAx%s", suffix);

      occa::properties dfloatKernelInfo = kernelInfo;
      occa::properties floatKernelInfo = kernelInfo;
      floatKernelInfo["defines/" "pfloat"]= "float";
      dfloatKernelInfo["defines/" "pfloat"]= dfloatString;

      elliptic->AxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);

      if(elliptic->elementType!=HEXAHEDRA){
        sprintf(kernelName, "ellipticPartialAx%s", suffix);
      }
      else{
        if(elliptic->options.compareArgs("ELEMENT MAP", "TRILINEAR")){
          sprintf(kernelName, "ellipticPartialAxTrilinear%s", suffix);
        }else{
          sprintf(kernelName, "ellipticPartialAx%s", suffix);
        }
      }

      elliptic->partialAxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);

      elliptic->partialFloatAxKernel = mesh->device.buildKernel(fileName,kernelName,floatKernelInfo);

      // Not implemented for Quad3D !!!!!
      if (options.compareArgs("BASIS","BERN")) {

        sprintf(fileName, DELLIPTIC "/okl/ellipticGradientBB%s.okl", suffix);
        sprintf(kernelName, "ellipticGradientBB%s", suffix);

        elliptic->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialGradientBB%s", suffix);
        elliptic->partialGradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdgBB%s.okl", suffix);
        sprintf(kernelName, "ellipticAxIpdgBB%s", suffix);
        elliptic->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialAxIpdgBB%s", suffix);
        elliptic->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      } else if (options.compareArgs("BASIS","NODAL")) {

        sprintf(fileName, DELLIPTIC "/okl/ellipticGradient%s.okl", suffix);
        sprintf(kernelName, "ellipticGradient%s", suffix);

        elliptic->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialGradient%s", suffix);
        elliptic->partialGradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdg%s.okl", suffix);
        sprintf(kernelName, "ellipticAxIpdg%s", suffix);
        elliptic->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialAxIpdg%s", suffix);
        elliptic->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
      }

      // Use the same kernel with quads for the following kenels
      if(elliptic->dim==3){
	if(elliptic->elementType==QUADRILATERALS)
	  suffix = strdup("Quad2D");
	else if(elliptic->elementType==TRIANGLES)
	  suffix = strdup("Tri2D");
      }
	
      sprintf(fileName, DELLIPTIC "/okl/ellipticPreconCoarsen%s.okl", suffix);
      sprintf(kernelName, "ellipticPreconCoarsen%s", suffix);
      elliptic->precon->coarsenKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(fileName, DELLIPTIC "/okl/ellipticPreconProlongate%s.okl", suffix);
      sprintf(kernelName, "ellipticPreconProlongate%s", suffix);
      elliptic->precon->prolongateKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      

      sprintf(fileName, DELLIPTIC "/okl/ellipticBlockJacobiPrecon.okl");
      sprintf(kernelName, "ellipticBlockJacobiPrecon");
      elliptic->precon->blockJacobiKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(kernelName, "ellipticPartialBlockJacobiPrecon");
      elliptic->precon->partialblockJacobiKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(fileName, DELLIPTIC "/okl/ellipticPatchSolver.okl");
      sprintf(kernelName, "ellipticApproxBlockJacobiSolver");
      elliptic->precon->approxBlockJacobiSolverKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      if (   elliptic->elementType == TRIANGLES
          || elliptic->elementType == TETRAHEDRA) {
        elliptic->precon->SEMFEMInterpKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMInterp.okl",
                     "ellipticSEMFEMInterp",
                     kernelInfo);

        elliptic->precon->SEMFEMAnterpKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMAnterp.okl",
                     "ellipticSEMFEMAnterp",
                     kernelInfo);
      }
    }
    MPI_Barrier(mesh->comm);
  }

  long long int pre = mesh->device.memoryAllocated();

  occaTimerTic(mesh->device,"PreconditionerSetup");
  ellipticPreconditionerSetup(elliptic, elliptic->ogs, lambda);
  occaTimerToc(mesh->device,"PreconditionerSetup");

  long long int usedBytes = mesh->device.memoryAllocated()-pre;

  elliptic->precon->preconBytes = usedBytes;
}

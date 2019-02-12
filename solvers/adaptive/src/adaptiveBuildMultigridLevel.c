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

// create adaptive and mesh structs for multigrid levels
adaptive_t *adaptiveBuildMultigridLevel(adaptive_t *baseAdaptive, int Nc, int Nf){

  adaptive_t *adaptive = (adaptive_t*) calloc(1, sizeof(adaptive_t));
  //adaptive_t *adaptive = new adaptive_t[1];

  adaptive->dim = baseAdaptive->dim;
  adaptive->elementType = baseAdaptive->elementType;
  adaptive->options = baseAdaptive->options;

  adaptive->tau = baseAdaptive->tau;
  adaptive->BCType = baseAdaptive->BCType;
  adaptive->allNeumann = baseAdaptive->allNeumann;
  adaptive->allNeumannPenalty = baseAdaptive->allNeumannPenalty;

  adaptive->sendBuffer = baseAdaptive->sendBuffer;
  adaptive->recvBuffer = baseAdaptive->recvBuffer;
  adaptive->gradSendBuffer = baseAdaptive->gradSendBuffer;
  adaptive->gradRecvBuffer = baseAdaptive->gradRecvBuffer;

  adaptive->defaultStream = baseAdaptive->defaultStream;
  adaptive->dataStream = baseAdaptive->dataStream;

  adaptive->o_EToB = baseAdaptive->o_EToB;

  adaptive->o_grad = baseAdaptive->o_grad;

  adaptive->o_EXYZ = baseAdaptive->o_EXYZ;

  adaptive->o_ggeoNoJW = baseAdaptive->o_ggeoNoJW; // ?
  
  adaptive->weightedInnerProduct1Kernel = baseAdaptive->weightedInnerProduct1Kernel;
  adaptive->weightedInnerProduct2Kernel = baseAdaptive->weightedInnerProduct2Kernel;
  adaptive->innerProductKernel = baseAdaptive->innerProductKernel;
  adaptive->weightedNorm2Kernel = baseAdaptive->weightedNorm2Kernel;
  adaptive->norm2Kernel = baseAdaptive->norm2Kernel;
  adaptive->scaledAddKernel = baseAdaptive->scaledAddKernel;
  adaptive->dotMultiplyKernel = baseAdaptive->dotMultiplyKernel;
  adaptive->dotDivideKernel = baseAdaptive->dotDivideKernel;

  //populate the mini-mesh using the mesh struct
  mesh_t *mesh = (mesh_t*) calloc(1,sizeof(mesh_t));
  //  mesh_t *mesh = new mesh_t[1];

#ifndef OCCA_VERSION_1_0
  memcpy(mesh,baseAdaptive->mesh,sizeof(mesh_t));
#else

  mesh->rank = baseAdaptive->mesh->rank;
  mesh->size = baseAdaptive->mesh->size;

  MPI_Comm_dup(baseAdaptive->mesh->comm, &(mesh->comm));

  mesh->dim = baseAdaptive->mesh->dim;
  mesh->Nverts        = baseAdaptive->mesh->Nverts;
  mesh->Nfaces        = baseAdaptive->mesh->Nfaces;
  mesh->NfaceVertices = baseAdaptive->mesh->NfaceVertices;

  mesh->Nfields = baseAdaptive->mesh->Nfields;
  
  mesh->Nnodes = baseAdaptive->mesh->Nnodes;
  mesh->EX = baseAdaptive->mesh->EX; // coordinates of vertices for each element
  mesh->EY = baseAdaptive->mesh->EY;
  mesh->EZ = baseAdaptive->mesh->EZ;

  mesh->Nelements = baseAdaptive->mesh->Nelements;
  mesh->EToV = baseAdaptive->mesh->EToV; // element-to-vertex connectivity
  mesh->EToE = baseAdaptive->mesh->EToE; // element-to-element connectivity
  mesh->EToF = baseAdaptive->mesh->EToF; // element-to-(local)face connectivity
  mesh->EToP = baseAdaptive->mesh->EToP; // element-to-partition/process connectivity
  mesh->EToB = baseAdaptive->mesh->EToB; // element-to-boundary condition type

  mesh->elementInfo = baseAdaptive->mesh->elementInfo;

  // boundary faces
  mesh->NboundaryFaces = baseAdaptive->mesh->NboundaryFaces;
  mesh->boundaryInfo = baseAdaptive->mesh->boundaryInfo;

  // MPI halo exchange info
  mesh->totalHaloPairs = baseAdaptive->mesh->totalHaloPairs;
  mesh->haloElementList = baseAdaptive->mesh->haloElementList;
  mesh->NhaloPairs = baseAdaptive->mesh->NhaloPairs;
  mesh->NhaloMessages = baseAdaptive->mesh->NhaloMessages;

  mesh->haloSendRequests = baseAdaptive->mesh->haloSendRequests;
  mesh->haloRecvRequests = baseAdaptive->mesh->haloRecvRequests;

  mesh->NinternalElements = baseAdaptive->mesh->NinternalElements;
  mesh->NnotInternalElements = baseAdaptive->mesh->NnotInternalElements;

  mesh->o_haloElementList = baseAdaptive->mesh->o_haloElementList;
  mesh->o_haloBuffer      = baseAdaptive->mesh->o_haloBuffer;
  mesh->o_internalElementIds    = baseAdaptive->mesh->o_internalElementIds;
  mesh->o_notInternalElementIds = baseAdaptive->mesh->o_notInternalElementIds;

  // volumeGeometricFactors;
  mesh->Nvgeo = baseAdaptive->mesh->Nvgeo;
  mesh->vgeo = baseAdaptive->mesh->vgeo;
  mesh->o_vgeo = baseAdaptive->mesh->o_vgeo;

  // second order volume geometric factors
  mesh->Nggeo = baseAdaptive->mesh->Nggeo;
  mesh->ggeo = baseAdaptive->mesh->ggeo;
  mesh->o_ggeo = baseAdaptive->mesh->o_ggeo;

  mesh->Nsgeo = baseAdaptive->mesh->Nsgeo;
  mesh->sgeo = baseAdaptive->mesh->sgeo;
  mesh->o_sgeo = baseAdaptive->mesh->o_sgeo;

  // occa stuff
  mesh->device = baseAdaptive->mesh->device;

#if USE_MASTER_NOEL==1
  //  mesh->device.UsePreCompiledKernels(mesh->rank!=0);

  int foo;
  // check to see if the options specify to use precompiled binaries
  if(adaptive->options.compareArgs("USE PRECOMPILED BINARIES", "TRUE")){
    mesh->device.UsePreCompiledKernels(1);
    occa::host().UsePreCompiledKernels(1);
  }
  else if(adaptive->options.compareArgs("USE PRECOMPILED BINARIES", "NONROOT")){
    mesh->device.UsePreCompiledKernels(mesh->rank!=0);
    occa::host().UsePreCompiledKernels(mesh->rank!=0);
  }else{
    mesh->device.UsePreCompiledKernels(0);
    occa::host().UsePreCompiledKernels(0);
  }
  
#endif
  
  mesh->defaultStream = baseAdaptive->mesh->defaultStream;
  mesh->dataStream = baseAdaptive->mesh->dataStream;

  mesh->haloExtractKernel = baseAdaptive->mesh->haloExtractKernel;
  mesh->addScalarKernel = baseAdaptive->mesh->addScalarKernel;
  mesh->maskKernel = baseAdaptive->mesh->maskKernel;
  mesh->sumKernel = baseAdaptive->mesh->sumKernel;
#endif

  adaptive->mesh = mesh;

  setupAide options = adaptive->options;

  meshLoadReferenceNodesHex3D(mesh, Nc);
  meshPhysicalNodesHex3D(mesh);
  meshGeometricFactorsHex3D(mesh);


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

  if (adaptive->dim==3) {
    mesh->z = (dfloat*) realloc(mesh->z, (localNodes+totalHaloNodes)*sizeof(dfloat));
    meshHaloExchange(mesh, mesh->Np*sizeof(dfloat), mesh->z, sendBuffer, mesh->z + localNodes);
  }

  if(!options.compareArgs("BOX DOMAIN", "TRUE"))
    meshConnectFaceNodes3D(mesh);
  else{
    dfloat XMIN = -1, XMAX = +1; // default bi-unit cube
    dfloat YMIN = -1, YMAX = +1;
    dfloat ZMIN = -1, ZMAX = +1;
    
    options.getArgs("BOX XMIN", XMIN);
    options.getArgs("BOX YMIN", YMIN);
    options.getArgs("BOX ZMIN", ZMIN);
    
    options.getArgs("BOX XMAX", XMAX);
    options.getArgs("BOX YMAX", YMAX);
    options.getArgs("BOX ZMAX", ZMAX);
    
    meshConnectPeriodicFaceNodes3D(mesh, XMAX-XMIN, YMAX-YMIN, ZMAX-ZMIN);
  }
  meshSurfaceGeometricFactorsHex3D(mesh);

  // global nodes
  meshParallelConnectNodes(mesh);

  //dont need these once vmap is made
  free(mesh->x);
  free(mesh->y);
  if (adaptive->dim==3) {
    free(mesh->z);
  }
  free(sendBuffer);

  dlong Ntotal = mesh->Np*mesh->Nelements;
  dlong Nblock = mymax(1,(Ntotal+blockSize-1)/blockSize);
  dlong Nblock2 = mymax(1,(Nblock+blockSize-1)/blockSize);

  adaptive->Nblock = Nblock;
  adaptive->Nblock2 = Nblock2;

  //lumped mass matrix
  mesh->MM = (dfloat *) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *DT = (dfloat*) calloc(mesh->Nq*mesh->Nq, sizeof(dfloat));
  
  for (int j=0;j<mesh->Nq;j++) {
    for (int i=0;i<mesh->Nq;i++) {
      DT[j*mesh->Nq+i] = mesh->D[i*mesh->Nq+j];
    }
  }
  
  for (int k=0;k<mesh->Nq;k++) {
    for (int j=0;j<mesh->Nq;j++) {
      for (int i=0;i<mesh->Nq;i++) {
	int n = i+j*mesh->Nq+k*mesh->Nq*mesh->Nq;
	mesh->MM[n+n*mesh->Np] = mesh->gllw[i]*mesh->gllw[j]*mesh->gllw[k];
      }
    }
  }
  
  mesh->o_D = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);
  mesh->o_Dmatrices = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);
  mesh->o_Smatrices = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), DT); // transpose(D)
  
  mesh->o_cubD = mesh->device.malloc(mesh->cubNq*mesh->cubNq*sizeof(dfloat), mesh->cubD);
  
  dfloat *cubInterpT = (dfloat*) calloc(mesh->cubNq*mesh->Nq, sizeof(dfloat));
  for(int n=0;n<mesh->Nq;++n){
    for(int m=0;m<mesh->cubNq;++m){        
        cubInterpT[m+n*mesh->cubNq] = mesh->cubInterp[m*mesh->Nq+n];
    }
  }
  
  mesh->o_cubInterpT = mesh->device.malloc(mesh->cubNq*mesh->Nq*sizeof(dfloat), cubInterpT);
  
  free(cubInterpT);
  
  mesh->o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*mesh->Np*sizeof(dfloat),
			mesh->vgeo);
  
  mesh->o_sgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nfp*mesh->Nsgeo*sizeof(dfloat),
			mesh->sgeo);
  
  mesh->o_ggeo =
    mesh->device.malloc(mesh->Nelements*mesh->Np*mesh->Nggeo*sizeof(dfloat),
			mesh->ggeo);
  
  mesh->o_cubggeo =
    mesh->device.malloc(mesh->Nelements*mesh->cubNp*mesh->Nggeo*sizeof(dfloat),
			mesh->cubggeo);
  
  
  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(dlong),
			mesh->vmapM);
  
  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(dlong),
			mesh->vmapP);
  
  mesh->LIFT = baseAdaptive->mesh->LIFT; //dummy buffer
  mesh->o_LIFTT = baseAdaptive->mesh->o_LIFTT; //dummy buffer

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

  mesh->o_MM =
    mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
                        mesh->MM);

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(int),
                        mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(int),
                        mesh->vmapP);


  //set the normalization constant for the allNeumann Poisson problem on this coarse mesh
  hlong localElements = (hlong) mesh->Nelements;
  hlong totalElements = 0;
  MPI_Allreduce(&localElements, &totalElements, 1, MPI_HLONG, MPI_SUM, mesh->comm);
  adaptive->allNeumannScale = 1.0/sqrt(mesh->Np*totalElements);

  adaptive->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));
  adaptive->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), adaptive->tmp);
  adaptive->o_tmp2 = mesh->device.malloc(Nblock2*sizeof(dfloat), adaptive->tmp);

  //tau
  adaptive->tau = 2.0*(mesh->N+1)*(mesh->N+3);

  //setup an unmasked gs handle
  int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;
  meshParallelGatherScatterSetup(mesh, Ntotal, mesh->globalIds, mesh->comm, verbose);

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  adaptive->mapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) adaptive->mapB[n+e*mesh->Np] = 1E9;
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
    if (adaptive->mapB[n] == 1E9) {
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

  // HERE
  occa::properties kernelInfo = adaptiveKernelInfo(mesh);
  
  //  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  // set kernel name suffix
  char *suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  for (int r=0;r<2;r++){

    MPI_Barrier(mesh->comm);

    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {
      
      kernelInfo["defines/" "p_blockSize"]= blockSize;
      
      // add custom defines
      kernelInfo["defines/" "p_NpP"]= (mesh->Np+mesh->Nfp*mesh->Nfaces);
      kernelInfo["defines/" "p_Nverts"]= mesh->Nverts;

      int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
      kernelInfo["defines/" "p_Nmax"]= Nmax;

      int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
      kernelInfo["defines/" "p_maxNodes"]= maxNodes;

      int NblockV = mymax(1,maxNthreads/mesh->Np); // works for CUDA
      kernelInfo["defines/" "p_NblockV"]= NblockV;

      int one = 1; //set to one for now. TODO: try optimizing over these
      kernelInfo["defines/" "p_NnodesV"]= one;

      int NblockS = mymax(1,maxNthreads/maxNodes); // works for CUDA
      kernelInfo["defines/" "p_NblockS"]= NblockS;

      int NblockP = mymax(1,maxNthreads/(4*mesh->Np)); // get close to maxNthreads threads
      kernelInfo["defines/" "p_NblockP"]= NblockP;

      int NblockG;
      if(mesh->Np<=32) NblockG = ( 32/mesh->Np );
      else NblockG = mymax(1,maxNthreads/mesh->Np);
      kernelInfo["defines/" "p_NblockG"]= NblockG;

      //add standard boundary functions
      char *boundaryHeaderFileName;
      if (adaptive->dim==2)
        boundaryHeaderFileName = strdup(DADAPTIVE "/data/adaptiveBoundary2D.h");
      else if (adaptive->dim==3)
        boundaryHeaderFileName = strdup(DADAPTIVE "/data/adaptiveBoundary3D.h");
      kernelInfo["includes"] += boundaryHeaderFileName;

      occa::properties dfloatKernelInfo = kernelInfo;
      occa::properties floatKernelInfo = kernelInfo;
      floatKernelInfo["defines/" "pfloat"]= "float";
      dfloatKernelInfo["defines/" "pfloat"]= dfloatString;

      printf("rank = %d, size = %d\n", mesh->rank, mesh->size);
      //      std::cout << dfloatKernelInfo << std::endl;
      
      sprintf(fileName, DADAPTIVE "/okl/adaptiveAx%s.okl", suffix);
      sprintf(kernelName, "adaptiveAx%s", suffix);
      adaptive->AxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo); 

      // check for trilinear
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

      //sprintf(kernelName, "adaptivePartialAx%s", suffix);

      adaptive->partialAxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);

      adaptive->partialFloatAxKernel = mesh->device.buildKernel(fileName,kernelName,floatKernelInfo);

      // only for Hex3D - cubature Ax
      if(adaptive->elementType==HEXAHEDRA){
	//	printf("BUILDING partialCubatureAxKernel\n");
	sprintf(fileName,  DADAPTIVE "/okl/adaptiveCubatureAx%s.okl", suffix);
	
	sprintf(kernelName, "adaptiveCubaturePartialAx%s", suffix);
	adaptive->partialCubatureAxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);
      }

      
      if (options.compareArgs("BASIS", "BERN")) {

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

      } else if (options.compareArgs("BASIS", "NODAL")) {

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
    }

    MPI_Barrier(mesh->comm);
  }

  //new precon struct
  adaptive->precon = (precon_t *) calloc(1,sizeof(precon_t));
  //  adaptive->precon = new precon_t[1];

  for (int r=0;r<2;r++){

    MPI_Barrier(mesh->comm);
    
    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {

      sprintf(fileName, DADAPTIVE "/okl/adaptiveBlockJacobiPrecon.okl");
      sprintf(kernelName, "adaptiveBlockJacobiPrecon");
      adaptive->precon->blockJacobiKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(kernelName, "adaptivePartialBlockJacobiPrecon");
      adaptive->precon->partialblockJacobiKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      //sizes for the coarsen and prolongation kernels. degree NFine to degree N
      int NqFine   = (Nf+1);
      int NqCoarse = (Nc+1);
      kernelInfo["defines/" "p_NqFine"]= Nf+1;
      kernelInfo["defines/" "p_NqCoarse"]= Nc+1;

      int NpFine, NpCoarse;
      NpFine   = (Nf+1)*(Nf+1)*(Nf+1);
      NpCoarse = (Nc+1)*(Nc+1)*(Nc+1);

      kernelInfo["defines/" "p_NpFine"]= NpFine;
      kernelInfo["defines/" "p_NpCoarse"]= NpCoarse;

      int NblockVFine = maxNthreads/NpFine;
      int NblockVCoarse = maxNthreads/NpCoarse;
      kernelInfo["defines/" "p_NblockVFine"]= NblockVFine;
      kernelInfo["defines/" "p_NblockVCoarse"]= NblockVCoarse;

      sprintf(fileName, DADAPTIVE "/okl/adaptivePreconCoarsen%s.okl", suffix);
      sprintf(kernelName, "adaptivePreconCoarsen%s", suffix);
      adaptive->precon->coarsenKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(fileName, DADAPTIVE "/okl/adaptivePreconProlongate%s.okl", suffix);
      sprintf(kernelName, "adaptivePreconProlongate%s", suffix);
      adaptive->precon->prolongateKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }
  
  if(adaptive->elementType==HEXAHEDRA){
    if(options.compareArgs("DISCRETIZATION","CONTINUOUS")){
      if(options.compareArgs("ELEMENT MAP", "TRILINEAR")){

        // pack gllz, gllw, and elementwise EXYZ
        dfloat *gllzw = (dfloat*) calloc(2*mesh->Nq, sizeof(dfloat));

        int sk = 0;
        for(int n=0;n<mesh->Nq;++n)
          gllzw[sk++] = mesh->gllz[n];
        for(int n=0;n<mesh->Nq;++n)
          gllzw[sk++] = mesh->gllw[n];

        adaptive->o_gllzw = mesh->device.malloc(2*mesh->Nq*sizeof(dfloat), gllzw);
        free(gllzw);
      }
    }
  }


  return adaptive;
}

#include "ellipticBenchmarkTet3D.h"

solver_t *ellipticSetupTet3D(mesh_t *mesh, dfloat tau, dfloat lambda, int*BCType,
			     occa::kernelInfo &kernelInfo, const char *options, const char *parAlmondOptions,
			     int Nblocks, int Nnodes){
  
  mesh->Nfields = 1;

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", rank%2);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");

  meshOccaSetup3D(mesh, deviceConfig, kernelInfo);

  int Ntotal = mesh->Np*mesh->Nelements;
  int Nblock = (Ntotal+blockSize-1)/blockSize;
  int Nhalo = mesh->Np*mesh->totalHaloPairs;
  int Nall   = Ntotal + Nhalo;

  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));

  solver->tau = tau;

  solver->mesh = mesh;

  solver->p   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->z   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->Ax  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->Ap  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  solver->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));

  solver->grad = (dfloat*) calloc(Nall*4, sizeof(dfloat));

  solver->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);

  solver->o_res = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
  solver->o_Sres = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
  solver->o_Ax  = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), solver->Ap);
  solver->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), solver->tmp);

  solver->o_grad  = mesh->device.malloc(Nall*4*sizeof(dfloat), solver->grad);

  //setup async halo stream
  solver->defaultStream = mesh->device.getStream();
  solver->dataStream = mesh->device.createStream();
  mesh->device.setStream(solver->defaultStream);

  int Nbytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  if(Nbytes>0){
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(Nbytes, NULL);

    solver->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    solver->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();

    occa::memory o_gradSendBuffer = mesh->device.mappedAlloc(2*Nbytes, NULL);
    occa::memory o_gradRecvBuffer = mesh->device.mappedAlloc(2*Nbytes, NULL);

    solver->gradSendBuffer = (dfloat*) o_gradSendBuffer.getMappedPointer();
    solver->gradRecvBuffer = (dfloat*) o_gradRecvBuffer.getMappedPointer();
  }else{
    solver->sendBuffer = NULL;
    solver->recvBuffer = NULL;
  }
  mesh->device.setStream(solver->defaultStream);

  solver->type = strdup(dfloatString);

  solver->Nblock = Nblock;

  //fill geometric factors in halo
  if(mesh->totalHaloPairs){
    int Nlocal = mesh->Nelements*mesh->Np;
    int Nhalo  = mesh->totalHaloPairs*mesh->Np;
    dfloat *vgeoSendBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Nvgeo, sizeof(dfloat));

    // import geometric factors from halo elements
    mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (mesh->Nelements+mesh->totalHaloPairs)*mesh->Nvgeo*sizeof(dfloat));

    meshHaloExchange(mesh,
        mesh->Nvgeo*sizeof(dfloat),
        mesh->vgeo,
        vgeoSendBuffer,
        mesh->vgeo + mesh->Nelements*mesh->Nvgeo);

    mesh->o_vgeo =
      mesh->device.malloc((mesh->Nelements + mesh->totalHaloPairs)*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);
  }

  //check all the bounaries for a Dirichlet
  bool allNeumann = (lambda==0) ? true :false;
  solver->allNeumannPenalty = 1;
  int totalElements = 0;
  MPI_Allreduce(&(mesh->Nelements), &totalElements, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  solver->allNeumannScale = 1.0/sqrt(mesh->Np*totalElements);

  solver->EToB = (int *) calloc(mesh->Nelements*mesh->Nfaces,sizeof(int));
  for (int e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if (bc>0) {
        int BC = BCType[bc]; //get the type of boundary
        solver->EToB[e*mesh->Nfaces+f] = BC; //record it
        if (BC!=2) allNeumann = false; //check if its a Dirchlet
      }
    }
  }
  MPI_Allreduce(&allNeumann, &(solver->allNeumann), 1, MPI::BOOL, MPI_LAND, MPI_COMM_WORLD);
  printf("allNeumann = %d \n", solver->allNeumann);

  solver->o_EToB = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int), solver->EToB);

  //add standard boundary functions
  char *boundaryHeaderFileName;
  boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTet3D/ellipticBoundary3D.h");
  kernelInfo.addInclude(boundaryHeaderFileName);

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");
  }

  kernelInfo.addDefine("p_blockSize", blockSize);

  // add custom defines
  kernelInfo.addDefine("p_NpP", (mesh->Np+mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);
  //sizes for the coarsen and prolongation kernels. degree N to degree 1
  kernelInfo.addDefine("p_NpFine", mesh->Np);
  kernelInfo.addDefine("p_NpCoarse", mesh->Nverts);

  int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nmax", Nmax);

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  kernelInfo.addDefine("p_NblockV", Nblocks);
  kernelInfo.addDefine("p_NnodesV", Nnodes);

  kernelInfo.addDefine("p_NblockS", Nblocks);
  kernelInfo.addDefine("p_NblockP", Nblocks);
  kernelInfo.addDefine("p_NblockG", Nblocks);
  int *globalGatherElementList    = (int*) calloc(mesh->Nelements, sizeof(int));
  int *localGatherElementList = (int*) calloc(mesh->Nelements, sizeof(int));
  int globalCount = 0;
  int localCount =0;

  for(int e=0;e<mesh->Nelements;++e){
    globalGatherElementList[globalCount++] = e;
    localGatherElementList[localCount++] = e;
  }

  solver->NglobalGatherElements = mesh->Nelements;
  solver->NlocalGatherElements = mesh->Nelements;

  solver->o_globalGatherElementList =
    mesh->device.malloc(globalCount*sizeof(int), globalGatherElementList);

  solver->o_localGatherElementList =
    mesh->device.malloc(localCount*sizeof(int), localGatherElementList);
  mesh->invSparseV = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  for (int n=0;n<mesh->Np*mesh->Np;n++)
    mesh->invSparseV[n] = mesh->sparseV[n];

  // TW: what is all this ? ----------------------------------------------------------------------------->
  //very simple, just look through the indices; there will be three nnz groups in worst case
  char *rowData = (char*) calloc(8*mesh->Np, sizeof(char));
  for (int i=0; i<mesh->Np; ++i){
    rowData[8*i] = (char)mesh->sparseStackedNZ[0*mesh->Np+i];    
    int nogroups = 1; 
    int j=0;   
    while (nogroups != 4){
      if (j<mesh->SparseNnzPerRow-1){      
        //printf("comparing %d and %d \n", mesh->sparseStackedNZ[j*mesh->Np+i], mesh->sparseStackedNZ[(j+1)*mesh->Np+i]);        
        if (mesh->sparseStackedNZ[j*mesh->Np+i]+1 == mesh->sparseStackedNZ[(j+1)*mesh->Np+i]){
          j++;
        } 
        else{
          if (mesh->sparseStackedNZ[(j+1)*mesh->Np+i] == 0){
            rowData[8*i + 2*(nogroups-1)+1] =  (char) mesh->sparseStackedNZ[j*mesh->Np+i];
            nogroups = 4;
          }
          else{
            //printf("ending %d starting at %d \n " , mesh->sparseStackedNZ[j*mesh->Np+i], mesh->sparseStackedNZ[(j+1)*mesh->Np+i]);
            rowData[8*i + 2*(nogroups-1)+1] = (char) mesh->sparseStackedNZ[j*mesh->Np+i];
            rowData[8*i+ 2*(nogroups-1)+2] = (char) mesh->sparseStackedNZ[(j+1)*mesh->Np+i];
            nogroups = nogroups+1;
            j++;
          }
        }
      }//if smaller then row lenght
      else{
        rowData[8*i +  nogroups] = (char) mesh->sparseStackedNZ[j*mesh->Np+i];
        nogroups = 4;
      }
    }


    printf("THIS IS ROW %d GROUP DATA\n", i);
    for (int j=0; j<8; j++){
      printf(" %d ", (int) rowData[8*i+j]);

    }
    if (rowData[8*i+2] == rowData[8*i+3]){
      rowData[8*i+2] = (char)1;
    }
    if (rowData[8*i+4] == rowData[8*i+5]){
      rowData[8*i+4] = (char)1;
    }   

    printf("\n");
  }

  //transpose
  char *rowDataTransposed = (char*) calloc(8*mesh->Np, sizeof(char));
  for (int i=0; i<mesh->Np; ++i){
    for (int j=0; j<4; j++){
      rowDataTransposed[4*i+j] = (char) rowData[8*i+j];    
    }
    for (int j=0; j<4; j++){
      rowDataTransposed[mesh->Np*4 + 4*i+j] = (char) rowData[8*i+j+4];
    }
  }
  for (int j=0; j<8*mesh->Np; ++j){
    if ((j)%4 == 0){
      printf(" || ");
    }    
    if ((j)%(4*mesh->Np) == 0){
      printf("\n\n");
    }
    printf(" %d ",  (int)rowDataTransposed[j]  );
  }
  // TW: <---------------------------------------------------------------------what is all this ? 
  printf("\n MAX IS %d \n", mesh->Np*8);
  mesh->o_rowData = mesh->device.malloc(mesh->Np*8*sizeof(char), rowDataTransposed);

  mesh->o_Srr = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSrrT);
  mesh->o_Srs = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSrsT);
  mesh->o_Srt = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSrtT);
  mesh->o_Sss = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSssT);  
  mesh->o_Sst = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSstT);  
  mesh->o_Stt = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSttT);  

  kernelInfo.addDefine("p_SparseNnzPerRow", mesh->SparseNnzPerRow);


  return solver;

}



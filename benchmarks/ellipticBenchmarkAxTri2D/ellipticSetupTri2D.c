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

#include "ellipticBenchmarkTri2D.h"

solver_t *ellipticSetupTri2D(mesh_t *mesh, dfloat tau, dfloat lambda, int*BCType,
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

  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

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
  boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTri2D/ellipticBoundary2D.h");
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
#if 0
  //int paddedRowSize = 4*((mesh->SparseNnzPerRow+3)/4); //make the nnz per row a multiple of 4
  int paddedRowSize = mesh->SparseNnzPerRow + (16-mesh->SparseNnzPerRow%16);
  printf("row size %d padded %d \n ", mesh->SparseNnzPerRow, paddedRowSize);
  char* IndTchar = (char*) calloc(paddedRowSize*mesh->Np,sizeof(char));

  for (int m=0;m<paddedRowSize/4;m++) {
    for (int n=0;n<mesh->Np;n++) {
      for (int k=0;k<4;k++) {
        if (k+4*m < mesh->SparseNnzPerRow) {
          IndTchar[k+4*n+m*4*mesh->Np] = mesh->sparseStackedNZ[n+(k+4*m)*mesh->Np];
          //printf("putting %d in indT %d \n", k+4*n+m*4*mesh->Np,mesh->sparseStackedNZ[n+(k+4*m)*mesh->Np]);  
        } else {
          IndTchar[k+4*n+m*4*mesh->Np] = 0;
        }
      }
    }
  }

  dfloat * sparseSrrTResized = (dfloat *) calloc(paddedRowSize*mesh->Np,sizeof(dfloat));
  dfloat * sparseSrsTResized = (dfloat *) calloc(paddedRowSize*mesh->Np,sizeof(dfloat));
  dfloat * sparseSssTResized = (dfloat *) calloc(paddedRowSize*mesh->Np,sizeof(dfloat));
  for (int m=0;m<paddedRowSize;m++) {
    for (int n=0;n<mesh->Np;n++) {
      if (m<mesh->SparseNnzPerRow){
        sparseSrrTResized[m*mesh->Np + n] = mesh->sparseSrrT[m*mesh->Np + n];
        sparseSrsTResized[m*mesh->Np + n] = mesh->sparseSrsT[m*mesh->Np + n];
        sparseSssTResized[m*mesh->Np + n] = mesh->sparseSssT[m*mesh->Np + n];

      }

    }
  }




  kernelInfo.addDefine("p_SparseNnzPerRow", paddedRowSize);
  printf("sparse NNZ %d \n", paddedRowSize);  
  mesh->SparseNnzPerRowNonPadded = mesh->SparseNnzPerRow;
  mesh->SparseNnzPerRow = paddedRowSize;


  mesh->o_sparseSrrT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), sparseSrrTResized);
  mesh->o_sparseSrsT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), sparseSrsTResized);
  mesh->o_sparseSssT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), sparseSssTResized);
  //  mesh->o_sparseSrrT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSrrT);
  //  mesh->o_sparseSrsT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSrsT);
  //  mesh->o_sparseSssT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSssT);

  mesh->o_sparseStackedNZ = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(char), IndTchar);

  free(IndTchar);
  free(sparseSrrTResized);
  free(sparseSrsTResized);
  free(sparseSssTResized);
  return solver;

  //special case 

  //if (mesh->SparseNnzPerRow>=mesh->Np){
  //pad it
  mesh->SparseNnzPerRowNonPadded = mesh->SparseNnzPerRow;
  mesh->SparseNnzPerRow = mesh->SparseNnzPerRow + 4-(mesh->SparseNnzPerRow%4);
  //}

  int * India = (int*) calloc((mesh->Np)+1, sizeof(int));

  char * Indja = (char*) calloc(mesh->Np*mesh->SparseNnzPerRow, sizeof(char));
  dfloat * Srr = (dfloat *) calloc(mesh->SparseNnzPerRow*mesh->Np,sizeof(dfloat));
  dfloat * Srs = (dfloat *) calloc(mesh->SparseNnzPerRow*mesh->Np,sizeof(dfloat));
  dfloat * Sss = (dfloat *) calloc(mesh->SparseNnzPerRow*mesh->Np,sizeof(dfloat));
  printf("arraySize %d\n", mesh->SparseNnzPerRow*mesh->Np);

  //CSR setup(brute force, no padding)
  /*
     int count = -1;
     India[0] = 1;
     for (int m=0;m<mesh->Np;m++) {
     int countRow = 0;

     for (int n=0;n<mesh->SparseNnzPerRow;n++) {
     if (mesh->sparseStackedNZ[m+mesh->Np*n]){
     count++;
     countRow++;
     Indja[count] = (char)mesh->sparseStackedNZ[m+mesh->Np*n];
     Srr[count] =  mesh->sparseSrrT[m+mesh->Np*n];
     Srs[count] =  mesh->sparseSrsT[m+mesh->Np*n];
     Sss[count] =  mesh->sparseSssT[m+mesh->Np*n];

     }
     }
     India[m+1] = India[m]+countRow;
     }*/
  //CSR padded
  int count = -1;
  India[0] = 1;
  for (int m=0;m<mesh->Np;m++) {
    int countRow = 0;
    for (int n=0;n<mesh->SparseNnzPerRowNonPadded;n++) {
      if (mesh->sparseStackedNZ[m+mesh->Np*n]){
        count++;
        countRow++;

        Indja[count] = (char)mesh->sparseStackedNZ[m+mesh->Np*n];
        Srr[count] =  mesh->sparseSrrT[m+mesh->Np*n];
        Srs[count] =  mesh->sparseSrsT[m+mesh->Np*n];
        Sss[count] =  mesh->sparseSssT[m+mesh->Np*n];

      }
    }
    printf("this is row %d, I counted %d \n", m, countRow);
    while (countRow%4){
      //add some extra zeros
      count++;
      countRow++; 
      Srr[count] = 0.0f;
      Srs[count] = 0.0f;
      Sss[count] = 0.0f;
      printf("putting 0 in place %d\n", count);
      Indja[count] =(char) 0;
    }
    India[m+1] = India[m]+countRow;
  }

  for (int m=0; m<mesh->Np+1; ++m){
    printf( "%d, ", India[m]);
    if (m<mesh->Np){
      printf("\n Row %d \n", m);
      for (int k=India[m]-1; k<India[m+1]-1; ++k){
        printf(" %d, ", Indja[k]);      
      }
      printf("\n");

      for (int k=India[m]-1; k<India[m+1]-1; ++k){
        printf(" %.16f, ", Srr[k]);      
      }
      printf("\n"); 
      for (int k=India[m]-1; k<India[m+1]-1; ++k){
        printf(" %.16f, ", Srs[k]);      
      }
      printf("\n");
      for (int k=India[m]-1; k<India[m+1]-1; ++k){
        printf(" %.16f, ", Sss[k]);      
      }


      printf("\n\n");
    }
  }


  printf("trying to allocate, size %d, nnz %d, max Nz per row %d \n", mesh->Np+1, India[mesh->Np], mesh->SparseNnzPerRow);
  //comment
  //comment
  //
  //  mesh->o_Srr = mesh->device.malloc((India[mesh->Np])*sizeof(dfloat), Srr);
  //  mesh->o_Srs = mesh->device.malloc((India[mesh->Np])*sizeof(dfloat), Srs);
  //  mesh->o_Sss = mesh->device.malloc((India[mesh->Np])*sizeof(dfloat), Sss);
  mesh->o_Srr = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), Srr);
  mesh->o_Srs = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), Srs);
  mesh->o_Sss = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), Sss);
  printf(" allocated 0\n");

  //  mesh->o_Indja = mesh->device.malloc((India[mesh->Np])*sizeof(char), mesh->SparseNnzPerRow*mesh->Np);
  mesh->o_Indja = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(char), Indja);
  printf(" allocated 1\n");

  mesh->o_India = mesh->device.malloc(((mesh->Np)+1)*sizeof(int), India);
  printf(" allocated 2\n");

  // mesh->SparseNnzPerRowNonPadded = mesh->SparseNnzPerRow;

  kernelInfo.addDefine("p_SparseNnzPerRow", mesh->SparseNnzPerRow);
  kernelInfo.addDefine("p_MaxChar4", mesh->SparseNnzPerRow/4);
  printf("maxnnz = %d maxchar4 = %d \n", mesh->SparseNnzPerRow, mesh->SparseNnzPerRow/4);
  kernelInfo.addDefine("p_NnzTotal", India[mesh->Np]-1);
  //printf("info: Nnodes = %d Nblocks = %d Nngeo = %d ", );

#else
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
            //  rowData[6*i + 2*(nogroups-1)+2] = 1;
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
      //    rowDataTransposed[mesh->Np*j+i]  = rowData[8*i+j];
       rowDataTransposed[4*i+j] = (char) rowData[8*i+j];    
  /*    if (j==0){
        rowDataTransposed[4*i+0] = (char) 1;
        rowDataTransposed[4*i+1] = (char) 4;
        rowDataTransposed[4*i+2] = (char) 1;
        rowDataTransposed[4*i+3] = (char) 0;  
      }*/
      //      printf("putting %d in place %d \n", rowData[8*i+j],4*i+j );  
    }
    for (int j=0; j<4; j++){
      rowDataTransposed[mesh->Np*4 + 4*i+j] = (char) rowData[8*i+j+4];
      /*if (j==0){
        rowDataTransposed[mesh->Np*4 + 4*i+0] = (char) 1;
        rowDataTransposed[mesh->Np*4 + 4*i+1] = (char) 0;
        rowDataTransposed[mesh->Np*4 + 4*i+2] = (char) 0;
        rowDataTransposed[mesh->Np*4 + 4*i+3] = (char) 0;  
      }*/

      //    printf("putting %d in place %d \n", rowData[8*i+j],mesh->Np*4 + 4*i+j );

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

  printf("\n MAX IS %d \n", mesh->Np*8);
  mesh->o_rowData = mesh->device.malloc(mesh->Np*8*sizeof(char), rowDataTransposed);
  printf("alloc!\n");  
  mesh->o_Srr = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat),mesh->sparseSrrT);
  mesh->o_Srs = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSrsT);
  mesh->o_Sss = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSssT);  
  kernelInfo.addDefine("p_SparseNnzPerRow", mesh->SparseNnzPerRow);


  return solver;
#endif
}



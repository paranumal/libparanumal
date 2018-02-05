#include "ellipticTri2D.h"

void ellipticSetupTri2D(mesh2D *mesh, occa::kernelInfo &kernelInfo, const char* options){

  mesh->Nfields = 1;

  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", rank%2);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");

  if (strstr(options,"SPARSE")) {//connectivity fix for the sparse basis 
    //connect the face modes between each element 
    meshConnectFaceModes2D(mesh,mesh->FaceModes,mesh->sparseV);
    //use the mmaps constructed and overwrite vmap and FaceNodes
    for (iint n=0;n<mesh->Nfp*mesh->Nfaces*mesh->Nelements;n++) {
      mesh->vmapM[n] = mesh->mmapM[n];
      mesh->vmapP[n] = mesh->mmapP[n];
    }
    for (int n=0;n<mesh->Nfaces*mesh->Nfp;n++) { //overwrite facenodes
      mesh->faceNodes[n] = mesh->FaceModes[n];
    }
    for (int n=0;n<mesh->Nverts;n++) { //overwrite vertex nodes (assumes their first in the list)
      mesh->vertexNodes[n] = n;
    }
    //free the old gather scatter arrays and re-run the connect nodes function using this updated connectivity
    free(mesh->gatherLocalIds );
    free(mesh->gatherBaseIds  );
    free(mesh->gatherBaseRanks);
    free(mesh->gatherHaloFlags);
    free(mesh->globalIds);
    free(mesh->globalOwners);
    free(mesh->globalHaloFlags);
    meshParallelConnectNodes(mesh);
  }

  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);
}

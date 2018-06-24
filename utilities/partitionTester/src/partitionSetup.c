#include "partition.h"

void partitionSetup(mesh_t *mesh){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  mesh->Nfields = 4;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  // fix this later (initial conditions)
  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];

      mesh->q[0+mesh->Nfields*(n+e*mesh->Np)] = rank;
      
    }
  }

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", 0);

  occa::kernelInfo kernelInfo;

  // generic occa device set up
  if(mesh->dim==2)
    meshOccaSetup2D(mesh, deviceConfig, kernelInfo);
  else
    meshOccaSetup3D(mesh, deviceConfig, kernelInfo);
  
}

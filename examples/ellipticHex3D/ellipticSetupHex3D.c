#include "ellipticHex3D.h"

void ellipticSetupHex3D(mesh3D *mesh, occa::kernelInfo &kernelInfo){

  mesh->Nfields = 1;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  iint rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);

  void meshOccaSetup3D(mesh3D *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo);
  meshOccaSetup3D(mesh, deviceConfig, kernelInfo);
}

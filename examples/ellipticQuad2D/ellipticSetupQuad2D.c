#include "ellipticQuad2D.h"

void ellipticSetupQuad2D(mesh2D *mesh, occa::kernelInfo &kernelInfo){

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

  void meshOccaSetup2D(mesh2D *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo);
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);
}

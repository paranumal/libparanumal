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

#include "massHex3D.h"

void massSetupHex3D(mesh3D *mesh, occa::kernelInfo &kernelInfo){

  mesh->Nfields = 1;

  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // use rank to choose DEVICE

  printf("rank %d hostid %d \n", rank, gethostid());

  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", rank%2);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = %d, platformID = 0", rank%2);
  //  sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");

  void meshOccaSetup3D(mesh3D *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo);
  meshOccaSetup3D(mesh, deviceConfig, kernelInfo);
}

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

#include "advection.hpp"

advection_t& advection_t::Setup(mesh_t& mesh, linAlg_t& linAlg){

  advection_t* advection = new advection_t(mesh, linAlg);

  settings_t& settings = advection->settings;

  //setup timeStepper
  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;
  advection->timeStepper = timeStepper_t::Setup(Nlocal, Nhalo, *advection);
  advection->timeStepper->Init();

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat cfl = 0.5; // depends on the stability region size

  dfloat dt = cfl*hmin/((mesh.N+1.)*(mesh.N+1.));
  advection->timeStepper->SetTimeStep(dt);

  //setup linear algebra module
  advection->linAlg.InitKernels({"innerProd"}, mesh.comm);

  /*setup trace halo exchange */
  advection->traceHalo = mesh.HaloTraceSetup(1); //one field

  // compute samples of q at interpolation nodes
  advection->q = (dfloat*) calloc(Nlocal+Nhalo, sizeof(dfloat));
  advection->o_q = mesh.device.malloc((Nlocal+Nhalo)*sizeof(dfloat), advection->q);

  //storage for M*q during reporting
  advection->o_Mq = mesh.device.malloc((Nlocal+Nhalo)*sizeof(dfloat), advection->q);

  // OCCA build stuff
  occa::properties kernelInfo = advection->props; //copy base occa properties

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= 1;

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = 1024/mesh.Np; // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int cubMaxNodes = mymax(mesh.Np, (mesh.intNfp*mesh.Nfaces));
  kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
  int cubMaxNodes1 = mymax(mesh.Np, (mesh.intNfp));
  kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  // set kernel name suffix
  char *suffix;
  if(mesh.elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(mesh.elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(mesh.elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(mesh.elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  // kernels from volume file
  sprintf(fileName, DADVECTION "/okl/advectionVolume%s.okl", suffix);
  sprintf(kernelName, "advectionVolume%s", suffix);

  advection->volumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);
  // kernels from surface file
  sprintf(fileName, DADVECTION "/okl/advectionSurface%s.okl", suffix);
  sprintf(kernelName, "advectionSurface%s", suffix);

  advection->surfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);

  // mass matrix operator
  sprintf(fileName, LIBP_DIR "/core/okl/MassMatrixOperator%s.okl", suffix);
  sprintf(kernelName, "MassMatrixOperator%s", suffix);

  advection->MassMatrixKernel = buildKernel(mesh.device, fileName, kernelName,
                                            kernelInfo, mesh.comm);


  if (mesh.dim==2) {
    sprintf(fileName, DADVECTION "/okl/advectionInitialCondition2D.okl");
    sprintf(kernelName, "advectionInitialCondition2D");
  } else {
    sprintf(fileName, DADVECTION "/okl/advectionInitialCondition3D.okl");
    sprintf(kernelName, "advectionInitialCondition3D");
  }

  advection->initialConditionKernel = buildKernel(mesh.device, fileName, kernelName,
                                                  kernelInfo, mesh.comm);

  return *advection;
}

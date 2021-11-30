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

#include "hjs.hpp"

hjs_t& hjs_t::Setup(platform_t& platform, mesh_t& mesh,
                                 hjsSettings_t& settings){

  hjs_t* hjs = new hjs_t(platform, mesh, settings);

  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;

  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    hjs->timeStepper = new TimeStepper::ab3(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, 1, *hjs);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    hjs->timeStepper = new TimeStepper::lserk4(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, 1, *hjs);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    hjs->timeStepper = new TimeStepper::dopri5(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, 1, *hjs, mesh.comm);
  }

  //setup linear algebra module
  platform.linAlg.InitKernels({"innerProd", "max"});

  /*setup trace halo exchange */
  hjs->traceHalo = mesh.HaloTraceSetup(1); //one field

  // compute samples of q at interpolation nodes
  hjs->q = (dfloat*) calloc(Nlocal+Nhalo, sizeof(dfloat));
  hjs->o_q = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), hjs->q);

  // // storage for the pi, qi 
  // hjs->gradq = (dfloat*) calloc((Nlocal+Nhalo)*mesh.dim, sizeof(dfloat));
  // hjs->o_gradq = platform.malloc((Nlocal+Nhalo)*mesh.dim*sizeof(dfloat), hjs->gradq);

  // storage for the pi, qi 
  hjs->gradq = (dfloat*) calloc(2*Nlocal*mesh.dim, sizeof(dfloat));
  hjs->o_gradq = platform.malloc(2*Nlocal*mesh.dim*sizeof(dfloat), hjs->gradq);



  //storage for M*q during reporting
  hjs->o_Mq = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), hjs->q);
  mesh.MassMatrixKernelSetup(1); // mass matrix operator

  // OCCA build stuff
  occa::properties kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= 1;

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  // int NblockV = mymax(1, blockMax/mesh.Np);
  // kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockV = 1;
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

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
  sprintf(fileName, DHJS "/okl/hjsVolume%s.okl", suffix);
  sprintf(kernelName, "hjsVolume%s", suffix);

  hjs->volumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);
  // kernels from surface file
  sprintf(fileName, DHJS "/okl/hjsSurface%s.okl", suffix);
  sprintf(kernelName, "hjsSurface%s", suffix);

  hjs->surfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);


  // kernels from surface file
  sprintf(fileName, DHJS "/okl/hjsHamiltonian%s.okl", suffix);
  sprintf(kernelName, "hjsHamiltonian%s", suffix);

  hjs->hamiltonianKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  if (mesh.dim==2) {
    sprintf(fileName, DHJS "/okl/hjsInitialCondition2D.okl");
    sprintf(kernelName, "hjsInitialCondition2D");
  } else {
    sprintf(fileName, DHJS "/okl/hjsInitialCondition3D.okl");
    sprintf(kernelName, "hjsInitialCondition3D");
  }

  hjs->initialConditionKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  // sprintf(fileName, DADVECTION "/okl/advectionMaxWaveSpeed%s.okl", suffix);
  // sprintf(kernelName, "advectionMaxWaveSpeed%s", suffix);

  // advection->maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  return *hjs;
}

hjs_t::~hjs_t() {
  volumeKernel.free();
  surfaceKernel.free();
  initialConditionKernel.free();
  // maxWaveSpeedKernel.free();

  if (timeStepper) delete timeStepper;
  if (traceHalo) traceHalo->Free();
}

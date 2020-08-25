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

#include "acoustics.hpp"

acoustics_t& acoustics_t::Setup(platform_t& platform, mesh_t& mesh,
                                acousticsSettings_t& settings){

  acoustics_t* acoustics = new acoustics_t(platform, mesh, settings);

  acoustics->Nfields = (mesh.dim==3) ? 4:3;

  dlong Nlocal = mesh.Nelements*mesh.Np*acoustics->Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*acoustics->Nfields;

  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    acoustics->timeStepper = new TimeStepper::ab3(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, acoustics->Nfields, *acoustics);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    acoustics->timeStepper = new TimeStepper::lserk4(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, acoustics->Nfields, *acoustics);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    acoustics->timeStepper = new TimeStepper::dopri5(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, acoustics->Nfields, *acoustics);
  }

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat cfl = 0.5; // depends on the stability region size

  dfloat dt = cfl*hmin/((mesh.N+1.)*(mesh.N+1.));
  acoustics->timeStepper->SetTimeStep(dt);

  //setup linear algebra module
  platform.linAlg.InitKernels({"innerProd"});

  // set penalty parameter
  dfloat Lambda2 = 0.5;

  /*setup trace halo exchange */
  acoustics->traceHalo = mesh.HaloTraceSetup(acoustics->Nfields);

  // compute samples of q at interpolation nodes
  acoustics->q = (dfloat*) calloc(Nlocal+Nhalo, sizeof(dfloat));
  acoustics->o_q = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), acoustics->q);

  //storage for M*q during reporting
  acoustics->o_Mq = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), acoustics->q);
  mesh.MassMatrixKernelSetup(acoustics->Nfields); // mass matrix operator

  // OCCA build stuff
  occa::properties kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;


  kernelInfo["defines/" "p_Nfields"]= acoustics->Nfields;

  const dfloat p_half = 1./2.;
  kernelInfo["defines/" "p_half"]= p_half;

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = mymax(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  kernelInfo["defines/" "p_Lambda2"]= Lambda2;

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
  sprintf(fileName, DACOUSTICS "/okl/acousticsVolume%s.okl", suffix);
  sprintf(kernelName, "acousticsVolume%s", suffix);

  acoustics->volumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  // kernels from surface file
  sprintf(fileName, DACOUSTICS "/okl/acousticsSurface%s.okl", suffix);
  sprintf(kernelName, "acousticsSurface%s", suffix);

  acoustics->surfaceKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  if (mesh.dim==2) {
    sprintf(fileName, DACOUSTICS "/okl/acousticsInitialCondition2D.okl");
    sprintf(kernelName, "acousticsInitialCondition2D");
  } else {
    sprintf(fileName, DACOUSTICS "/okl/acousticsInitialCondition3D.okl");
    sprintf(kernelName, "acousticsInitialCondition3D");
  }

  acoustics->initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);

  return *acoustics;
}

acoustics_t::~acoustics_t() {
  volumeKernel.free();
  surfaceKernel.free();
  initialConditionKernel.free();

  if (timeStepper) delete timeStepper;
  if (traceHalo) traceHalo->Free();
}

/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

void acoustics_t::Setup(platform_t& _platform, mesh_t& _mesh,
                        acousticsSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  Nfields = (mesh.dim==3) ? 4:3;

  dlong Nlocal = mesh.Nelements*mesh.Np*Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*Nfields;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(Nfields);

  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
                                        mesh.totalHaloPairs,
                                        mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  }

  // set penalty parameter
  dfloat Lambda2 = 0.5;

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal+Nhalo);
  o_q = platform.malloc<dfloat>(q);

  //storage for M*q during reporting
  o_Mq = platform.malloc<dfloat>(q);
  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= Nfields;

  const dfloat p_half = 1./2.;
  kernelInfo["defines/" "p_half"]= p_half;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  kernelInfo["defines/" "p_Lambda2"]= Lambda2;

  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES)
    suffix = "Tri2D";
  if(mesh.elementType==Mesh::QUADRILATERALS)
    suffix = "Quad2D";
  if(mesh.elementType==Mesh::TETRAHEDRA)
    suffix = "Tet3D";
  if(mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "Hex3D";

  std::string oklFilePrefix = DACOUSTICS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // kernels from volume file
  fileName   = oklFilePrefix + "acousticsVolume" + suffix + oklFileSuffix;
  kernelName = "acousticsVolume" + suffix;

  volumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  // kernels from surface file
  fileName   = oklFilePrefix + "acousticsSurface" + suffix + oklFileSuffix;
  kernelName = "acousticsSurface" + suffix;

  surfaceKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "acousticsInitialCondition2D" + oklFileSuffix;
    kernelName = "acousticsInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "acousticsInitialCondition3D" + oklFileSuffix;
    kernelName = "acousticsInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);
}

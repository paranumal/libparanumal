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

#include "advection.hpp"

void advection_t::Setup(platform_t& _platform, mesh_t& _mesh,
                         advectionSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = mesh.comm;
  settings = _settings;

  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd", "max"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(1); //one field

  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
                                        mesh.totalHaloPairs,
                                        mesh.Np, 1, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, 1, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, 1, platform, comm);
  }

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal+Nhalo);
  o_q = platform.malloc<dfloat>(q);

  //storage for M*q during reporting
  o_Mq = platform.malloc<dfloat>(q);
  mesh.MassMatrixKernelSetup(1); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

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

  std::string oklFilePrefix = DADVECTION "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // kernels from volume file
  fileName   = oklFilePrefix + "advectionVolume" + suffix + oklFileSuffix;
  kernelName = "advectionVolume" + suffix;

  volumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

  // kernels from surface file
  fileName   = oklFilePrefix + "advectionSurface" + suffix + oklFileSuffix;
  kernelName = "advectionSurface" + suffix;

  surfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "advectionInitialCondition2D" + oklFileSuffix;
    kernelName = "advectionInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "advectionInitialCondition3D" + oklFileSuffix;
    kernelName = "advectionInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "advectionMaxWaveSpeed" + suffix + oklFileSuffix;
  kernelName = "advectionMaxWaveSpeed" + suffix;

  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
}

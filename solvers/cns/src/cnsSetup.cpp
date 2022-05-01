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

#include "cns.hpp"

void cns_t::Setup(platform_t& _platform, mesh_t& _mesh,
                  cnsSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //get physical paramters
  settings.getSetting("VISCOSITY", mu);
  settings.getSetting("GAMMA", gamma);

  cubature   = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;
  isothermal = (settings.compareSetting("ISOTHERMAL", "TRUE")) ? 1:0;

  //setup cubature
  if (cubature) {
    mesh.CubatureSetup();
    mesh.CubaturePhysicalNodes();
  }

  Nfields   = (mesh.dim==3) ? 4:3;
  Ngrads = mesh.dim*mesh.dim;

  if (!isothermal) Nfields++; //include energy equation

  dlong NlocalFields = mesh.Nelements*mesh.Np*Nfields;
  dlong NhaloFields  = mesh.totalHaloPairs*mesh.Np*Nfields;
  dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
  dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;

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

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd", "max"});

  /*setup trace halo exchange */
  fieldTraceHalo = mesh.HaloTraceSetup(Nfields);
  gradTraceHalo  = mesh.HaloTraceSetup(Ngrads);

  // compute samples of q at interpolation nodes
  q.malloc(NlocalFields+NhaloFields);
  o_q = platform.malloc<dfloat>(q);

  gradq.malloc(NlocalGrads+NhaloGrads);
  o_gradq = platform.malloc<dfloat>(gradq);

  Vort.malloc(mesh.dim*mesh.Nelements*mesh.Np);
  o_Vort = platform.malloc<dfloat>(Vort);

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
  kernelInfo["defines/" "p_Ngrads"]= Ngrads;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  if (cubature) {
    int cubMaxNodes = std::max(mesh.Np, (mesh.intNfp*mesh.Nfaces));
    kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
    int cubMaxNodes1 = std::max(mesh.Np, (mesh.intNfp));
    kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

    int cubNblockV = std::max(1, blockMax/mesh.cubNp);
    kernelInfo["defines/" "p_cubNblockV"]= cubNblockV;

    int cubNblockS = std::max(1, blockMax/cubMaxNodes);
    kernelInfo["defines/" "p_cubNblockS"]= cubNblockS;
  }

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

  std::string oklFilePrefix = DCNS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  if (isothermal) {
    if (cubature) {
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsIsothermalCubatureVolume" + suffix + oklFileSuffix;
      kernelName = "cnsIsothermalCubatureVolume" + suffix;

      cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
      // kernels from surface file
      fileName   = oklFilePrefix + "cnsIsothermalCubatureSurface" + suffix + oklFileSuffix;
      kernelName = "cnsIsothermalCubatureSurface" + suffix;

      cubatureSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
    } else {
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsIsothermalVolume" + suffix + oklFileSuffix;
      kernelName = "cnsIsothermalVolume" + suffix;

      volumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      // kernels from surface file
      fileName   = oklFilePrefix + "cnsIsothermalSurface" + suffix + oklFileSuffix;
      kernelName = "cnsIsothermalSurface" + suffix;

      surfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    }
  } else {
    if (cubature) {
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsCubatureVolume" + suffix + oklFileSuffix;
      kernelName = "cnsCubatureVolume" + suffix;

      cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
      // kernels from surface file
      fileName   = oklFilePrefix + "cnsCubatureSurface" + suffix + oklFileSuffix;
      kernelName = "cnsCubatureSurface" + suffix;

      cubatureSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
    } else {
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsVolume" + suffix + oklFileSuffix;
      kernelName = "cnsVolume" + suffix;

      volumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      // kernels from surface file
      fileName   = oklFilePrefix + "cnsSurface" + suffix + oklFileSuffix;
      kernelName = "cnsSurface" + suffix;

      surfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    }
  }

  // kernels from volume file
  fileName   = oklFilePrefix + "cnsGradVolume" + suffix + oklFileSuffix;
  kernelName = "cnsGradVolume" + suffix;

  gradVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  // kernels from surface file
  fileName   = oklFilePrefix + "cnsGradSurface" + suffix + oklFileSuffix;
  kernelName = "cnsGradSurface" + suffix;

  gradSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

  // vorticity calculation
  fileName   = oklFilePrefix + "cnsVorticity" + suffix + oklFileSuffix;
  kernelName = "cnsVorticity" + suffix;

  vorticityKernel = platform.buildKernel(fileName, kernelName,
                                     kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "cnsInitialCondition2D" + oklFileSuffix;
    if (isothermal)
      kernelName = "cnsIsothermalInitialCondition2D";
    else
      kernelName = "cnsInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "cnsInitialCondition3D" + oklFileSuffix;
    if (isothermal)
      kernelName = "cnsIsothermalInitialCondition3D";
    else
      kernelName = "cnsInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                            kernelInfo);

  fileName   = oklFilePrefix + "cnsMaxWaveSpeed" + suffix + oklFileSuffix;
  if (isothermal) {
    kernelName = "cnsIsothermalMaxWaveSpeed" + suffix;
  } else {
    kernelName = "cnsMaxWaveSpeed" + suffix;
  }

  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName,
                                            kernelInfo);
}

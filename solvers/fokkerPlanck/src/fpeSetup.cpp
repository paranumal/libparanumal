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

#include "fpe.hpp"

void fpe_t::Setup(platform_t& _platform, mesh_t& _mesh,
                  fpeSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  settings.getSetting("VISCOSITY", mu);

  cubature = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;

  //setup cubature
  if (cubature) {
    mesh.CubatureSetup();
    mesh.CubaturePhysicalNodes();
  }

  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;

  //setup timeStepper
  dfloat gamma = 0.0;
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
  } else if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")){
    timeStepper.Setup<TimeStepper::extbdf3>(mesh.Nelements,
                                            mesh.totalHaloPairs,
                                            mesh.Np, 1, platform, comm);
    gamma = timeStepper.GetGamma();
  } else if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")){
    timeStepper.Setup<TimeStepper::ssbdf3>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, 1, platform, comm);
    gamma = timeStepper.GetGamma();
  }

  Nsubcycles=1;
  if (settings.compareSetting("TIME INTEGRATOR","SSBDF3"))
    settings.getSetting("NUMBER OF SUBCYCLES", Nsubcycles);

  //Setup Elliptic solver
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    ||settings.compareSetting("TIME INTEGRATOR","SSBDF3")){

    int NBCTypes = 7;
    memory<int> BCType(NBCTypes);
    // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
    BCType[0] = 0;
    BCType[1] = 1;
    BCType[2] = 1;
    BCType[3] = 2;
    BCType[4] = 1;
    BCType[5] = 1;
    BCType[6] = 1;

    ellipticSettings = _settings.extractEllipticSettings();

    //make a guess at dt for the lambda value
    //TODO: we should allow preconditioners to be re-setup if lambda is updated
    dfloat hmin = mesh.MinCharacteristicLength();
    dfloat dtAdvc = Nsubcycles*hmin/((mesh.N+1.)*(mesh.N+1.));
    dfloat lambda = gamma/(dtAdvc*mu);

    elliptic.Setup(platform, mesh, ellipticSettings,
                   lambda, NBCTypes, BCType);
    tau = elliptic.tau;

    if (ellipticSettings.compareSetting("LINEAR SOLVER","NBPCG")){
      linearSolver.Setup<LinearSolver::nbpcg>(elliptic.Ndofs, elliptic.Nhalo,
                                              platform, ellipticSettings, comm);
    } else if (ellipticSettings.compareSetting("LINEAR SOLVER","NBFPCG")){
      linearSolver.Setup<LinearSolver::nbfpcg>(elliptic.Ndofs, elliptic.Nhalo,
                                              platform, ellipticSettings, comm);
    } else if (ellipticSettings.compareSetting("LINEAR SOLVER","PCG")){
      linearSolver.Setup<LinearSolver::pcg>(elliptic.Ndofs, elliptic.Nhalo,
                                              platform, ellipticSettings, comm);
    } else if (ellipticSettings.compareSetting("LINEAR SOLVER","PGMRES")){
      linearSolver.Setup<LinearSolver::pgmres>(elliptic.Ndofs, elliptic.Nhalo,
                                              platform, ellipticSettings, comm);
    } else if (ellipticSettings.compareSetting("LINEAR SOLVER","PMINRES")){
      linearSolver.Setup<LinearSolver::pminres>(elliptic.Ndofs, elliptic.Nhalo,
                                              platform, ellipticSettings, comm);
    }
  } else {
    //set penalty
    if (mesh.elementType==Mesh::TRIANGLES ||
        mesh.elementType==Mesh::QUADRILATERALS){
      tau = 2.0*(mesh.N+1)*(mesh.N+2)/2.0;
      if(mesh.dim==3)
        tau *= 1.5;
    } else
      tau = 2.0*(mesh.N+1)*(mesh.N+3);
  }

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd", "axpy", "max"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(1); //one field

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal+Nhalo, 0.0);
  o_q = platform.malloc<dfloat>(q);

  //storage for M*q during reporting
  o_Mq = platform.malloc<dfloat>(q);
  mesh.MassMatrixKernelSetup(1); // mass matrix operator

  grad.malloc((Nlocal+Nhalo)*4, 0.0);
  o_grad  = platform.malloc<dfloat>(grad);

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= 1;

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

  std::string oklFilePrefix = DFPE "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // advection kernels
  if (cubature) {
    fileName   = oklFilePrefix + "fpeCubatureAdvection" + suffix + oklFileSuffix;
    kernelName = "fpeAdvectionCubatureVolume" + suffix;
    advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
    kernelName = "fpeAdvectionCubatureSurface" + suffix;
    advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    fileName   = oklFilePrefix + "fpeAdvection" + suffix + oklFileSuffix;
    kernelName = "fpeAdvectionVolume" + suffix;
    advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
    kernelName = "fpeAdvectionSurface" + suffix;
    advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }


  // diffusion kernels
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    ||settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    fileName   = oklFilePrefix + "fpeDiffusionRhs" + suffix + oklFileSuffix;
    kernelName = "fpeDiffusionRhs" + suffix;
    diffusionRhsKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    // gradient kernel
    fileName   = oklFilePrefix + "fpeGradient" + suffix + oklFileSuffix;
    kernelName = "fpeGradient" + suffix;
    gradientKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    fileName   = oklFilePrefix + "fpeDiffusion" + suffix + oklFileSuffix;
    kernelName = "fpeDiffusion" + suffix;
    diffusionKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "fpeInitialCondition2D" + oklFileSuffix;
    kernelName = "fpeInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "fpeInitialCondition3D" + oklFileSuffix;
    kernelName = "fpeInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);

  fileName   = oklFilePrefix + "fpeMaxWaveSpeed" + suffix + oklFileSuffix;
  kernelName = "fpeMaxWaveSpeed" + suffix;

  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  //build subcycler
  if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    subcycler.platform = platform;
    subcycler.mesh = mesh;
    subcycler.comm = comm;
    subcycler.settings = settings;

    subcycler.cubature = cubature;
    subcycler.traceHalo = traceHalo;
    subcycler.advectionVolumeKernel = advectionVolumeKernel;
    subcycler.advectionSurfaceKernel = advectionSurfaceKernel;

    if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","AB3")){
      subStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
                                         mesh.totalHaloPairs,
                                         mesh.Np, 1, platform, comm);
    } else if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","LSERK4")){
      subStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                            mesh.totalHaloPairs,
                                            mesh.Np, 1, platform, comm);
    } else if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","DOPRI5")){
      subStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
                                            mesh.totalHaloPairs,
                                            mesh.Np, 1, platform, comm);
    }
  }
}

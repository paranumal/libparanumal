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

#include "bns.hpp"

void bns_t::Setup(platform_t& _platform, mesh_t& _mesh,
                  bnsSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  //get physical paramters
  settings.getSetting("SPEED OF SOUND", c);
  settings.getSetting("VISCOSITY", nu);
  RT     = c*c;
  tauInv = RT/nu;

  Nfields    = (mesh.dim==3) ? 10:6;
  Npmlfields = mesh.dim*Nfields;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //setup cubature
  mesh.CubatureSetup();

  //Setup PML
  PmlSetup();

  //setup timeStepper
  dlong Nlocal = mesh.Nelements*mesh.Np*Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*Nfields;

  semiAnalytic = 0;
  if (settings.compareSetting("TIME INTEGRATOR","SARK4")
    ||settings.compareSetting("TIME INTEGRATOR","SARK5")
    ||settings.compareSetting("TIME INTEGRATOR","SAAB3")
    ||settings.compareSetting("TIME INTEGRATOR","MRSAAB3"))
    semiAnalytic = 1;

  //semi-analytic exponential coefficients
  memory<dfloat> lambda(Nfields);
  for (int i=0;i<mesh.dim+1;i++) lambda[i] = 0.0;
  for (int i=mesh.dim+1;i<Nfields;i++) lambda[i] = -tauInv;

  //make array of time step estimates for each element
  memory<dfloat> EtoDT(mesh.Nelements);
  dfloat vmax = MaxWaveSpeed();
  for(dlong e=0;e<mesh.Nelements;++e){
    dfloat h = mesh.ElementCharacteristicLength(e);
    dfloat dtAdv  = h/(vmax*(mesh.N+1.)*(mesh.N+1.));
    dfloat dtVisc = 1.0/tauInv;

    if (semiAnalytic)
      EtoDT[e] = dtAdv;
    else
      EtoDT[e] = std::min(dtAdv, dtVisc);

    /*
    Artificial warping of time step size for multirate testing
    */
#if 0
    dfloat x=0., y=0, z=0;
    for (int v=0;v<mesh.Nverts;v++){
      x += mesh.EX[e*mesh.Nverts+v];
      y += mesh.EY[e*mesh.Nverts+v];
      if (mesh.dim==3)
        z += mesh.EZ[e*mesh.Nverts+v];
    }
    x /=mesh.Nverts;
    y /=mesh.Nverts;
    z /=mesh.Nverts;

    dfloat c = mymin(fabs(x),fabs(y));
    if (mesh.dim==3)
      c = mymin(c,fabs(z));

    c = std::max(0.5, c);
    EtoDT[e] *= c;
#endif
  }

  mesh.mrNlevels=0;
  if (settings.compareSetting("TIME INTEGRATOR","MRAB3") ||
      settings.compareSetting("TIME INTEGRATOR","MRSAAB3")) {
    mesh.MultiRateSetup(EtoDT);
    mesh.MultiRatePmlSetup();
    multirateTraceHalo = mesh.MultiRateHaloTraceSetup(Nfields);
  }

  if (settings.compareSetting("TIME INTEGRATOR","MRAB3")){
    timeStepper.Setup<TimeStepper::mrab3_pml>(mesh.Nelements, mesh.NpmlElements,
                                              mesh.totalHaloPairs,
                                              mesh.Np, Nfields, Npmlfields,
                                              platform, mesh);
  } else if (settings.compareSetting("TIME INTEGRATOR","MRSAAB3")){
    timeStepper.Setup<TimeStepper::mrsaab3_pml>(mesh.Nelements, mesh.NpmlElements,
                                                mesh.totalHaloPairs,
                                                mesh.Np, Nfields, Npmlfields,
                                                lambda, platform, mesh);
  } else if (settings.compareSetting("TIME INTEGRATOR","SAAB3")) {
    timeStepper.Setup<TimeStepper::saab3_pml>(mesh.Nelements, mesh.NpmlElements,
                                              mesh.totalHaloPairs,
                                              mesh.Np, Nfields, Npmlfields,
                                              lambda, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3_pml>(mesh.Nelements, mesh.NpmlElements,
                                            mesh.totalHaloPairs,
                                            mesh.Np, Nfields, Npmlfields,
                                            platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4_pml>(mesh.Nelements, mesh.NpmlElements,
                                               mesh.totalHaloPairs,
                                               mesh.Np, Nfields, Npmlfields,
                                               platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5_pml>(mesh.Nelements, mesh.NpmlElements,
                                               mesh.totalHaloPairs,
                                               mesh.Np, Nfields, Npmlfields,
                                               platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","SARK4")) {
    timeStepper.Setup<TimeStepper::sark4_pml>(mesh.Nelements, mesh.NpmlElements,
                                              mesh.totalHaloPairs,
                                              mesh.Np, Nfields, Npmlfields,
                                              lambda, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","SARK5")) {
    timeStepper.Setup<TimeStepper::sark5_pml>(mesh.Nelements, mesh.NpmlElements,
                                              mesh.totalHaloPairs,
                                              mesh.Np, Nfields, Npmlfields,
                                              lambda, platform, comm);
  } else {
    LIBP_FORCE_ABORT("Requested TIME INTEGRATOR not found.");
  }

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(Nfields);

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal+Nhalo, 0.0);
  o_q = platform.malloc<dfloat>(q);

  Vort.malloc(mesh.dim*mesh.Nelements*mesh.Np, 0.0);
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
  kernelInfo["defines/" "p_Npmlfields"]= Npmlfields;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode()=="CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int NblockCub = std::max(1, blockMax/mesh.cubNp);
  kernelInfo["defines/" "p_NblockCub"]= NblockCub;

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

  std::string oklFilePrefix = DBNS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // kernels from volume file
  fileName   = oklFilePrefix + "bnsVolume" + suffix + oklFileSuffix;
  kernelName = "bnsVolume" + suffix;
  volumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  if (pmlcubature) {
    kernelName = "bnsPmlVolumeCub" + suffix;
    pmlVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  } else {
    kernelName = "bnsPmlVolume" + suffix;
    pmlVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  }

  // kernels from relaxation file
  fileName   = oklFilePrefix + "bnsRelaxation" + suffix + oklFileSuffix;
  kernelName = "bnsRelaxation" + suffix;
  relaxationKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  if (pmlcubature) {
    kernelName = "bnsPmlRelaxationCub" + suffix;
    pmlRelaxationKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    pmlRelaxationKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }


  // kernels from surface file
  fileName   = oklFilePrefix + "bnsSurface" + suffix + oklFileSuffix;
  if (settings.compareSetting("TIME INTEGRATOR","MRAB3") ||
      settings.compareSetting("TIME INTEGRATOR","MRSAAB3")) {
    kernelName = "bnsMRSurface" + suffix;
    surfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    kernelName = "bnsMRPmlSurface" + suffix;
    pmlSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    kernelName = "bnsSurface" + suffix;
    surfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    kernelName = "bnsPmlSurface" + suffix;
    pmlSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }

  // vorticity calculation
  fileName   = oklFilePrefix + "bnsVorticity" + suffix + oklFileSuffix;
  kernelName = "bnsVorticity" + suffix;

  vorticityKernel = platform.buildKernel(fileName, kernelName,
                                     kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "bnsInitialCondition2D" + oklFileSuffix;
    kernelName = "bnsInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "bnsInitialCondition3D" + oklFileSuffix;
    kernelName = "bnsInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                            kernelInfo);
}

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

#include "lbs.hpp"
#define D2Q9 1
#define D3Q15 2

void lbs_t::Setup(platform_t& _platform, mesh_t& _mesh,
                  lbsSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  // Set reference lattice-Boltzmann data  
  latticeSetup();
  
  Npmlfields = mesh.dim*Nfields;

  // AK: not in use yet ... Setup PML
  // PmlSetup();

  //setup timeStepper
  dlong Nlocal = mesh.Nelements*mesh.Np*Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*Nfields;

  //make array of time step estimates for each element
  memory<dfloat> EtoDT(mesh.Nelements);
  dfloat vmax = MaxWaveSpeed();
  for(dlong e=0;e<mesh.Nelements;++e){
    dfloat h = mesh.ElementCharacteristicLength(e);
    dfloat dtAdv  = h/(vmax*(mesh.N+1.)*(mesh.N+1.));
    EtoDT[e] = dtAdv;
  }

  if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields,
                                           platform, comm);
  }else {
    LIBP_FORCE_ABORT("Requested TIME INTEGRATOR not found.");
  }

  
  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(Nfields);

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal+Nhalo, 0.0);
  o_q = platform.malloc<dfloat>(q);

  F.malloc(Nlocal+Nhalo, 0.0);
  o_F = platform.malloc<dfloat>(F);


  Vort.malloc(mesh.dim*mesh.Nelements*mesh.Np, 0.0);
  o_Vort = platform.malloc<dfloat>(Vort);

  // Hold macro quantites i.e. density + velocities
  U.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*Nmacro, 0.0);
  o_U = platform.malloc<dfloat>(U);

  //storage for M*q during reporting
  o_Mq = platform.malloc<dfloat>(U);
  mesh.MassMatrixKernelSetup(Nmacro); // mass matrix operator

  // // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= Nfields;
  // kernelInfo["defines/" "p_Npmlfields"]= Npmlfields;
  kernelInfo["defines/" "p_Nmacro"] = Nmacro;

  kernelInfo["defines/" "p_c"] = c;
  kernelInfo["defines/" "p_ic2"] = 1.0/ pow(c,2);
  kernelInfo["defines/" "p_ic4"] = 1.0/ pow(c,4);

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode()=="CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  // set kernel name suffix
  char *suffix;
  if(mesh.elementType==mesh_t::TRIANGLES)
    suffix = strdup("Tri2D");
  if(mesh.elementType==mesh_t::QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(mesh.elementType==mesh_t::TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(mesh.elementType==mesh_t::HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  if (mesh.dim==2) {
    sprintf(fileName, DLBS "/okl/lbsInitialCondition2D.okl");
    sprintf(kernelName, "lbsInitialCondition2D");
  } else {
    sprintf(fileName, DLBS "/okl/lbsInitialCondition3D.okl");
    sprintf(kernelName, "lbsInitialCondition3D");
  }
  initialConditionKernel = platform.buildKernel(fileName, kernelName,
						     kernelInfo);
  
  // kernels from volume file
  sprintf(fileName, DLBS "/okl/lbsCollision%s.okl", suffix);  

  sprintf(kernelName, "lbsCollision%s", suffix);
  collisionKernel =  platform.buildKernel(fileName, kernelName,
					       kernelInfo);

  sprintf(kernelName, "lbsForcing%s", suffix);
  forcingKernel =  platform.buildKernel(fileName, kernelName,
					     kernelInfo);

  sprintf(kernelName, "lbsMoments%s", suffix);
  momentsKernel =  platform.buildKernel(fileName, kernelName,
					     kernelInfo);

  sprintf(kernelName, "lbsPhaseField%s", suffix);
  phaseFieldKernel =  platform.buildKernel(fileName, kernelName,
						kernelInfo);

  // kernels from volume file
  sprintf(fileName, DLBS "/okl/lbsVolume%s.okl", suffix);
  sprintf(kernelName, "lbsVolume%s", suffix);
  volumeKernel =  platform.buildKernel(fileName, kernelName,
					    kernelInfo);

  // kernels from surface file
  sprintf(fileName, DLBS "/okl/lbsSurface%s.okl", suffix);
  
  sprintf(kernelName, "lbsSurface%s", suffix);
  surfaceKernel = platform.buildKernel(fileName, kernelName,
					    kernelInfo);

  // vorticity calculation
  sprintf(fileName, DLBS "/okl/lbsVorticity%s.okl", suffix);
  sprintf(kernelName, "lbsVorticity%s", suffix);

  vorticityKernel = platform.buildKernel(fileName, kernelName,
					      kernelInfo);
}

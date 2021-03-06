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

#include "cns.hpp"

cns_t& cns_t::Setup(platform_t& platform, mesh_t& mesh,
                    cnsSettings_t& settings){

  cns_t* cns = new cns_t(platform, mesh, settings);

  //get physical paramters
  settings.getSetting("VISCOSITY", cns->mu);
  settings.getSetting("GAMMA", cns->gamma);

  cns->cubature   = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;
  cns->isothermal = (settings.compareSetting("ISOTHERMAL", "TRUE")) ? 1:0;

  //setup cubature
  if (cns->cubature) {
    mesh.CubatureSetup();
    mesh.CubatureNodes();
  }

  cns->Nfields   = (mesh.dim==3) ? 4:3;
  cns->Ngrads = mesh.dim*mesh.dim;

  if (!cns->isothermal) cns->Nfields++; //include energy equation

  dlong NlocalFields = mesh.Nelements*mesh.Np*cns->Nfields;
  dlong NhaloFields  = mesh.totalHaloPairs*mesh.Np*cns->Nfields;
  dlong NlocalGrads = mesh.Nelements*mesh.Np*cns->Ngrads;
  dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*cns->Ngrads;

  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    cns->timeStepper = new TimeStepper::ab3(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, cns->Nfields, *cns);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    cns->timeStepper = new TimeStepper::lserk4(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, cns->Nfields, *cns);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    cns->timeStepper = new TimeStepper::dopri5(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, cns->Nfields, *cns, mesh.comm);
  }

  //setup linear algebra module
  platform.linAlg.InitKernels({"innerProd", "max"});

  /*setup trace halo exchange */
  cns->fieldTraceHalo = mesh.HaloTraceSetup(cns->Nfields);
  cns->gradTraceHalo  = mesh.HaloTraceSetup(cns->Ngrads);

  // compute samples of q at interpolation nodes
  cns->q = (dfloat*) calloc(NlocalFields+NhaloFields, sizeof(dfloat));
  cns->o_q = platform.malloc((NlocalFields+NhaloFields)*sizeof(dfloat),
                                              cns->q);

  cns->gradq = (dfloat*) calloc(NlocalGrads+NhaloGrads, sizeof(dfloat));
  cns->o_gradq = platform.malloc((NlocalGrads+NhaloGrads)*sizeof(dfloat),
                                              cns->gradq);

  cns->Vort = (dfloat*) calloc(mesh.dim*mesh.Nelements*mesh.Np, sizeof(dfloat));
  cns->o_Vort = platform.malloc((mesh.dim*mesh.Nelements*mesh.Np)*sizeof(dfloat),
                                              cns->Vort);

  //storage for M*q during reporting
  cns->o_Mq = platform.malloc((NlocalFields+NhaloFields)*sizeof(dfloat), cns->q);
  mesh.MassMatrixKernelSetup(cns->Nfields); // mass matrix operator

  // OCCA build stuff
  occa::properties kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= cns->Nfields;
  kernelInfo["defines/" "p_Ngrads"]= cns->Ngrads;

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = mymax(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  if (cns->cubature) {
    int cubMaxNodes = mymax(mesh.Np, (mesh.intNfp*mesh.Nfaces));
    kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
    int cubMaxNodes1 = mymax(mesh.Np, (mesh.intNfp));
    kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

    int cubNblockV = mymax(1, blockMax/mesh.cubNp);
    kernelInfo["defines/" "p_cubNblockV"]= cubNblockV;

    int cubNblockS = mymax(1, blockMax/cubMaxNodes);
    kernelInfo["defines/" "p_cubNblockS"]= cubNblockS;
  }

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

  if (cns->isothermal) {
    if (cns->cubature) {
      // kernels from volume file
      sprintf(fileName, DCNS "/okl/cnsIsothermalCubatureVolume%s.okl", suffix);
      sprintf(kernelName, "cnsIsothermalCubatureVolume%s", suffix);

      cns->cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
      // kernels from surface file
      sprintf(fileName, DCNS "/okl/cnsIsothermalCubatureSurface%s.okl", suffix);
      sprintf(kernelName, "cnsIsothermalCubatureSurface%s", suffix);

      cns->cubatureSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
    } else {
      // kernels from volume file
      sprintf(fileName, DCNS "/okl/cnsIsothermalVolume%s.okl", suffix);
      sprintf(kernelName, "cnsIsothermalVolume%s", suffix);

      cns->volumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      // kernels from surface file
      sprintf(fileName, DCNS "/okl/cnsIsothermalSurface%s.okl", suffix);
      sprintf(kernelName, "cnsIsothermalSurface%s", suffix);

      cns->surfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    }
  } else {
    if (cns->cubature) {
      // kernels from volume file
      sprintf(fileName, DCNS "/okl/cnsCubatureVolume%s.okl", suffix);
      sprintf(kernelName, "cnsCubatureVolume%s", suffix);

      cns->cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
      // kernels from surface file
      sprintf(fileName, DCNS "/okl/cnsCubatureSurface%s.okl", suffix);
      sprintf(kernelName, "cnsCubatureSurface%s", suffix);

      cns->cubatureSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
    } else {
      // kernels from volume file
      sprintf(fileName, DCNS "/okl/cnsVolume%s.okl", suffix);
      sprintf(kernelName, "cnsVolume%s", suffix);

      cns->volumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      // kernels from surface file
      sprintf(fileName, DCNS "/okl/cnsSurface%s.okl", suffix);
      sprintf(kernelName, "cnsSurface%s", suffix);

      cns->surfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    }
  }

  // kernels from volume file
  sprintf(fileName, DCNS "/okl/cnsGradVolume%s.okl", suffix);
  sprintf(kernelName, "cnsGradVolume%s", suffix);

  cns->gradVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  // kernels from surface file
  sprintf(fileName, DCNS "/okl/cnsGradSurface%s.okl", suffix);
  sprintf(kernelName, "cnsGradSurface%s", suffix);

  cns->gradSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

  // vorticity calculation
  sprintf(fileName, DCNS "/okl/cnsVorticity%s.okl", suffix);
  sprintf(kernelName, "cnsVorticity%s", suffix);

  cns->vorticityKernel = platform.buildKernel(fileName, kernelName,
                                     kernelInfo);

  if (mesh.dim==2) {
    sprintf(fileName, DCNS "/okl/cnsInitialCondition2D.okl");
    if (cns->isothermal)
      sprintf(kernelName, "cnsIsothermalInitialCondition2D");
    else
      sprintf(kernelName, "cnsInitialCondition2D");
  } else {
    sprintf(fileName, DCNS "/okl/cnsInitialCondition3D.okl");
    if (cns->isothermal)
      sprintf(kernelName, "cnsIsothermalInitialCondition3D");
    else
      sprintf(kernelName, "cnsInitialCondition3D");
  }


  cns->initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                            kernelInfo);

  sprintf(fileName, DCNS "/okl/cnsMaxWaveSpeed%s.okl", suffix);
  if (cns->isothermal) {
    sprintf(kernelName, "cnsIsothermalMaxWaveSpeed%s", suffix);
  } else {
    sprintf(kernelName, "cnsMaxWaveSpeed%s", suffix);
  }

  cns->maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName,
                                            kernelInfo);

  return *cns;
}

cns_t::~cns_t() {
  volumeKernel.free();
  surfaceKernel.free();
  cubatureVolumeKernel.free();
  cubatureSurfaceKernel.free();
  gradVolumeKernel.free();
  gradSurfaceKernel.free();
  vorticityKernel.free();
  constrainKernel.free();
  initialConditionKernel.free();
  maxWaveSpeedKernel.free();

  if (timeStepper) delete timeStepper;
  if (fieldTraceHalo) fieldTraceHalo->Free();
  if (gradTraceHalo) gradTraceHalo->Free();
}

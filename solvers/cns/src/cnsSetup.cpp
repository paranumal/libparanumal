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

cns_t& cns_t::Setup(mesh_t& mesh, linAlg_t& linAlg){

  cns_t* cns = new cns_t(mesh, linAlg);

  settings_t& settings = cns->settings;

  //get physical paramters
  settings.getSetting("VISCOSITY", cns->mu);
  settings.getSetting("GAMMA", cns->gamma);

  cns->cubature   = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;
  cns->isothermal = (settings.compareSetting("ISOTHERMAL", "TRUE")) ? 1:0;

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
                                              mesh.Np, cns->Nfields, *cns);
  }

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat cfl = 0.5; // depends on the stability region size

  dfloat dtAdv  = hmin/((mesh.N+1.)*(mesh.N+1.)*cns->gamma); //just an estimate
  dfloat dtVisc = pow(hmin, 2)/(pow(mesh.N+1,4)*cns->mu);

  dfloat dt = cfl*mymin(dtAdv, dtVisc);
  cns->timeStepper->SetTimeStep(dt);

  //setup linear algebra module
  cns->linAlg.InitKernels({"innerProd"}, mesh.comm);

  /*setup trace halo exchange */
  cns->fieldTraceHalo = mesh.HaloTraceSetup(cns->Nfields);
  cns->gradTraceHalo  = mesh.HaloTraceSetup(cns->Ngrads);

  // compute samples of q at interpolation nodes
  cns->q = (dfloat*) calloc(NlocalFields+NhaloFields, sizeof(dfloat));
  cns->o_q = mesh.device.malloc((NlocalFields+NhaloFields)*sizeof(dfloat),
                                              cns->q);

  cns->gradq = (dfloat*) calloc(NlocalGrads+NhaloGrads, sizeof(dfloat));
  cns->o_gradq = mesh.device.malloc((NlocalGrads+NhaloGrads)*sizeof(dfloat),
                                              cns->gradq);

  cns->Vort = (dfloat*) calloc(mesh.dim*mesh.Nelements*mesh.Np, sizeof(dfloat));
  cns->o_Vort = mesh.device.malloc((mesh.dim*mesh.Nelements*mesh.Np)*sizeof(dfloat),
                                              cns->Vort);

  //storage for M*q during reporting
  cns->o_Mq = mesh.device.malloc((NlocalFields+NhaloFields)*sizeof(dfloat), cns->q);

  // OCCA build stuff
  occa::properties kernelInfo = cns->props; //copy base occa properties

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= cns->Nfields;
  kernelInfo["defines/" "p_Ngrads"]= cns->Ngrads;

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = 512/mesh.Np; // works for CUDA
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

  if (cns->isothermal) {
    if (cns->cubature) {
      // kernels from volume file
      sprintf(fileName, DCNS "/okl/cnsIsothermalCubatureVolume%s.okl", suffix);
      sprintf(kernelName, "cnsIsothermalCubatureVolume%s", suffix);

      cns->cubatureVolumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                               kernelInfo, mesh.comm);
      // kernels from surface file
      sprintf(fileName, DCNS "/okl/cnsIsothermalCubatureSurface%s.okl", suffix);
      sprintf(kernelName, "cnsIsothermalCubatureSurface%s", suffix);

      cns->cubatureSurfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                               kernelInfo, mesh.comm);
    } else {
      // kernels from volume file
      sprintf(fileName, DCNS "/okl/cnsIsothermalVolume%s.okl", suffix);
      sprintf(kernelName, "cnsIsothermalVolume%s", suffix);

      cns->volumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                             kernelInfo, mesh.comm);
      // kernels from surface file
      sprintf(fileName, DCNS "/okl/cnsIsothermalSurface%s.okl", suffix);
      sprintf(kernelName, "cnsIsothermalSurface%s", suffix);

      cns->surfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                             kernelInfo, mesh.comm);
    }
  } else {
    if (cns->cubature) {
      // kernels from volume file
      sprintf(fileName, DCNS "/okl/cnsCubatureVolume%s.okl", suffix);
      sprintf(kernelName, "cnsCubatureVolume%s", suffix);

      cns->cubatureVolumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                               kernelInfo, mesh.comm);
      // kernels from surface file
      sprintf(fileName, DCNS "/okl/cnsCubatureSurface%s.okl", suffix);
      sprintf(kernelName, "cnsCubatureSurface%s", suffix);

      cns->cubatureSurfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                               kernelInfo, mesh.comm);
    } else {
      // kernels from volume file
      sprintf(fileName, DCNS "/okl/cnsVolume%s.okl", suffix);
      sprintf(kernelName, "cnsVolume%s", suffix);

      cns->volumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                             kernelInfo, mesh.comm);
      // kernels from surface file
      sprintf(fileName, DCNS "/okl/cnsSurface%s.okl", suffix);
      sprintf(kernelName, "cnsSurface%s", suffix);

      cns->surfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                             kernelInfo, mesh.comm);
    }
  }

  // kernels from volume file
  sprintf(fileName, DCNS "/okl/cnsGradVolume%s.okl", suffix);
  sprintf(kernelName, "cnsGradVolume%s", suffix);

  cns->gradVolumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                           kernelInfo, mesh.comm);
  // kernels from surface file
  sprintf(fileName, DCNS "/okl/cnsGradSurface%s.okl", suffix);
  sprintf(kernelName, "cnsGradSurface%s", suffix);

  cns->gradSurfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                           kernelInfo, mesh.comm);

  // mass matrix operator
  sprintf(fileName, LIBP_DIR "/core/okl/MassMatrixOperator%s.okl", suffix);
  sprintf(kernelName, "MassMatrixOperator%s", suffix);

  cns->MassMatrixKernel = buildKernel(mesh.device, fileName, kernelName,
                                      kernelInfo, mesh.comm);

  // vorticity calculation
  sprintf(fileName, DCNS "/okl/cnsVorticity%s.okl", suffix);
  sprintf(kernelName, "cnsVorticity%s", suffix);

  cns->vorticityKernel = buildKernel(mesh.device, fileName, kernelName,
                                     kernelInfo, mesh.comm);

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

  cns->initialConditionKernel = buildKernel(mesh.device, fileName, kernelName,
                                            kernelInfo, mesh.comm);

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
  MassMatrixKernel.free();
  initialConditionKernel.free();

  if (timeStepper) delete timeStepper;
  if (fieldTraceHalo) fieldTraceHalo->Free();
  if (gradTraceHalo) gradTraceHalo->Free();
}
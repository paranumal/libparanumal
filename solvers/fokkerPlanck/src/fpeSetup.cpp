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

#include "fpe.hpp"

fpe_t& fpe_t::Setup(platform_t& platform, mesh_t& mesh,
                    fpeSettings_t& settings){

  fpe_t* fpe = new fpe_t(platform, mesh, settings);

  settings.getSetting("VISCOSITY", fpe->mu);

  fpe->cubature = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;

  //setup cubature
  if (fpe->cubature) {
    mesh.CubatureSetup();
    mesh.CubatureNodes();
  }

  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;

  //setup timeStepper
  dfloat gamma = 0.0;
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    fpe->timeStepper = new TimeStepper::ab3(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, 1, *fpe);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    fpe->timeStepper = new TimeStepper::lserk4(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, 1, *fpe);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    fpe->timeStepper = new TimeStepper::dopri5(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, 1, *fpe);
  } else if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")){
    fpe->timeStepper = new TimeStepper::extbdf3(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, 1, *fpe);
    gamma = ((TimeStepper::extbdf3*) fpe->timeStepper)->getGamma();
  } else if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")){
    fpe->timeStepper = new TimeStepper::ssbdf3(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, 1, *fpe);
    gamma = ((TimeStepper::ssbdf3*) fpe->timeStepper)->getGamma();
  }

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat cfl = 0.5; // depends on the stability region size

  dfloat dtAdvc = cfl*hmin/((mesh.N+1.)*(mesh.N+1.));
  dfloat dtDiff = fpe->mu>0.0 ? pow(hmin, 2)/(pow(mesh.N+1,4)*fpe->mu) : 1.0e9;
  dfloat dt = mymin(dtAdvc, dtDiff);

  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3"))
    fpe->timeStepper->SetTimeStep(dtAdvc);
  else if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    fpe->Nsubcycles=1;
    settings.getSetting("NUMBER OF SUBCYCLES", fpe->Nsubcycles);
    fpe->timeStepper->SetTimeStep(fpe->Nsubcycles*dtAdvc);
  } else
    fpe->timeStepper->SetTimeStep(dt);

  //Setup Elliptic solver
  fpe->elliptic=NULL;
  fpe->linearSolver=NULL;
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    ||settings.compareSetting("TIME INTEGRATOR","SSBDF3")){

    int NBCTypes = 7;
    int BCType[NBCTypes] = {0,1,1,2,1,1,1}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.

    fpe->ellipticSettings = settings.extractEllipticSettings();
    dfloat lambda = gamma/(dtAdvc*fpe->mu);
    fpe->elliptic = &(elliptic_t::Setup(platform, mesh, *(fpe->ellipticSettings),
                                             lambda, NBCTypes, BCType));
    fpe->tau = fpe->elliptic->tau;

    int weighted = settings.compareSetting("ELLIPTIC DISCRETIZATION", "CONTINUOUS") ? 1 : 0;
    fpe->linearSolver = linearSolver_t::Setup(Nlocal, Nhalo,
                                              platform, *(fpe->ellipticSettings), mesh.comm,
                                              weighted, fpe->elliptic->o_weight);
  } else {
    //set penalty
    if (mesh.elementType==TRIANGLES ||
        mesh.elementType==QUADRILATERALS){
      fpe->tau = 2.0*(mesh.N+1)*(mesh.N+2)/2.0;
      if(mesh.dim==3)
        fpe->tau *= 1.5;
    } else
      fpe->tau = 2.0*(mesh.N+1)*(mesh.N+3);
  }

  //setup linear algebra module
  platform.linAlg.InitKernels({"innerProd", "axpy"});

  /*setup trace halo exchange */
  fpe->traceHalo = mesh.HaloTraceSetup(1); //one field

  // compute samples of q at interpolation nodes
  fpe->q = (dfloat*) calloc(Nlocal+Nhalo, sizeof(dfloat));
  fpe->o_q = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), fpe->q);

  //storage for M*q during reporting
  fpe->o_Mq = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), fpe->q);
  mesh.MassMatrixKernelSetup(1); // mass matrix operator

  fpe->grad = (dfloat*) calloc((Nlocal+Nhalo)*4, sizeof(dfloat));
  fpe->o_grad  = platform.malloc((Nlocal+Nhalo)*4*sizeof(dfloat), fpe->grad);

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

  int NblockV = mymax(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  if (fpe->cubature) {
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

  // advection kernels
  if (fpe->cubature) {
    sprintf(fileName, DFPE "/okl/fpeCubatureAdvection%s.okl", suffix);
    sprintf(kernelName, "fpeAdvectionCubatureVolume%s", suffix);
    fpe->advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
    sprintf(kernelName, "fpeAdvectionCubatureSurface%s", suffix);
    fpe->advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    sprintf(fileName, DFPE "/okl/fpeAdvection%s.okl", suffix);
    sprintf(kernelName, "fpeAdvectionVolume%s", suffix);
    fpe->advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
    sprintf(kernelName, "fpeAdvectionSurface%s", suffix);
    fpe->advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }


  // diffusion kernels
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    ||settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    sprintf(fileName, DFPE "/okl/fpeDiffusionRhs%s.okl", suffix);
    sprintf(kernelName, "fpeDiffusionRhs%s", suffix);
    fpe->diffusionRhsKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    // gradient kernel
    sprintf(fileName, DFPE "/okl/fpeGradient%s.okl", suffix);
    sprintf(kernelName, "fpeGradient%s", suffix);
    fpe->gradientKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    sprintf(fileName, DFPE "/okl/fpeDiffusion%s.okl", suffix);
    sprintf(kernelName, "fpeDiffusion%s", suffix);
    fpe->diffusionKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }

  if (mesh.dim==2) {
    sprintf(fileName, DFPE "/okl/fpeInitialCondition2D.okl");
    sprintf(kernelName, "fpeInitialCondition2D");
  } else {
    sprintf(fileName, DFPE "/okl/fpeInitialCondition3D.okl");
    sprintf(kernelName, "fpeInitialCondition3D");
  }

  fpe->initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);

  //build subcycler
  fpe->subcycler=NULL;
  fpe->subStepper=NULL;
  if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    fpe->subcycler  = new subcycler_t(*fpe);
    if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","AB3")){
      fpe->subStepper = new TimeStepper::ab3(mesh.Nelements, mesh.totalHaloPairs,
                                                mesh.Np, 1, *(fpe->subcycler));
    } else if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","LSERK4")){
      fpe->subStepper = new TimeStepper::lserk4(mesh.Nelements, mesh.totalHaloPairs,
                                                mesh.Np, 1, *(fpe->subcycler));
    } else if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","DOPRI5")){
      fpe->subStepper = new TimeStepper::dopri5(mesh.Nelements, mesh.totalHaloPairs,
                                                mesh.Np, 1, *(fpe->subcycler));
    }
    fpe->subStepper->SetTimeStep(dtAdvc);
  }

  return *fpe;
}

fpe_t::~fpe_t() {
  advectionVolumeKernel.free();
  advectionSurfaceKernel.free();
  gradientKernel.free();
  diffusionKernel.free();
  diffusionRhsKernel.free();
  initialConditionKernel.free();

  if (elliptic) delete elliptic;
  if (timeStepper) delete timeStepper;
  if (linearSolver) delete linearSolver;
  if (subStepper) delete subStepper;
  if (subcycler) delete subcycler;
  if (traceHalo) traceHalo->Free();
}

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

#include "ins.hpp"

ins_t& ins_t::Setup(platform_t& platform, mesh_t& mesh,
                    insSettings_t& settings){

  ins_t* ins = new ins_t(platform, mesh, settings);

  ins->NVfields = (mesh.dim==3) ? 3:2; // Total Number of Velocity Fields
  ins->NTfields = (mesh.dim==3) ? 4:3; // Total Velocity + Pressure

  settings.getSetting("VISCOSITY", ins->nu);

  ins->cubature = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;
  ins->pressureIncrement = (settings.compareSetting("PRESSURE INCREMENT", "TRUE")) ? 1:0;

  //setup cubature
  if (ins->cubature) {
    mesh.CubatureSetup();
    mesh.CubatureNodes();
  }

  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;

  //setup timeStepper
  dfloat gamma = 0.0;
  int Nstages=1;
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")){
    ins->timeStepper = new TimeStepper::extbdf3(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, ins->NVfields, *ins);
    gamma = ((TimeStepper::extbdf3*) ins->timeStepper)->getGamma();
    Nstages = ((TimeStepper::extbdf3*) ins->timeStepper)->Nstages;
  } else if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")){
    ins->timeStepper = new TimeStepper::ssbdf3(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, ins->NVfields, *ins);
    gamma = ((TimeStepper::ssbdf3*) ins->timeStepper)->getGamma();
    Nstages = ((TimeStepper::ssbdf3*) ins->timeStepper)->Nstages;
  }

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat cfl = 0.25; // depends on the stability region size

  dfloat dtAdvc = cfl*hmin/((mesh.N+1.)*(mesh.N+1.));
  dfloat dtDiff = ins->nu>0.0 ? pow(hmin, 2)/(pow(mesh.N+1,4)*ins->nu) : 1.0e9;
  dfloat dt = mymin(dtAdvc, dtDiff);

  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3"))
    ins->timeStepper->SetTimeStep(dtAdvc);
  else if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    ins->Nsubcycles=1;
    settings.getSetting("NUMBER OF SUBCYCLES", ins->Nsubcycles);
    ins->timeStepper->SetTimeStep(ins->Nsubcycles*dtAdvc);
  } else
    ins->timeStepper->SetTimeStep(dt);

  //Setup velocity Elliptic solvers
  ins->uSolver=NULL;
  ins->vSolver=NULL;
  ins->wSolver=NULL;
  ins->uLinearSolver=NULL;
  ins->vLinearSolver=NULL;
  ins->wLinearSolver=NULL;
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    ||settings.compareSetting("TIME INTEGRATOR","SSBDF3")){

    // SetUp Boundary Flags types for Elliptic Solve
    // bc = 1 -> wall
    // bc = 2 -> inflow
    // bc = 3 -> outflow
    // bc = 4 -> x-aligned slip
    // bc = 5 -> y-aligned slip
    // bc = 6 -> z-aligned slip
    int NBCTypes = 7;
    int uBCType[NBCTypes] = {0,1,1,2,1,2,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
    int vBCType[NBCTypes] = {0,1,1,2,2,1,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
    int wBCType[NBCTypes] = {0,1,1,2,2,2,1}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.

    ins->vSettings = settings.extractVelocitySettings();
    dfloat lambda = gamma/(dtAdvc*ins->nu);
    ins->uSolver = &(elliptic_t::Setup(platform, mesh, *(ins->vSettings),
                                             lambda, NBCTypes, uBCType));
    ins->vSolver = &(elliptic_t::Setup(platform, mesh, *(ins->vSettings),
                                             lambda, NBCTypes, vBCType));
    ins->wSolver = &(elliptic_t::Setup(platform, mesh, *(ins->vSettings),
                                             lambda, NBCTypes, wBCType));
    ins->vTau = ins->uSolver->tau;

    ins->vDisc_c0 = settings.compareSetting("VELOCITY DISCRETIZATION", "CONTINUOUS") ? 1 : 0;

    //ins->vLinearSolver = linearSolver_t::Setup(Nlocal, Nhalo,
    //                                           platform, *(ins->vSettings), mesh.comm,
    //                                           ins->vDisc_c0, ins->vSolver->o_weight);

    ins->uLinearSolver = initialGuessSolver_t::Setup(Nlocal, Nhalo,
                                                    platform, *(ins->vSettings), mesh.comm,
                                                    ins->vDisc_c0, ins->uSolver->o_weight);

    ins->vLinearSolver = initialGuessSolver_t::Setup(Nlocal, Nhalo,
                                                    platform, *(ins->vSettings), mesh.comm,
                                                    ins->vDisc_c0, ins->vSolver->o_weight);
    if (mesh.dim == 3) {
      ins->wLinearSolver = initialGuessSolver_t::Setup(Nlocal, Nhalo,
                                                       platform, *(ins->vSettings), mesh.comm,
                                                       ins->vDisc_c0, ins->wSolver->o_weight);
    }

  } else {
    ins->vDisc_c0 = 0;

    //set penalty
    if (mesh.elementType==TRIANGLES ||
        mesh.elementType==QUADRILATERALS){
      ins->vTau = 2.0*(mesh.N+1)*(mesh.N+2)/2.0;
      if(mesh.dim==3)
        ins->vTau *= 1.5;
    } else
      ins->vTau = 2.0*(mesh.N+1)*(mesh.N+3);
  }

  //Setup pressure Elliptic solver
  {
    int NBCTypes = 7;
    int pBCType[NBCTypes] = {0,2,2,1,2,2,2}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.

    ins->pSettings = settings.extractPressureSettings();
    ins->pSolver = &(elliptic_t::Setup(platform, mesh, *(ins->pSettings),
                                             0.0, NBCTypes, pBCType));
    ins->pTau = ins->pSolver->tau;

    ins->pDisc_c0 = settings.compareSetting("PRESSURE DISCRETIZATION", "CONTINUOUS") ? 1 : 0;

    //ins->pLinearSolver = linearSolver_t::Setup(Nlocal, Nhalo,
    //                                           platform, *(ins->pSettings), mesh.comm,
    //                                           ins->pDisc_c0, ins->pSolver->o_weight);

    ins->pLinearSolver = initialGuessSolver_t::Setup(Nlocal, Nhalo,
                                                     platform, *(ins->pSettings), mesh.comm,
                                                     ins->pDisc_c0, ins->pSolver->o_weight);
  }

  //Solver tolerances
  ins->presTOL = 1E-8;
  ins->velTOL  = 1E-8;

  //build node-wise boundary flag
  ins->BoundarySetup();

  //setup linear algebra module
  platform.linAlg.InitKernels({"innerProd", "axpy"});

  /*setup trace halo exchange */
  ins->pTraceHalo = mesh.HaloTraceSetup(1); //one field
  ins->vTraceHalo = mesh.HaloTraceSetup(ins->NVfields); //one field

  // u and p at interpolation nodes
  ins->u = (dfloat*) calloc((Nlocal+Nhalo)*ins->NVfields, sizeof(dfloat));
  ins->o_u = platform.malloc((Nlocal+Nhalo)*ins->NVfields*sizeof(dfloat), ins->u);

  ins->p = (dfloat*) calloc(Nlocal+Nhalo, sizeof(dfloat));
  ins->o_p = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), ins->p);

  //storage for velocity gradient
  if ( !settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    && !settings.compareSetting("TIME INTEGRATOR","SSBDF3"))
    ins->o_GU = platform.malloc((Nlocal+Nhalo)*4*sizeof(dfloat));

  //extra buffers for solvers
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    ||settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    ins->o_UH = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat));
    ins->o_VH = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat));
    if (mesh.dim==3)
      ins->o_WH = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat));
    else
      ins->o_WH = platform.malloc((1)*sizeof(dfloat));

    ins->o_rhsU = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat));
    ins->o_rhsV = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat));
    if (mesh.dim==3)
      ins->o_rhsW = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat));
    else
      ins->o_rhsW = platform.malloc((1)*sizeof(dfloat));
  }

  if (ins->pressureIncrement)
    ins->o_PI = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), ins->p);

  ins->o_rhsP = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat));

  //storage for M*u during reporting
  ins->o_MU = platform.malloc((Nlocal+Nhalo)*ins->NVfields*sizeof(dfloat), ins->u);
  mesh.MassMatrixKernelSetup(ins->NVfields); // mass matrix operator

  if (mesh.dim==2) {
    ins->Vort = (dfloat*) calloc((Nlocal+Nhalo), sizeof(dfloat));
    ins->o_Vort = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), ins->Vort);
  } else {
    ins->Vort = (dfloat*) calloc((Nlocal+Nhalo)*ins->NVfields, sizeof(dfloat));
    ins->o_Vort = platform.malloc((Nlocal+Nhalo)*ins->NVfields*sizeof(dfloat), ins->Vort);
  }

  // OCCA build stuff
  occa::properties kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"] = ins->NVfields;
  kernelInfo["defines/" "p_NVfields"]= ins->NVfields;
  kernelInfo["defines/" "p_NTfields"]= ins->NTfields;

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;

  int NblockV = mymax(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1,blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  if (ins->cubature) {
    int cubMaxNodes = mymax(mesh.Np, (mesh.intNfp*mesh.Nfaces));
    kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
    int cubMaxNodes1 = mymax(mesh.Np, (mesh.intNfp));
    kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

    int cubNblockV = mymax(1,blockMax/mesh.cubNp);
    kernelInfo["defines/" "p_cubNblockV"]= cubNblockV;

    int cubNblockS = mymax(1,blockMax/cubMaxNodes);
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
  ins->subcycler=NULL;
  ins->subStepper=NULL;
  if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    //subcycle kernels
    if (ins->cubature) {
      sprintf(fileName, DINS "/okl/insSubcycleCubatureAdvection%s.okl", suffix);
      sprintf(kernelName, "insSubcycleAdvectionCubatureVolume%s", suffix);
      ins->advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      sprintf(kernelName, "insSubcycleAdvectionCubatureSurface%s", suffix);
      ins->advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      sprintf(kernelName, "insSubcycleAdvectionCubatureInterpolation%s", suffix);
      ins->advectionInterpolationKernel = platform.buildKernel(fileName, kernelName,
							       kernelInfo);

    } else {
      sprintf(fileName, DINS "/okl/insSubcycleAdvection%s.okl", suffix);
      sprintf(kernelName, "insSubcycleAdvectionVolume%s", suffix);
      ins->advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      sprintf(kernelName, "insSubcycleAdvectionSurface%s", suffix);
      ins->advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    }

    //build subcycler
    ins->subcycler  = new subcycler_t(*ins);
    if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","AB3")){
      ins->subStepper = new TimeStepper::ab3(mesh.Nelements, mesh.totalHaloPairs,
                                                mesh.Np, ins->NVfields, *(ins->subcycler));
    } else if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","LSERK4")){
      ins->subStepper = new TimeStepper::lserk4(mesh.Nelements, mesh.totalHaloPairs,
                                                mesh.Np, ins->NVfields, *(ins->subcycler));
    } else if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","DOPRI5")){
      ins->subStepper = new TimeStepper::dopri5(mesh.Nelements, mesh.totalHaloPairs,
                                                mesh.Np, ins->NVfields, *(ins->subcycler), mesh.comm);
    }
    ins->subStepper->SetTimeStep(dtAdvc);

    sprintf(fileName, DINS "/okl/insSubcycleAdvection.okl");
    sprintf(kernelName, "insSubcycleAdvectionKernel");
    ins->subcycler->subCycleAdvectionKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);

    if (ins->cubature) {
      dlong cubNlocal = mesh.Nelements*mesh.cubNp;
      dlong cubNhalo  = mesh.totalHaloPairs*mesh.cubNp;

      //history space of U interpolated to cubature nodes
      ins->subcycler->o_cUh = platform.malloc(Nstages*cubNlocal*NVfields*sizeof(dfloat));

      ins->subcycler->o_Ue = platform.malloc((cubNlocal+cubNhalo)*ins->NVfields*sizeof(dfloat));
    } else {
      ins->subcycler->o_Ue = platform.malloc((Nlocal+Nhalo)*ins->NVfields*sizeof(dfloat), ins->u);
    }

  } else {
    //regular advection kernels
    if (ins->cubature) {
      sprintf(fileName, DINS "/okl/insCubatureAdvection%s.okl", suffix);
      sprintf(kernelName, "insAdvectionCubatureVolume%s", suffix);
      ins->advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      sprintf(kernelName, "insAdvectionCubatureSurface%s", suffix);
      ins->advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    } else {
      sprintf(fileName, DINS "/okl/insAdvection%s.okl", suffix);
      sprintf(kernelName, "insAdvectionVolume%s", suffix);
      ins->advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      sprintf(kernelName, "insAdvectionSurface%s", suffix);
      ins->advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    }
  }

  // diffusion kernels
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    ||settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    sprintf(fileName, DINS "/okl/insVelocityRhs%s.okl", suffix);

    if (ins->vDisc_c0)
      sprintf(kernelName, "insVelocityRhs%s", suffix);
    else
      sprintf(kernelName, "insVelocityIpdgRhs%s", suffix);
    ins->velocityRhsKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    sprintf(kernelName, "insVelocityBC%s", suffix);
    ins->velocityBCKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    // gradient kernel
    sprintf(fileName, DINS "/okl/insVelocityGradient%s.okl", suffix);
    sprintf(kernelName, "insVelocityGradient%s", suffix);
    ins->velocityGradientKernel =  platform.buildKernel(fileName, kernelName,
                                               kernelInfo);

    sprintf(fileName, DINS "/okl/insDiffusion%s.okl", suffix);
    sprintf(kernelName, "insDiffusion%s", suffix);
    ins->diffusionKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }

  //pressure gradient kernels
  sprintf(fileName, DINS "/okl/insGradient%s.okl", suffix);
  sprintf(kernelName, "insGradientVolume%s", suffix);
  ins->gradientVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  sprintf(kernelName, "insGradientSurface%s", suffix);
  ins->gradientSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  //velocity divergence kernels
  sprintf(fileName, DINS "/okl/insDivergence%s.okl", suffix);
  sprintf(kernelName, "insDivergenceVolume%s", suffix);
  ins->divergenceVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  sprintf(kernelName, "insDivergenceSurface%s", suffix);
  ins->divergenceSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  //pressure solver kernels
  if (ins->pressureIncrement) {
    sprintf(fileName, DINS "/okl/insPressureIncrementRhs%s.okl", suffix);

    if (ins->pDisc_c0)
      sprintf(kernelName, "insPressureIncrementRhs%s", suffix);
    else
      sprintf(kernelName, "insPressureIncrementIpdgRhs%s", suffix);
    ins->pressureIncrementRhsKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    sprintf(kernelName, "insPressureIncrementBC%s", suffix);
    ins->pressureIncrementBCKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    sprintf(fileName, DINS "/okl/insPressureRhs%s.okl", suffix);
    if (ins->pDisc_c0)
      sprintf(kernelName, "insPressureRhs%s", suffix);
    else
      sprintf(kernelName, "insPressureIpdgRhs%s", suffix);
    ins->pressureRhsKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    sprintf(kernelName, "insPressureBC%s", suffix);
    ins->pressureBCKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }

  sprintf(fileName, DINS "/okl/insVorticity%s.okl", suffix);
  sprintf(kernelName, "insVorticity%s", suffix);
  ins->vorticityKernel =  platform.buildKernel(fileName, kernelName,
                                            kernelInfo);

  if (mesh.dim==2) {
    sprintf(fileName, DINS "/okl/insInitialCondition2D.okl");
    sprintf(kernelName, "insInitialCondition2D");
  } else {
    sprintf(fileName, DINS "/okl/insInitialCondition3D.okl");
    sprintf(kernelName, "insInitialCondition3D");
  }

  ins->initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);



  return *ins;
}

ins_t::~ins_t() {
  advectionVolumeKernel.free();
  advectionSurfaceKernel.free();
  divergenceVolumeKernel.free();
  divergenceSurfaceKernel.free();
  gradientVolumeKernel.free();
  gradientSurfaceKernel.free();
  velocityGradientKernel.free();
  diffusionKernel.free();
  velocityRhsKernel.free();
  velocityBCKernel.free();
  pressureRhsKernel.free();
  pressureBCKernel.free();
  vorticityKernel.free();
  initialConditionKernel.free();

  if (pSolver) delete pSolver;
  if (uSolver) delete uSolver;
  if (vSolver) delete vSolver;
  if (wSolver) delete wSolver;
  if (timeStepper) delete timeStepper;
  if (pLinearSolver) delete pLinearSolver;
  if (uLinearSolver) delete uLinearSolver;
  if (vLinearSolver) delete vLinearSolver;
  if (wLinearSolver) delete wLinearSolver;
  if (subStepper) delete subStepper;
  if (subcycler) {
    subcycler->subCycleAdvectionKernel.free();
    delete subcycler;
  }

  if (vTraceHalo) vTraceHalo->Free();
  if (pTraceHalo) pTraceHalo->Free();
}

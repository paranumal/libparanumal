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

#include "ins.hpp"

void ins_t::Setup(platform_t& _platform, mesh_t& _mesh,
                  insSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  NVfields = (mesh.dim==3) ? 3:2; // Total Number of Velocity Fields
  NTfields = (mesh.dim==3) ? 4:3; // Total Velocity + Pressure

  settings.getSetting("VISCOSITY", nu);

  cubature = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;
  pressureIncrement = (settings.compareSetting("PRESSURE INCREMENT", "TRUE")) ? 1:0;

  //setup cubature
  if (cubature) {
    mesh.CubatureSetup();
    mesh.CubaturePhysicalNodes();
  }

  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;

  //setup timeStepper
  dfloat gamma = 0.0;
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")){
    timeStepper.Setup<TimeStepper::extbdf3>(mesh.Nelements,
                                            mesh.totalHaloPairs,
                                            mesh.Np, NVfields, platform, comm);
    gamma = timeStepper.GetGamma();
  } else if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")){
    timeStepper.Setup<TimeStepper::ssbdf3>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, NVfields, platform, comm);
    gamma = timeStepper.GetGamma();
  }

  Nsubcycles=1;
  if (settings.compareSetting("TIME INTEGRATOR","SSBDF3"))
    settings.getSetting("NUMBER OF SUBCYCLES", Nsubcycles);

  //Setup velocity Elliptic solvers
  dlong uNlocal=0, vNlocal=0, wNlocal=0;
  dlong uNhalo=0, vNhalo=0, wNhalo=0;
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
    memory<int> uBCType(NBCTypes);
    // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
    uBCType[0] = 0;
    uBCType[1] = 1;
    uBCType[2] = 1;
    uBCType[3] = 2;
    uBCType[4] = 1;
    uBCType[5] = 2;
    uBCType[6] = 2;

    memory<int> vBCType(NBCTypes);
    // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
    vBCType[0] = 0;
    vBCType[1] = 1;
    vBCType[2] = 1;
    vBCType[3] = 2;
    vBCType[4] = 2;
    vBCType[5] = 1;
    vBCType[6] = 2;

    memory<int> wBCType(NBCTypes);
    // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
    wBCType[0] = 0;
    wBCType[1] = 1;
    wBCType[2] = 1;
    wBCType[3] = 2;
    wBCType[4] = 2;
    wBCType[5] = 2;
    wBCType[6] = 1;

    vSettings = _settings.extractVelocitySettings();

    //make a guess at dt for the lambda value
    //TODO: we should allow preconditioners to be re-setup if lambda is updated
    dfloat hmin = mesh.MinCharacteristicLength();
    dfloat dtAdvc = Nsubcycles*hmin/((mesh.N+1.)*(mesh.N+1.));
    dfloat lambda = gamma/(dtAdvc*nu);
    uSolver.Setup(platform, mesh, vSettings,
                  lambda, NBCTypes, uBCType);
    vSolver.Setup(platform, mesh, vSettings,
                  lambda, NBCTypes, vBCType);
    if (mesh.dim == 3)
      wSolver.Setup(platform, mesh, vSettings,
                    lambda, NBCTypes, wBCType);

    vTau = uSolver.tau;

    vDisc_c0 = settings.compareSetting("VELOCITY DISCRETIZATION", "CONTINUOUS") ? 1 : 0;

    uNlocal = uSolver.Ndofs;
    vNlocal = vSolver.Ndofs;
    if (mesh.dim == 3) wNlocal = wSolver.Ndofs;

    uNhalo = uSolver.Nhalo;
    vNhalo = vSolver.Nhalo;
    if (mesh.dim == 3) wNhalo = wSolver.Nhalo;

    if (vSettings.compareSetting("LINEAR SOLVER","NBPCG")){

      uLinearSolver.Setup<LinearSolver::nbpcg>(uNlocal, uNhalo, platform, vSettings, comm);
      vLinearSolver.Setup<LinearSolver::nbpcg>(vNlocal, vNhalo, platform, vSettings, comm);
      if (mesh.dim==3)
        wLinearSolver.Setup<LinearSolver::nbpcg>(wNlocal, wNhalo, platform, vSettings, comm);

    } else if (vSettings.compareSetting("LINEAR SOLVER","NBFPCG")){

      uLinearSolver.Setup<LinearSolver::nbfpcg>(uNlocal, uNhalo, platform, vSettings, comm);
      vLinearSolver.Setup<LinearSolver::nbfpcg>(vNlocal, vNhalo, platform, vSettings, comm);
      if (mesh.dim==3)
        wLinearSolver.Setup<LinearSolver::nbfpcg>(wNlocal, wNhalo, platform, vSettings, comm);

    } else if (vSettings.compareSetting("LINEAR SOLVER","PCG")){

      uLinearSolver.Setup<LinearSolver::pcg>(uNlocal, uNhalo, platform, vSettings, comm);
      vLinearSolver.Setup<LinearSolver::pcg>(vNlocal, vNhalo, platform, vSettings, comm);
      if (mesh.dim==3)
        wLinearSolver.Setup<LinearSolver::pcg>(wNlocal, wNhalo, platform, vSettings, comm);

    } else if (vSettings.compareSetting("LINEAR SOLVER","PGMRES")){

      uLinearSolver.Setup<LinearSolver::pgmres>(uNlocal, uNhalo, platform, vSettings, comm);
      vLinearSolver.Setup<LinearSolver::pgmres>(vNlocal, vNhalo, platform, vSettings, comm);
      if (mesh.dim==3)
        wLinearSolver.Setup<LinearSolver::pgmres>(wNlocal, wNhalo, platform, vSettings, comm);

    } else if (vSettings.compareSetting("LINEAR SOLVER","PMINRES")){

      uLinearSolver.Setup<LinearSolver::pminres>(uNlocal, uNhalo, platform, vSettings, comm);
      vLinearSolver.Setup<LinearSolver::pminres>(vNlocal, vNhalo, platform, vSettings, comm);
      if (mesh.dim==3)
        wLinearSolver.Setup<LinearSolver::pminres>(wNlocal, wNhalo, platform, vSettings, comm);
    }

    if (vSettings.compareSetting("INITIAL GUESS STRATEGY", "NONE")) {

      uLinearSolver.SetupInitialGuess<InitialGuess::Default>(uNlocal, platform, vSettings, comm);
      vLinearSolver.SetupInitialGuess<InitialGuess::Default>(vNlocal, platform, vSettings, comm);
      if (mesh.dim==3)
        wLinearSolver.SetupInitialGuess<InitialGuess::Default>(wNlocal, platform, vSettings, comm);

    } else if (vSettings.compareSetting("INITIAL GUESS STRATEGY", "ZERO")) {

      uLinearSolver.SetupInitialGuess<InitialGuess::Zero>(uNlocal, platform, vSettings, comm);
      vLinearSolver.SetupInitialGuess<InitialGuess::Zero>(vNlocal, platform, vSettings, comm);
      if (mesh.dim==3)
        wLinearSolver.SetupInitialGuess<InitialGuess::Zero>(wNlocal, platform, vSettings, comm);

    } else if (vSettings.compareSetting("INITIAL GUESS STRATEGY", "CLASSIC")) {

      uLinearSolver.SetupInitialGuess<InitialGuess::ClassicProjection>(uNlocal, platform, vSettings, comm);
      vLinearSolver.SetupInitialGuess<InitialGuess::ClassicProjection>(vNlocal, platform, vSettings, comm);
      if (mesh.dim==3)
        wLinearSolver.SetupInitialGuess<InitialGuess::ClassicProjection>(wNlocal, platform, vSettings, comm);

    } else if (vSettings.compareSetting("INITIAL GUESS STRATEGY", "QR")) {

      uLinearSolver.SetupInitialGuess<InitialGuess::RollingQRProjection>(uNlocal, platform, vSettings, comm);
      vLinearSolver.SetupInitialGuess<InitialGuess::RollingQRProjection>(vNlocal, platform, vSettings, comm);
      if (mesh.dim==3)
        wLinearSolver.SetupInitialGuess<InitialGuess::RollingQRProjection>(wNlocal, platform, vSettings, comm);

    } else if (vSettings.compareSetting("INITIAL GUESS STRATEGY", "EXTRAP")) {

      uLinearSolver.SetupInitialGuess<InitialGuess::Extrap>(uNlocal, platform, vSettings, comm);
      vLinearSolver.SetupInitialGuess<InitialGuess::Extrap>(vNlocal, platform, vSettings, comm);
      if (mesh.dim==3)
        wLinearSolver.SetupInitialGuess<InitialGuess::Extrap>(wNlocal, platform, vSettings, comm);

    }

  } else {
    vDisc_c0 = 0;

    //set penalty
    if (mesh.elementType==Mesh::TRIANGLES ||
        mesh.elementType==Mesh::QUADRILATERALS){
      vTau = 2.0*(mesh.N+1)*(mesh.N+2)/2.0;
      if(mesh.dim==3)
        vTau *= 1.5;
    } else
      vTau = 2.0*(mesh.N+1)*(mesh.N+3);
  }

  //Setup pressure Elliptic solver
  dlong pNlocal=0, pNhalo=0;
  {
    int NBCTypes = 7;
    memory<int> pBCType(NBCTypes);
    // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.
    pBCType[0] = 0;
    pBCType[1] = 2;
    pBCType[2] = 2;
    pBCType[3] = 1;
    pBCType[4] = 2;
    pBCType[5] = 2;
    pBCType[6] = 2;

    pSettings = _settings.extractPressureSettings();
    pSolver.Setup(platform, mesh, pSettings,
                  0.0, NBCTypes, pBCType);
    pTau = pSolver.tau;

    pDisc_c0 = settings.compareSetting("PRESSURE DISCRETIZATION", "CONTINUOUS") ? 1 : 0;

    if (pDisc_c0) {
      pNlocal = pSolver.ogsMasked.Ngather;
      pNhalo  = pSolver.gHalo.Nhalo;
    } else {
      pNlocal = mesh.Nelements*mesh.Np;
      pNhalo  = mesh.totalHaloPairs*mesh.Np;
    }

    if (vSettings.compareSetting("LINEAR SOLVER","NBPCG")){
      pLinearSolver.Setup<LinearSolver::nbpcg>(pNlocal, pNhalo, platform, pSettings, comm);
    } else if (pSettings.compareSetting("LINEAR SOLVER","NBFPCG")){
      pLinearSolver.Setup<LinearSolver::nbfpcg>(pNlocal, pNhalo, platform, pSettings, comm);
    } else if (pSettings.compareSetting("LINEAR SOLVER","PCG")){
      pLinearSolver.Setup<LinearSolver::pcg>(pNlocal, pNhalo, platform, pSettings, comm);
    } else if (pSettings.compareSetting("LINEAR SOLVER","PGMRES")){
      pLinearSolver.Setup<LinearSolver::pgmres>(pNlocal, pNhalo, platform, pSettings, comm);
    } else if (pSettings.compareSetting("LINEAR SOLVER","PMINRES")){
      pLinearSolver.Setup<LinearSolver::pminres>(pNlocal, pNhalo, platform, pSettings, comm);
    }

    if (pSettings.compareSetting("INITIAL GUESS STRATEGY", "NONE")) {
      pLinearSolver.SetupInitialGuess<InitialGuess::Default>(pNlocal, platform, pSettings, comm);
    } else if (pSettings.compareSetting("INITIAL GUESS STRATEGY", "ZERO")) {
      pLinearSolver.SetupInitialGuess<InitialGuess::Zero>(pNlocal, platform, pSettings, comm);
    } else if (pSettings.compareSetting("INITIAL GUESS STRATEGY", "CLASSIC")) {
      pLinearSolver.SetupInitialGuess<InitialGuess::ClassicProjection>(pNlocal, platform, pSettings, comm);
    } else if (pSettings.compareSetting("INITIAL GUESS STRATEGY", "QR")) {
      pLinearSolver.SetupInitialGuess<InitialGuess::RollingQRProjection>(pNlocal, platform, pSettings, comm);
    } else if (pSettings.compareSetting("INITIAL GUESS STRATEGY", "EXTRAP")) {
      pLinearSolver.SetupInitialGuess<InitialGuess::Extrap>(pNlocal, platform, pSettings, comm);
    }
  }

  //Solver tolerances
  if (sizeof(dfloat)==sizeof(double)) {
    presTOL = 1.0E-8;
    velTOL  = 1.0E-8;
  } else {
    presTOL = 1.0E-5;
    velTOL  = 1.0E-5;
  }

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd", "axpy", "max"});

  /*setup trace halo exchange */
  pTraceHalo = mesh.HaloTraceSetup(1); //one field
  vTraceHalo = mesh.HaloTraceSetup(NVfields); //one field

  // u and p at interpolation nodes
  u.malloc((Nlocal+Nhalo)*NVfields, 0.0);
  o_u = platform.malloc<dfloat>(u);

  p.malloc(Nlocal+Nhalo, 0.0);
  o_p = platform.malloc<dfloat>(p);

  //storage for velocity gradient
  if ( !settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    && !settings.compareSetting("TIME INTEGRATOR","SSBDF3"))
    o_GU = platform.malloc<dfloat>((Nlocal+Nhalo)*4);

  //extra buffers for solvers
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    ||settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    o_UH = platform.malloc<dfloat>(Nlocal+Nhalo, u);
    o_VH = platform.malloc<dfloat>(Nlocal+Nhalo, u);
    if (mesh.dim==3)
      o_WH = platform.malloc<dfloat>(Nlocal+Nhalo, u);

    o_rhsU = platform.malloc<dfloat>(Nlocal+Nhalo, u);
    o_rhsV = platform.malloc<dfloat>(Nlocal+Nhalo, u);
    if (mesh.dim==3)
      o_rhsW = platform.malloc<dfloat>(Nlocal+Nhalo, u);

    if (vDisc_c0) {
      o_GUH = platform.malloc<dfloat>(uNlocal+uNhalo, u);
      o_GVH = platform.malloc<dfloat>(vNlocal+vNhalo, u);
      if (mesh.dim==3)
        o_GWH = platform.malloc<dfloat>(wNlocal+wNhalo, u);

      o_GrhsU = platform.malloc<dfloat>(uNlocal+uNhalo, u);
      o_GrhsV = platform.malloc<dfloat>(vNlocal+vNhalo, u);
      if (mesh.dim==3)
        o_GrhsW = platform.malloc<dfloat>(wNlocal+wNhalo, u);
    }
  }

  if (pressureIncrement) {
    o_PI  = platform.malloc<dfloat>(p);
    o_GPI = platform.malloc<dfloat>(p);
  }

  o_rhsP = platform.malloc<dfloat>(p);
  if (pDisc_c0) {
    o_GP    = platform.malloc<dfloat>(p);
    o_GrhsP = platform.malloc<dfloat>(p);
  }

  //storage for M*u during reporting
  o_MU = platform.malloc<dfloat>(u);
  mesh.MassMatrixKernelSetup(NVfields); // mass matrix operator

  if (mesh.dim==2) {
    Vort.malloc(Nlocal+Nhalo, 0.0);
    o_Vort = platform.malloc<dfloat>(Vort);
  } else {
    Vort.malloc((Nlocal+Nhalo)*NVfields, 0.0);
    o_Vort = platform.malloc<dfloat>(Vort);
  }

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"] = NVfields;
  kernelInfo["defines/" "p_NVfields"]= NVfields;
  kernelInfo["defines/" "p_NTfields"]= NTfields;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;

  int NblockV = std::max(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1,blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  if (cubature) {
    int cubMaxNodes = std::max(mesh.Np, (mesh.intNfp*mesh.Nfaces));
    kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
    int cubMaxNodes1 = std::max(mesh.Np, (mesh.intNfp));
    kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

    int cubNblockV = std::max(1,blockMax/mesh.cubNp);
    kernelInfo["defines/" "p_cubNblockV"]= cubNblockV;

    int cubNblockS = std::max(1,blockMax/cubMaxNodes);
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

  std::string oklFilePrefix = DINS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // advection kernels
  if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    //subcycle kernels
    if (cubature) {
      fileName   = oklFilePrefix + "insSubcycleCubatureAdvection" + suffix + oklFileSuffix;
      kernelName = "insSubcycleAdvectionCubatureVolume" + suffix;
      advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      kernelName = "insSubcycleAdvectionCubatureSurface" + suffix;
      advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    } else {
      fileName   = oklFilePrefix + "insSubcycleAdvection" + suffix + oklFileSuffix;
      kernelName = "insSubcycleAdvectionVolume" + suffix;
      advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      kernelName = "insSubcycleAdvectionSurface" + suffix;
      advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    }

    //build subcycler
    subcycler.platform = platform;
    subcycler.mesh = mesh;
    subcycler.comm = comm;
    subcycler.settings = settings;

    subcycler.NVfields = NVfields;
    subcycler.nu = nu;
    subcycler.cubature = cubature;
    subcycler.vTraceHalo = vTraceHalo;
    subcycler.advectionVolumeKernel = advectionVolumeKernel;
    subcycler.advectionSurfaceKernel = advectionSurfaceKernel;

    if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","AB3")){
      subStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
                                         mesh.totalHaloPairs,
                                         mesh.Np, NVfields, platform, comm);
    } else if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","LSERK4")){
      subStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                            mesh.totalHaloPairs,
                                            mesh.Np, NVfields, platform, comm);
    } else if (settings.compareSetting("SUBCYCLING TIME INTEGRATOR","DOPRI5")){
      subStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
                                            mesh.totalHaloPairs,
                                            mesh.Np, NVfields, platform, comm);
    }

    fileName   = oklFilePrefix + "insSubcycleAdvection" + oklFileSuffix;
    kernelName = "insSubcycleAdvectionKernel";
    subcycler.subCycleAdvectionKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);

    subcycler.o_Ue = platform.malloc<dfloat>(u);

  } else {
    //regular advection kernels
    if (cubature) {
      fileName   = oklFilePrefix + "insCubatureAdvection" + suffix + oklFileSuffix;
      kernelName = "insAdvectionCubatureVolume" + suffix;
      advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      kernelName = "insAdvectionCubatureSurface" + suffix;
      advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    } else {
      fileName   = oklFilePrefix + "insAdvection" + suffix + oklFileSuffix;
      kernelName = "insAdvectionVolume" + suffix;
      advectionVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      kernelName = "insAdvectionSurface" + suffix;
      advectionSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    }
  }

  // diffusion kernels
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")
    ||settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    fileName   = oklFilePrefix + "insVelocityRhs" + suffix + oklFileSuffix;

    if (vDisc_c0)
      kernelName = "insVelocityRhs" + suffix;
    else
      kernelName = "insVelocityIpdgRhs" + suffix;
    velocityRhsKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    kernelName = "insVelocityBC" + suffix;
    velocityBCKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    // gradient kernel
    fileName   = oklFilePrefix + "insVelocityGradient" + suffix + oklFileSuffix;
    kernelName = "insVelocityGradient" + suffix;
    velocityGradientKernel =  platform.buildKernel(fileName, kernelName,
                                               kernelInfo);

    fileName   = oklFilePrefix + "insDiffusion" + suffix + oklFileSuffix;
    kernelName = "insDiffusion" + suffix;
    diffusionKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }

  //pressure gradient kernels
  fileName   = oklFilePrefix + "insGradient" + suffix + oklFileSuffix;
  kernelName = "insGradientVolume" + suffix;
  gradientVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  kernelName = "insGradientSurface" + suffix;
  gradientSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  //velocity divergence kernels
  fileName   = oklFilePrefix + "insDivergence" + suffix + oklFileSuffix;
  kernelName = "insDivergenceVolume" + suffix;
  divergenceVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  kernelName = "insDivergenceSurface" + suffix;
  divergenceSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  //pressure solver kernels
  if (pressureIncrement) {
    fileName   = oklFilePrefix + "insPressureIncrementRhs" + suffix + oklFileSuffix;

    if (pDisc_c0)
      kernelName = "insPressureIncrementRhs" + suffix;
    else
      kernelName = "insPressureIncrementIpdgRhs" + suffix;
    pressureIncrementRhsKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    kernelName = "insPressureIncrementBC" + suffix;
    pressureIncrementBCKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  } else {
    fileName   = oklFilePrefix + "insPressureRhs" + suffix + oklFileSuffix;
    if (pDisc_c0)
      kernelName = "insPressureRhs" + suffix;
    else
      kernelName = "insPressureIpdgRhs" + suffix;
    pressureRhsKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    kernelName = "insPressureBC" + suffix;
    pressureBCKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  }

  fileName   = oklFilePrefix + "insVorticity" + suffix + oklFileSuffix;
  kernelName = "insVorticity" + suffix;
  vorticityKernel =  platform.buildKernel(fileName, kernelName,
                                            kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "insInitialCondition2D" + oklFileSuffix;
    kernelName = "insInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "insInitialCondition3D" + oklFileSuffix;
    kernelName = "insInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);

  fileName   = oklFilePrefix + "insMaxWaveSpeed" + suffix + oklFileSuffix;
  kernelName = "insMaxWaveSpeed" + suffix;

  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
}

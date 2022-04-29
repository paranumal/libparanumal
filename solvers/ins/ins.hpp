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

#ifndef INS_HPP
#define INS_HPP

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"
#include "elliptic.hpp"
#include "initialGuess.hpp"

#define DINS LIBP_DIR"/solvers/ins/"

using namespace libp;

class insSettings_t: public settings_t {
public:
  insSettings_t(comm_t& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);

  ellipticSettings_t extractVelocitySettings();
  ellipticSettings_t extractPressureSettings();
};

class ins_t;

class subcycler_t: public solver_t {
public:
  mesh_t mesh;

  int cubature;
  ogs::halo_t vTraceHalo;
  kernel_t advectionVolumeKernel;
  kernel_t advectionSurfaceKernel;

  kernel_t subCycleAdvectionKernel;

  int NVfields;
  int order, maxOrder, shiftIndex;
  dfloat nu, T0, dt;

  deviceMemory<dfloat> o_Ue, o_Uh;

  subcycler_t() = default;

  void Report(dfloat time, int tstep){};

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
};

class ins_t: public solver_t {
public:
  mesh_t mesh;
  timeStepper_t timeStepper;

  ogs::halo_t vTraceHalo;
  ogs::halo_t pTraceHalo;

  ellipticSettings_t vSettings, pSettings;
  elliptic_t uSolver, vSolver, wSolver;
  elliptic_t pSolver;

  linearSolver_t uLinearSolver;
  linearSolver_t vLinearSolver;
  linearSolver_t wLinearSolver;
  linearSolver_t pLinearSolver;

  int NVfields, NTfields;

  int NiterU, NiterV, NiterW, NiterP;

  int cubature, pressureIncrement;
  int vDisc_c0, pDisc_c0;
  dfloat velTOL, presTOL;

  dfloat nu;
  dfloat vTau, pTau;

  memory<dfloat> u, p;
  deviceMemory<dfloat> o_u, o_p;

  deviceMemory<dfloat> o_GU;

  deviceMemory<dfloat> o_MU;

  memory<dfloat> Vort;
  deviceMemory<dfloat> o_Vort;

  //extra buffers for solvers
  deviceMemory<dfloat> o_UH, o_VH, o_WH;
  deviceMemory<dfloat> o_rhsU, o_rhsV, o_rhsW;
  deviceMemory<dfloat> o_rhsP, o_PI;

  deviceMemory<dfloat> o_GUH, o_GVH, o_GWH;
  deviceMemory<dfloat> o_GrhsU, o_GrhsV, o_GrhsW;
  deviceMemory<dfloat> o_GrhsP, o_GP, o_GPI;

  //subcycling
  int Nsubcycles;
  timeStepper_t subStepper;
  subcycler_t subcycler;

  kernel_t advectionVolumeKernel;
  kernel_t advectionSurfaceKernel;

  kernel_t divergenceVolumeKernel;
  kernel_t divergenceSurfaceKernel;

  kernel_t gradientVolumeKernel;
  kernel_t gradientSurfaceKernel;

  kernel_t velocityGradientKernel;
  kernel_t diffusionKernel;

  kernel_t velocityRhsKernel;
  kernel_t velocityBCKernel;

  kernel_t pressureRhsKernel;
  kernel_t pressureBCKernel;

  kernel_t pressureIncrementRhsKernel;
  kernel_t pressureIncrementBCKernel;

  kernel_t vorticityKernel;

  kernel_t initialConditionKernel;
  kernel_t maxWaveSpeedKernel;

  ins_t() = default;
  ins_t(platform_t &_platform, mesh_t &_mesh,
        insSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             insSettings_t& _settings);

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat>& U, memory<dfloat>& P, memory<dfloat>& V, std::string fileName);

  dfloat MaxWaveSpeed(deviceMemory<dfloat>& o_U, const dfloat T);

  // void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  void rhs_imex_f(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
  // void rhs_imex_g(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  void rhs_imex_invg(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat gamma, const dfloat time);

  void rhs_subcycle_f(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_QHAT,
                      const dfloat T, const dfloat dt, const memory<dfloat> B,
                      const int order, const int shiftIndex, const int maxOrder);

  void Advection(const dfloat alpha, deviceMemory<dfloat>& o_U,
                 const dfloat beta,  deviceMemory<dfloat>& o_RHS,
                 const dfloat T);
  void Diffusion(const dfloat alpha, deviceMemory<dfloat>& o_U,
                 const dfloat beta,  deviceMemory<dfloat>& o_RHS,
                 const dfloat T);
  void Divergence(const dfloat alpha, deviceMemory<dfloat>& o_U,
                 const dfloat beta,  deviceMemory<dfloat>& o_RHS,
                 const dfloat T);
  void Gradient(const dfloat alpha, deviceMemory<dfloat>& o_P,
                 const dfloat beta,  deviceMemory<dfloat>& o_RHS,
                 const dfloat T);

  void VelocitySolve(deviceMemory<dfloat>& o_U, deviceMemory<dfloat>& o_RHS,
                     const dfloat gamma, const dfloat T);
  void PressureSolve(deviceMemory<dfloat>& o_P, deviceMemory<dfloat>& o_RHS,
                     const dfloat gamma, const dfloat T);
  void PressureIncrementSolve(deviceMemory<dfloat>& o_P, deviceMemory<dfloat>& o_RHS,
                     const dfloat gamma, const dfloat T, const dfloat dt);
};

#endif

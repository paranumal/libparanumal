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

#ifndef FPE_HPP
#define FPE_HPP

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"
#include "elliptic.hpp"

#define DFPE LIBP_DIR"/solvers/fokkerPlanck/"

using namespace libp;

class fpeSettings_t: public settings_t {
public:
  fpeSettings_t(comm_t& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);

  ellipticSettings_t extractEllipticSettings();
};

class fpe_t;

class subcycler_t: public solver_t {
public:
  mesh_t mesh;

  int cubature;
  ogs::halo_t traceHalo;
  kernel_t advectionVolumeKernel;
  kernel_t advectionSurfaceKernel;

  subcycler_t() = default;

  void Report(dfloat time, int tstep){};

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
};

class fpe_t: public solver_t {
public:
  mesh_t mesh;
  timeStepper_t timeStepper;

  ogs::halo_t traceHalo;

  ellipticSettings_t ellipticSettings;
  elliptic_t elliptic;
  linearSolver_t linearSolver;

  int Nfields;

  int cubature;

  dfloat mu;
  dfloat tau;

  memory<dfloat> q;
  deviceMemory<dfloat> o_q;

  deviceMemory<dfloat> o_Mq;

  memory<dfloat> grad;
  deviceMemory<dfloat> o_grad;

  //subcycling
  int Nsubcycles;
  timeStepper_t subStepper;
  subcycler_t subcycler;

  kernel_t advectionVolumeKernel;
  kernel_t advectionSurfaceKernel;
  kernel_t gradientKernel;
  kernel_t diffusionKernel;
  kernel_t diffusionRhsKernel;

  kernel_t initialConditionKernel;
  kernel_t maxWaveSpeedKernel;

  fpe_t() = default;
  fpe_t(platform_t &_platform, mesh_t &_mesh,
        fpeSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             fpeSettings_t& _settings);

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat>& Q, std::string fileName);

  dfloat MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  void rhs_imex_f(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
  void rhs_imex_g(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  void rhs_imex_invg(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat gamma, const dfloat time);

  void rhs_subcycle_f(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_QHAT,
                      const dfloat T, const dfloat dt, const memory<dfloat> B,
                      const int order, const int shiftIndex, const int maxOrder);

  void Advection(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T);
  void Diffusion(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T);
};

#endif

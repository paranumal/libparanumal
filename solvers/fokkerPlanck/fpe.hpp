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

class fpeSettings_t: public settings_t {
public:
  fpeSettings_t(MPI_Comm& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const string filename);

  ellipticSettings_t* extractEllipticSettings();
};

class fpe_t;

class subcycler_t: public solver_t {
public:
  mesh_t& mesh;

  int cubature;
  halo_t* traceHalo;
  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;

  subcycler_t() = delete;
  subcycler_t(fpe_t& fpe);

  ~subcycler_t(){};

  void Report(dfloat time, int tstep){};

  void rhsf(occa::memory& o_q, occa::memory& o_rhs, const dfloat time);
};

class fpe_t: public solver_t {
public:
  mesh_t& mesh;
  TimeStepper::timeStepper_t* timeStepper;

  halo_t* traceHalo;

  ellipticSettings_t *ellipticSettings;
  elliptic_t *elliptic;
  linearSolver_t *linearSolver;

  int Nfields;

  int cubature;

  dfloat mu;
  dfloat tau;

  dfloat *q;
  occa::memory o_q;

  occa::memory o_Mq;

  dfloat *grad;
  occa::memory o_grad;

  //subcycling
  int Nsubcycles;
  TimeStepper::timeStepper_t* subStepper;
  subcycler_t *subcycler;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel gradientKernel;
  occa::kernel diffusionKernel;
  occa::kernel diffusionRhsKernel;

  occa::kernel MassMatrixKernel;

  occa::kernel initialConditionKernel;

  fpe_t() = delete;
  fpe_t(platform_t &_platform, mesh_t &_mesh,
        settings_t& _settings):
    solver_t(_platform, _settings), mesh(_mesh) {}

  ~fpe_t();

  //setup
  static fpe_t& Setup(platform_t& platform, mesh_t& mesh,
                      fpeSettings_t& settings);

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(dfloat* Q, char *fileName);

  void rhsf(occa::memory& o_q, occa::memory& o_rhs, const dfloat time);

  void rhs_imex_f(occa::memory& o_q, occa::memory& o_rhs, const dfloat time);
  void rhs_imex_g(occa::memory& o_q, occa::memory& o_rhs, const dfloat time);

  void rhs_imex_invg(occa::memory& o_q, occa::memory& o_rhs, const dfloat gamma, const dfloat time);

  void rhs_subcycle_f(occa::memory& o_Q, occa::memory& o_QHAT,
                           const dfloat T, const dfloat dt, const dfloat* B,
                           const int order, const int shiftIndex, const int maxOrder);

  void Advection(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T);
  void Diffusion(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T);
};

#endif

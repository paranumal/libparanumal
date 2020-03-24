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

#ifndef INS_HPP
#define INS_HPP

#include "core.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"
#include "elliptic.hpp"

#define DINS LIBP_DIR"/solvers/ins/"

class insSettings_t: public settings_t {
public:
  insSettings_t(MPI_Comm& _comm);
  void report();
  void parseFromFile(occaSettings_t& occaSettings,
                     meshSettings_t& meshSettings,
                     const string filename);

  ellipticSettings_t* extractVelocitySettings();
  ellipticSettings_t* extractPressureSettings();
};

class ins_t;

class subcycler_t: public solver_t {
public:
  int cubature;
  halo_t* vTraceHalo;
  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;

  occa::kernel subCycleAdvectionKernel;

  int NVfields;
  int order, maxOrder, shiftIndex;
  dfloat nu, T0, dt;

  occa::memory o_Ue, o_Uh;

  subcycler_t() = delete;
  subcycler_t(ins_t& ins);

  ~subcycler_t(){};

  void Report(dfloat time, int tstep){};

  void rhsf(occa::memory& o_q, occa::memory& o_rhs, const dfloat time);
};

class ins_t: public solver_t {
public:
  TimeStepper::timeStepper_t* timeStepper;

  halo_t* vTraceHalo;
  halo_t* pTraceHalo;

  ellipticSettings_t *vSettings, *pSettings;
  elliptic_t *uSolver, *vSolver, *wSolver;
  elliptic_t *pSolver;

  linearSolver_t *vLinearSolver;
  linearSolver_t *pLinearSolver;

  int NVfields, NTfields;

  int NiterU, NiterV, NiterW, NiterP;

  int cubature, pressureIncrement;
  int vDisc_c0, pDisc_c0;
  dfloat velTOL, presTOL;

  dfloat nu;
  dfloat vTau, pTau;

  dfloat *u, *p;
  occa::memory o_u, o_p;

  occa::memory o_GU;

  occa::memory o_MU;

  dfloat *Vort;
  occa::memory o_Vort;

  //extra buffers for solvers
  occa::memory o_UH, o_VH, o_WH;
  occa::memory o_rhsU, o_rhsV, o_rhsW;
  occa::memory o_rhsP, o_PI;

  int *mapB; //node-wise boundary flag
  occa::memory o_mapB;

  //subcycling
  int Nsubcycles;
  TimeStepper::timeStepper_t* subStepper;
  subcycler_t *subcycler;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;

  occa::kernel divergenceVolumeKernel;
  occa::kernel divergenceSurfaceKernel;

  occa::kernel gradientVolumeKernel;
  occa::kernel gradientSurfaceKernel;

  occa::kernel velocityGradientKernel;
  occa::kernel diffusionKernel;

  occa::kernel velocityRhsKernel;
  occa::kernel velocityBCKernel;

  occa::kernel pressureRhsKernel;
  occa::kernel pressureBCKernel;

  occa::kernel pressureIncrementRhsKernel;
  occa::kernel pressureIncrementBCKernel;

  occa::kernel MassMatrixKernel;
  occa::kernel vorticityKernel;

  occa::kernel initialConditionKernel;

  ins_t() = delete;
  ins_t(mesh_t& _mesh, linAlg_t& _linAlg, settings_t& _settings):
    solver_t(_mesh, _linAlg, _settings) {}

  ~ins_t();

  //setup
  static ins_t& Setup(mesh_t& mesh, linAlg_t& linAlg,
                      insSettings_t& settings);

  void BoundarySetup();

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(dfloat* U, dfloat* P, dfloat *V, char *fileName);

  // void rhsf(occa::memory& o_q, occa::memory& o_rhs, const dfloat time);

  void rhs_imex_f(occa::memory& o_q, occa::memory& o_rhs, const dfloat time);
  // void rhs_imex_g(occa::memory& o_q, occa::memory& o_rhs, const dfloat time);

  void rhs_imex_invg(occa::memory& o_q, occa::memory& o_rhs, const dfloat gamma, const dfloat time);

  void rhs_subcycle_f(occa::memory& o_Q, occa::memory& o_QHAT,
                           const dfloat T, const dfloat dt, const dfloat* B,
                           const int order, const int shiftIndex, const int maxOrder);

  void Advection(const dfloat alpha, occa::memory& o_U,
                 const dfloat beta,  occa::memory& o_RHS,
                 const dfloat T);
  void Diffusion(const dfloat alpha, occa::memory& o_U,
                 const dfloat beta,  occa::memory& o_RHS,
                 const dfloat T);
  void Divergence(const dfloat alpha, occa::memory& o_U,
                 const dfloat beta,  occa::memory& o_RHS,
                 const dfloat T);
  void Gradient(const dfloat alpha, occa::memory& o_P,
                 const dfloat beta,  occa::memory& o_RHS,
                 const dfloat T);

  void VelocitySolve(occa::memory& o_U, occa::memory& o_RHS,
                     const dfloat gamma, const dfloat T);
  void PressureSolve(occa::memory& o_P, occa::memory& o_RHS,
                     const dfloat gamma, const dfloat T);
  void PressureIncrementSolve(occa::memory& o_P, occa::memory& o_RHS,
                     const dfloat gamma, const dfloat T, const dfloat dt);
};

#endif

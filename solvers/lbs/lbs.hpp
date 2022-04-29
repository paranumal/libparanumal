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

#ifndef LBS_HPP
#define LBS_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#define DLBS LIBP_DIR"/solvers/lbs/"

using namespace libp;

class lbsSettings_t: public settings_t {
public:
  lbsSettings_t(comm_t& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);
};

class lbs_t: public solver_t {
public:
  mesh_t mesh;

  int Nfields;
  int Nmacro;
  int Npmlfields;
  int velModel;

  timeStepper_t timeStepper;

  ogs::halo_t traceHalo;
  memory<ogs::halo_t> multirateTraceHalo;

  // dfloat RT, c, tauInv, Ma, Re, nu; // Flow parameters
  dfloat RT, c, tauInv, Re, nu, alpha; // Flow parameters

  // Pml
  int pmlOrder;
  dfloat  sigmaXmax, sigmaYmax, sigmaZmax;
  memory<dfloat> pmlSigma;
  dfloat pmlAlpha;

  // Flag for using cubature integration for sigma terms in pml
  int pmlcubature;

  // Flag for semi-analytic timestepping
  int semiAnalytic;

  memory<dfloat> q;
  deviceMemory<dfloat> o_q;

  // external forcing in velocity space
  memory<dfloat> F;
  deviceMemory<dfloat> o_F;

  // Macro quantities i.e. density + velocity
  memory<dfloat> U;
  deviceMemory<dfloat> o_U;

  memory<dfloat> LBM;
  deviceMemory<dfloat> o_LBM;

  memory<int> LMAP;
  deviceMemory<int> o_LMAP;

  deviceMemory<dfloat> o_Mq;

  memory<dfloat> Vort, VortMag;
  deviceMemory<dfloat> o_Vort, o_VortMag;

  deviceMemory<dfloat> o_pmlSigma;

  kernel_t collisionKernel;
  kernel_t forcingKernel;
  kernel_t momentsKernel;
  kernel_t phaseFieldKernel;

  kernel_t volumeKernel;
  kernel_t surfaceKernel;
  kernel_t relaxationKernel;

  kernel_t pmlVolumeKernel;
  kernel_t pmlSurfaceKernel;
  kernel_t pmlRelaxationKernel;

  kernel_t vorticityKernel;

  kernel_t initialConditionKernel;

  lbs_t() = default;
  lbs_t(platform_t &_platform, mesh_t &_mesh,
        lbsSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             lbsSettings_t& _settings);

  void PmlSetup();

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat>& Q, memory<dfloat>& V, std::string fileName);

  dfloat MaxWaveSpeed();

  void rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T);

  void rhsVolume(dlong N, deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T);

  void rhsSurface(dlong N, deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T);

  void latticeSetup(); 
};
#endif





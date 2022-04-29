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

#ifndef BNS_HPP
#define BNS_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#define DBNS LIBP_DIR"/solvers/bns/"

using namespace libp;

class bnsSettings_t: public settings_t {
public:
  bnsSettings_t(comm_t& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);
};

class bns_t: public solver_t {
public:
  mesh_t mesh;

  int Nfields;
  int Npmlfields;

  timeStepper_t timeStepper;

  ogs::halo_t traceHalo;
  memory<ogs::halo_t> multirateTraceHalo;

  dfloat RT, c, tauInv, Ma, Re, nu; // Flow parameters

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

  deviceMemory<dfloat> o_Mq;

  memory<dfloat> Vort, VortMag;
  deviceMemory<dfloat> o_Vort, o_VortMag;

  deviceMemory<dfloat> o_pmlSigma;

  kernel_t volumeKernel;
  kernel_t surfaceKernel;
  kernel_t relaxationKernel;

  kernel_t pmlVolumeKernel;
  kernel_t pmlSurfaceKernel;
  kernel_t pmlRelaxationKernel;

  kernel_t vorticityKernel;

  kernel_t initialConditionKernel;

  bns_t() = default;
  bns_t(platform_t &_platform, mesh_t &_mesh,
              bnsSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             bnsSettings_t& _settings);

  void PmlSetup();

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat>& Q, memory<dfloat>& V, std::string fileName);

  dfloat MaxWaveSpeed();

  void rhsf_pml(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
                deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS, const dfloat T);

  void rhsf_MR_pml(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
                   deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS,
                   deviceMemory<dfloat>& o_fQM, const dfloat T, const int lev);

  //seperate components of rhs evaluation
  void rhsVolume(dlong N, deviceMemory<dlong>& o_ids,
                 deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T);
  void rhsPmlVolume(dlong N, deviceMemory<dlong>& o_ids, deviceMemory<dlong>& o_pmlids,
                    deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
                    deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS, const dfloat T);
  void rhsRelaxation(dlong N, deviceMemory<dlong>& o_ids,
                     deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS);
  void rhsPmlRelaxation(dlong N, deviceMemory<dlong>& o_ids, deviceMemory<dlong>& o_pmlids,
                        deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
                        deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS);
  void rhsSurface(dlong N, deviceMemory<dlong>& o_ids,
                  deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T);
  void rhsPmlSurface(dlong N, deviceMemory<dlong>& o_ids, deviceMemory<dlong>& o_pmlids,
                     deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
                     deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS, const dfloat T);
  void rhsSurfaceMR(dlong N, deviceMemory<dlong>& o_ids,
                    deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS,
                    deviceMemory<dfloat>& o_fQM, const dfloat T);
  void rhsPmlSurfaceMR(dlong N, deviceMemory<dlong>& o_ids, deviceMemory<dlong>& o_pmlids,
                       deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
                       deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS,
                       deviceMemory<dfloat>& o_fQM, const dfloat T);
};
#endif


/*
  // Some Iso-surfacing variables
  int isoField, isoColorField, isoNfields, isoNlevels, isoMaxNtris, *isoNtris;
  dfloat isoMinVal, isoMaxVal, *isoLevels, *isoq;
  size_t isoMax;

  deviceMemory<dfloat> o_isoLevels, o_isoq, o_isoNtris;

  // MRSAAB Coefficients
  dfloat *MRSAAB_A, *MRSAAB_B, *MRSAAB_C, *MRAB_A, *MRAB_B, *MRAB_C;
  // SARK and RK3 Coefficients
  dfloat RK_A[5][5], RK_B[5], RK_C[5], SARK_A[5][5], SARK_B[5], SARK_C[5];

  int *isoGNlevels, isoGNgroups;
  dfloat **isoGLvalues;

  deviceMemory<dfloat> *o_isoGLvalues;

  // NBN: add storage for compacted isosurf data for gmsh write
  std::vector<dfloat> iso_nodes;
  std::vector<int> iso_tris;
*/



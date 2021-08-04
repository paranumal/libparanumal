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

#ifndef LBS_HPP
#define LBS_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#define DLBS LIBP_DIR"/solvers/lbs/"

class lbsSettings_t: public settings_t {
public:
  lbsSettings_t(MPI_Comm& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const string filename);
};

class lbs_t: public solver_t {
public:
  mesh_t& mesh;

  int Nfields;
  int Nmacro;
  int Npmlfields;
  int velModel;

  TimeStepper::timeStepper_t* timeStepper;

  halo_t* traceHalo;
  halo_t** multirateTraceHalo;

  // dfloat RT, c, tauInv, Ma, Re, nu; // Flow parameters
  dfloat RT, c, tauInv, Re, nu, alpha; // Flow parameters

  // Pml
  int pmlOrder;
  dfloat  sigmaXmax, sigmaYmax, sigmaZmax;
  dfloat *pmlSigma;
  dfloat pmlAlpha;

  // Flag for using cubature integration for sigma terms in pml
  int pmlcubature;

  // Flag for semi-analytic timestepping
  int semiAnalytic;

  dfloat *q;
  occa::memory o_q;
  
  // Macro quantities i.e. density + velocity
  dfloat *U; 
  occa::memory o_U; 

  dfloat *LBM; 
  occa::memory o_LBM;

  int *LMAP; 
  occa::memory o_LMAP; 

  occa::memory o_Mq;

  dfloat *Vort, *VortMag;
  occa::memory o_Vort, o_VortMag;

  occa::memory o_pmlSigma;

  occa::kernel collisionKernel; 
  occa::kernel momentsKernel; 
  occa::kernel phaseFieldKernel; 

  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel relaxationKernel;

  occa::kernel pmlVolumeKernel;
  occa::kernel pmlSurfaceKernel;
  occa::kernel pmlRelaxationKernel;

  occa::kernel vorticityKernel;

  occa::kernel initialConditionKernel;

  lbs_t() = delete;
  lbs_t(platform_t &_platform, mesh_t &_mesh,
              lbsSettings_t& _settings):
    solver_t(_platform, _settings), mesh(_mesh) {}

  ~lbs_t();

  //setup
  static lbs_t& Setup(platform_t& platform, mesh_t& mesh,
                      lbsSettings_t& settings);

  void PmlSetup();

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(dfloat* Q, dfloat* V, char *fileName);

  dfloat MaxWaveSpeed();

  void rhsf_pml(occa::memory& o_Q, occa::memory& o_pmlQ,
                occa::memory& o_RHS, occa::memory& o_pmlRHS, const dfloat T);

  // void rhsf_MR_pml(occa::memory& o_Q, occa::memory& o_pmlQ,
  //                  occa::memory& o_RHS, occa::memory& o_pmlRHS,
  //                  occa::memory& o_fQM, const dfloat T, const int lev);

  //seperate components of rhs evaluation
  void rhsVolume(dlong N, occa::memory& o_ids,
                 occa::memory& o_Q, occa::memory& o_RHS, const dfloat T);
  // void rhsPmlVolume(dlong N, occa::memory& o_ids, occa::memory& o_pmlids,
                    // occa::memory& o_Q, occa::memory& o_pmlQ,
                    // occa::memory& o_RHS, occa::memory& o_pmlRHS, const dfloat T);
  // void rhsRelaxation(dlong N, occa::memory& o_ids,
                     // occa::memory& o_Q, occa::memory& o_RHS);
  // void rhsPmlRelaxation(dlong N, occa::memory& o_ids, occa::memory& o_pmlids,
                        // occa::memory& o_Q, occa::memory& o_pmlQ,
                        // occa::memory& o_RHS, occa::memory& o_pmlRHS);
  void rhsSurface(dlong N, occa::memory& o_ids,
                  occa::memory& o_Q, occa::memory& o_RHS, const dfloat T);
  // void rhsPmlSurface(dlong N, occa::memory& o_ids, occa::memory& o_pmlids,
  //                    occa::memory& o_Q, occa::memory& o_pmlQ,
  //                    occa::memory& o_RHS, occa::memory& o_pmlRHS, const dfloat T);
  // void rhsSurfaceMR(dlong N, occa::memory& o_ids,
  //                   occa::memory& o_Q, occa::memory& o_RHS,
  //                   occa::memory& o_fQM, const dfloat T);
  // void rhsPmlSurfaceMR(dlong N, occa::memory& o_ids, occa::memory& o_pmlids,
  //                      occa::memory& o_Q, occa::memory& o_pmlQ,
  //                      occa::memory& o_RHS, occa::memory& o_pmlRHS,
  //                      occa::memory& o_fQM, const dfloat T);

  void latticeSetup(); 
};
#endif





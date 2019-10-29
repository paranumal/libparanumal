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

#ifndef CNS_HPP
#define CNS_HPP 1

#include "core.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#define DCNS LIBP_DIR"/solvers/cns/"

class cnsSettings_t: public settings_t {
public:
  cnsSettings_t(MPI_Comm& _comm);
  void report();
};

class cns_t: public solver_t {
public:
  int Nfields;
  int Ngrads;

  dfloat mu;
  dfloat gamma;

  int cubature;
  int isothermal;

  timeStepper_t* timeStepper;

  halo_t* fieldTraceHalo;
  halo_t* gradTraceHalo;

  dfloat *q;
  occa::memory o_q;

  dfloat *gradq;
  occa::memory o_gradq;

  dfloat *Vort;
  occa::memory o_Vort;

  occa::memory o_Mq;

  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel cubatureVolumeKernel;
  occa::kernel cubatureSurfaceKernel;

  occa::kernel gradVolumeKernel;
  occa::kernel gradSurfaceKernel;

  occa::kernel vorticityKernel;

  occa::kernel constrainKernel;

  occa::kernel MassMatrixKernel;

  occa::kernel initialConditionKernel;

  cns_t() = delete;
  cns_t(mesh_t& _mesh, linAlg_t& _linAlg):
    solver_t(_mesh, _linAlg) {}

  ~cns_t();

  //setup
  static cns_t& Setup(mesh_t& mesh, linAlg_t& linAlg);

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(dfloat* Q, dfloat *V, char *fileName);

  void rhsf(occa::memory& o_q, occa::memory& o_rhs, const dfloat time);
};

#endif

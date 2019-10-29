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

#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

#include <occa.hpp>
#include "types.h"
#include "utils.hpp"
#include "solver.hpp"
#include "settings.hpp"

//virtual base time stepper class
class timeStepper_t {
public:
  dlong N;
  dlong Nhalo;

  solver_t& solver;

  MPI_Comm& comm;
  occa::device& device;
  settings_t& settings;
  occa::properties& props;

  dfloat dt;

  timeStepper_t(dlong _N, dlong _Nhalo, solver_t& _solver):
    N(_N), Nhalo(_Nhalo),
    solver(_solver),
    comm(_solver.comm),
    device(_solver.device),
    settings(_solver.settings),
    props(_solver.props) {}

  static timeStepper_t* Setup(dlong N, dlong Nhalo, solver_t& solver);

  virtual ~timeStepper_t() {};

  virtual void Init()=0;
  virtual void Run(occa::memory& o_q, dfloat start, dfloat end)=0;

  void SetTimeStep(dfloat dt_) {dt = dt_;};
};

class lserk4: public timeStepper_t {
private:
  int Nrk;
  dfloat *rka, *rkb, *rkc;

  occa::memory o_rhsq;
  occa::memory o_resq;

  occa::kernel updateKernel;

  void Step(occa::memory& o_q, dfloat time, dfloat dt);

public:
  lserk4(dlong N, dlong _Nhalo, solver_t& solver);
  ~lserk4();

  void Init();
  void Run(occa::memory& o_q, dfloat start, dfloat end);
};

class dopri5: public timeStepper_t {
private:
  int Nrk;

  dlong Nblock;

  dfloat *rkC, *rkA, *rkE;
  occa::memory o_rkA, o_rkE;

  dfloat *errtmp;

  occa::memory o_rhsq;
  occa::memory o_rkq;
  occa::memory o_rkrhsq;
  occa::memory o_rkerr;

  occa::memory o_saveq;

  occa::memory o_errtmp;

  occa::kernel rkUpdateKernel;
  occa::kernel rkStageKernel;
  occa::kernel rkErrorEstimateKernel;

  dfloat dtMIN; //minumum allowed timestep
  dfloat ATOL;  //absolute error tolerance
  dfloat RTOL;  //relative error tolerance
  dfloat safe;   //safety factor

  //error control parameters
  dfloat beta;
  dfloat factor1;
  dfloat factor2;

  dfloat exp1;
  dfloat invfactor1;
  dfloat invfactor2;
  dfloat facold;
  dfloat sqrtinvNtotal;

  void Step(occa::memory& o_q, dfloat time, dfloat dt);

  dfloat Estimater(occa::memory& o_q);

public:
  dopri5(dlong N, dlong _Nhalo, solver_t& solver);
  ~dopri5();

  void Init();
  void Run(occa::memory& o_q, dfloat start, dfloat end);
};



#endif

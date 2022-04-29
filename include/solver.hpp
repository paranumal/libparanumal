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

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "settings.hpp"
#include "platform.hpp"
#include "operator.hpp"

namespace libp {

class solver_t: public operator_t {
public:
  platform_t platform;
  settings_t settings;
  comm_t comm;

  solver_t() = default;

  solver_t(platform_t& _platform, settings_t& _settings, comm_t _comm):
    platform(_platform),
    settings(_settings),
    comm(_comm) {};

  virtual void Run() {
    LIBP_FORCE_ABORT("Run not implemented in this solver");
  };
  virtual void Report(dfloat time, int tstep) {
    LIBP_FORCE_ABORT("Report not implemented in this solver");
  }

  //Full rhs evaluation of solver in form dq/dt = rhsf(q,t)
  virtual void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time) {
    LIBP_FORCE_ABORT("rhsf not implemented in this solver");
  }

  // Partial rhs evaluation of f with solver in form dq/dt = f(q,t) + g(q,t)
  virtual void rhs_imex_f(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time) {
    LIBP_FORCE_ABORT("rhs_imex_f not implemented in this solver");
  }

  // Partial rhs evaluation of g with solver in form dq/dt = f(q,t) + g(q,t)
  virtual void rhs_imex_g(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time) {
    LIBP_FORCE_ABORT("rhs_imex_g not implemented in this solver");
  }

  // Inversion of g function with solver in form dq/dt = f(q,t) + g(q,t)
  //  Solves gamma*q - g(q,t) = rhs for q
  virtual void rhs_imex_invg(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_q, const dfloat gamma, const dfloat time) {
    LIBP_FORCE_ABORT("rhs_imex_invg not implemented in this solver");
  }

  // Evolve rhs f function via a sub-timestepper
  virtual void rhs_subcycle_f(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_QHAT,
                           const dfloat T, const dfloat dt, const memory<dfloat> B,
                           const int order, const int shiftIndex, const int maxOrder) {
    LIBP_FORCE_ABORT("Subcycling not implemented in this solver");
  }

  //Full rhs evaluation of solver in form dq/dt = rhsf(q,t) for multi-rate timestepping
  virtual void rhsf_MR(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_fQM, const dfloat time, const int level) {
    LIBP_FORCE_ABORT("rhsf_MR not implemented in this solver");
  }

  //Full rhs evaluation of solver in form dq/dt = rhsf(q,t) with a perfectly matched layer (PML)
  virtual void rhsf_pml(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_pmlq,
                        deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_pmlrhs, const dfloat time) {
    LIBP_FORCE_ABORT("rhsf_pml not implemented in this solver");
  }

  //Full rhs evaluation of solver in form dq/dt = rhsf(q,t) for multi-rate timestepping with a PML
  virtual void rhsf_MR_pml(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_pmlq,
                           deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_pmlrhs,
                           deviceMemory<dfloat>& o_fQM, const dfloat time, const int level) {
    LIBP_FORCE_ABORT("rhsf_MR_pml not implemented in this solver");
  }

  //Evaluation of solver as a operator in the form A(q)
  virtual void Operator(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_Aq) {
    LIBP_FORCE_ABORT("Operator not implemented in this solver");
  }
};

} //namespace libp

#endif

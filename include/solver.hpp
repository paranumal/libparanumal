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

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "mesh.hpp"
#include "linAlg.hpp"

class solver_t {
public:
  mesh_t& mesh;

  MPI_Comm& comm;
  occa::device& device;
  settings_t& settings;
  occa::properties& props;
  linAlg_t& linAlg;

  solver_t() = delete;

  solver_t(mesh_t& _mesh, linAlg_t& _linAlg,
           settings_t& _settings):
    mesh(_mesh),
    comm(_mesh.comm),
    device(_mesh.device),
    settings(_settings),
    props(_mesh.props),
    linAlg(_linAlg) {};

  virtual ~solver_t(){}

  virtual void Run() {
    LIBP_ABORT(string("Run not implemented in this solver"))
  };
  virtual void Report(dfloat time=0.0, int tstep=0) {
    LIBP_ABORT(string("Report not implemented in this solver"))
  }

  //Full rhs evaluation of solver in form dq/dt = rhsf(q,t)
  virtual void rhsf(occa::memory& o_q, occa::memory& o_rhs, const dfloat time) {
    LIBP_ABORT(string("rhsf not implemented in this solver"))
  }

  // Partial rhs evaluation of f with solver in form dq/dt = f(q,t) + g(q,t)
  virtual void rhs_imex_f(occa::memory& o_q, occa::memory& o_rhs, const dfloat time) {
    LIBP_ABORT(string("rhs_imex_f not implemented in this solver"))
  }

  // Partial rhs evaluation of g with solver in form dq/dt = f(q,t) + g(q,t)
  virtual void rhs_imex_g(occa::memory& o_q, occa::memory& o_rhs, const dfloat time) {
    LIBP_ABORT(string("rhs_imex_g not implemented in this solver"))
  }

  // Inversion of g function with solver in form dq/dt = f(q,t) + g(q,t)
  //  Solves gamma*q - g(q,t) = rhs for q
  virtual int rhs_imex_invg(occa::memory& o_rhs, occa::memory& o_q, const dfloat gamma, const dfloat time) {
    LIBP_ABORT(string("rhs_imex_invg not implemented in this solver"))
    return 0;
  }

  // Evolve rhs f function via a sub-timestepper
  virtual void rhs_subcycle_f(occa::memory& o_Q, occa::memory& o_QHAT,
                           const dfloat T, const dfloat dt, const dfloat* B,
                           const int order, const int shiftIndex, const int maxOrder) {
    LIBP_ABORT(string("Subcycling not implemented in this solver"))
  }

  //Full rhs evaluation of solver in form dq/dt = rhsf(q,t) for multi-rate timestepping
  virtual void rhsf_MR(occa::memory& o_q, occa::memory& o_rhs, occa::memory& o_fQM, const dfloat time, const int level) {
    LIBP_ABORT(string("rhsf_MR not implemented in this solver"))
  }

  //Full rhs evaluation of solver in form dq/dt = rhsf(q,t) with a perfectly matched layer (PML)
  virtual void rhsf_pml(occa::memory& o_q, occa::memory& o_pmlq,
                        occa::memory& o_rhs, occa::memory& o_pmlrhs, const dfloat time) {
    LIBP_ABORT(string("rhsf_pml not implemented in this solver"))
  }

  //Full rhs evaluation of solver in form dq/dt = rhsf(q,t) for multi-rate timestepping with a PML
  virtual void rhsf_MR_pml(occa::memory& o_q, occa::memory& o_pmlq,
                           occa::memory& o_rhs, occa::memory& o_pmlrhs,
                           occa::memory& o_fQM, const dfloat time, const int level) {
    LIBP_ABORT(string("rhsf_MR_pml not implemented in this solver"))
  }

  //Evaluation of solver as a operator in the form A(q)
  virtual void Operator(occa::memory& o_q, occa::memory& o_Aq) {
    LIBP_ABORT(string("Operator not implemented in this solver"))
  }
};

#endif
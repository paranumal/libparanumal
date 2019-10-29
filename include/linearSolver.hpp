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

#ifndef LINEARSOLVER_HPP
#define LINEARSOLVER_HPP

#include <occa.hpp>
#include "types.h"
#include "utils.hpp"
#include "core.hpp"
#include "settings.hpp"
#include "solver.hpp"
#include "precon.hpp"
#include "linAlg.hpp"

//virtual base linear solver class
class linearSolver_t {
public:
  MPI_Comm& comm;
  occa::device& device;
  settings_t& settings;
  occa::properties& props;
  mesh_t& mesh;
  linAlg_t& linAlg;

  dlong N;

  linearSolver_t(solver_t& solver):
    comm(solver.comm),
    device(solver.device),
    settings(solver.settings),
    props(solver.props),
    mesh(solver.mesh),
    linAlg(solver.linAlg) {}

  static linearSolver_t* Setup(solver_t& solver);

  virtual void Init(int _weighted, occa::memory& o_weight)=0;
  virtual int Solve(solver_t& solver, precon_t& precon,
                    occa::memory& o_x, occa::memory& o_rhs,
                    const dfloat tol, const int MAXIT, const int verbose)=0;

  virtual ~linearSolver_t(){}
};

//Preconditioned Conjugate Gradient
class pcg: public linearSolver_t {
private:
  occa::memory o_p, o_Ap, o_z, o_Ax, o_w;

  dfloat* tmprdotr;
  occa::memory h_tmprdotr;
  occa::memory o_tmprdotr;

  int flexible, weighted;

  occa::kernel updatePCGKernel;

  dfloat UpdatePCG(const dfloat alpha, occa::memory &o_x, occa::memory &o_r);

public:
  pcg(solver_t& solver);
  ~pcg();

  void Init(int _weighted, occa::memory& o_weight);
  int Solve(solver_t& solver, precon_t& precon,
            occa::memory& o_x, occa::memory& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);
};

//Non-Blocking Preconditioned Conjugate Gradient
class nbpcg: public linearSolver_t {
private:
  occa::memory o_p, o_s, o_S, o_z, o_Z, o_Ax, o_weight;

  dfloat* tmpdots;
  occa::memory h_tmpdots;
  occa::memory o_tmpdots;

  int weighted;

  occa::kernel update1NBPCGKernel;
  occa::kernel update2NBPCGKernel;

  dfloat *localdots, *globaldots;

  MPI_Request request;
  MPI_Status  status;

  void Update1NBPCG(const dfloat beta);
  void Update2NBPCG(const dfloat alpha, occa::memory &o_r);

public:
  nbpcg(solver_t& solver);
  ~nbpcg();

  void Init(int _weighted, occa::memory& o_weight_);
  int Solve(solver_t& solver, precon_t& precon,
            occa::memory& o_x, occa::memory& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);
};

//Non-Blocking Flexible Preconditioned Conjugate Gradient
class nbfpcg: public linearSolver_t {
private:
  occa::memory o_u, o_p, o_w, o_n, o_m, o_s, o_z, o_q, o_Ax, o_weight;

  dfloat* tmpdots;
  occa::memory h_tmpdots;
  occa::memory o_tmpdots;

  int weighted;

  occa::kernel update0NBFPCGKernel;
  occa::kernel update1NBFPCGKernel;

  dfloat *localdots, *globaldots;

  MPI_Request request;
  MPI_Status  status;

  void Update0NBFPCG(occa::memory &o_r);
  void Update1NBFPCG(const dfloat alpha, occa::memory &o_x, occa::memory &o_r);

public:
  nbfpcg(solver_t& solver);
  ~nbfpcg();

  void Init(int _weighted, occa::memory& o_weight_);
  int Solve(solver_t& solver, precon_t& precon,
            occa::memory& o_x, occa::memory& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);
};

#endif

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

#include "platform.hpp"
#include "solver.hpp"
#include "precon.hpp"

//virtual base linear solver class
class linearSolver_t {
public:
  platform_t& platform;
  settings_t& settings;

  dlong N;
  dlong Nhalo;

  linearSolver_t(dlong _N, dlong _Nhalo,
                 platform_t& _platform, settings_t& _settings):
    platform(_platform), settings(_settings), N(_N), Nhalo(_Nhalo) {}

  static linearSolver_t* Setup(dlong _N, dlong _Nhalo,
                               platform_t& platform, settings_t& settings,
                               int _weighted, occa::memory& _o_weight);

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
  pcg(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings,
       int _weighted, occa::memory& _o_weight);
  ~pcg();

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
  nbpcg(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings,
       int _weighted, occa::memory& _o_weight);
  ~nbpcg();

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
  nbfpcg(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings,
       int _weighted, occa::memory& _o_weight);
  ~nbfpcg();

  int Solve(solver_t& solver, precon_t& precon,
            occa::memory& o_x, occa::memory& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);
};

#endif

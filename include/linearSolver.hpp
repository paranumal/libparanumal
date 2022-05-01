/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Anthony Austin

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

#include "core.hpp"
#include "platform.hpp"
#include "solver.hpp"
#include "precon.hpp"
#include "initialGuess.hpp"

namespace libp {

namespace LinearSolver { class linearSolverBase_t; }

/* General LinearSolver object*/
class linearSolver_t {
 public:
  linearSolver_t() = default;

  /*Generic setup. Create a Solver object and wrap it in a shared_ptr*/
  template<class Solver, class... Args>
  void Setup(Args&& ... args) {
    ls = std::make_shared<Solver>(args...);

    /*Make an initial guess strategy if we dont have one setup yet*/
    if (ig==nullptr) {
      MakeDefaultInitialGuessStrategy();
    }
  }

  /*Generic setup. Create a InitialGuess object and wrap it in a shared_ptr*/
  template<class InitialGuess, class... Args>
  void SetupInitialGuess(Args&& ... args) {
    ig = std::make_shared<InitialGuess>(args...);
  }

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);

 private:
  std::shared_ptr<LinearSolver::linearSolverBase_t> ls=nullptr;
  std::shared_ptr<InitialGuess::initialGuessStrategy_t> ig=nullptr;

  void MakeDefaultInitialGuessStrategy();

  void assertInitialized();
};


namespace LinearSolver {

//virtual base linear solver class
class linearSolverBase_t {
public:
  platform_t platform;
  settings_t settings;
  comm_t comm;

  dlong N;
  dlong Nhalo;

  linearSolverBase_t(dlong _N, dlong _Nhalo,
                 platform_t& _platform, settings_t& _settings, comm_t _comm):
    platform(_platform), settings(_settings), comm(_comm),
    N(_N), Nhalo(_Nhalo) {}

  virtual int Solve(operator_t& linearOperator, operator_t& precon,
                    deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs,
                    const dfloat tol, const int MAXIT, const int verbose)=0;
};

//Preconditioned Conjugate Gradient
class pcg: public linearSolverBase_t {
private:
  deviceMemory<dfloat> o_p, o_Ap, o_z, o_Ax;

  pinnedMemory<dfloat> rdotr;
  deviceMemory<dfloat> o_rdotr;

  int flexible;

  kernel_t updatePCGKernel;

  dfloat UpdatePCG(const dfloat alpha, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_r);

public:
  pcg(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings, comm_t _comm);

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);
};

//Preconditioned GMRES
class pgmres: public linearSolverBase_t {
private:
  deviceMemory<dfloat> o_Ax, o_z, o_r;
  memory<deviceMemory<dfloat>> o_V;

  int restart;

  memory<dfloat> H, sn, cs, s, y;

  void UpdateGMRES(deviceMemory<dfloat>& o_x, const int I);

public:
  pgmres(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings, comm_t _comm);

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);
};

// Preconditioned MINRES
class pminres : public linearSolverBase_t {
private:
  deviceMemory<dfloat> o_p;
  deviceMemory<dfloat> o_z;
  deviceMemory<dfloat> o_r;
  deviceMemory<dfloat> o_r_old;
  deviceMemory<dfloat> o_q;
  deviceMemory<dfloat> o_q_old;

  kernel_t updateMINRESKernel;

  dfloat innerProd(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_y);
  void UpdateMINRES(const dfloat ma2, const dfloat ma3, const dfloat alpha, const dfloat beta);

public:
  pminres(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, comm_t _comm);

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);
};

//Non-Blocking Preconditioned Conjugate Gradient
class nbpcg: public linearSolverBase_t {
private:
  deviceMemory<dfloat> o_p, o_s, o_S, o_z, o_Z, o_Ax;

  pinnedMemory<dfloat> dots;
  deviceMemory<dfloat> o_dots;

  kernel_t update1NBPCGKernel;
  kernel_t update2NBPCGKernel;

  Comm::request_t request;

  void Update1NBPCG(const dfloat beta);
  void Update2NBPCG(const dfloat alpha, deviceMemory<dfloat>& o_r);

public:
  nbpcg(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings, comm_t _comm);

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);
};

//Non-Blocking Flexible Preconditioned Conjugate Gradient
class nbfpcg: public linearSolverBase_t {
private:
  deviceMemory<dfloat> o_u, o_p, o_w, o_n, o_m, o_s, o_z, o_q, o_Ax;

  pinnedMemory<dfloat> dots;
  deviceMemory<dfloat> o_dots;

  kernel_t update0NBFPCGKernel;
  kernel_t update1NBFPCGKernel;

  Comm::request_t request;

  void Update0NBFPCG(deviceMemory<dfloat>& o_r);
  void Update1NBFPCG(const dfloat alpha, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_r);

public:
  nbfpcg(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings, comm_t _comm);

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);
};

} //namespace LinearSolver

} //namespace libp

#endif

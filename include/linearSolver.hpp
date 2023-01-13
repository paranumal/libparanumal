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
template <typename T>
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
            deviceMemory<T>& o_x, deviceMemory<T>& o_rhs,
            const T tol, const int MAXIT, const int verbose);

  bool isInitialized();

 private:
  std::shared_ptr<LinearSolver::linearSolverBase_t> ls=nullptr;
  std::shared_ptr<InitialGuess::initialGuessStrategy_t> ig=nullptr;

  void MakeDefaultInitialGuessStrategy();

  void assertInitialized();
};

  template class linearSolver_t<double>;
  template class linearSolver_t<float>;

  
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
	    deviceMemory<float>& o_x, deviceMemory<float>& o_rhs,
	    const float tol, const int MAXIT, const int verbose){ printf("sjfosef\n"); return 0;}

  virtual int Solve(operator_t& linearOperator, operator_t& precon,
	    deviceMemory<double>& o_x, deviceMemory<double>& o_rhs,
	    const double tol, const int MAXIT, const int verbose){ printf("sjfosef\n"); return 0;}

};

  
//Preconditioned Conjugate Gradient
template<typename T>
class pcg: public linearSolverBase_t {
private:
  int flexible;

  kernel_t updatePCGKernel;

  T UpdatePCG(const T alpha,
                   deviceMemory<T>& o_p,
                   deviceMemory<T>& o_Ap,
                   deviceMemory<T>& o_x,
                   deviceMemory<T>& o_r);

public:
  pcg(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings, comm_t _comm);

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<T>& o_x, deviceMemory<T>& o_rhs,
            const T tol, const int MAXIT, const int verbose);
};

template class pcg<float>;
template class pcg<double>;

  
//Preconditioned GMRES
template <typename T>
class pgmres: public linearSolverBase_t {
private:
  int restart;

  memory<T> H, sn, cs, s, y;

  void UpdateGMRES(memory<deviceMemory<T>>& o_V,
                   deviceMemory<T>& o_x,
                   const int I);

public:
  pgmres(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings, comm_t _comm);

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<T>& o_x, deviceMemory<T>& o_rhs,
            const T tol, const int MAXIT, const int verbose);
};

  template class pgmres<float>;
template class pgmres<double>;

  
// Preconditioned MINRES
template <typename T>
class pminres : public linearSolverBase_t {  
private:
  kernel_t updateMINRESKernel;

  void UpdateMINRES(const T ma2,
                    const T ma3,
                    const T alpha,
                    const T beta,
                    deviceMemory<T>& o_z,
                    deviceMemory<T>& o_q_old,
                    deviceMemory<T>& o_q,
                    deviceMemory<T>& o_r_old,
                    deviceMemory<T>& o_r,
                    deviceMemory<T>& o_p);

public:
  pminres(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, comm_t _comm);

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<T>& o_x, deviceMemory<T>& o_rhs,
            const T tol, const int MAXIT, const int verbose);
};

template class pminres<float>;
template class pminres<double>;
  
//Non-Blocking Preconditioned Conjugate Gradient
template <typename T>
class nbpcg: public linearSolverBase_t {
private:
  pinnedMemory<T> dots;

  kernel_t update1NBPCGKernel;
  kernel_t update2NBPCGKernel;

  Comm::request_t request;

  void Update1NBPCG(const T beta,
                    deviceMemory<T>& o_z,
                    deviceMemory<T>& o_Z,
                    deviceMemory<T>& o_p,
                    deviceMemory<T>& o_s);
  void Update2NBPCG(const T alpha,
                    deviceMemory<T>& o_s,
                    deviceMemory<T>& o_S,
                    deviceMemory<T>& o_r,
                    deviceMemory<T>& o_z);

public:
  nbpcg(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings, comm_t _comm);

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<T>& o_x, deviceMemory<T>& o_rhs,
            const T tol, const int MAXIT, const int verbose);
};

template class nbpcg<float>;
template class nbpcg<double>;
  
//Non-Blocking Flexible Preconditioned Conjugate Gradient
template <typename T>
class nbfpcg: public linearSolverBase_t {
private:
  pinnedMemory<T> dots;

  kernel_t update0NBFPCGKernel;
  kernel_t update1NBFPCGKernel;

  Comm::request_t request;

  void Update0NBFPCG(deviceMemory<T>& o_u,
                     deviceMemory<T>& o_r,
                     deviceMemory<T>& o_w);
  void Update1NBFPCG(const T alpha,
                     deviceMemory<T>& o_p,
                     deviceMemory<T>& o_s,
                     deviceMemory<T>& o_q,
                     deviceMemory<T>& o_z,
                     deviceMemory<T>& o_x,
                     deviceMemory<T>& o_r,
                     deviceMemory<T>& o_u,
                     deviceMemory<T>& o_w);

public:
  nbfpcg(dlong _N, dlong _Nhalo,
       platform_t& _platform, settings_t& _settings, comm_t _comm);

  int Solve(operator_t& linearOperator, operator_t& precon,
            deviceMemory<T>& o_x, deviceMemory<T>& o_rhs,
            const T tol, const int MAXIT, const int verbose);
};

template class nbfpcg<float>;
template class nbfpcg<double>;
  
} //namespace LinearSolver

} //namespace libp

#endif

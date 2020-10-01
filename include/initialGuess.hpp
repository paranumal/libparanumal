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

#ifndef INITIALGUESS_HPP
#define INITIALGUESS_HPP

#include "linearSolver.hpp"

// Abstract base class for different initial guess strategies.
class initialGuessStrategy_t {
protected:
  platform_t& platform;
  settings_t& settings;
  MPI_Comm   comm;

  dlong Ntotal;     // Degrees of freedom

public:
  initialGuessStrategy_t(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm _comm, int _weighted, occa::memory& _o_weight);
  virtual ~initialGuessStrategy_t();

  virtual void FormInitialGuess(occa::memory& o_x, occa::memory& o_rhs) = 0;
  virtual void Update(solver_t& solver, occa::memory& o_x, occa::memory& o_rhs) = 0;
};

// Default initial guess strategy:  use whatever the user gave us.
class igDefaultStrategy : public initialGuessStrategy_t {
public:
  igDefaultStrategy(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm _comm, int _weighted, occa::memory& _o_weight);

  void FormInitialGuess(occa::memory& o_x, occa::memory& o_rhs);
  void Update(solver_t &solver, occa::memory& o_x, occa::memory& o_rhs);
};

// Zero initial guess strategy:  use a zero initial guess.
class igZeroStrategy : public initialGuessStrategy_t {
public:
  igZeroStrategy(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm _comm, int _weighted, occa::memory& _o_weight);

  void FormInitialGuess(occa::memory& o_x, occa::memory& o_rhs);
  void Update(solver_t &solver, occa::memory& o_x, occa::memory& o_rhs);
};

// Linear solver with successive-RHS initial-guess generation.
class initialGuessSolver_t : public linearSolver_t {
protected:
  initialGuessStrategy_t* igStrategy;   // The initial guess strategy.
  linearSolver_t*         linearSolver; // The linearSolver_t that does the solve.

public:
  initialGuessSolver_t(dlong _N, dlong _Nhalo, platform_t& _platform, settings_t& _settings, MPI_Comm _comm);
  ~initialGuessSolver_t();

  static initialGuessSolver_t* Setup(dlong _N, dlong _Nhalo,
                                     platform_t& platform, settings_t& settings, MPI_Comm _comm,
                                     int _weighted, occa::memory& _o_weight);

  int Solve(solver_t& solver, precon_t& precon,
            occa::memory& o_x, occa::memory& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);
};

#endif /* INITIALGUESS_HPP */

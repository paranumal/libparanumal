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

#include "initialGuess.hpp"

initialGuessSolver_t* initialGuessSolver_t::Setup(dlong N, dlong Nhalo, platform_t& platform, settings_t& settings, MPI_Comm comm, int weighted, occa::memory& o_weight)
{
  initialGuessSolver_t* initialGuessSolver = new initialGuessSolver_t(N, Nhalo, platform, settings, comm);
  initialGuessSolver->linearSolver = linearSolver_t::Setup(N, Nhalo, platform, settings, comm, weighted, o_weight);
  initialGuessSolver->igStrategy = nullptr;

  if (settings.compareSetting("INITIAL GUESS STRATEGY", "NONE")) {
    initialGuessSolver->igStrategy = new igDefaultStrategy(N, platform, settings, comm, weighted, o_weight);
  } else {
    LIBP_ABORT("Requested INITIAL GUESS STRATEGY not found.");
  }

  return initialGuessSolver;
}

initialGuessSolver_t::initialGuessSolver_t(dlong _N, dlong _Nhalo, platform_t& _platform, settings_t& _settings, MPI_Comm _comm):
  linearSolver_t(_N, _Nhalo, _platform, _settings, _comm),
  igStrategy(nullptr),
  linearSolver(nullptr)
{
  return;
}

initialGuessSolver_t::~initialGuessSolver_t()
{
  delete igStrategy;
  delete linearSolver;
}

int initialGuessSolver_t::Solve(solver_t& solver, precon_t& precon, occa::memory& o_x, occa::memory& o_rhs, const dfloat tol, const int MAXIT, const int verbose)
{
  int iter = 0;

  igStrategy->FormInitialGuess(o_x, o_rhs);
  iter = linearSolver->Solve(solver, precon, o_x, o_rhs, tol, MAXIT, verbose);
  igStrategy->Update(solver, o_x, o_rhs);

  return iter;
}

/*****************************************************************************/

initialGuessStrategy_t::initialGuessStrategy_t(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm& _comm, int _weighted, occa::memory& _o_weight):
  platform(_platform),
  settings(_settings),
  comm(_comm),
  Ntotal(_N)
{
  return;
}

initialGuessStrategy_t::~initialGuessStrategy_t()
{
  return;
}

/*****************************************************************************/

igDefaultStrategy::igDefaultStrategy(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm& _comm, int _weighted, occa::memory& _o_weight):
  initialGuessStrategy_t(_N, _platform, _settings, _comm, _weighted, _o_weight)
{
  return;
}

void igDefaultStrategy::FormInitialGuess(occa::memory& o_x, occa::memory& o_rhs)
{
  return;
}

void igDefaultStrategy::Update(solver_t &solver, occa::memory& o_x, occa::memory& o_rhs)
{
  return;
}

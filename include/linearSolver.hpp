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
#include "linAlg.hpp"

// #include "../solvers/elliptic/elliptic.hpp"
class elliptic_t: public solver_t {
public:

  linAlg_t* linAlg;

  occa::memory o_invDegree;

  void Preconditioner(occa::memory &o_r, occa::memory &o_z);
  void Operator(occa::memory &o_x, occa::memory &o_Ax);

};

//virtual base linear solver class
class linearSolver_t {
public:
  elliptic_t& elliptic; //must be made with an elliptic_t solver for now

  MPI_Comm& comm;
  occa::device& device;
  settings_t& settings;
  occa::properties& props;

  dlong N;

  linearSolver_t(elliptic_t& _elliptic);

  static linearSolver_t* Setup(elliptic_t& elliptic);

  virtual void Init()=0;
  virtual int Solve(occa::memory& o_x, occa::memory& o_rhs,
                    const dfloat tol, const int MAXIT, const int verbose)=0;
};

class pcg: public linearSolver_t {
private:
  occa::memory o_p, o_Ap, o_z, o_Ax, o_w;

  dfloat* tmprdotr;
  occa::memory h_tmprdotr;
  occa::memory o_tmprdotr;

  int flexible, weighted;

  occa::kernel updatePCGKernel;
  occa::kernel weightedUpdatePCGKernel;

public:
  pcg(elliptic_t& elliptic);
  ~pcg();

  void Init();
  int Solve(occa::memory& o_x, occa::memory& o_rhs,
            const dfloat tol, const int MAXIT, const int verbose);

  dfloat UpdatePCG(const dfloat alpha, occa::memory &o_x, occa::memory &o_r);
};

#endif

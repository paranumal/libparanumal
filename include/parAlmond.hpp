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

#ifndef PARALMOND_HPP
#define PARALMOND_HPP

#include "core.hpp"
#include "settings.hpp"
#include "platform.hpp"
#include "ogs.hpp"
#include "solver.hpp"
#include "precon.hpp"
#include "linearSolver.hpp"

namespace parAlmond {

#define PARALMOND_MAX_LEVELS 100

typedef enum {VCYCLE=0,KCYCLE=1,EXACT=3} CycleType;
typedef enum {PCG=0,GMRES=1} KrylovType;
typedef enum {DAMPED_JACOBI=0,CHEBYSHEV=1} SmoothType;
typedef enum {RUGESTUBEN=0,SYMMETRIC=1} StrengthType;

void AddSettings(settings_t& settings, const string prefix="");
void ReportSettings(settings_t& settings);

extern MPI_Datatype MPI_NONZERO_T;

//distributed matrix class passed to AMG setup
class parCOO {
public:
  platform_t &platform;
  MPI_Comm comm;

  dlong nnz=0;
  hlong *globalRowStarts=nullptr;
  hlong *globalColStarts=nullptr;

  //non-zero matrix entries
  struct nonZero_t {
    hlong row;
    hlong col;
    dfloat val;
  };
  nonZero_t *entries=nullptr;

  parCOO(platform_t &_platform, MPI_Comm _comm):
    platform(_platform), comm(_comm) {};

  ~parCOO() {
    if(entries) free(entries);
    if(globalRowStarts) free(globalRowStarts);
    if(globalColStarts) free(globalColStarts);
  }
};

//abstract multigrid level
// Class is derived from solver, and must have Operator defined
class multigridLevel: public solver_t  {
public:
  dlong Nrows=0, Ncols=0;

  //switch for weighted inner products
  bool weighted=false;
  dfloat *weight=nullptr;
  occa::memory o_weight;

  occa::memory o_scratch;

  multigridLevel(dlong N, dlong M, platform_t& _platform,
                 settings_t& _settings):
    solver_t(_platform, _settings), Nrows(N), Ncols(M) {}
  virtual ~multigridLevel() {};

  virtual void smooth(occa::memory& o_rhs, occa::memory& o_x, bool x_is_zero)=0;
  virtual void residual(occa::memory& o_rhs, occa::memory& o_x, occa::memory& o_res)=0;
  virtual void coarsen(occa::memory& o_x, occa::memory& o_Cx)=0;
  virtual void prolongate(occa::memory& o_x, occa::memory& o_Px)=0;
  virtual void Report()=0;
};

//forward declaration
//multigrid preconditioner
class multigrid_t;

class parAlmond_t: public precon_t {
public:
  parAlmond_t(platform_t& _platform, settings_t& settings_, MPI_Comm comm);
  ~parAlmond_t();

  //Add level to multigrid heirarchy
  void AddLevel(multigridLevel* level);

  // Setup AMG
  //-- Local A matrix data must be globally indexed & row sorted
  void AMGSetup(parCOO& A,
               bool nullSpace,
               dfloat *nullVector,
               dfloat nullSpacePenalty);

  void Operator(occa::memory& o_rhs, occa::memory& o_x);

  void Report();

  dlong getNumCols(int k);
  dlong getNumRows(int k);
private:
  platform_t& platform;
  settings_t& settings;

  multigrid_t *multigrid=nullptr;
};

} //namespace parAlmond

#endif

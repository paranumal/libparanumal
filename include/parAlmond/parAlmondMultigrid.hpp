/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#ifndef PARALMOND_MULTIGRID_HPP
#define PARALMOND_MULTIGRID_HPP

#include "settings.hpp"
#include "platform.hpp"
#include "solver.hpp"
#include "precon.hpp"
#include "parAlmond.hpp"
#include "parAlmond/parAlmondDefines.hpp"
#include "parAlmond/parAlmondCoarseSolver.hpp"

namespace parAlmond {

#define PARALMOND_MAX_LEVELS 100

typedef enum {VCYCLE=0,KCYCLE=1,EXACT=3} CycleType;
typedef enum {SMOOTHED=0,UNSMOOTHED=1} AggType;
typedef enum {PCG=0,GMRES=1} KrylovType;
typedef enum {DAMPED_JACOBI=0,CHEBYSHEV=1} SmoothType;
typedef enum {RUGESTUBEN=0,SYMMETRIC=1} StrengthType;

//multigrid preconditioner
class multigrid_t: public precon_t {
public:
  platform_t& platform;
  settings_t& settings;
  MPI_Comm comm;

  bool exact;
  linearSolver_t *linearSolver=nullptr;

  CycleType ctype;
  AggType aggtype;
  StrengthType strtype;

  int numLevels=0;
  int baseLevel=0;
  multigridLevel *levels[PARALMOND_MAX_LEVELS];

  occa::memory o_rhs[PARALMOND_MAX_LEVELS];
  occa::memory o_x[PARALMOND_MAX_LEVELS];

  coarseSolver_t *coarseSolver;

  //scratch space for smoothing and temporary residual vector
  size_t scratchSpaceBytes=0;
  occa::memory o_scratch;

  KrylovType ktype;

  occa::memory o_ck[PARALMOND_MAX_LEVELS];
  occa::memory o_vk[PARALMOND_MAX_LEVELS];
  occa::memory o_wk[PARALMOND_MAX_LEVELS];

  //scratch space
  size_t reductionScratchBytes=0;
  void *reductionScratch=nullptr;
  occa::memory h_reductionScratch;
  occa::memory o_reductionScratch;

  multigrid_t(platform_t& _platform, settings_t& _settings, MPI_Comm _comm);
  ~multigrid_t();

  void AddLevel(multigridLevel* level);

  void Operator(occa::memory& o_RHS, occa::memory& o_X);

  void vcycle(const int k, occa::memory& o_RHS, occa::memory& o_X);
  void kcycle(const int k, occa::memory& o_RHS, occa::memory& o_X);

private:
  void kcycleOp1(multigridLevel* level,
                 occa::memory& o_X,  occa::memory& o_RHS,
                 occa::memory& o_CK, occa::memory& o_VK,
                 dfloat *alpha1, dfloat *rho1,
                 dfloat *norm_rhs, dfloat *norm_rhstilde);

  void kcycleOp2(multigridLevel* level,
                occa::memory& o_X,  occa::memory& o_RHS,
                occa::memory& o_CK, occa::memory& o_VK, occa::memory& o_WK,
                const dfloat alpha1, const dfloat rho1);

  void kcycleCombinedOp1(multigridLevel* level, dfloat *aDotbc, occa::memory& o_a,
                       occa::memory& o_b, occa::memory& o_c);
  void kcycleCombinedOp2(multigridLevel* level, dfloat *aDotbcd,
                        occa::memory& o_a, occa::memory& o_b,
                        occa::memory& o_c, occa::memory& o_d);
  dfloat vectorAddInnerProd(multigridLevel* level,
                          const dfloat alpha, occa::memory& o_x,
                          const dfloat beta,  occa::memory& o_y);
};

}

#endif
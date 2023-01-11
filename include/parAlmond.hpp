/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

namespace libp {

namespace parAlmond {

void AddSettings(settings_t& settings, const std::string prefix="");
void ReportSettings(settings_t& settings);

//distributed matrix class passed to AMG setup
class parCOO {
public:
  platform_t platform;
  comm_t comm;

  dlong nnz=0;
  memory<hlong> globalRowStarts;
  memory<hlong> globalColStarts;

  //non-zero matrix entries
  struct nonZero_t {
    hlong row;
    hlong col;
    pfloat val;
  };
  memory<nonZero_t> entries;

  parCOO() = default;
  parCOO(platform_t &_platform, comm_t _comm):
    platform(_platform), comm(_comm) {};
};

//abstract multigrid level
// Class is derived from solver, and must have Operator defined
class multigridLevel: public operator_t  {
public:
  platform_t platform;
  settings_t settings;
  comm_t comm;

  dlong Nrows=0, Ncols=0;

  multigridLevel() = default;
  multigridLevel(dlong N, dlong M, platform_t& _platform,
                 settings_t& _settings, comm_t _comm):
    platform(_platform), settings(_settings),
    comm(_comm), Nrows(N), Ncols(M) {}

  virtual void smooth(deviceMemory<pfloat>& o_rhs, deviceMemory<pfloat>& o_x, bool x_is_zero)=0;
  virtual void residual(deviceMemory<pfloat>& o_rhs, deviceMemory<pfloat>& o_x, deviceMemory<pfloat>& o_res)=0;
  virtual void coarsen(deviceMemory<pfloat>& o_x, deviceMemory<pfloat>& o_Cx)=0;
  virtual void prolongate(deviceMemory<pfloat>& o_x, deviceMemory<pfloat>& o_Px)=0;
  virtual void Report()=0;
  virtual size_t SmootherScratchSize()=0;
};

typedef enum {VCYCLE=0,KCYCLE=1,EXACT=3} CycleType;
typedef enum {SMOOTHED=0,UNSMOOTHED=1} AggType;
typedef enum {PCG=0,GMRES=1} KrylovType;
typedef enum {DAMPED_JACOBI=0,CHEBYSHEV=1} SmoothType;
typedef enum {RUGESTUBEN=0,SYMMETRIC=1} StrengthType;
typedef enum {COARSEEXACT=0,COARSEOAS=1} CoarseType;

class coarseSolver_t;

//multigrid preconditioner
class multigrid_t: public operator_t {
public:
  platform_t platform;
  settings_t settings;
  comm_t comm;

  bool exact=false;
  linearSolver_t linearSolver;

  CycleType ctype;
  AggType aggtype;
  StrengthType strtype;
  CoarseType coarsetype;

  int numLevels=0;
  int baseLevel=0;
  static constexpr int PARALMOND_MAX_LEVELS=100;
  std::shared_ptr<multigridLevel> levels[PARALMOND_MAX_LEVELS];

  std::shared_ptr<coarseSolver_t> coarseSolver;

  size_t scratchEntries=0;

  KrylovType ktype;

  multigrid_t() = default;
  multigrid_t(platform_t& _platform, settings_t& _settings,
              comm_t _comm);

  template<class Level, class... Args>
  Level& AddLevel(Args&& ... args) {
    levels[numLevels++] = std::make_shared<Level>(args...);
    return dynamic_cast<Level&>(*levels[numLevels-1]);
  }
  template<class Level>
  Level& GetLevel(const int l) {
    return dynamic_cast<Level&>(*levels[l]);
  }

  void EstimateScratchSpace();
  size_t LevelScratchSpace(const int k);

  void Operator(deviceMemory<pfloat>& o_rhs, deviceMemory<pfloat>& o_x);

  void vcycle(const int k, deviceMemory<pfloat>& o_rhs, deviceMemory<pfloat>& o_x);
  void kcycle(const int k, deviceMemory<pfloat>& o_rhs, deviceMemory<pfloat>& o_x);

private:
  void kcycleOp1(multigridLevel& level,
                 deviceMemory<pfloat>& o_x,  deviceMemory<pfloat>& o_rhs,
                 deviceMemory<pfloat>& o_ck, deviceMemory<pfloat>& o_vk,
                 pfloat& alpha1, pfloat& rho1,
                 pfloat& norm_rhs, pfloat& norm_rhstilde);

  void kcycleOp2(multigridLevel& level,
                deviceMemory<pfloat>& o_x,  deviceMemory<pfloat>& o_rhs,
                deviceMemory<pfloat>& o_ck, deviceMemory<pfloat>& o_vk,
                deviceMemory<pfloat>& o_wk,
                const pfloat alpha1, const pfloat rho1);

  void kcycleCombinedOp1(multigridLevel& level,
                        deviceMemory<pfloat>& o_a,
                        deviceMemory<pfloat>& o_b,
                        deviceMemory<pfloat>& o_c,
                        pfloat& aDotb,
                        pfloat& aDotc,
                        pfloat& bDotb);
  void kcycleCombinedOp2(multigridLevel& level,
                        deviceMemory<pfloat>& o_a,
                        deviceMemory<pfloat>& o_b,
                        deviceMemory<pfloat>& o_c,
                        deviceMemory<pfloat>& o_d,
                        pfloat& aDotb,
                        pfloat& aDotc,
                        pfloat& aDotd);
  pfloat vectorAddInnerProd(multigridLevel& level,
                          const pfloat alpha, deviceMemory<pfloat>& o_x,
                          const pfloat beta,  deviceMemory<pfloat>& o_y);
};

class parAlmond_t: public operator_t {
public:
  parAlmond_t() = default;
  parAlmond_t(platform_t& _platform, settings_t& _settings, comm_t _comm) {
    Setup(_platform, _settings, _comm);
  }

  void Setup(platform_t& _platform, settings_t& _settings, comm_t _comm);

  template<class Level, class... Args>
  Level& AddLevel(Args&& ... args) {
    return multigrid->AddLevel<Level, Args...>(args...);
  }
  template<class Level>
  Level& GetLevel(const int l) {
    return multigrid->GetLevel<Level>(l);
  }

  int NumLevels();

  // Setup AMG
  //-- Local A matrix data must be globally indexed & row sorted
  void AMGSetup(parCOO& A,
               bool nullSpace,
               memory<pfloat> nullVector,
               pfloat nullSpacePenalty);

  void Operator(deviceMemory<pfloat>& o_rhs, deviceMemory<pfloat>& o_x);

  void Report();

  dlong getNumCols(int k);
  dlong getNumRows(int k);
private:
  platform_t platform;
  settings_t settings;

  std::shared_ptr<multigrid_t> multigrid=nullptr;
};

} //namespace parAlmond

} //namespace libp

#endif

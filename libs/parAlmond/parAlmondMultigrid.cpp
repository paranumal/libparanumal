/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"
#include "parAlmond/parAlmondDefines.hpp"
#include "parAlmond/parAlmondCoarseSolver.hpp"

namespace libp {

namespace parAlmond {

void multigrid_t::Operator(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_x) {

  /*Pre-reserve memory pool space to avoid some unnecessary re-sizing*/
  platform.reserve<dfloat>(scratchEntries);

  if (ctype == KCYCLE) {
    kcycle(0, o_rhs, o_x);
  } else {
    vcycle(0, o_rhs, o_x);
  }
}

multigrid_t::multigrid_t(platform_t& _platform, settings_t& _settings,
                         comm_t _comm):
    platform(_platform), settings(_settings), comm(_comm) {

  //determine what sort of multigrid cycle to construct
  if(settings.compareSetting("PARALMOND CYCLE", "KCYCLE")) {
    ctype = KCYCLE;
  } else {
    ctype = VCYCLE;
  }

  //strength type
  if(settings.compareSetting("PARALMOND STRENGTH", "SYMMETRIC")) {
    strtype = SYMMETRIC;
  } else {
    strtype = RUGESTUBEN;
  }

  //aggregation type
  if(settings.compareSetting("PARALMOND AGGREGATION", "UNSMOOTHED")) {
    aggtype = UNSMOOTHED;
  } else {
    aggtype = SMOOTHED;
  }

  if (settings.compareSetting("PARALMOND CYCLE", "NONSYM")) {
    ktype = GMRES;
  } else {
    ktype = PCG;
  }

  //determine whether parAlmond should do an exact solve when called
  if(settings.compareSetting("PARALMOND CYCLE", "EXACT"))
    exact = true;
  else
    exact = false;

  //Hard code for now
  coarsetype=COARSEEXACT;

  if (coarsetype==COARSEEXACT) {
    coarseSolver = std::make_shared<exactSolver_t>(_platform, _settings, _comm);
  } else {
    coarseSolver = std::make_shared<oasSolver_t>(_platform, _settings, _comm);
  }
}

void multigrid_t::EstimateScratchSpace() {
  scratchEntries = LevelScratchSpace(0);
}

size_t multigrid_t::LevelScratchSpace(const int k) {
  multigridLevel& level = *levels[k];

  // space for rhs and x
  size_t baseEntries = 0;
  if (k>0) {
    baseEntries += 2*level.Ncols + 2 * platform.memPoolAlignment<dfloat>();
  }

  // space for residual vector
  size_t residualEntries = 0;
  if (k<numLevels-1) {
    residualEntries = level.Ncols + platform.memPoolAlignment<dfloat>();
  }

  //Scratch space needed for smoothing
  size_t smootherEntries = 0;
  if (k<numLevels-1) {
    smootherEntries += level.SmootherScratchSize();
  }

  //residual can re-use space for smoothing
  smootherEntries = std::max(smootherEntries, residualEntries);

  size_t cycleEntries = 0;
  //extra stroage for kcycle vectors ck, vk, wk and reduction scratch
  if (ctype==KCYCLE && k>0 && k<NUMKCYCLES+1) {
    cycleEntries += level.Ncols + 2*level.Nrows + 3 * PARALMOND_NBLOCKS
                  + 4 * platform.memPoolAlignment<dfloat>();
  }

  if (k==numLevels-1) {
    //base level
    cycleEntries += coarseSolver->scratchSize();
  } else {
    //add usage from coarse levels
    cycleEntries += LevelScratchSpace(k+1);
  }

  return baseEntries + std::max(smootherEntries, cycleEntries);
}

} //namespace parAlmond

} //namespace libp

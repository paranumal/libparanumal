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

void multigrid_t::Operator(deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_X) {
  if (ctype == KCYCLE) {
    kcycle(0, o_RHS, o_X);
  } else {
    vcycle(0, o_RHS, o_X);
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

void multigrid_t::AllocateLevelWorkSpace(const int k){

  multigridLevel& level = *levels[k];

  //If using an exact solver and this is the first level, setup a linearSovler
  if (exact && k==0) {
    if (settings.compareSetting("PARALMOND CYCLE", "NONSYM"))
      linearSolver.Setup<LinearSolver::pgmres>(level.Nrows,
                                               level.Ncols - level.Nrows,
                                               platform, settings, comm);
    else
      linearSolver.Setup<LinearSolver::pcg>(level.Nrows,
                                            level.Ncols - level.Nrows,
                                            platform, settings, comm);
  }

  if (ctype==KCYCLE) {
    //first level
    if (NreductionScratch==0) {
      NreductionScratch = 3*PARALMOND_NBLOCKS;
      memory<dfloat> dummy(3*PARALMOND_NBLOCKS, 0.0);
      reductionScratch = platform.hostMalloc<dfloat>(NreductionScratch, dummy);
      o_reductionScratch = platform.malloc<dfloat>(NreductionScratch, dummy);
    }

    //extra stroage for kcycle vectors
    if (k>0 && k<NUMKCYCLES+1) {
      memory<dfloat> dummy(level.Ncols,0.0);
      o_ck[k] = platform.malloc<dfloat>(level.Ncols,dummy);
      o_vk[k] = platform.malloc<dfloat>(level.Nrows,dummy);
      o_wk[k] = platform.malloc<dfloat>(level.Nrows,dummy);
    }
  }

  //allocate space for coarse rhs and x
  if (k>0) {
    memory<dfloat> dummy(level.Ncols,0.0);
    o_x[k]   = platform.malloc<dfloat>(level.Ncols,dummy);
    o_rhs[k] = platform.malloc<dfloat>(level.Ncols,dummy);
  }

  //scratch space includes space for residual and 2 vectors used in Chebyshev smoothing
  size_t Nrequired = 2*level.Ncols;
  if (Nrequired>NscratchSpace) {
    NscratchSpace = Nrequired;
    memory<dfloat> dummy(2*level.Ncols,0.0);
    o_scratch = platform.malloc<dfloat>(Nrequired, dummy);
  }

  level.o_scratch = o_scratch;
}

} //namespace parAlmond

} //namespace libp

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

#include "parAlmond.hpp"
#include "parAlmond/parAlmondDefines.hpp"
#include "parAlmond/parAlmondMultigrid.hpp"
#include "parAlmond/parAlmondCoarseSolver.hpp"

namespace parAlmond {

void multigrid_t::Operator(occa::memory& o_RHS, occa::memory& o_X) {
  if (ctype == KCYCLE) {
    kcycle(0, o_RHS, o_X);
  } else {
    vcycle(0, o_RHS, o_X);
  }
}

multigrid_t::multigrid_t(platform_t& _platform, settings_t& _settings,
                         MPI_Comm _comm):
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
    coarseSolver = new exactSolver_t(_platform, _settings, _comm);
  } else {
    coarseSolver = new oasSolver_t(_platform, _settings, _comm);
  }
}

multigrid_t::~multigrid_t() {
  if (linearSolver) delete linearSolver;
  if (coarseSolver) delete coarseSolver;
  for (int n=0;n<numLevels;n++) delete levels[n];
}

void multigrid_t::AddLevel(multigridLevel* level){

  //If using an exact solver and this is the first level, setup a linearSovler
  if (exact && numLevels==0) {
    if (settings.compareSetting("PARALMOND CYCLE", "NONSYM"))
      linearSolver = new pgmres(level->Nrows, level->Ncols - level->Nrows,
                             platform, settings, comm);
    else
      linearSolver = new pcg(level->Nrows, level->Ncols - level->Nrows,
                             platform, settings, comm);
  }

  if (ctype==KCYCLE) {
    //first level
    if (reductionScratchBytes==0) {
      reductionScratchBytes = 3*PARALMOND_NBLOCKS*sizeof(dfloat);
      dfloat *dummy = (dfloat *) calloc(3*PARALMOND_NBLOCKS,sizeof(dfloat));
      o_reductionScratch = platform.malloc(reductionScratchBytes, dummy);
      reductionScratch = platform.hostMalloc(reductionScratchBytes, NULL, h_reductionScratch);
      free(dummy);
    }

    //extra stroage for kcycle vectors
    if (numLevels>0 && numLevels<NUMKCYCLES+1) {
      dfloat *dummy = (dfloat *) calloc(level->Ncols,sizeof(dfloat));
      o_ck[numLevels] = platform.malloc(level->Ncols*sizeof(dfloat),dummy);
      o_vk[numLevels] = platform.malloc(level->Nrows*sizeof(dfloat),dummy);
      o_wk[numLevels] = platform.malloc(level->Nrows*sizeof(dfloat),dummy);
      free(dummy);
    }
  }

  //allocate space for coarse rhs and x
  if (numLevels>0) {
    dfloat *dummy = (dfloat *) calloc(level->Ncols,sizeof(dfloat));
    o_x[numLevels]   = platform.malloc(level->Ncols*sizeof(dfloat),dummy);
    o_rhs[numLevels] = platform.malloc(level->Ncols*sizeof(dfloat),dummy);
    free(dummy);
  }

  //scratch space includes space for residual and 2 vectors used in Chebyshev smoothing
  size_t requiredBytes = 2*level->Ncols*sizeof(dfloat);
  if (requiredBytes>scratchSpaceBytes) {
    scratchSpaceBytes = requiredBytes;
    dfloat *dummy = (dfloat *) calloc(2*level->Ncols,sizeof(dfloat));
    o_scratch = platform.malloc(requiredBytes, dummy);
    free(dummy);
  }

  level->o_scratch = o_scratch;

  levels[numLevels++] = level;
}

} //namespace parAlmond

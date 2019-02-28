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

#include "adaptive.h"

int adaptiveSolve(adaptive_t *adaptive, dfloat lambda, dfloat tol,
                  occa::memory &o_b, occa::memory &o_x){
  
  setupAide options = adaptive->options;

  int Niter = 0;
  int maxIter = 1000; 

  level_t *level = adaptive->lvl;

  options.getArgs("MAXIMUM ITERATIONS", maxIter);
  options.getArgs("SOLVER TOLERANCE", tol);

  // (Snc*S*G*Gnc)*AL*xL = Snc*S*G*Gnc*ML*fL
  // xL = Snc*S*xg
#if USE_GASPAR==1
  adaptiveGatherScatter(level, o_b);
#else
  // gather over noncon faces to coarse side dofs
  level->gather_noncon(level->Klocal, level->o_EToC, level->o_Ib, level->o_It, o_b);
  
  // add noncon gs around this
  ogsGatherScatter(o_b, ogsDfloat, ogsAdd, level->ogs);
#endif

#if USE_NULL_PROJECTION==1
  if(adaptive->allNeumann) // zero mean of RHS
    adaptiveZeroMean(adaptive, adaptive->lvl, o_b);
#endif

  Niter = pcg (adaptive, lambda, o_b, o_x, tol, maxIter);
  
#if USE_NULL_PROJECTION==1
  if(adaptive->allNeumann) // zero mean of RHS
    adaptiveZeroMean(adaptive, adaptive->lvl, o_x);
#endif

#if USE_GASPAR==0
  level->scatter_noncon(level->Klocal, level->o_EToC, level->o_Ib, level->o_It, o_x);
#endif
  
  //  adaptiveRankOneProjection(level, level->o_filtU, level->o_filtV, o_x, o_x);

#if 0
  adaptive->dotMultiplyKernel(level->Np*level->Klocal, level->o_invDegree, o_x, o_x);

  adaptiveGatherScatter(level, o_x);
#endif
  
  return Niter;
}

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

#include "elliptic.hpp"

int elliptic_t::Solve(occa::memory &o_x, occa::memory &o_r,
                      const dfloat tol, const int MAXIT, const int verbose){


#if USE_NULL_PROJECTION==1
  if(allNeumann) // zero mean of RHS
    ZeroMean(o_r);
#endif

  int Niter = linearSolver->Solve(o_x, o_r, tol, MAXIT, verbose);

  // if(!options.compareArgs("KRYLOV SOLVER", "NONBLOCKING"))
  //   Niter = pcg (elliptic, lambda, o_r, o_x, tol, maxIter, verbose);
  // else{
  //   if(!options.compareArgs("KRYLOV SOLVER", "FLEXIBLE")){
  //     Niter = nbpcg (elliptic, lambda, o_r, o_x, tol, maxIter);
  //   }
  //   else{
  //     Niter = nbfpcg (elliptic, lambda, o_r, o_x, tol, maxIter);
  //   }
  // }

#if USE_NULL_PROJECTION==1
  if(allNeumann) // zero mean of RHS
    ZeroMean(o_x);
#endif

  return Niter;
}

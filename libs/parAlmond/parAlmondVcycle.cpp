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
#include "parAlmond/parAlmondCoarseSolver.hpp"

namespace libp {

namespace parAlmond {

void multigrid_t::vcycle(const int k, deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_X){

  //check for base level
  if(k==baseLevel) {
    coarseSolver->solve(o_RHS, o_X);
    return;
  }

  multigridLevel& level = *levels[k];
  deviceMemory<dfloat>& o_RHSC = o_rhs[k+1];
  deviceMemory<dfloat>& o_XC   = o_x[k+1];
  deviceMemory<dfloat>& o_RES  = o_scratch;

  //apply smoother to x and then compute res = rhs-Ax
  level.smooth(o_RHS, o_X, true);
  level.residual(o_RHS, o_X, o_RES);

  // rhsC = P^T res
  level.coarsen(o_RES, o_RHSC);

  vcycle(k+1, o_RHSC, o_XC);

  // x = x + P xC
  level.prolongate(o_XC, o_X);

  level.smooth(o_RHS, o_X, false);
}

} //namespace parAlmond

} //namespace libp

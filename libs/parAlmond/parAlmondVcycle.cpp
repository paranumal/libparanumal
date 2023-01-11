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

void multigrid_t::vcycle(const int k, deviceMemory<pfloat>& o_rhs, deviceMemory<pfloat>& o_x){

  //check for base level
  if(k==baseLevel) {
    coarseSolver->solve(o_rhs, o_x);
    return;
  }

  multigridLevel& level = *levels[k];
  multigridLevel& levelC = *levels[k+1];

  //apply smoother to x and then compute res = rhs-Ax
  level.smooth(o_rhs, o_x, true);

  deviceMemory<pfloat> o_rhsC = platform.reserve<pfloat>(levelC.Ncols);
  deviceMemory<pfloat> o_xC   = platform.reserve<pfloat>(levelC.Ncols);
  deviceMemory<pfloat> o_res  = platform.reserve<pfloat>(level.Ncols);

  level.residual(o_rhs, o_x, o_res);

  // rhsC = P^T res
  level.coarsen(o_res, o_rhsC);
  o_res.free();

  vcycle(k+1, o_rhsC, o_xC);

  // x = x + P xC
  level.prolongate(o_xC, o_x);
  o_xC.free(); o_rhsC.free();

  level.smooth(o_rhs, o_x, false);
}

} //namespace parAlmond

} //namespace libp

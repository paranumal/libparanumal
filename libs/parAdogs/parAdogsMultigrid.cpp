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

#include "parAdogs.hpp"
#include "parAdogs/parAdogsGraph.hpp"
#include "parAdogs/parAdogsMultigrid.hpp"

namespace libp {

namespace paradogs {

/****************************************/
/* Multigrid vcycle                     */
/****************************************/
void graph_t::MultigridVcycle(const int l,
                              memory<dfloat>& r,
                              memory<dfloat>& x) {

  //check for base level
  if(l==Nlevels-1) {
    coarseSolver.Solve(r, x);
    return;
  }

  mgLevel_t& Lf = L[l];
  memory<dfloat>& res = Lf.RES;

  mgLevel_t& Lc = L[l+1];
  memory<dfloat>& rC = Lc.RHS;
  memory<dfloat>& xC = Lc.X;

  //Pre smooth and then compute res = rhs-Ax
  Lf.Smooth(r, x, true);
  Lf.Residual(r, x, res);

  // rhsC = P^T res
  Lf.Coarsen(res, rC);

  // Recursive call
  MultigridVcycle(l+1, rC, xC);
  // for (int n=0;n<Lc.Nrows;++n) xC[n] = rC[n];

  // x = x + P xC
  Lf.Prolongate(xC, x);

  // Post smooth
  Lf.Smooth(r, x, false);
}

void mgLevel_t::Residual(memory<dfloat>& r, memory<dfloat>& x, memory<dfloat>& res) {
  A.SpMV(-1.0, x, 1.0, r, res);
}

void mgLevel_t::Coarsen(memory<dfloat>& x, memory<dfloat>& xC) {
  R.SpMV(1.0, x, 0.0, xC);
}

void mgLevel_t::Prolongate(memory<dfloat>& xC, memory<dfloat>& x) {
  P.SpMV(1.0, xC, 1.0, x);
}

void mgLevel_t::Smooth(memory<dfloat>& r, memory<dfloat>& x, const bool xIsZero) {
  const int ChebyshevIterations=2;
  A.SmoothChebyshev(r, x, lambda0, lambda1,
                     xIsZero, scratch,
                     ChebyshevIterations);
}

} //namespace paradogs

} //namespace libp

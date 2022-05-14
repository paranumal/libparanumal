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
#include "parAlmond/parAlmondAMGLevel.hpp"
#include "parAlmond/parAlmondKernels.hpp"

namespace libp {

namespace parAlmond {

void parCSR::smoothDampedJacobi(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_x,
                                const dfloat lambda, bool x_is_zero,
                                deviceMemory<dfloat>& o_scratch){

  if(x_is_zero){
    // x = lambda*inv(D)*r
    platform.linAlg().amxpy(Nrows, lambda, o_diagInv, o_r, 0.0, o_x);
    return;
  }

  deviceMemory<dfloat> o_d = o_scratch;

  halo.ExchangeStart(o_x, 1);

  // d = lambda*inv(D)*(r-A*x)
  if (diag.NrowBlocks)
    SmoothJacobiCSRKernel(diag.NrowBlocks,
                         diag.o_blockRowStarts, diag.o_rowStarts,
                         diag.o_cols, diag.o_vals,
                         lambda, o_diagInv,
                         o_r, o_x, o_d);

  halo.ExchangeFinish(o_x, 1);

  if (offd.NrowBlocks)
    SmoothJacobiMCSRKernel(offd.NrowBlocks,
                           offd.o_blockRowStarts, offd.o_mRowStarts,
                           offd.o_rows, offd.o_cols, offd.o_vals,
                           lambda, o_diagInv, o_x, o_d);

  platform.linAlg().axpy(Nrows, 1.0, o_d, 1.0, o_x);
}

void parCSR::smoothChebyshev(deviceMemory<dfloat>& o_b, deviceMemory<dfloat>& o_x,
                             const dfloat lambda0, const dfloat lambda1,
                             bool x_is_zero, deviceMemory<dfloat>& o_scratch,
                             const int ChebyshevIterations) {

  const dfloat theta = 0.5*(lambda1+lambda0);
  const dfloat delta = 0.5*(lambda1-lambda0);
  const dfloat invTheta = 1.0/theta;
  const dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  deviceMemory<dfloat> o_d = o_scratch + 0*Ncols;
  deviceMemory<dfloat> o_r = o_scratch + 1*Ncols;


  if(x_is_zero){ //skip the Ax if x is zero
    //r = D^{-1}b
    //d = invTheta*r
    //x = d
    if (Nrows)
      SmoothChebyshevStartKernel(Nrows, invTheta, o_diagInv,
                                 o_b, o_r, o_d, o_x);
  } else {
    //r = D^{-1}(b-A*x)
    halo.ExchangeStart(o_x, 1);

    const dfloat alpha = 0.0;
    const dfloat beta = 1.0;

    if (diag.NrowBlocks)
      SmoothChebyshevCSRKernel(diag.NrowBlocks,
                               diag.o_blockRowStarts, diag.o_rowStarts,
                               diag.o_cols, diag.o_vals,
                               alpha, beta, o_diagInv,
                               o_b, o_x, o_r);

    halo.ExchangeFinish(o_x, 1);

    if (offd.NrowBlocks)
      SmoothChebyshevMCSRKernel(offd.NrowBlocks,
                               offd.o_blockRowStarts, offd.o_mRowStarts,
                               offd.o_rows, offd.o_cols, offd.o_vals,
                               o_diagInv, o_x, o_r);

    const int last_it = (ChebyshevIterations==0) ? 1 : 0;

    //d = invTheta*r
    //x = x + d
    if (Nrows)
      SmoothChebyshevUpdateKernel(Nrows, alpha, invTheta,
                                  last_it, o_r, o_d, o_x);
  }

  for (int k=0;k<ChebyshevIterations;k++) {

    const dfloat alpha = 1.0;
    const dfloat beta = 0.0;

    //r_k+1 = r_k - D^{-1}Ad_k
    halo.ExchangeStart(o_d, 1);

    if (diag.NrowBlocks)
      SmoothChebyshevCSRKernel(diag.NrowBlocks,
                               diag.o_blockRowStarts, diag.o_rowStarts,
                               diag.o_cols, diag.o_vals,
                               alpha, beta, o_diagInv,
                               o_b, o_d, o_r);

    halo.ExchangeFinish(o_d, 1);

    if (offd.NrowBlocks)
      SmoothChebyshevMCSRKernel(offd.NrowBlocks,
                               offd.o_blockRowStarts, offd.o_mRowStarts,
                               offd.o_rows, offd.o_cols, offd.o_vals,
                               o_diagInv, o_d, o_r);

    const int last_it = (k==ChebyshevIterations-1) ? 1 : 0;

    rho_np1 = 1.0/(2.*sigma-rho_n);

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    //x_k+1 = x_k + d_k+1
    if (Nrows)
      SmoothChebyshevUpdateKernel(Nrows,
                                  rho_np1*rho_n,
                                  dfloat(2.0)*rho_np1/delta,
                                  last_it,
                                  o_r, o_d, o_x);

    rho_n = rho_np1;
  }
}

} //namespace parAlmond

} //namespace libp

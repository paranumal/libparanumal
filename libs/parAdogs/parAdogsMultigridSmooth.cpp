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
#include "parAdogs/parAdogsMatrix.hpp"
#include "parAdogs/parAdogsPartition.hpp"

namespace libp {

namespace paradogs {

void parCSR::SmoothChebyshev(memory<dfloat>& b, memory<dfloat>& x,
                             const dfloat lambda0, const dfloat lambda1,
                             const bool xIsZero, memory<dfloat>& scratch,
                             const int ChebyshevIterations) {

  const dfloat theta = 0.5*(lambda1+lambda0);
  const dfloat delta = 0.5*(lambda1-lambda0);
  const dfloat invTheta = 1.0/theta;
  const dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  memory<dfloat> d = scratch + 0*Ncols;
  memory<dfloat> r = scratch + 1*Ncols;

  if(xIsZero){ //skip the Ax if x is zero
    // r = D^{-1}b
    // d = invTheta*r
    // x = d
    #pragma omp parallel for
    for (dlong n=0;n<Nrows;++n) {
      const dfloat r_r = diagInv[n]*b[n];
      r[n] = r_r;
      d[n] = invTheta*r_r;
      x[n] = invTheta*r_r;
    }

  } else {
    halo.ExchangeStart(x, 1);

    //r = D^{-1}(b-A*x)
    #pragma omp parallel for
    for (dlong n=0;n<Nrows;++n) {

      dfloat rn = b[n];

      const dlong start = diag.rowStarts[n];
      const dlong end   = diag.rowStarts[n+1];
      for (dlong j=start;j<end;++j) {
        rn -= diag.vals[j]*x[diag.cols[j]];
      }

      r[n] = diagInv[n]*rn;
    }

    halo.ExchangeFinish(x, 1);

    #pragma omp parallel for
    for(dlong n=0; n<offd.nzRows; n++){ //local

      dfloat rn = 0.0;

      const dlong row = offd.rows[n];
      const dlong start = offd.mRowStarts[n];
      const dlong end   = offd.mRowStarts[n+1];
      for(dlong j=start; j<end; ++j) {
        rn -= offd.vals[j]*x[offd.cols[j]];
      }

      r[row] += diagInv[row]*rn;
    }

    const int last_it = (ChebyshevIterations==0) ? 1 : 0;

    //d = invTheta*r
    //x = x + d
    if (last_it) {
      #pragma omp parallel for
      for (dlong n=0;n<Nrows;++n) {
        x[n] += invTheta*r[n];
      }
    } else {
      #pragma omp parallel for
      for (dlong n=0;n<Nrows;++n) {
        d[n] = invTheta*r[n];
        x[n] += d[n];
      }
    }
  }

  for (int k=0;k<ChebyshevIterations;k++) {

    halo.ExchangeStart(d, 1);

    //r_k+1 = r_k - D^{-1}Ad_k
    #pragma omp parallel for
    for (dlong n=0;n<Nrows;++n) {

      dfloat rn = 0.0;

      const dlong start = diag.rowStarts[n];
      const dlong end   = diag.rowStarts[n+1];
      for (dlong j=start;j<end;++j) {
        rn -= diag.vals[j]*d[diag.cols[j]];
      }

      r[n] += diagInv[n]*rn;
    }

    halo.ExchangeFinish(d, 1);

    #pragma omp parallel for
    for(dlong n=0; n<offd.nzRows; n++){ //local

      dfloat rn = 0.0;

      const dlong row = offd.rows[n];
      const dlong start = offd.mRowStarts[n];
      const dlong end   = offd.mRowStarts[n+1];
      for(dlong j=start; j<end; ++j) {
        rn -= offd.vals[j]*d[offd.cols[j]];
      }

      r[row] += diagInv[row]*rn;
    }

    const int last_it = (k==ChebyshevIterations-1) ? 1 : 0;

    rho_np1 = 1.0/(2.*sigma-rho_n);

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    //x_k+1 = x_k + d_k+1
    if (last_it) {
      #pragma omp parallel for
      for (dlong n=0;n<Nrows;++n) {
        const dfloat d_np1 = (rho_np1*rho_n)*d[n] + (2.0*rho_np1/delta)*r[n];
        x[n] += d_np1;
      }
    } else {
      #pragma omp parallel for
      for (dlong n=0;n<Nrows;++n) {
        d[n] = (rho_np1*rho_n)*d[n] + (2.0*rho_np1/delta)*r[n];
        x[n] += d[n];
      }
    }

    rho_n = rho_np1;
  }
}

} //namespace paradogs

} //namespace libp

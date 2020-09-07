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
#include "parAlmond/parAlmondAMGLevel.hpp"

namespace parAlmond {

void amgLevel::smoothDampedJacobi(occa::memory& o_r, occa::memory& o_X,
                                   bool x_is_zero){

  if(x_is_zero){
    platform.linAlg.amxpy(Nrows, lambda, A->o_diagInv, o_r, 0.0, o_X);
    return;
  }

  occa::memory& o_RES = o_scratch;

  // RES = r-A*x
  A->SpMV(-1.0, o_X, 1.0, o_r, o_RES);

  // x = x + alpha*inv(D)*RES
  platform.linAlg.amxpy(Nrows, lambda, A->o_diagInv, o_RES, 1.0, o_X);
}

void amgLevel::smoothChebyshev(occa::memory& o_r, occa::memory& o_X,
                                bool x_is_zero) {

  const dfloat theta = 0.5*(lambda1+lambda0);
  const dfloat delta = 0.5*(lambda1-lambda0);
  const dfloat invTheta = 1.0/theta;
  const dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  occa::memory o_RES = o_scratch + 0*Ncols*sizeof(dfloat);
  occa::memory o_Ad  = o_scratch + 1*Ncols*sizeof(dfloat);
  occa::memory o_d   = o_scratch + 2*Ncols*sizeof(dfloat);


  if(x_is_zero){ //skip the Ax if x is zero
    //RES = D^{-1}r
    platform.linAlg.amxpy(Nrows, 1.0, A->o_diagInv, o_r, 0.0, o_RES);
    platform.linAlg.set(Nrows, 0.0, o_X);
    //d = invTheta*RES
    platform.linAlg.axpy(Nrows, invTheta, o_RES, 0.0, o_d);
  } else {
    //RES = D^{-1}(r-Ax)
    A->SpMV(-1.0, o_X, 1.0, o_r, o_RES);
    platform.linAlg.amx(Nrows, 1.0, A->o_diagInv, o_RES);

    //d = invTheta*RES
    platform.linAlg.axpy(Nrows, invTheta, o_RES, 0.0, o_d);
  }

  for (int k=0;k<ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    platform.linAlg.axpy(Nrows, 1.0, o_d, 1.0, o_X);

    //r_k+1 = r_k - D^{-1}Ad_k
    A->SpMV(1.0, o_d, 0.0, o_Ad);
    platform.linAlg.amxpy(Nrows, -1.0, A->o_diagInv, o_Ad, 1.0, o_RES);

    rho_np1 = 1.0/(2.*sigma-rho_n);

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    platform.linAlg.axpy(Nrows, 2.0*rho_np1/delta, o_RES, rho_np1*rho_n, o_d);
    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  platform.linAlg.axpy(Nrows, 1.0, o_d, 1.0, o_X);
}

} //namespace parAlmond

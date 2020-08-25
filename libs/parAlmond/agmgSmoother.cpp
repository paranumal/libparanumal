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

namespace parAlmond {

void agmgLevel::smoothJacobi(dfloat *r, dfloat *X,
                             const bool x_is_zero) {

  // X = X + inv(D)*(b-A*X)
  if(x_is_zero){
    vectorDotStar(Nrows,1.0,A->diagInv,r,0.0,X);
    return;
  }

  dfloat *RES = (dfloat *) scratch;

  A->SpMV(-1.0, X, 1.0, r, RES);
  vectorDotStar(Nrows, 1.0, A->diagInv, RES, 1.0, X);
}


void agmgLevel::smoothDampedJacobi(dfloat *r, dfloat *X,
                                   const bool x_is_zero) {

  // X = X + alpha*inv(D)*(b-A*X)
  if(x_is_zero){
    vectorDotStar(Nrows,lambda,A->diagInv,r,0.0,X);
    return;
  }

  dfloat *RES = (dfloat *) scratch;

  A->SpMV(-1.0, X, 1.0, r, RES);
  vectorDotStar(Nrows, lambda, A->diagInv, RES, 1.0, X);
}

void agmgLevel::smoothChebyshev(dfloat *r, dfloat *X,
                                const bool x_is_zero) {

  const dfloat theta = 0.5*(lambda1+lambda0);
  const dfloat delta = 0.5*(lambda1-lambda0);
  const dfloat invTheta = 1.0/theta;
  const dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  dfloat *RES = ((dfloat*) scratch) + 0*Ncols;
  dfloat *Ad  = ((dfloat*) scratch) + 1*Ncols;
  dfloat *d   = ((dfloat*) scratch) + 2*Ncols;

  if(x_is_zero){ //skip the Ax if X is zero
    //RES = D^{-1}r
    vectorDotStar(Nrows, 1.0, A->diagInv, r, 0.0, RES);
    vectorSet(Nrows, 0.0, X);
    //d = invTheta*RES
    vectorAdd(Nrows, invTheta, RES, 0.0, d);
  } else {
    //RES = D^{-1}(r-Ax)
    A->SpMV(-1.0, X, 1.0, r, RES);
    vectorDotStar(Nrows, A->diagInv, RES);

    //d = invTheta*RES
    vectorAdd(Nrows, invTheta, RES, 0.0, d);
  }

  for (int k=0;k<ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    vectorAdd(Nrows, 1.0, d, 1.0, X);

    //r_k+1 = r_k - D^{-1}Ad_k
    A->SpMV(1.0, d, 0.0, Ad);
    vectorDotStar(Nrows, -1.0, A->diagInv, Ad, 1.0, RES);

    rho_np1 = 1.0/(2.*sigma-rho_n);

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    vectorAdd(Nrows, 2.0*rho_np1/delta, RES, rho_np1*rho_n, d);
    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  vectorAdd(Nrows, 1.0, d, 1.0, X);
}

void agmgLevel::smoothJacobi(occa::memory& o_r, occa::memory& o_X,
                             bool x_is_zero) {

  // occaTimerTic(parAlmond->device,"device smoothJacobi");
  if(x_is_zero){
    vectorDotStar(Nrows, 1.0, o_A->o_diagInv, o_r, 0.0, o_X);
    // occaTimerToc(parAlmond->device,"device smoothJacobi");
    return;
  }

  occa::memory& o_RES = o_scratch;

  // RES = r-A*x
  o_A->SpMV(-1.0, o_X, 1.0, o_r, o_RES);

  // x = x + alpha*inv(D)*RES
  vectorDotStar(Nrows, 1.0, o_A->o_diagInv, o_RES, 1.0, o_X);
  // occaTimerToc(parAlmond->device,"hyb smoothJacobi");
}

void agmgLevel::smoothDampedJacobi(occa::memory& o_r, occa::memory& o_X,
                                   bool x_is_zero){

  // occaTimerTic(parAlmond->device,"device smoothDampedJacobi");
  if(x_is_zero){
    vectorDotStar(Nrows, lambda, o_A->o_diagInv, o_r, 0.0, o_X);
    // occaTimerToc(parAlmond->device,"device smoothDampedJacobi");
    return;
  }

  occa::memory& o_RES = o_scratch;

  // RES = r-A*x
  o_A->SpMV(-1.0, o_X, 1.0, o_r, o_RES);

  // x = x + alpha*inv(D)*RES
  vectorDotStar(Nrows, lambda, o_A->o_diagInv, o_RES, 1.0, o_X);
  // occaTimerToc(parAlmond->device,"device smoothDampedJacobi");
}

void agmgLevel::smoothChebyshev(occa::memory& o_r, occa::memory& o_X,
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

  // occaTimerTic(parAlmond->device,"device smoothChebyshev");

  if(x_is_zero){ //skip the Ax if x is zero
    //RES = D^{-1}r
    vectorDotStar(Nrows, 1.0, o_A->o_diagInv, o_r, 0.0, o_RES);
    vectorSet(Nrows, 0.0, o_X);
    //d = invTheta*RES
    vectorAdd(Nrows, invTheta, o_RES, 0.0, o_d);
  } else {
    //RES = D^{-1}(r-Ax)
    o_A->SpMV(-1.0, o_X, 1.0, o_r, o_RES);
    vectorDotStar(Nrows, o_A->o_diagInv, o_RES);

    //d = invTheta*RES
    vectorAdd(Nrows, invTheta, o_RES, 0.0, o_d);
  }

  for (int k=0;k<ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    vectorAdd(Nrows, 1.0, o_d, 1.0, o_X);

    //r_k+1 = r_k - D^{-1}Ad_k
    o_A->SpMV(1.0, o_d, 0.0, o_Ad);
    vectorDotStar(Nrows, -1.0, o_A->o_diagInv, o_Ad, 1.0, o_RES);

    rho_np1 = 1.0/(2.*sigma-rho_n);

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    vectorAdd(Nrows, 2.0*rho_np1/delta, o_RES, rho_np1*rho_n, o_d);
    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  vectorAdd(Nrows, 1.0, o_d, 1.0, o_X);

  // occaTimerToc(parAlmond->device,"device smoothChebyshev");
}

} //namespace parAlmond

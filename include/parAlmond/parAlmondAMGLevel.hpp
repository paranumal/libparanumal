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

#ifndef PARALMOND_AMGLEVEL_HPP
#define PARALMOND_AMGLEVEL_HPP

#include "parAlmond.hpp"
#include "parAlmond/parAlmondparCSR.hpp"


namespace parAlmond {

class amgLevel: public multigridLevel {

public:
  parCSR *A=nullptr, *P=nullptr, *R=nullptr;

  SmoothType stype;
  dfloat lambda, lambda1, lambda0; //smoothing params

  int ChebyshevIterations=2;

  bool gatherLevel=false;
  ogs_t *ogs=nullptr;
  occa::memory o_gatherWeight;
  occa::memory o_Sx, o_Gx;

  amgLevel(parCSR *AA, settings_t& _settings);
  amgLevel(parCSR *AA, parCSR *PP, parCSR *RR, settings_t& _settings);
  ~amgLevel();

  void Operator(occa::memory& o_x, occa::memory& o_Ax);
  void residual(occa::memory& o_rhs, occa::memory& o_x, occa::memory& o_res);
  void coarsen(occa::memory& o_x, occa::memory& o_Cx);
  void prolongate(occa::memory& o_x, occa::memory& o_Px);

  void smooth(occa::memory& o_rhs, occa::memory& o_x, bool x_is_zero);
  void smoothDampedJacobi(occa::memory& o_r, occa::memory& o_x, bool x_is_zero);
  void smoothChebyshev(occa::memory& o_r, occa::memory& o_x, bool x_is_zero);

  void Report();

  /*   Setup routines */
  void setupSmoother();
  void syncToDevice();
};

}

#endif

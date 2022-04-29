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

#ifndef PARALMOND_AMGLEVEL_HPP
#define PARALMOND_AMGLEVEL_HPP

#include "parAlmond.hpp"
#include "parAlmond/parAlmondparCSR.hpp"

namespace libp {

namespace parAlmond {

class amgLevel: public multigridLevel {

public:
  parCSR A, P, R;

  SmoothType stype;
  dfloat lambda, lambda1, lambda0; //smoothing params

  int ChebyshevIterations=2;

  amgLevel() = default;
  amgLevel(parCSR& AA, settings_t& _settings);

  void Operator(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_Ax);
  void residual(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_res);
  void coarsen(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_Cx);
  void prolongate(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_Px);

  void smooth(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_x, bool x_is_zero);
  void smoothDampedJacobi(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_x, bool x_is_zero);
  void smoothChebyshev(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_x, bool x_is_zero);

  void Report();

  /*   Setup routines */
  void setupSmoother();
  void syncToDevice();
};

} //namespace parAlmond

} //namespace libp

#endif

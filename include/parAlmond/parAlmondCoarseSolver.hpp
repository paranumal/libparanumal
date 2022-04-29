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

#ifndef PARALMOND_COARSESOLVE_HPP
#define PARALMOND_COARSESOLVE_HPP

#include "settings.hpp"
#include "platform.hpp"
#include "solver.hpp"
#include "parAlmond.hpp"
#include "parAlmond/parAlmondDefines.hpp"
#include "parAlmond/parAlmondparCSR.hpp"

namespace libp {

namespace parAlmond {

class coarseSolver_t: public operator_t {

public:
  platform_t platform;
  settings_t settings;
  comm_t comm;

  int Nrows;
  int Ncols;

  int rank, size;

  coarseSolver_t(platform_t& _platform, settings_t& _settings,
                 comm_t _comm):
    platform(_platform), settings(_settings),
    comm(_comm) {}

  virtual int getTargetSize()=0;

  virtual void setup(parCSR& A, bool nullSpace,
                     memory<dfloat> nullVector, dfloat nullSpacePenalty)=0;

  virtual void syncToDevice()=0;

  virtual void Report(int lev)=0;

  virtual void solve(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_x)=0;
};

class exactSolver_t: public coarseSolver_t {

public:
  parCSR A;

  int coarseTotal;
  int coarseOffset;
  memory<int> coarseOffsets;
  memory<int> coarseCounts;
  memory<int> sendOffsets;
  memory<int> sendCounts;

  int N;
  int offdTotal=0;

  memory<dfloat> diagInvAT, offdInvAT;
  deviceMemory<dfloat> o_diagInvAT, o_offdInvAT;

  memory<dfloat> diagRhs, offdRhs;
  deviceMemory<dfloat> o_offdRhs;

  exactSolver_t(platform_t& _platform, settings_t& _settings,
                comm_t _comm):
    coarseSolver_t(_platform, _settings, _comm) {}

  int getTargetSize();

  void setup(parCSR& A, bool nullSpace,
             memory<dfloat> nullVector, dfloat nullSpacePenalty);

  void syncToDevice();

  void Report(int lev);

  void solve(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_x);
};

class oasSolver_t: public coarseSolver_t {

public:
  parCSR A;

  int N;
  int diagTotal=0, offdTotal=0;

  memory<dfloat> diagInvAT, offdInvAT;
  deviceMemory<dfloat> o_diagInvAT, o_offdInvAT;

  oasSolver_t(platform_t& _platform, settings_t& _settings,
              comm_t _comm):
    coarseSolver_t(_platform, _settings, _comm) {}

  int getTargetSize();

  void setup(parCSR& A, bool nullSpace,
             memory<dfloat> nullVector, dfloat nullSpacePenalty);

  void syncToDevice();

  void Report(int lev);

  void solve(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_x);
};

} //namespace parAlmond

} //namespace libp

#endif

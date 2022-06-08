/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Anthony Austin

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

#ifndef INITIALGUESS_HPP
#define INITIALGUESS_HPP

#include "core.hpp"
#include "platform.hpp"
#include "solver.hpp"

namespace libp {

namespace InitialGuess {


void AddSettings(settings_t& settings, const std::string prefix = "");

// Abstract base class for different initial guess strategies.
class initialGuessStrategy_t {
 protected:
  platform_t platform;
  settings_t settings;
  comm_t   comm;

  dlong Ntotal;     // Degrees of freedom

 public:
  initialGuessStrategy_t(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
    platform(_platform), settings(_settings), comm(_comm), Ntotal(_N) {}

  virtual void FormInitialGuess(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs) = 0;
  virtual void Update(operator_t& linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs) = 0;
};

// Default initial guess strategy:  use whatever the user gave us.
class Default : public initialGuessStrategy_t {
public:
  Default(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);

  void FormInitialGuess(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs);
  void Update(operator_t &linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs);
};

// Zero initial guess strategy:  use a zero initial guess.
class Zero : public initialGuessStrategy_t {
public:
  Zero(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);

  void FormInitialGuess(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs);
  void Update(operator_t &linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs);
};

// Initial guess strategies based on RHS projection.
class Projection : public initialGuessStrategy_t {
protected:
  dlong curDim;           // Current dimension of the initial guess space
  dlong maxDim;           // Maximum dimension of the initial guess space

  deviceMemory<dfloat> o_btilde;  //  vector (e.g., to be added to space)
  deviceMemory<dfloat> o_xtilde;  // Solution vector corresponding to o_btilde
  deviceMemory<dfloat> o_Btilde;  //  space (orthogonalized)
  deviceMemory<dfloat> o_Xtilde;  // Solution space corresponding to  space

  // temporary buffer for basis inner product output
  dlong        ctmpNblocks;
  pinnedMemory<dfloat> ctmp;
  deviceMemory<dfloat> o_ctmp;

  pinnedMemory<dfloat> alphas;    // Buffers for storing inner products.
  deviceMemory<dfloat> o_alphas;

  kernel_t igBasisInnerProductsKernel;
  kernel_t igReconstructKernel;
  kernel_t igScaleKernel;
  kernel_t igUpdateKernel;

  void igBasisInnerProducts(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_c, pinnedMemory<dfloat>& c);
  void igReconstruct(deviceMemory<dfloat>& o_u, dfloat a, deviceMemory<dfloat>& o_c, deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_unew);

public:
  Projection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);

  virtual void FormInitialGuess(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs);
  virtual void Update(operator_t& linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs) = 0;
};

// "Classic" initial guess strategy from Fischer's 1998 paper.
class ClassicProjection : public Projection {
public:
  ClassicProjection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);

  void Update(operator_t &linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs);
};

// Rolling QR update for projection history space a la Christensen's thesis.
class RollingQRProjection : public Projection {
private:
  pinnedMemory<dfloat>   R;   // R factor in QR decomposition (row major)
  deviceMemory<dfloat> o_R;

  kernel_t igDropQRFirstColumnKernel;

	void givensRotation(dfloat a, dfloat b, dfloat& c, dfloat& s);

public:
  RollingQRProjection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);

  void Update(operator_t &linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs);
};

// Extrapolation initial guess strategy.
class Extrap : public initialGuessStrategy_t {
private:
  int Nhistory;
  int shift;
  int entry;
  deviceMemory<dfloat> o_xh;
  deviceMemory<dfloat> o_coeffs;
  kernel_t igExtrapKernel;
  kernel_t igExtrapSparseKernel;

  int Nsparse;
  deviceMemory<int> o_sparseIds;
  deviceMemory<dfloat> o_sparseCoeffs;

  void extrapCoeffs(int m, int M, memory<dfloat> c);

public:
  Extrap(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);

  void FormInitialGuess(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs);
  void Update(operator_t &linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs);
};

} //namespace InitialGuess

} //namespace libp

#endif /* INITIALGUESS_HPP */

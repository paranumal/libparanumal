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

#include "initialGuess.hpp"
#include "mesh.hpp"

namespace libp {

namespace InitialGuess {

#define IG_BLOCKSIZE 256

void AddSettings(settings_t& settings, const std::string prefix)
{
  settings.newSetting(prefix + "INITIAL GUESS STRATEGY",
                      "NONE",
                      "Strategy for selecting initial guess for linear solver",
                      {"NONE", "ZERO", "CLASSIC", "QR", "EXTRAP"});

  settings.newSetting(prefix + "INITIAL GUESS HISTORY SPACE DIMENSION",
                      "-1",
                      "Dimension of the initial guess space");

  settings.newSetting(prefix + "INITIAL GUESS EXTRAP DEGREE",
                      "-1",
                      "Degree used for EXTRAP initial guess schemes.");

  settings.newSetting(prefix + "INITIAL GUESS EXTRAP COEFFS METHOD",
                      "MINNORM",
                      "Method for selecting coefficients with EXTRAP initial guess schemes.",
                      {"MINNORM", "CPQR"});
}

/*****************************************************************************/

Default::Default(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  initialGuessStrategy_t(_N, _platform, _settings, _comm)
{}

void Default::FormInitialGuess(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs)
{}

void Default::Update(operator_t &linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs)
{}

/*****************************************************************************/

Zero::Zero(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  initialGuessStrategy_t(_N, _platform, _settings, _comm)
{
  platform.linAlg().InitKernels({"set"});
}

void Zero::FormInitialGuess(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs)
{
  platform.linAlg().set(Ntotal, 0.0, o_x);
}

void Zero::Update(operator_t &linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs)
{}

/*****************************************************************************/

Projection::Projection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  initialGuessStrategy_t(_N, _platform, _settings, _comm)
{
  curDim = 0;
  settings.getSetting("INITIAL GUESS HISTORY SPACE DIMENSION", maxDim);

  o_btilde = platform.malloc<dfloat>(Ntotal);
  o_xtilde = platform.malloc<dfloat>(Ntotal);
  o_Btilde = platform.malloc<dfloat>(Ntotal*maxDim);
  o_Xtilde = platform.malloc<dfloat>(Ntotal*maxDim);

  alphas = platform.hostMalloc<dfloat>(maxDim);
  o_alphas = platform.malloc<dfloat>(maxDim);

  ctmpNblocks = (Ntotal + IG_BLOCKSIZE - 1)/IG_BLOCKSIZE;
  ctmp = platform.hostMalloc<dfloat>(ctmpNblocks*maxDim);
  o_ctmp = platform.malloc<dfloat>(ctmpNblocks*maxDim);

  // Build kernels.
  platform.linAlg().InitKernels({"set"});

  properties_t kernelInfo = platform.props();
  kernelInfo["defines/" "p_igNhist"] = maxDim;

  igBasisInnerProductsKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/igBasisInnerProducts.okl", "igBasisInnerProducts", kernelInfo);
  igReconstructKernel        = platform.buildKernel(LINEARSOLVER_DIR "/okl/igReconstruct.okl",        "igReconstruct",        kernelInfo);
  igScaleKernel              = platform.buildKernel(LINEARSOLVER_DIR "/okl/igScale.okl",              "igScale",              kernelInfo);
  igUpdateKernel             = platform.buildKernel(LINEARSOLVER_DIR "/okl/igUpdate.okl",             "igUpdate",             kernelInfo);
}

void Projection::FormInitialGuess(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs)
{
  if (curDim > 0) {
    igBasisInnerProducts(o_rhs, o_Btilde, o_alphas, alphas);
    platform.linAlg().set(Ntotal, 0.0, o_x);
    igReconstruct(o_x, 1.0, o_alphas, o_Xtilde, o_x);
  }
}

void Projection::igBasisInnerProducts(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_c, pinnedMemory<dfloat>& c)
{
  igBasisInnerProductsKernel(Ntotal, ctmpNblocks, curDim, o_x, o_Q, o_ctmp);

  ctmp.copyFrom(o_ctmp, ctmpNblocks*curDim);

  for (int m = 0; m < curDim; ++m) {
    c[m] = 0;
    for (int n = 0; n < ctmpNblocks; ++n) {
      c[m] += ctmp[m*ctmpNblocks + n];
    }
  }

  comm.Allreduce(c, Comm::Sum, curDim);
  c.copyTo(o_c, curDim);
}

void Projection::igReconstruct(deviceMemory<dfloat>& o_u, dfloat a, deviceMemory<dfloat>& o_c, deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_unew)
{
  igReconstructKernel(Ntotal, curDim, o_u, a, o_c, o_Q, o_unew);
}


/*****************************************************************************/

ClassicProjection::ClassicProjection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  Projection(_N, _platform, _settings, _comm)
{}

void ClassicProjection::Update(operator_t &linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs)
{
  // Compute RHS corresponding to the approximate solution obtained.
  linearOperator.Operator(o_x, o_btilde);

  // Insert new solution into the initial guess space.
  if ((curDim >= maxDim) || (curDim == 0)) {
    dfloat normbtilde = 0.0;

    normbtilde = platform.linAlg().norm2(Ntotal, o_btilde, comm);

    if (normbtilde > 0) {
      igScaleKernel(Ntotal, dfloat(1.0)/normbtilde, o_btilde, o_Btilde);
      igScaleKernel(Ntotal, dfloat(1.0)/normbtilde, o_x,      o_Xtilde);

      curDim = 1;
    }
  } else {
    dfloat    invnormbtilde = 0.0;
    const int Nreorth = 2;

    o_x.copyTo(o_xtilde, Ntotal);

    // Orthogonalize new RHS against previous ones.
    for (int n = 0; n < Nreorth; n++) {
      igBasisInnerProducts(o_btilde, o_Btilde, o_alphas, alphas);
      igReconstruct(o_btilde, -1.0, o_alphas, o_Btilde, o_btilde);
      igReconstruct(o_xtilde, -1.0, o_alphas, o_Xtilde, o_xtilde);
    }

    // Normalize.
    invnormbtilde = platform.linAlg().norm2(Ntotal, o_btilde, comm);
    invnormbtilde = 1.0/invnormbtilde;

#if 0
    igScaleKernel(Ntotal, invnormbtilde, o_btilde, o_btilde);
    igScaleKernel(Ntotal, invnormbtilde, o_xtilde, o_xtilde);

    // Store.
    o_btilde.copyTo(o_Btilde + curDim*Ntotal*sizeof(dfloat));
    o_xtilde.copyTo(o_Xtilde + curDim*Ntotal*sizeof(dfloat));
#else
    igUpdateKernel(Ntotal, curDim, invnormbtilde, o_btilde, 1, o_Btilde, o_xtilde, 1, o_Xtilde);
#endif

    curDim++;
  }
}

/*****************************************************************************/

RollingQRProjection::RollingQRProjection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  Projection(_N, _platform, _settings, _comm)
{
  R = platform.hostMalloc<dfloat>(maxDim*maxDim);
  o_R = platform.malloc<dfloat>(maxDim*maxDim);

  properties_t kernelInfo = platform.props();
  kernelInfo["defines/" "p_igNhist"] = maxDim;

  igDropQRFirstColumnKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/igDropQRFirstColumn.okl", "igDropQRFirstColumn", kernelInfo);
}

void RollingQRProjection::Update(operator_t &linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs)
{
  // Compute RHS corresponding to the approximate solution obtained.
  linearOperator.Operator(o_x, o_btilde);

  // Rotate the history space (QR update).
  if (curDim == maxDim) {
    // Drop the first column in the QR factorization:  R = R(:, 2:end).
    for (int j = 0; j < maxDim; j++) {
      for (int i = 0; i < maxDim - 1; i++)
        R[j*maxDim + i] = R[j*maxDim + (i + 1)];
      R[j*maxDim + (maxDim - 1)] = 0.0;
    }

    R.copyTo(o_R);

    // Update the RHS and solution spaces.
		igDropQRFirstColumnKernel(Ntotal, o_Btilde, o_Xtilde, o_R);

    // Restore R to triangular form (overlapped with Q update).
    for (int j = 0; j < maxDim - 1 ; j++) {
      dfloat c = 0.0, s = 0.0;
      dfloat Rjj   = R[j*maxDim + j];
      dfloat Rjp1j = R[(j + 1)*maxDim + j];

      givensRotation(Rjj, Rjp1j, c, s);

      for (int i = j; i < maxDim; i++) {
        dfloat Rji   = R[j*maxDim + i];
        dfloat Rjp1i = R[(j + 1)*maxDim + i];

        R[j*maxDim + i]       =  c*Rji + s*Rjp1i;
        R[(j + 1)*maxDim + i] = -s*Rji + c*Rjp1i;
      }
    }

    // Copy the updated R back to the device.
    platform.finish();
    R.copyTo(o_R);

    curDim--;
  }

  // Orthogonalize and tack on the new column.
  if (curDim == 0) {
    dfloat normbtilde = 0.0;

    normbtilde = platform.linAlg().norm2(Ntotal, o_btilde, comm);

    if (normbtilde > 0) {
#if 0
      igScaleKernel(Ntotal, 1.0/normbtilde, o_btilde, o_Btilde);
      igScaleKernel(Ntotal, 1.0/normbtilde, o_x,      o_Xtilde);
#else
      dfloat invnormbtilde = 1.0/normbtilde;
      igUpdateKernel(Ntotal, 0, invnormbtilde, o_btilde, 0, o_Btilde, o_x, 0, o_Xtilde);
#endif

      R[0] = normbtilde;

      curDim = 1;
    }
  } else {
    dfloat    normbtilde = 0.0, normbtildeproj = 0.0;;
    const int Nreorth = 2;

    o_x.copyTo(o_xtilde, Ntotal);

    // Compute the initial norm of the new vector.
    normbtilde = platform.linAlg().norm2(Ntotal, o_btilde, comm);

    // Zero the entries above/on the diagonal of the column of R into which we want to write.
    for (int i = 0; i < curDim; i++)
      R[i*maxDim + curDim] = 0.0;

    // Orthogonalize new RHS against previous ones.
    for (int n = 0; n < Nreorth; n++) {
      igBasisInnerProducts(o_btilde, o_Btilde, o_alphas, alphas);
      igReconstruct(o_btilde, (dfloat)(-1.0), o_alphas, o_Btilde, o_btilde);
      igReconstruct(o_xtilde, (dfloat)(-1.0), o_alphas, o_Xtilde, o_xtilde);

      for (int i = 0; i < curDim; i++)
        R[i*maxDim + curDim] += alphas[i];
    }

    // Normalize.
    normbtildeproj = platform.linAlg().norm2(Ntotal, o_btilde, comm);

    // Only add if the remainder after projection is large enough.
    //
    // TODO:  What is the appropriate criterion here?
    if (normbtildeproj/normbtilde > 1.0e-10) {
#if 0
      igScaleKernel(Ntotal, 1.0/normbtildeproj, o_btilde, o_btilde);
      igScaleKernel(Ntotal, 1.0/normbtildeproj, o_xtilde, o_xtilde);

      // Store.
      o_btilde.copyTo(o_Btilde + curDim*Ntotal*sizeof(dfloat));
      o_xtilde.copyTo(o_Xtilde + curDim*Ntotal*sizeof(dfloat));
#else
      dfloat invnormbtildeproj = 1.0/normbtildeproj;
      igUpdateKernel(Ntotal, curDim, invnormbtildeproj, o_btilde, 1, o_Btilde, o_xtilde, 1, o_Xtilde);
#endif

      R[curDim*maxDim + curDim] = normbtildeproj;

      curDim++;
    }
  }

  R.copyTo(o_R);
}

void RollingQRProjection::givensRotation(dfloat a, dfloat b, dfloat& c, dfloat& s)
{
	// Compute a Givens rotation that zeros the bottom component of [a ; b].
  if (b != 0) {
    dfloat h = hypot(a, b);
    dfloat d = 1.0/h;
    c = std::abs(a)*d;
    s = std::copysign(d, a)*b;
  } else {
    c = 1.0;
    s = 0.0;
  }
}

/*****************************************************************************/

Extrap::Extrap(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  initialGuessStrategy_t(_N, _platform, _settings, _comm)
{
  int M, m;
  settings.getSetting("INITIAL GUESS HISTORY SPACE DIMENSION", M);
  settings.getSetting("INITIAL GUESS EXTRAP DEGREE", m);

  memory<dfloat> c(M);
  extrapCoeffs(m, M, c);

  Nhistory = M;

  entry = 0;

  o_coeffs = platform.malloc<dfloat>(Nhistory, c);

  shift = 0;

  o_xh = platform.malloc<dfloat>(Nhistory*Ntotal);

  platform.linAlg().InitKernels({"set"});

  properties_t kernelInfo = platform.props();
  kernelInfo["defines/" "p_igNhist"] = Nhistory;

  igExtrapKernel       = platform.buildKernel(LINEARSOLVER_DIR "/okl/igExtrap.okl",       "igExtrap",   kernelInfo);
  igExtrapSparseKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/igExtrap.okl", "igExtrapSparse",   kernelInfo);

  platform.linAlg().set(Nhistory*Ntotal, 0.0, o_xh);
}

void Extrap::FormInitialGuess(deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs)
{
  if (entry < Nhistory) {
    int M, m;
    if (entry == Nhistory - 1) {
      settings.getSetting("INITIAL GUESS HISTORY SPACE DIMENSION", M);
      settings.getSetting("INITIAL GUESS EXTRAP DEGREE", m);
    } else {
      M = std::max(1, entry + 1);
      m = sqrt(static_cast<double>(M));
    }

    // Construct the extrapolation coefficients.
    memory<dfloat> c(Nhistory);
    memory<dfloat> d(Nhistory);
    memory<dfloat> sparseCoeffs(Nhistory);
    for (int n = 0; n < Nhistory; ++n) {
      c[n] = 0;
      d[n] = 0;
      sparseCoeffs[n] = 0;
    }

    if (M == 1) {
      d[Nhistory - 1] = 1.0;
    } else {
      extrapCoeffs(m, M, c);

      // need d[0:M-1] = {0, 0, 0, .., c[0], c[1], .., c[M-1]}
      for (int i = 0; i < M; i++)
        d[Nhistory - M + i] = c[i];
    }

    memory<int> sparseIds(Nhistory);
    Nsparse = 0;
    for (int n = 0; n < Nhistory; ++n) {
      if (std::abs(d[n]) > 1e-14) { // hmm
        sparseIds[Nsparse] = n;
        sparseCoeffs[Nsparse] = d[n];
        ++Nsparse;
      }
    }

    o_coeffs = platform.malloc<dfloat>(d);
    o_sparseIds = platform.malloc<int>(sparseIds);
    o_sparseCoeffs = platform.malloc<dfloat>(sparseCoeffs);

    ++entry;
  }

  if (settings.compareSetting("INITIAL GUESS EXTRAP COEFFS METHOD", "MINNORM"))
    igExtrapKernel(Ntotal, Nhistory, shift, o_coeffs, o_xh, o_x);
  else {
    igExtrapSparseKernel(Ntotal, Nhistory, shift, Nsparse, o_sparseIds, o_sparseCoeffs, o_xh, o_x);
  }
}

void Extrap::Update(operator_t &linearOperator, deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_rhs)
{
  deviceMemory<dfloat> o_tmp = o_xh + Ntotal*shift;
  o_x.copyTo(o_tmp, Ntotal);
  shift = (shift + 1) % Nhistory;
}

void Extrap::extrapCoeffs(int m, int M, memory<dfloat> c)
{
  LIBP_ABORT("Extrapolation space dimension (" << M << ") too low for degree (" << m << ").",
             M < m + 1);

  const dfloat h = 2.0/(M - 1);
  memory<dfloat> r(M);
  for (int i = 0; i < M; i++)
    r[i] = -1.0 + i*h;

  memory<dfloat> ro(1);
  ro[0] = 1.0 + h;  // Evaluation point.

  memory<dfloat> V;
  mesh_t::Vandermonde1D(m, r, V);

  memory<dfloat> b;
  mesh_t::Vandermonde1D(m, ro, b);

  if (settings.compareSetting("INITIAL GUESS EXTRAP COEFFS METHOD", "MINNORM")) {
    linAlg_t::matrixUnderdeterminedRightSolveMinNorm(M, m + 1, V, b, c);
  } else if (settings.compareSetting("INITIAL GUESS EXTRAP COEFFS METHOD", "CPQR")) {
    linAlg_t::matrixUnderdeterminedRightSolveCPQR(M, m + 1, V, b, c);
  }
}

} //namespace InitialGuess

} //namespace libp

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
#include "core.hpp"
#include "platform.hpp"
#include "solver.hpp"
#include "initialGuess.hpp"
#include "mesh.hpp"

namespace libp {

namespace InitialGuess {

#define IG_BLOCKSIZE 256

void AddSettings(settings_t& settings, const std::string prefix)
{
  settings.newSetting(prefix + "INITIAL GUESS STRATEGY",
                      "LAST",
                      "Strategy for selecting initial guess for linear solver",
                      {"LAST", "ZERO", "CLASSIC", "QR", "EXTRAP"});

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

template<typename T>  
Last::Last(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
initialGuessStrategy_t(_N, _platform, _settings, _comm)
{
  this->platform.linAlg().InitKernels({"set"});
  o_xLast = platform.malloc<T>(Ntotal);
  this->platform.linAlg().set(Ntotal, (T)0.0, o_xLast);
}

template<typename T>
void Last::FormInitialGuess(deviceMemory<T>& o_x, deviceMemory<T>& o_rhs)
{
  o_x.copyFrom(o_xLast, Ntotal, 0, properties_t("async", true));
}

template<typename T>
void Last::Update(operator_t &linearOperator, deviceMemory<T>& o_x, deviceMemory<T>& o_rhs)
{
  o_xLast.copyFrom(o_x, Ntotal, 0, properties_t("async", true));
}

/*****************************************************************************/

template<typename T>  
Zero::Zero(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  initialGuessStrategy_t(_N, _platform, _settings, _comm)
{
  platform.linAlg().InitKernels({"set"});
}

template<typename T>
void Zero::FormInitialGuess(deviceMemory<T>& o_x, deviceMemory<T>& o_rhs)
{
  platform.linAlg().set(Ntotal, (T)0.0, o_x);
}

template<typename T>
void Zero::Update(operator_t &linearOperator, deviceMemory<T>& o_x, deviceMemory<T>& o_rhs)
{}

/*****************************************************************************/

template<typename T> 
Projection::Projection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  initialGuessStrategy_t(_N, _platform, _settings, _comm)
{
  curDim = 0;
  settings.getSetting("INITIAL GUESS HISTORY SPACE DIMENSION", maxDim);

  o_Btilde = platform.malloc<T>(Ntotal*maxDim);
  o_Xtilde = platform.malloc<T>(Ntotal*maxDim);

  // Build kernels.
  platform.linAlg().InitKernels({"set", "axpy"});

  properties_t kernelInfo = platform.props();
  kernelInfo["defines/" "p_igNhist"] = maxDim;

  igBasisInnerProductsKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/igBasisInnerProducts.okl", "igBasisInnerProducts", kernelInfo);
  igReconstructKernel        = platform.buildKernel(LINEARSOLVER_DIR "/okl/igReconstruct.okl",        "igReconstruct",        kernelInfo);
  igUpdateKernel             = platform.buildKernel(LINEARSOLVER_DIR "/okl/igUpdate.okl",             "igUpdate",             kernelInfo);
}

template<typename T> 
void Projection::FormInitialGuess(deviceMemory<T>& o_x,
                                  deviceMemory<T>& o_rhs)
{
  if (curDim > 0) {
    deviceMemory<T> o_alphas = platform.reserve<T>(maxDim);
    pinnedMemory<T> h_alphas = platform.hostReserve<T>(maxDim);
    igBasisInnerProducts(o_rhs, o_Btilde, o_alphas, h_alphas);
    igReconstruct(0.0, o_x, 1.0, o_alphas, o_Xtilde, o_x);
  } else {
    platform.linAlg().set(Ntotal, (T)0.0, o_x);
  }
}

template<typename T> 
void Projection::igBasisInnerProducts(deviceMemory<T>& o_x,
                                      deviceMemory<T>& o_Q,
                                      deviceMemory<T>& o_alphas,
                                      pinnedMemory<T>& h_alphas)
{
  int Nblocks = (Ntotal+IG_BLOCKSIZE-1)/IG_BLOCKSIZE;
  Nblocks = std::min(Nblocks, IG_BLOCKSIZE); //limit to IG_BLOCKSIZE entries

  //pinned tmp buffer for reductions
  deviceMemory<T> o_scratch = platform.reserve<T>(maxDim*IG_BLOCKSIZE);

  igBasisInnerProductsKernel(Ntotal, Nblocks, curDim, o_x, o_Q, o_scratch, o_alphas);

  if (Nblocks>0) {
    h_alphas.copyFrom(o_alphas, curDim, 0, properties_t("async", true));
    platform.finish();
  } else {
    for (int m = 0; m < curDim; ++m) {
      h_alphas[m] = 0.0;
    }
  }

  comm.Allreduce(h_alphas, Comm::Sum, curDim);
  h_alphas.copyTo(o_alphas, curDim, 0, properties_t("async", true));
}

template<typename T> 
void Projection::igReconstruct(const T a,
                               deviceMemory<T>& o_u,
                               const T b,
                               deviceMemory<T>& o_alphas,
                               deviceMemory<T>& o_Q,
                               deviceMemory<T>& o_unew)
{
  igReconstructKernel(Ntotal, curDim, a, o_u, b, o_alphas, o_Q, o_unew);
}


/*****************************************************************************/
template<typename T>
ClassicProjection::ClassicProjection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  Projection(_N, _platform, _settings, _comm)
{}

template<typename T>
void ClassicProjection::Update(operator_t &linearOperator, deviceMemory<T>& o_x, deviceMemory<T>& o_rhs)
{
  // Compute RHS corresponding to the approximate solution obtained.
  deviceMemory<T> o_btilde = platform.reserve<T>(Ntotal);
  linearOperator.Operator(o_x, o_btilde);

  // Insert new solution into the initial guess space.
  if ((curDim >= maxDim) || (curDim == 0)) {
    T normbtilde = platform.linAlg().norm2(Ntotal, o_btilde, comm);

    if (normbtilde > 0) {
      platform.linAlg().axpy(Ntotal, (T)1.0/normbtilde, o_btilde, (T)0.0, o_Btilde);
      platform.linAlg().axpy(Ntotal, (T)1.0/normbtilde, o_x,      (T)0.0, o_Xtilde);

      curDim = 1;
    }
  } else {
    const int Nreorth = 2;

    deviceMemory<T> o_xtilde = platform.reserve<T>(Ntotal);
    deviceMemory<T> o_alphas = platform.reserve<T>(maxDim);
    pinnedMemory<T> h_alphas = platform.hostReserve<T>(maxDim);

    // Orthogonalize new RHS against previous ones.
    igBasisInnerProducts(o_btilde, o_Btilde, o_alphas, h_alphas);
    igReconstruct(1.0, o_btilde, -1.0, o_alphas, o_Btilde, o_btilde);
    igReconstruct(1.0,      o_x, -1.0, o_alphas, o_Xtilde, o_xtilde);

    for (int n = 1; n < Nreorth; n++) {
      igBasisInnerProducts(o_btilde, o_Btilde, o_alphas, h_alphas);
      igReconstruct(1.0, o_btilde, -1.0, o_alphas, o_Btilde, o_btilde);
      igReconstruct(1.0, o_xtilde, -1.0, o_alphas, o_Xtilde, o_xtilde);
    }

    // Normalize.
    T invnormbtilde = 1.0/platform.linAlg().norm2(Ntotal, o_btilde, comm);
    igUpdateKernel(Ntotal, curDim, invnormbtilde, o_btilde, o_Btilde, o_xtilde, o_Xtilde);

    curDim++;
  }
}

/*****************************************************************************/

template<typename T> 
RollingQRProjection::RollingQRProjection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  Projection(_N, _platform, _settings, _comm)
{
  R.malloc(maxDim*maxDim);

  h_c = platform.hostMalloc<T>(maxDim);
  h_s = platform.hostMalloc<T>(maxDim);
  o_c = platform.malloc<T>(maxDim);
  o_s = platform.malloc<T>(maxDim);

  properties_t kernelInfo = platform.props();
  kernelInfo["defines/" "p_igNhist"] = maxDim;

  igDropQRFirstColumnKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/igDropQRFirstColumn.okl", "igDropQRFirstColumn", kernelInfo);
}

template<typename T>
void RollingQRProjection::Update(operator_t &linearOperator, deviceMemory<T>& o_x, deviceMemory<T>& o_rhs)
{
  // Rotate the history space (QR update).
  if (curDim == maxDim) {
    // Drop the first column in the QR factorization:  R = R(:, 2:end).
    for (int j = 0; j < maxDim; j++) {
      for (int i = 0; i < maxDim - 1; i++)
        R[j*maxDim + i] = R[j*maxDim + (i + 1)];
      R[j*maxDim + (maxDim - 1)] = 0.0;
    }

    // Restore R to triangular form (overlapped with Q update).
    for (int j = 0; j < maxDim - 1 ; j++) {
      T Rjj   = R[j*maxDim + j];
      T Rjp1j = R[(j + 1)*maxDim + j];

      h_c[j] = 0.0, h_s[j] = 0.0;
      givensRotation(Rjj, Rjp1j, h_c[j], h_s[j]);

      for (int i = j; i < maxDim; i++) {
        T Rji   = R[j*maxDim + i];
        T Rjp1i = R[(j + 1)*maxDim + i];

        R[j*maxDim + i]       =  h_c[j]*Rji + h_s[j]*Rjp1i;
        R[(j + 1)*maxDim + i] = -h_s[j]*Rji + h_c[j]*Rjp1i;
      }
    }

    // Copy c and s to device
    h_c.copyTo(o_c, maxDim, 0, properties_t("async", true));
    h_s.copyTo(o_s, maxDim, 0, properties_t("async", true));

    // Update the RHS and solution spaces.
    igDropQRFirstColumnKernel(Ntotal, o_c, o_s, o_Btilde, o_Xtilde);

    curDim--;
  }

  // Compute RHS corresponding to the approximate solution obtained.
  deviceMemory<T> o_btilde = platform.reserve<T>(Ntotal);
  linearOperator.Operator(o_x, o_btilde);

  // Zero the column of R into which we want to write.
  for (int i = 0; i < maxDim; i++)
    R[i*maxDim + curDim] = 0.0;

  // Compute the initial norm of the new vector.
  T normbtilde = platform.linAlg().norm2(Ntotal, o_btilde, comm);

  // Orthogonalize and tack on the new column.
  if (curDim == 0) {
    if (normbtilde > 0) {
      T invnormbtilde = 1.0/normbtilde;
      igUpdateKernel(Ntotal, 0, invnormbtilde, o_btilde, o_Btilde, o_x, o_Xtilde);

      R[0] = normbtilde;

      curDim = 1;
    }
  } else {
    const int Nreorth = 2;

    deviceMemory<T> o_xtilde = platform.reserve<T>(Ntotal);
    deviceMemory<T> o_alphas = platform.reserve<T>(maxDim);
    pinnedMemory<T> h_alphas = platform.hostReserve<T>(maxDim);

    // Orthogonalize new RHS against previous ones.
    igBasisInnerProducts(o_btilde, o_Btilde, o_alphas, h_alphas);
    igReconstruct(1.0, o_btilde, -1.0, o_alphas, o_Btilde, o_btilde);
    igReconstruct(1.0,      o_x, -1.0, o_alphas, o_Xtilde, o_xtilde);

    for (int i = 0; i < curDim; i++)
      R[i*maxDim + curDim] += h_alphas[i];

    for (int n = 1; n < Nreorth; n++) {
      igBasisInnerProducts(o_btilde, o_Btilde, o_alphas, h_alphas);
      igReconstruct(1.0, o_btilde, -1.0, o_alphas, o_Btilde, o_btilde);
      igReconstruct(1.0, o_xtilde, -1.0, o_alphas, o_Xtilde, o_xtilde);

      for (int i = 0; i < curDim; i++)
        R[i*maxDim + curDim] += h_alphas[i];
    }

    // Normalize.
    T normbtildeproj = platform.linAlg().norm2(Ntotal, o_btilde, comm);

    // Only add if the remainder after projection is large enough.
    //
    // TODO:  What is the appropriate criterion here?
    if (normbtildeproj/normbtilde > 1.0e-10) {
      T invnormbtildeproj = 1.0/normbtildeproj;
      igUpdateKernel(Ntotal, curDim, invnormbtildeproj, o_btilde, o_Btilde, o_xtilde, o_Xtilde);

      R[curDim*maxDim + curDim] = normbtildeproj;

      curDim++;
    }
  }
}

template<typename T>
void RollingQRProjection::givensRotation(T a, T b, T& c, T& s)
{
  // Compute a Givens rotation that zeros the bottom component of [a ; b].
  if (b != 0) {
    T h = hypot(a, b);
    T d = 1.0/h;
    c = std::abs(a)*d;
    s = std::copysign(d, a)*b;
  } else {
    c = 1.0;
    s = 0.0;
  }
}

/*****************************************************************************/
template<typename T>
Extrap::Extrap(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
  initialGuessStrategy_t(_N, _platform, _settings, _comm)
{
  platform.linAlg().InitKernels({"set"});

  settings.getSetting("INITIAL GUESS HISTORY SPACE DIMENSION", Nhistory);
  settings.getSetting("INITIAL GUESS EXTRAP DEGREE", ExtrapDegree);

  entry = 0;
  shift = 0;

  h_coeffs = platform.hostMalloc<T>(Nhistory);
  h_sparseIds = platform.hostMalloc<int>(Nhistory);
  h_sparseCoeffs = platform.hostMalloc<T>(Nhistory);

  o_coeffs = platform.malloc<T>(Nhistory);
  o_sparseIds = platform.malloc<int>(Nhistory);
  o_sparseCoeffs = platform.malloc<T>(Nhistory);

  o_xh = platform.malloc<T>(Nhistory*Ntotal);

  properties_t kernelInfo = platform.props();
  kernelInfo["defines/" "p_igNhist"] = Nhistory;

  igExtrapKernel       = platform.buildKernel(LINEARSOLVER_DIR "/okl/igExtrap.okl",       "igExtrap",   kernelInfo);
  igExtrapSparseKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/igExtrap.okl", "igExtrapSparse",   kernelInfo);
}

template<typename T>
void Extrap::FormInitialGuess(deviceMemory<T>& o_x, deviceMemory<T>& o_rhs)
{
  if (entry == 0) {
    platform.linAlg().set(Ntotal, (T)0.0, o_x);
    return;
  }

  if (entry <= Nhistory) {
    int M = entry;
    int m;
    if (entry == Nhistory) {
      m = ExtrapDegree;
    } else {
      m = sqrt(static_cast<double>(M));
    }

    // Construct the extrapolation coefficients.
    for (int n = 0; n < Nhistory; ++n) {
      h_coeffs[n] = 0;
      h_sparseIds[n] = 0;
      h_sparseCoeffs[n] = 0;
    }

    if (M == 1) {
      h_coeffs[Nhistory - 1] = 1.0;
    } else {
      memory<T> c(Nhistory);
      extrapCoeffs(m, M, c);

      // need d[0:M-1] = {0, 0, 0, .., c[0], c[1], .., c[M-1]}
      for (int i = 0; i < M; i++)
        h_coeffs[Nhistory - M + i] = c[i];
    }


    Nsparse = 0;
    for (int n = 0; n < Nhistory; ++n) {
      if (std::abs(h_coeffs[n]) > 1e-14) { // hmm
        h_sparseIds[Nsparse] = n;
        h_sparseCoeffs[Nsparse] = h_coeffs[n];
        ++Nsparse;
      }
    }

    h_coeffs.copyTo(o_coeffs, properties_t("async", true));
    h_sparseIds.copyTo(o_sparseIds, properties_t("async", true));
    h_sparseCoeffs.copyTo(o_sparseCoeffs, properties_t("async", true));
  }

  if (settings.compareSetting("INITIAL GUESS EXTRAP COEFFS METHOD", "MINNORM"))
    igExtrapKernel(Ntotal, Nhistory, shift, o_coeffs, o_xh, o_x);
  else {
    igExtrapSparseKernel(Ntotal, Nhistory, shift, Nsparse, o_sparseIds, o_sparseCoeffs, o_xh, o_x);
  }
}

template<typename T>
void Extrap::Update(operator_t &linearOperator, deviceMemory<T>& o_x, deviceMemory<T>& o_rhs)
{
  deviceMemory<T> o_tmp = o_xh + Ntotal*shift;
  o_x.copyTo(o_tmp, Ntotal, 0, properties_t("async", true));
  shift = (shift + 1) % Nhistory;
  entry = std::min(Nhistory+1, entry+1);
}

template<typename T>
void Extrap::extrapCoeffs(int m, int M, memory<T> c)
{
  LIBP_ABORT("Extrapolation space dimension (" << M << ") too low for degree (" << m << ").",
             M < m + 1);

  const T h = 2.0/(M - 1);
  memory<T> r(M);
  for (int i = 0; i < M; i++)
    r[i] = -1.0 + i*h;

  memory<T> ro(1);
  ro[0] = 1.0 + h;  // Evaluation point.

  memory<T> V;
  mesh_t::Vandermonde1D(m, r, V);

  memory<T> b;
  mesh_t::Vandermonde1D(m, ro, b);

  if (settings.compareSetting("INITIAL GUESS EXTRAP COEFFS METHOD", "MINNORM")) {
    linAlg_t::matrixUnderdeterminedRightSolveMinNorm(M, m + 1, V, b, c);
  } else if (settings.compareSetting("INITIAL GUESS EXTRAP COEFFS METHOD", "CPQR")) {
    linAlg_t::matrixUnderdeterminedRightSolveCPQR(M, m + 1, V, b, c);
  }
}

} //namespace InitialGuess

} //namespace libp

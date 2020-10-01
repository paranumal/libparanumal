/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

initialGuessSolver_t* initialGuessSolver_t::Setup(dlong N, dlong Nhalo, platform_t& platform, settings_t& settings, MPI_Comm comm, int weighted, occa::memory& o_weight)
{
  initialGuessSolver_t* initialGuessSolver = new initialGuessSolver_t(N, Nhalo, platform, settings, comm);
  initialGuessSolver->linearSolver = linearSolver_t::Setup(N, Nhalo, platform, settings, comm, weighted, o_weight);
  initialGuessSolver->igStrategy = nullptr;

  if (settings.compareSetting("INITIAL GUESS STRATEGY", "NONE")) {
    initialGuessSolver->igStrategy = new igDefaultStrategy(N, platform, settings, comm, weighted, o_weight);
  } else if (settings.compareSetting("INITIAL GUESS STRATEGY", "ZERO")) {
    initialGuessSolver->igStrategy = new igZeroStrategy(N, platform, settings, comm, weighted, o_weight);
  } else if (settings.compareSetting("INITIAL GUESS STRATEGY", "CLASSIC")) {
    initialGuessSolver->igStrategy = new igClassicProjectionStrategy(N, platform, settings, comm, weighted, o_weight);
  } else if (settings.compareSetting("INITIAL GUESS STRATEGY", "QR")) {
    initialGuessSolver->igStrategy = new igRollingQRProjectionStrategy(N, platform, settings, comm, weighted, o_weight);
  } else {
    LIBP_ABORT("Requested INITIAL GUESS STRATEGY not found.");
  }

  return initialGuessSolver;
}

initialGuessSolver_t::initialGuessSolver_t(dlong _N, dlong _Nhalo, platform_t& _platform, settings_t& _settings, MPI_Comm _comm):
  linearSolver_t(_N, _Nhalo, _platform, _settings, _comm),
  igStrategy(nullptr),
  linearSolver(nullptr)
{
  return;
}

initialGuessSolver_t::~initialGuessSolver_t()
{
  delete igStrategy;
  delete linearSolver;
}

int initialGuessSolver_t::Solve(solver_t& solver, precon_t& precon, occa::memory& o_x, occa::memory& o_rhs, const dfloat tol, const int MAXIT, const int verbose)
{
  int iter = 0;

  igStrategy->FormInitialGuess(o_x, o_rhs);
  iter = linearSolver->Solve(solver, precon, o_x, o_rhs, tol, MAXIT, verbose);
  igStrategy->Update(solver, o_x, o_rhs);

  return iter;
}

/*****************************************************************************/

initialGuessStrategy_t::initialGuessStrategy_t(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm _comm, int _weighted, occa::memory& _o_weight):
  platform(_platform),
  settings(_settings),
  comm(_comm),
  Ntotal(_N)
{
  return;
}

initialGuessStrategy_t::~initialGuessStrategy_t()
{
  return;
}

/*****************************************************************************/

igDefaultStrategy::igDefaultStrategy(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm _comm, int _weighted, occa::memory& _o_weight):
  initialGuessStrategy_t(_N, _platform, _settings, _comm, _weighted, _o_weight)
{
  return;
}

void igDefaultStrategy::FormInitialGuess(occa::memory& o_x, occa::memory& o_rhs)
{
  return;
}

void igDefaultStrategy::Update(solver_t &solver, occa::memory& o_x, occa::memory& o_rhs)
{
  return;
}

/*****************************************************************************/

igZeroStrategy::igZeroStrategy(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm _comm, int _weighted, occa::memory& _o_weight):
  initialGuessStrategy_t(_N, _platform, _settings, _comm, _weighted, _o_weight)
{
  return;
}

void igZeroStrategy::FormInitialGuess(occa::memory& o_x, occa::memory& o_rhs)
{
  platform.linAlg.set(Ntotal, 0.0, o_x);
  return;
}

void igZeroStrategy::Update(solver_t &solver, occa::memory& o_x, occa::memory& o_rhs)
{
  return;
}

/*****************************************************************************/

igProjectionStrategy::igProjectionStrategy(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm _comm, int _weighted, occa::memory& _o_weight):
  initialGuessStrategy_t(_N, _platform, _settings, _comm, _weighted, _o_weight)
{
  curDim = 0;
  settings.getSetting("INITIAL GUESS HISTORY SPACE DIMENSION", maxDim);

  o_btilde = platform.malloc(Ntotal*sizeof(dfloat));
  o_xtilde = platform.malloc(Ntotal*sizeof(dfloat));
  o_Btilde = platform.malloc(Ntotal*maxDim*sizeof(dfloat));
  o_Xtilde = platform.malloc(Ntotal*maxDim*sizeof(dfloat));

  o_w = _o_weight;

  alphas = new dfloat[maxDim]();
  alphasThisRank = new dfloat[maxDim]();
  o_alphas = platform.malloc(maxDim*sizeof(dfloat));

  ctmpNblocks = (Ntotal + BLOCKSIZE - 1)/BLOCKSIZE;
  ctmp = (dfloat*)calloc(ctmpNblocks*maxDim, sizeof(dfloat));
  o_ctmp = platform.malloc(ctmpNblocks*maxDim*sizeof(dfloat), ctmp);

  // Build kernels.
  occa::properties kernelInfo = platform.props;
  kernelInfo["defines/" "p_igNhist"] = maxDim;

  igBasisInnerProductsKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/igBasisInnerProducts.okl", "igBasisInnerProducts", kernelInfo);
  igReconstructKernel        = platform.buildKernel(LINEARSOLVER_DIR "/okl/igReconstruct.okl",        "igReconstruct",        kernelInfo);
  igScaleKernel              = platform.buildKernel(LINEARSOLVER_DIR "/okl/igScale.okl",              "igScale",              kernelInfo);
  igUpdateKernel             = platform.buildKernel(LINEARSOLVER_DIR "/okl/igUpdate.okl",             "igUpdate",             kernelInfo);

  return;
}

igProjectionStrategy::~igProjectionStrategy()
{
  if (ctmp)
    delete[] ctmp;
  if (alphas)
    delete[] alphas;
  if (alphasThisRank)
    delete[] alphasThisRank;

  return;
}

void igProjectionStrategy::FormInitialGuess(occa::memory& o_x, occa::memory& o_rhs)
{
  if (curDim > 0) {
    igBasisInnerProducts(o_rhs, o_Btilde, o_alphas, alphas, alphasThisRank);
    platform.linAlg.set(Ntotal, 0.0, o_x);
    igReconstruct(o_x, 1.0, o_alphas, o_Xtilde, o_x);
  }

  return;
}

void igProjectionStrategy::igBasisInnerProducts(occa::memory& o_x, occa::memory& o_Q, occa::memory& o_c, dfloat *c, dfloat *cThisRank)
{
  igBasisInnerProductsKernel(Ntotal, ctmpNblocks, o_w, curDim, o_x, o_Q, o_ctmp);

  o_ctmp.copyTo(ctmp, ctmpNblocks*curDim*sizeof(dfloat));

  dlong cnt = 0;
  for (int m = 0; m < curDim; ++m) {
    cThisRank[m] = 0;
    for (int n = 0; n < ctmpNblocks; ++n) {
      cThisRank[m] += ctmp[cnt];
      ++cnt;
    }
  }

  MPI_Allreduce(cThisRank, c, curDim, MPI_DFLOAT, MPI_SUM, comm);
  o_c.copyFrom(c, curDim*sizeof(dfloat));

  return;
}

void igProjectionStrategy::igReconstruct(occa::memory& o_u, dfloat a, occa::memory& o_c, occa::memory& o_Q, occa::memory& o_unew)
{
  igReconstructKernel(Ntotal, curDim, o_u, a, o_c, o_Q, o_unew);
  return;
}


/*****************************************************************************/

igClassicProjectionStrategy::igClassicProjectionStrategy(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm _comm, int _weighted, occa::memory& _o_weight):
  igProjectionStrategy(_N, _platform, _settings, _comm, _weighted, _o_weight)
{
  return;
}

void igClassicProjectionStrategy::Update(solver_t &solver, occa::memory& o_x, occa::memory& o_rhs)
{
  // Compute RHS corresponding to the approximate solution obtained.
  solver.Operator(o_x, o_btilde);

  // Insert new solution into the initial guess space.
  if ((curDim >= maxDim) || (curDim == 0)) {
    dfloat normbtilde = 0.0;

    normbtilde = platform.linAlg.weightedNorm2(Ntotal,  o_w, o_btilde, comm);

    if (normbtilde > 0) {
      igScaleKernel(Ntotal, 1.0/normbtilde, o_btilde, o_Btilde);
      igScaleKernel(Ntotal, 1.0/normbtilde, o_x,      o_Xtilde);

      curDim = 1;
    }
  } else {
    dfloat    invnormbtilde = 0.0;
    const int Nreorth = 2;

    o_x.copyTo(o_xtilde, Ntotal*sizeof(dfloat));

    // Orthogonalize new RHS against previous ones.
    for (int n = 0; n < Nreorth; n++) {
      igBasisInnerProducts(o_btilde, o_Btilde, o_alphas, alphas, alphasThisRank);
      igReconstruct(o_btilde, (dfloat)(-1.0), o_alphas, o_Btilde, o_btilde);
      igReconstruct(o_xtilde, (dfloat)(-1.0), o_alphas, o_Xtilde, o_xtilde);
    }

    // Normalize.
    invnormbtilde = platform.linAlg.weightedNorm2(Ntotal, o_w, o_btilde, comm);
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

  return;
}

/*****************************************************************************/

igRollingQRProjectionStrategy::igRollingQRProjectionStrategy(dlong _N, platform_t& _platform, settings_t& _settings, MPI_Comm _comm, int _weighted, occa::memory& _o_weight):
  igProjectionStrategy(_N, _platform, _settings, _comm, _weighted, _o_weight)
{
  R = new dfloat[maxDim*maxDim]();
  o_R = platform.malloc(maxDim*maxDim*sizeof(dfloat));

  occa::properties kernelInfo = platform.props;
  kernelInfo["defines/" "p_igNhist"] = maxDim;

  igDropQRFirstColumnKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/igDropQRFirstColumn.okl", "igDropQRFirstColumn", kernelInfo);

  return;
}

igRollingQRProjectionStrategy::~igRollingQRProjectionStrategy()
{
  if (R)
    delete[] R;

  return;
}

void igRollingQRProjectionStrategy::Update(solver_t &solver, occa::memory& o_x, occa::memory& o_rhs)
{
  // Compute RHS corresponding to the approximate solution obtained.
  solver.Operator(o_x, o_btilde);

  // Rotate the history space (QR update).
  if (curDim == maxDim) {
    // Drop the first column in the QR factorization:  R = R(:, 2:end).
    for (int j = 0; j < maxDim; j++) {
      for (int i = 0; i < maxDim - 1; i++)
        R[j*maxDim + i] = R[j*maxDim + (i + 1)];
      R[j*maxDim + (maxDim - 1)] = 0.0;
    }

    o_R.copyFrom(R);

    // Update the RHS and solution spaces.
		igDropQRFirstColumnKernel(Ntotal, o_Btilde, o_Xtilde, o_R);

    // Restore R to triangular form (overlapped with Q update).
    for (int j = 0; j < maxDim - 1 ; j++) {
      dfloat c = 0.0, s = 0.0;
      dfloat Rjj   = R[j*maxDim + j];
      dfloat Rjp1j = R[(j + 1)*maxDim + j];

      givensRotation(Rjj, Rjp1j, &c, &s);

      for (int i = j; i < maxDim; i++) {
        dfloat Rji   = R[j*maxDim + i];
        dfloat Rjp1i = R[(j + 1)*maxDim + i];

        R[j*maxDim + i]       =  c*Rji + s*Rjp1i;
        R[(j + 1)*maxDim + i] = -s*Rji + c*Rjp1i;
      }
    }

    // Copy the updated R back to the device.
    platform.device.finish();
    o_R.copyFrom(R);

    curDim--;
  }

  // Orthogonalize and tack on the new column.
  if (curDim == 0) {
    dfloat normbtilde = 0.0;

    normbtilde = platform.linAlg.weightedNorm2(Ntotal, o_w, o_btilde, comm);

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

    o_x.copyTo(o_xtilde, Ntotal*sizeof(dfloat));

    // Compute the initial norm of the new vector.
    normbtilde = platform.linAlg.weightedNorm2(Ntotal, o_w, o_btilde, comm);

    // Zero the entries above/on the diagonal of the column of R into which we want to write.
    for (int i = 0; i < curDim; i++)
      R[i*maxDim + curDim] = 0.0;

    // Orthogonalize new RHS against previous ones.
    for (int n = 0; n < Nreorth; n++) {
      igBasisInnerProducts(o_btilde, o_Btilde, o_alphas, alphas, alphasThisRank);
      igReconstruct(o_btilde, (dfloat)(-1.0), o_alphas, o_Btilde, o_btilde);
      igReconstruct(o_xtilde, (dfloat)(-1.0), o_alphas, o_Xtilde, o_xtilde);

      for (int i = 0; i < curDim; i++)
        R[i*maxDim + curDim] += alphas[i];
    }

    // Normalize.
    normbtildeproj = platform.linAlg.weightedNorm2(Ntotal, o_w, o_btilde, comm);

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

  o_R.copyFrom(R);
}

void igRollingQRProjectionStrategy::givensRotation(dfloat a, dfloat b, dfloat *c, dfloat *s)
{
	// Compute a Givens rotation that zeros the bottom component of [a ; b].
  if (b != 0) {
    dfloat h = hypot(a, b);
    dfloat d = 1.0/h;
    *c = fabs(a)*d;
    *s = copysign(d, a)*b;
  } else {
    *c = 1.0;
    *s = 0.0;
  }

  return;
}

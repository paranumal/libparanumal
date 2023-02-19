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

#include "linearSolver.hpp"

namespace libp {

namespace LinearSolver {

template<typename T>
pminres<T>::pminres(dlong _N, dlong _Nhalo,
                    platform_t& _platform, settings_t& _settings, comm_t _comm):
  linearSolverBase_t(_N, _Nhalo, _platform, _settings, _comm)
{
  platform.linAlg().InitKernels({"axpy", "scale", "innerProd"});

  properties_t kernelInfo = platform.props();
  updateMINRESKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/linearSolverUpdateMINRES.okl", "updateMINRES", kernelInfo);
}

template<typename T>
int pminres<T>::Solve(operator_t& linearOperator, operator_t& precon,
                      deviceMemory<T>& o_x, deviceMemory<T>& o_b,
                      const T tol, const int MAXIT, const int verbose)
{
  int iter;
  T a0, a1, a2, a3, del, gam, gamp, c, cp, s, sp, eta;
  T TOL;

  int rank = comm.rank();
  linAlg_t &linAlg = platform.linAlg();

  deviceMemory<pfloat> o_pfloat_r;
  deviceMemory<pfloat> o_pfloat_z;

  /*Pre-reserve memory pool space to avoid some unnecessary re-sizing*/
  dlong Ntotal = N + Nhalo;
  if constexpr (sizeof(pfloat)==sizeof(T)) {
    platform.reserve<T>(6*Ntotal
                        + 6 * platform.memPoolAlignment<T>());
  } else {
    platform.reserve<std::byte>(sizeof(T) * (6*Ntotal)
                              + sizeof(pfloat) * 2 * Ntotal
                              + 8 * platform.memPoolAlignment());

    o_pfloat_r  = platform.reserve<pfloat>(Ntotal);
    o_pfloat_z  = platform.reserve<pfloat>(Ntotal);
  }

  deviceMemory<T> o_p     = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_z     = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_r     = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_r_old = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_q     = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_q_old = platform.reserve<T>(Ntotal);


  linearOperator.Operator(o_x, o_r);            // r = b - A*x
  linAlg.axpy(N, (T)1.0, o_b, (T)-1.0, o_r);


  // z = M\r
  if constexpr(sizeof(pfloat)==sizeof(T)){
    precon.Operator(o_r, o_z);
  } else {
    linAlg.d2p(N, o_r, o_pfloat_r);
    precon.Operator(o_pfloat_r, o_pfloat_z);
    linAlg.p2d(N, o_pfloat_z, o_z);
  }


  gamp = 0.0;
  gam  = sqrt(linAlg.innerProd(N, o_z, o_r, comm)); // gam = sqrt(z . r);
  eta  = gam;
  sp   = 0.0;
  s    = 0.0;
  cp   = 1.0;
  c    = 1.0;

  TOL = std::max(tol*std::abs(eta), tol);
  if (verbose && (rank == 0)) {
    printf("PMINRES:  initial eta = % .15e, target %.15e\n", eta, tol);
  }

  // MINRES iteration loop.
  iter = 0;
  while (iter < MAXIT) {
    if (verbose && (rank == 0)) {
      printf("PMINRES:  it %3d  eta = % .15e, gamma = %.15e\n", iter, eta, gam);
    }

    if ((std::abs(eta) < TOL) && (iter >= 1)) {
      if (verbose && (rank == 0)) {
        printf("PMINRES converged in %d iterations (eta = % .15e).\n", iter, eta);
      }
      return iter;
    }

    linAlg.scale(N, (T)1.0/gam, o_z);                    // z = z/gam
    linearOperator.Operator(o_z, o_p);                // p = A*z
    del = linAlg.innerProd(N, o_p, o_z, comm);        // del = p . z
    a0 = c*del - cp*s*gam;
    a2 = s*del + cp*c*gam;
    a3 = sp*gam;

    if (iter == 0) {
      o_q.copyFrom(o_z, properties_t("async", true));     // q = z
      o_r_old.copyFrom(o_r, properties_t("async", true)); // r_old = r
      linAlg.axpy(N, (T)1.0, o_p, (T)-(del/gam), o_r);          // r = p - (del/gam)*r
    }
    else if (iter == 1) {
      o_q_old.copyFrom(o_q, properties_t("async", true)); // q_old = q
      linAlg.axpy(N, (T)1.0, o_z, -a2, o_q);                 // q = z - a2*q
      o_z.copyFrom(o_r, properties_t("async", true));     // z = r (save r in z)
      linAlg.axpy(N, (T)1.0, o_p, -(del/gam), o_r);          // r = p - (del/gam)*r
      linAlg.axpy(N, (T)-(gam/gamp), o_r_old, (T)1.0, o_r);     // r = r - (gam/gamp)*r_old
      o_r_old.copyFrom(o_z, properties_t("async", true)); // r_old = z (i.e. r saved in z)
    }
    else {
#if 0
    linAlg.axpy(N, -a2, o_q, 1.0, o_z);               // z = z - a2*q - a3*q_old
    linAlg.axpy(N, -a3, o_q_old, 1.0, o_z);
    o_q_old.copyFrom(o_q);                            // q_old = q
    o_q.copyFrom(o_z);                                // q = z
    o_z.copyFrom(o_r);                                // z = r
    linAlg.axpy(N, 1.0, o_p, -(del/gam), o_r);        // r = p - (del/gam)*r
    if (iter > 0)
      linAlg.axpy(N, -(gam/gamp), o_r_old, 1.0, o_r); // r = r - (gam/gamp)*r_old
    o_r_old.copyFrom(o_z);                            // r_old = z
#else
      UpdateMINRES(-a2, -a3, -del/gam, -gam/gamp,
                   o_z, o_q_old, o_q, o_r_old, o_r, o_p);
#endif
    }

    // z = M\r
    if constexpr(sizeof(pfloat)==sizeof(T)){
      precon.Operator(o_r, o_z);
    } else {
      linAlg.d2p(N, o_r, o_pfloat_r);
      precon.Operator(o_pfloat_r, o_pfloat_z);
      linAlg.p2d(N, o_pfloat_z, o_z);
    }

    gamp = gam;
    gam  = sqrt(linAlg.innerProd(N, o_z, o_r, comm)); // gam = sqrt(z . r)
    a1   = sqrt(a0*a0 + gam*gam);
    cp   = c;
    c    = a0/a1;
    sp   = s;
    s    = gam/a1;
    linAlg.scale(N, (T)1.0/a1, o_q);                     // q = q/a1
    linAlg.axpy(N, c*eta, o_q, (T)1.0, o_x);             // x = x + c*eta*q
    eta = -s*eta;

    iter++;
  }

  return iter;
}

template<typename T>
void pminres<T>::UpdateMINRES(const T ma2,
                              const T ma3,
                              const T alpha,
                              const T beta,
                              deviceMemory<T>& o_z,
                              deviceMemory<T>& o_q_old,
                              deviceMemory<T>& o_q,
                              deviceMemory<T>& o_r_old,
                              deviceMemory<T>& o_r,
                              deviceMemory<T>& o_p)
{
  updateMINRESKernel(N, ma2, ma3, alpha, beta,
                     o_z, o_q_old, o_q, o_r_old, o_r, o_p);
}

template class pminres<float>;
template class pminres<double>;

} //namespace LinearSolver

} //namespace libp

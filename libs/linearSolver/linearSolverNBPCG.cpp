/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#define NBPCG_BLOCKSIZE 512

template <typename T>
nbpcg<T>::nbpcg(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, comm_t _comm):
  linearSolverBase_t(_N, _Nhalo, _platform, _settings, _comm) {

  platform.linAlg().InitKernels({"axpy"});

  //pinned buffer for reductions
  dots = platform.hostMalloc<T>(3);

  /* build kernels */
  properties_t kernelInfo = platform.props(); //copy base properties

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)NBPCG_BLOCKSIZE;

  // combined NBPCG update kernels
  update1NBPCGKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/linearSolverUpdateNBPCG.okl",
                                "update1NBPCG", kernelInfo);
  update2NBPCGKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/linearSolverUpdateNBPCG.okl",
                                "update2NBPCG", kernelInfo);
}

template <typename T>
int nbpcg<T>::Solve(operator_t& linearOperator, operator_t& precon,
                    deviceMemory<T>& o_x, deviceMemory<T>& o_r,
                    const T tol, const int MAXIT, const int verbose,
                    stoppingCriteria_t<T> *stoppingCriteria,
                    std::shared_ptr<InitialGuess::initialGuessStrategy_t> ig){

  int rank = comm.rank();
  linAlg_t &linAlg = platform.linAlg();

  deviceMemory<pfloat> o_pfloat_s;
  deviceMemory<pfloat> o_pfloat_S;

  /*Pre-reserve memory pool space to avoid some unnecessary re-sizing*/
  dlong Ntotal = N + Nhalo;
  if constexpr (sizeof(pfloat)==sizeof(T)) {
    platform.reserve<T>(6*Ntotal + 3*NBPCG_BLOCKSIZE
                        + 7 * platform.memPoolAlignment<T>());
  } else {
    platform.reserve<std::byte>(sizeof(T) * (6*Ntotal + 3*NBPCG_BLOCKSIZE)
                              + sizeof(pfloat) * 2 * Ntotal
                              + 9 * platform.memPoolAlignment());

    o_pfloat_s  = platform.reserve<pfloat>(Ntotal);
    o_pfloat_S  = platform.reserve<pfloat>(Ntotal);
  }

  /*aux variables */
  deviceMemory<T> o_p  = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_s  = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_S  = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_z  = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_Z  = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_Ax = platform.reserve<T>(Ntotal);

  // register scalars
  T zdotz0 = 0;
  T rdotr0 = 0;

  T alpha0 = 0;
  T beta0  = 0;
  T gamma0 = 0;
  T delta0 = 0;

  T gamma1 = 0; // history gamma

  // compute A*x
  linearOperator.Operator(o_x, o_Ax);

  // subtract r = r - A*x
  linAlg.axpy(N, (T) -1.f, o_Ax,  (T) 1.f, o_r);

   // z = M*r [ Gropp notation ]
  if constexpr (sizeof(pfloat)==sizeof(T)) {
    precon.Operator(o_r, o_z);
  } else {
    linAlg.d2p(N, o_r, o_pfloat_s);
    precon.Operator(o_pfloat_s, o_pfloat_S);
    linAlg.p2d(N, o_pfloat_S, o_z);
  }

  // set alpha = 0 to get
  // r.z and z.z
  alpha0 = 0;
  Update2NBPCG(alpha0, o_s, o_S, o_r, o_z);

  linearOperator.Operator(o_z, o_Z);

  comm.Wait(request);
  gamma0 = dots[0]; // rdotz
  zdotz0 = dots[1];
  rdotr0 = dots[2];

  T TOL = std::max(tol*tol*rdotr0,tol*tol);

  if (verbose&&(rank==0))
    printf("NBPCG: initial res norm %12.12f \n", sqrt(rdotr0));

  int iter;
  beta0 = 0;
  for(iter=0;iter<MAXIT;++iter){

    //exit if tolerance is reached
    if(rdotr0<=TOL) break;

    // p <= z + beta*p
    // s <= Z + beta*s
    // delta <= pdots
    Update1NBPCG(beta0, o_z, o_Z, o_p, o_s);

    // z = Precon^{-1} r
    if constexpr(sizeof(pfloat)==sizeof(T)) {
      precon.Operator(o_s, o_S);
    } else {
      linAlg.d2p(N, o_s, o_pfloat_s);
      precon.Operator(o_pfloat_s, o_pfloat_S);
      linAlg.p2d(N, o_pfloat_S, o_S);
    }

    // block for delta
    comm.Wait(request);
    delta0 = dots[0];

    // alpha = gamma/delta
    alpha0 = gamma0/delta0;

    // r <= r - alpha*s
    // z <= z - alpha*S
    // r.z
    // z.z
    // r.r
    Update2NBPCG(alpha0, o_s, o_S, o_r, o_z);

    // x <= x + alpha*p (delayed)
    linAlg.axpy(N, alpha0, o_p,  (T) 1.0, o_x);

    // Z = A*z
    linearOperator.Operator(o_z, o_Z);

    // block for delta
    comm.Wait(request);
    gamma1 = gamma0;
    gamma0 = dots[0]; // gamma = r.z
    zdotz0 = dots[1]; //
    rdotr0 = dots[2]; //

    beta0 = gamma0/gamma1;


    if (verbose&&(rank==0)) {
      if(rdotr0<0)
        printf("WARNING NBPCG: rdotr = %17.15lf\n", rdotr0);

      printf("NBPCG: it %d, r norm %12.12le, gamma = %le zdotz = %le \n", iter+1, sqrt(rdotr0), gamma0, zdotz0);
    }
  }

  return iter;
}

template <typename T>
void nbpcg<T>::Update1NBPCG(const T beta,
                            deviceMemory<T>& o_z,
                            deviceMemory<T>& o_Z,
                            deviceMemory<T>& o_p,
                            deviceMemory<T>& o_s){

  // p <= z + beta*p
  // s <= Z + beta*s
  // dot(p,s)
  int Nblocks = (N+NBPCG_BLOCKSIZE-1)/NBPCG_BLOCKSIZE;
  Nblocks = std::min(Nblocks,NBPCG_BLOCKSIZE); //limit to NBPCG_BLOCKSIZE entries

  //buffer for reductions
  deviceMemory<T> o_dots = platform.reserve<T>(NBPCG_BLOCKSIZE);

  update1NBPCGKernel(N, Nblocks, o_z, o_Z, beta, o_p, o_s, o_dots);

  if (Nblocks>0) {
    dots.copyFrom(o_dots, 1, 0, properties_t("async", true));
    platform.finish();
  } else {
    dots[0] = 0.0;
  }

  comm.Iallreduce(dots, Comm::Sum, 1, request);
}

template <typename T>
void nbpcg<T>::Update2NBPCG(const T alpha,
                            deviceMemory<T>& o_s,
                            deviceMemory<T>& o_S,
                            deviceMemory<T>& o_r,
                            deviceMemory<T>& o_z){

  // r <= r - alpha*s
  // z <= z - alpha*S
  // dot(r,z)
  // dot(z,z)
  // dot(r,r)
  int Nblocks = (N+NBPCG_BLOCKSIZE-1)/NBPCG_BLOCKSIZE;
  Nblocks = (Nblocks>NBPCG_BLOCKSIZE) ? NBPCG_BLOCKSIZE : Nblocks; //limit to NBPCG_BLOCKSIZE entries

  //buffer for reductions
  deviceMemory<T> o_dots = platform.reserve<T>(3*NBPCG_BLOCKSIZE);

  update2NBPCGKernel(N, Nblocks, o_s, o_S, alpha, o_r, o_z, o_dots);

  if (Nblocks>0) {
    dots.copyFrom(o_dots, 3, 0, properties_t("async", true));
    platform.finish();
  } else {
    dots[0] = 0.0;
    dots[1] = 0.0;
    dots[2] = 0.0;
  }

  comm.Iallreduce(dots, Comm::Sum, 3, request);
}

template class nbpcg<float>;
template class nbpcg<double>;

} //namespace LinearSolver

} //namespace libp

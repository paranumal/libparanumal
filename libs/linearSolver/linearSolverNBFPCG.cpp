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

#define NBFPCG_BLOCKSIZE 512

nbfpcg::nbfpcg(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, comm_t _comm):
  linearSolverBase_t(_N, _Nhalo, _platform, _settings, _comm) {

  platform.linAlg().InitKernels({"axpy", "zaxpy"});

  //pinned buffer for reductions
  dots = platform.hostMalloc<dfloat>(4);

  /* build kernels */
  properties_t kernelInfo = platform.props(); //copy base properties

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)NBFPCG_BLOCKSIZE;

  // combined NBFPCG update kernels
  update0NBFPCGKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/linearSolverUpdateNBFPCG.okl",
                                "update0NBFPCG", kernelInfo);
  update1NBFPCGKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/linearSolverUpdateNBFPCG.okl",
                                "update1NBFPCG", kernelInfo);
}

int nbfpcg::Solve(operator_t& linearOperator, operator_t& precon,
                  deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_r,
                  const dfloat tol, const int MAXIT, const int verbose) {

  int rank = comm.rank();
  linAlg_t &linAlg = platform.linAlg();

  /*Pre-reserve memory pool space to avoid some unnecessary re-sizing*/
  dlong Ntotal = N + Nhalo;
  platform.reserve<dfloat>(8*Ntotal + 4*NBFPCG_BLOCKSIZE
                           + 9 * platform.memPoolAlignment<dfloat>());

  /*aux variables */
  deviceMemory<dfloat> o_u = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_p = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_w = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_n = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_m = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_s = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_z = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_q = platform.reserve<dfloat>(Ntotal);

  // register scalars
  dfloat alpha0 = 0;
  dfloat beta0  = 0;
  dfloat gamma0 = 0;
  dfloat delta0 = 0;
  dfloat eta0   = 0;
  dfloat rdotr0 = 0;

  // compute A*x
  linearOperator.Operator(o_x, o_u);

  // subtract r = r - A*x
  linAlg.axpy(N, -1.f, o_u, 1.f, o_r);

  // u = M*r [ Sanan notation ]
  precon.Operator(o_r, o_u);

  // p = u
  o_p.copyFrom(o_u, properties_t("async", true));

  // w = A*p
  linearOperator.Operator(o_p, o_w);

  // gamma = u.r
  // delta = u.w
  Update0NBFPCG(o_u, o_r, o_w);

  precon.Operator(o_w, o_m);

  linearOperator.Operator(o_m, o_n);
  o_s.copyFrom(o_w, properties_t("async", true));
  o_q.copyFrom(o_m, properties_t("async", true));
  o_z.copyFrom(o_n, properties_t("async", true));

  comm.Wait(request);
  gamma0 = dots[0]; // udotr
  delta0 = dots[1]; // udotw
  rdotr0 = dots[2]; // rdotr
  eta0   = delta0;
  alpha0 = gamma0/eta0;

  dfloat TOL = std::max(tol*tol*rdotr0,tol*tol);

  if (verbose&&(rank==0))
    printf("NBFPCG: initial res norm %12.12f \n", sqrt(rdotr0));

  int iter;
  beta0 = 0;
  for(iter=0;iter<MAXIT;++iter){

    //exit if tolerance is reached
    if(rdotr0<=TOL) break;

    // x <= x + alpha*p
    // r <= r - alpha*s
    // u <= u - alpha*q
    // w <= w - alpha*z
    // gamma <= u.r
    // beta  <= -u.s/eta
    // delta <= u.w
    Update1NBFPCG(alpha0, o_p, o_s, o_q, o_z, o_x, o_r, o_u, o_w);

    // n <= (w-r)
    linAlg.zaxpy(N, 1.0, o_w, -1.0, o_r, o_n);

    // m <= M*(w-r)
    precon.Operator(o_n, o_m);

    // m <= u + M*(w-r)
    linAlg.axpy(N, 1.0, o_u, 1.0, o_m);

    // n = A*m
    linearOperator.Operator(o_m, o_n);

    // block for delta
    comm.Wait(request);
    gamma0 = dots[0];       //  u.r
    beta0  = -dots[1]/eta0; // -u.s/eta
    delta0 = dots[2];       //  u.w
    rdotr0 = dots[3];       // r.r

    //  p <= u + beta*p
    linAlg.axpy(N, 1.0, o_u, beta0, o_p);

    //  s <= w + beta*s
    linAlg.axpy(N, 1.0, o_w, beta0, o_s);

    //  q <= m + beta*q
    linAlg.axpy(N, 1.0, o_m, beta0, o_q);

    //  z <= n + beta*z
    linAlg.axpy(N, 1.0, o_n, beta0, o_z);

    // eta = delta - beta^2*eta
    eta0 = delta0 - beta0*beta0*eta0;

    // alpha = gamma/eta
    alpha0 = gamma0/eta0;


    if (verbose&&(rank==0)) {
      if(rdotr0<0)
        printf("WARNING NBFPCG: rdotr = %17.15lf\n", rdotr0);

      printf("NBFPCG: it %d, r norm %12.12le, alpha = %12.12le beta = %le "
       "rdotr = %le gamma = %le delta = %le, eta = %le \n", iter+1, sqrt(rdotr0), alpha0, beta0, rdotr0, gamma0, delta0, eta0);
    }
  }

  return iter;
}

void nbfpcg::Update0NBFPCG(deviceMemory<dfloat>& o_u,
                           deviceMemory<dfloat>& o_r,
                           deviceMemory<dfloat>& o_w){

  // (u.r)
  // (u.w)
  // (r.r)
  int Nblocks = (N+NBFPCG_BLOCKSIZE-1)/NBFPCG_BLOCKSIZE;
  Nblocks = std::min(Nblocks, NBFPCG_BLOCKSIZE); //limit to NBFPCG_BLOCKSIZE entries

  deviceMemory<dfloat> o_dots = platform.reserve<dfloat>(3*NBFPCG_BLOCKSIZE);

  update0NBFPCGKernel(N, Nblocks, o_u, o_r, o_w, o_dots);

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

void nbfpcg::Update1NBFPCG(const dfloat alpha,
                           deviceMemory<dfloat>& o_p,
                           deviceMemory<dfloat>& o_s,
                           deviceMemory<dfloat>& o_q,
                           deviceMemory<dfloat>& o_z,
                           deviceMemory<dfloat>& o_x,
                           deviceMemory<dfloat>& o_r,
                           deviceMemory<dfloat>& o_u,
                           deviceMemory<dfloat>& o_w){

  // p <= z + beta*p
  // s <= Z + beta*s
  // dot(p,s)
  int Nblocks = (N+NBFPCG_BLOCKSIZE-1)/NBFPCG_BLOCKSIZE;
  Nblocks = std::min(Nblocks,NBFPCG_BLOCKSIZE); //limit to NBFPCG_BLOCKSIZE entries

  deviceMemory<dfloat> o_dots = platform.reserve<dfloat>(4*NBFPCG_BLOCKSIZE);

  update1NBFPCGKernel(N, Nblocks, o_p, o_s, o_q, o_z, alpha, o_x, o_r, o_u, o_w, o_dots);

  if (Nblocks>0) {
    dots.copyFrom(o_dots, 4, 0, properties_t("async", true));
    platform.finish();
  } else {
    dots[0] = 0.0;
    dots[1] = 0.0;
    dots[2] = 0.0;
    dots[3] = 0.0;
  }

  comm.Iallreduce(dots, Comm::Sum, 4, request);
}

} //namespace LinearSolver

} //namespace libp

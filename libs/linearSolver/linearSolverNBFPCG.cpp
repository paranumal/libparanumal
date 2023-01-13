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

template<typename T>
nbfpcg<T>::nbfpcg(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, comm_t _comm):
  linearSolverBase_t(_N, _Nhalo, _platform, _settings, _comm) {

  platform.linAlg().InitKernels({"axpy", "zaxpy"});

  //pinned buffer for reductions
  dots = platform.hostMalloc<T>(4);

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

  template<typename T>
int nbfpcg<T>::Solve(operator_t& linearOperator, operator_t& precon,
                  deviceMemory<T>& o_x, deviceMemory<T>& o_r,
                  const T tol, const int MAXIT, const int verbose) {

  int rank = comm.rank();
  linAlg_t &linAlg = platform.linAlg();

  /*Pre-reserve memory pool space to avoid some unnecessary re-sizing*/
  dlong Ntotal = N + Nhalo;
  platform.reserve<T>(8*Ntotal + 4*NBFPCG_BLOCKSIZE
                           + 9 * platform.memPoolAlignment<T>());

  /*aux variables */
  deviceMemory<T> o_u = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_p = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_w = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_n = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_m = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_s = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_z = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_q = platform.reserve<T>(Ntotal);

  platform.reserve<pfloat>(2*Ntotal +
                           + 4 * platform.memPoolAlignment<T>());
  
  deviceMemory<pfloat> o_pfloat_tmp  = platform.reserve<pfloat>(Ntotal);
  deviceMemory<pfloat> o_pfloat_Ptmp  = platform.reserve<pfloat>(Ntotal);

  
  // register scalars
  T alpha0 = 0;
  T beta0  = 0;
  T gamma0 = 0;
  T delta0 = 0;
  T eta0   = 0;
  T rdotr0 = 0;

  // compute A*x
  linearOperator.Operator(o_x, o_u);

  // subtract r = r - A*x
  linAlg.axpy(N, (T)-1.f, o_u, (T)1.f, o_r);

  // u = M*r [ Sanan notation ]
  //  precon.Operator(o_r, o_u);
  if(sizeof(pfloat)==sizeof(T)){
    precon.Operator(o_r, o_u);
  }
  else{
    linAlg.d2p(N, o_r, o_pfloat_tmp);
    precon.Operator(o_pfloat_tmp, o_pfloat_Ptmp);
    linAlg.p2d(N, o_pfloat_Ptmp, o_u);
  }
  

  // p = u
  o_p.copyFrom(o_u, properties_t("async", true));

  // w = A*p
  linearOperator.Operator(o_p, o_w);

  // gamma = u.r
  // delta = u.w
  Update0NBFPCG(o_u, o_r, o_w);

  //  precon.Operator(o_w, o_m);
  if(sizeof(pfloat)==sizeof(T)){
    precon.Operator(o_w, o_m);
  }
  else{
    linAlg.d2p(N, o_w, o_pfloat_tmp);
    precon.Operator(o_pfloat_tmp, o_pfloat_Ptmp);
    linAlg.p2d(N, o_pfloat_Ptmp, o_m);
  }


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

  T TOL = std::max(tol*tol*rdotr0,tol*tol);

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
    linAlg.zaxpy(N, (T)1.0, o_w, (T)-1.0, o_r, o_n);

    // m <= M*(w-r)
    //    precon.Operator(o_n, o_m);
    if(sizeof(pfloat)==sizeof(T)){
      precon.Operator(o_n, o_m);
    }
    else{
      linAlg.d2p(N, o_n, o_pfloat_tmp);
      precon.Operator(o_pfloat_tmp, o_pfloat_Ptmp);
      linAlg.p2d(N, o_pfloat_Ptmp, o_m);
    }
    

    // m <= u + M*(w-r)
    linAlg.axpy(N, (T)1.0, o_u, (T)1.0, o_m);

    // n = A*m
    linearOperator.Operator(o_m, o_n);

    // block for delta
    comm.Wait(request);
    gamma0 = dots[0];       //  u.r
    beta0  = -dots[1]/eta0; // -u.s/eta
    delta0 = dots[2];       //  u.w
    rdotr0 = dots[3];       // r.r

    T one = 1.0;
    
    //  p <= u + beta*p
    linAlg.axpy(N, one, o_u, beta0, o_p);

    //  s <= w + beta*s
    linAlg.axpy(N, one, o_w, beta0, o_s);

    //  q <= m + beta*q
    linAlg.axpy(N, one, o_m, beta0, o_q);

    //  z <= n + beta*z
    linAlg.axpy(N, one, o_n, beta0, o_z);

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

  template<typename T>
void nbfpcg<T>::Update0NBFPCG(deviceMemory<T>& o_u,
                           deviceMemory<T>& o_r,
                           deviceMemory<T>& o_w){

  // (u.r)
  // (u.w)
  // (r.r)
  int Nblocks = (N+NBFPCG_BLOCKSIZE-1)/NBFPCG_BLOCKSIZE;
  Nblocks = std::min(Nblocks, NBFPCG_BLOCKSIZE); //limit to NBFPCG_BLOCKSIZE entries

  deviceMemory<T> o_dots = platform.reserve<T>(3*NBFPCG_BLOCKSIZE);

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

  template<typename T>
void nbfpcg<T>::Update1NBFPCG(const T alpha,
                           deviceMemory<T>& o_p,
                           deviceMemory<T>& o_s,
                           deviceMemory<T>& o_q,
                           deviceMemory<T>& o_z,
                           deviceMemory<T>& o_x,
                           deviceMemory<T>& o_r,
                           deviceMemory<T>& o_u,
                           deviceMemory<T>& o_w){

  // p <= z + beta*p
  // s <= Z + beta*s
  // dot(p,s)
  int Nblocks = (N+NBFPCG_BLOCKSIZE-1)/NBFPCG_BLOCKSIZE;
  Nblocks = std::min(Nblocks,NBFPCG_BLOCKSIZE); //limit to NBFPCG_BLOCKSIZE entries

  deviceMemory<T> o_dots = platform.reserve<T>(4*NBFPCG_BLOCKSIZE);

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

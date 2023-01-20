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

#define PCG_BLOCKSIZE 512

template <typename T>
pcg<T>::pcg(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, comm_t _comm):
  linearSolverBase_t(_N, _Nhalo, _platform, _settings, _comm) {

  platform.linAlg().InitKernels({"axpy", "innerProd", "norm2", "d2p", "p2d"});

  flexible = settings.compareSetting("LINEAR SOLVER", "FPCG");

  /* build kernels */
  properties_t kernelInfo = platform.props(); //copy base properties

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)PCG_BLOCKSIZE;

  // combined PCG update and r.r kernel
  updatePCGKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/linearSolverUpdatePCG.okl",
					 "updatePCG", kernelInfo);
}

template <typename T>
int pcg<T>::Solve(operator_t& linearOperator, operator_t& precon,
                  deviceMemory<T>& o_x, deviceMemory<T>& o_r,
                  const T tol, const int MAXIT, const int verbose) {

  int rank = comm.rank();
  linAlg_t &linAlg = platform.linAlg();

  // register scalars
  T rdotz1 = 0.0;
  T rdotz2 = 0.0;
  T alpha = 0.0, beta = 0.0, pAp = 0.0;
  T rdotr0 = 0.0;
  T TOL = 0.0;

  /*Pre-reserve memory pool space to avoid some unnecessary re-sizing*/
  dlong Ntotal = N + Nhalo;
  // TW 3=>5
  platform.reserve<T>(3*Ntotal + PCG_BLOCKSIZE
                           + 4 * platform.memPoolAlignment<T>());

  platform.reserve<pfloat>(2*Ntotal +
                           + 4 * platform.memPoolAlignment<T>());
  /*aux variables */
  deviceMemory<T> o_p  = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_z  = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_Ap = platform.reserve<T>(Ntotal);

  deviceMemory<pfloat> o_pfloat_r  = platform.reserve<pfloat>(Ntotal);
  deviceMemory<pfloat> o_pfloat_z  = platform.reserve<pfloat>(Ntotal);

  // Comput norm of RHS (for stopping tolerance).
  if (settings.compareSetting("LINEAR SOLVER STOPPING CRITERION", "ABS/REL-RHS-2NORM")) {
    T normb = linAlg.norm2(N, o_r, comm);
    TOL = std::max(tol*tol*normb*normb, tol*tol);
  }

  // compute A*x
  linearOperator.Operator(o_x, o_Ap);

  // subtract r = r - A*x
  linAlg.axpy(N, (T)-1.f, o_Ap, (T)1.f, o_r);

  rdotr0 = linAlg.norm2(N, o_r, comm);
  rdotr0 = rdotr0*rdotr0;

  if (settings.compareSetting("LINEAR SOLVER STOPPING CRITERION", "ABS/REL-INITRESID")) {
    TOL = std::max(tol*tol*rdotr0,tol*tol);
  }

  if (verbose&&(rank==0))
    printf("PCG: initial res norm %12.12f \n", sqrt(rdotr0));

  int iter;
  for(iter=0;iter<MAXIT;++iter){

    // Exit if tolerance is reached, taking at least one step.
    if (((iter == 0) && (rdotr0 == 0.0)) ||
        ((iter > 0) && (rdotr0 <= TOL))) {
      break;
    }

    // z = Precon^{-1} r
    if(sizeof(pfloat)==sizeof(T)){
      precon.Operator(o_r, o_z);
    }
    else{
      linAlg.d2p(N, o_r, o_pfloat_r);
      precon.Operator(o_pfloat_r, o_pfloat_z);
      linAlg.p2d(N, o_pfloat_z, o_z);
    }

    // r.z
    rdotz2 = rdotz1;
    rdotz1 = linAlg.innerProd(N, o_r, o_z, comm);

    if(flexible){
      if (iter==0) {
        beta = 0.0;
      } else {
        T zdotAp = linAlg.innerProd(N, o_z, o_Ap, comm);
        beta = -alpha*zdotAp/rdotz2;
      }
    } else {
      beta = (iter==0) ? 0.0 : rdotz1/rdotz2;
    }

    // p = z + beta*p
    linAlg.axpy(N, (T)1.f, o_z, beta, o_p);

    // A*p
    linearOperator.Operator(o_p, o_Ap);

    // p.Ap
    pAp =  linAlg.innerProd(N, o_p, o_Ap, comm);

    alpha = rdotz1/pAp;

    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)
    rdotr0 = UpdatePCG(alpha, o_p, o_Ap, o_x, o_r);

    if (verbose&&(rank==0)) {
      if(rdotr0<0)
        printf("WARNING CG: rdotr = %17.15lf\n", rdotr0);

      printf("CG: it %d, r norm %12.12le, alpha = %le \n", iter+1, sqrt(rdotr0), alpha);
    }
  }

  return iter;
}

template <typename T>
T pcg<T>::UpdatePCG(const T alpha,
                    deviceMemory<T>& o_p,
                    deviceMemory<T>& o_Ap,
                    deviceMemory<T>& o_x,
                    deviceMemory<T>& o_r){

  // x <= x + alpha*p
  // r <= r - alpha*A*p
  // dot(r,r)
  int Nblocks = (N+PCG_BLOCKSIZE-1)/PCG_BLOCKSIZE;
  Nblocks = std::min(Nblocks, PCG_BLOCKSIZE); //limit to PCG_BLOCKSIZE entries

  //pinned tmp buffer for reductions
  pinnedMemory<T> h_rdotr = platform.hostReserve<T>(1);
  deviceMemory<T> o_rdotr = platform.reserve<T>(PCG_BLOCKSIZE);

  updatePCGKernel(N, Nblocks, o_Ap, alpha, o_r, o_rdotr);

  T rdotr;
  if (Nblocks>0) {
    h_rdotr.copyFrom(o_rdotr, 1, 0, properties_t("async", true));
    platform.finish();
    rdotr = h_rdotr[0];
  } else {
    rdotr = 0.0;
  }

  // x <= x + alpha*p
  platform.linAlg().axpy(N, alpha, o_p, (T)1.0, o_x);

  /*Compute allreduce while axpy is running*/
  comm.Allreduce(rdotr);
  return rdotr;
}

} //namespace LinearSolver

} //namespace libp

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

nbpcg::nbpcg(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, comm_t _comm):
  linearSolverBase_t(_N, _Nhalo, _platform, _settings, _comm) {

  platform.linAlg().InitKernels({"axpy"});

  dlong Ntotal = N + Nhalo;

  memory<dfloat> dummy(Ntotal, 0.0);

  /*aux variables */
  o_p  = platform.malloc<dfloat>(Ntotal, dummy);
  o_s  = platform.malloc<dfloat>(Ntotal, dummy);
  o_S  = platform.malloc<dfloat>(Ntotal, dummy);
  o_z  = platform.malloc<dfloat>(Ntotal);
  o_Z  = platform.malloc<dfloat>(Ntotal);
  o_Ax = platform.malloc<dfloat>(Ntotal);

  //pinned tmp buffer for reductions
  dots = platform.hostMalloc<dfloat>(3*NBPCG_BLOCKSIZE);
  o_dots = platform.malloc<dfloat>(3*NBPCG_BLOCKSIZE);

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

int nbpcg::Solve(operator_t& linearOperator, operator_t& precon,
                 deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_r,
                 const dfloat tol, const int MAXIT, const int verbose) {

  int rank = comm.rank();
  linAlg_t &linAlg = platform.linAlg();

  // register scalars
  dfloat zdotz0 = 0;
  dfloat rdotr0 = 0;

  dfloat alpha0 = 0;
  dfloat beta0  = 0;
  dfloat gamma0 = 0;
  dfloat delta0 = 0;

  dfloat gamma1 = 0; // history gamma

  // compute A*x
  linearOperator.Operator(o_x, o_Ax);

  // subtract r = r - A*x
  linAlg.axpy(N, -1.f, o_Ax, 1.f, o_r);

   // z = M*r [ Gropp notation ]
  precon.Operator(o_r, o_z);

  // set alpha = 0 to get
  // r.z and z.z
  alpha0 = 0;
  Update2NBPCG(alpha0, o_r);

  linearOperator.Operator(o_z, o_Z);

  comm.Wait(request);
  gamma0 = dots[0]; // rdotz
  zdotz0 = dots[1];
  rdotr0 = dots[2];

  dfloat TOL = std::max(tol*tol*rdotr0,tol*tol);

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
    Update1NBPCG(beta0);

    // z = Precon^{-1} r
    precon.Operator(o_s, o_S);

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
    Update2NBPCG(alpha0, o_r);

    // x <= x + alpha*p (delayed)
    linAlg.axpy(N, alpha0, o_p, 1.0, o_x);

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

void nbpcg::Update1NBPCG(const dfloat beta){

  // p <= z + beta*p
  // s <= Z + beta*s
  // dot(p,s)
  int Nblocks = (N+NBPCG_BLOCKSIZE-1)/NBPCG_BLOCKSIZE;
  Nblocks = std::min(Nblocks,NBPCG_BLOCKSIZE); //limit to NBPCG_BLOCKSIZE entries

  update1NBPCGKernel(N, Nblocks, o_z, o_Z, beta, o_p, o_s, o_dots);

  if (Nblocks>0) {
    dots.copyFrom(o_dots, Nblocks);
  } else {
    dots[0] = 0.0;
  }

  for(int n=1;n<Nblocks;++n)
    dots[0] += dots[n];

  comm.Iallreduce(dots, Comm::Sum, 1, request);
}

void nbpcg::Update2NBPCG(const dfloat alpha, deviceMemory<dfloat>& o_r){

  // r <= r - alpha*s
  // z <= z - alpha*S
  // dot(r,z)
  // dot(z,z)
  // dot(r,r)
  int Nblocks = (N+NBPCG_BLOCKSIZE-1)/NBPCG_BLOCKSIZE;
  Nblocks = (Nblocks>NBPCG_BLOCKSIZE) ? NBPCG_BLOCKSIZE : Nblocks; //limit to NBPCG_BLOCKSIZE entries

  update2NBPCGKernel(N, Nblocks, o_s, o_S, alpha, o_r, o_z, o_dots);

  if (Nblocks>0) {
    dots.copyFrom(o_dots, 3*Nblocks);
  } else {
    dots[0] = 0.0;
    dots[1] = 0.0;
    dots[2] = 0.0;
  }

  for(int n=1;n<Nblocks;++n) {
    dots[0] += dots[0+3*n];
    dots[1] += dots[1+3*n];
    dots[2] += dots[2+3*n];
  }
  comm.Iallreduce(dots, Comm::Sum, 3, request);
}

} //namespace LinearSolver

} //namespace libp

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

#include "linearSolver.hpp"

#define NBFPCG_BLOCKSIZE 512

nbfpcg::nbfpcg(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, MPI_Comm _comm):
  linearSolver_t(_N, _Nhalo, _platform, _settings, _comm) {

  dlong Ntotal = N + Nhalo;

  /*aux variables */
  o_u  = platform.malloc(Ntotal*sizeof(dfloat));
  o_p  = platform.malloc(Ntotal*sizeof(dfloat));
  o_w  = platform.malloc(Ntotal*sizeof(dfloat));
  o_n  = platform.malloc(Ntotal*sizeof(dfloat));
  o_m  = platform.malloc(Ntotal*sizeof(dfloat));
  o_s  = platform.malloc(Ntotal*sizeof(dfloat));
  o_z  = platform.malloc(Ntotal*sizeof(dfloat));
  o_q  = platform.malloc(Ntotal*sizeof(dfloat));
  o_Ax = platform.malloc(Ntotal*sizeof(dfloat));

  localdots  = (dfloat*) calloc(4, sizeof(dfloat));
  globaldots = (dfloat*) calloc(4, sizeof(dfloat));

  //pinned tmp buffer for reductions
  tmpdots = (dfloat*) platform.hostMalloc(4*NBFPCG_BLOCKSIZE*sizeof(dfloat),
                                          NULL, h_tmpdots);

  o_tmpdots = platform.malloc(4*NBFPCG_BLOCKSIZE*sizeof(dfloat));

  /* build kernels */
  occa::properties kernelInfo = platform.props; //copy base properties

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)NBFPCG_BLOCKSIZE;

  // combined NBFPCG update kernels
  update0NBFPCGKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/linearSolverUpdateNBFPCG.okl",
                                "update0NBFPCG", kernelInfo);
  update1NBFPCGKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/linearSolverUpdateNBFPCG.okl",
                                "update1NBFPCG", kernelInfo);
}

int nbfpcg::Solve(solver_t& solver, precon_t& precon,
                  occa::memory &o_x, occa::memory &o_r,
                  const dfloat tol, const int MAXIT, const int verbose) {

  int rank;
  MPI_Comm_rank(comm, &rank);
  linAlg_t<dfloat> &linAlg = platform.linAlg;

  // register scalars
  dfloat alpha0 = 0;
  dfloat beta0  = 0;
  dfloat gamma0 = 0;
  dfloat delta0 = 0;
  dfloat eta0   = 0;
  dfloat rdotr0 = 0;

  // compute A*x
  solver.Operator(o_x, o_Ax);

  // subtract r = r - A*x
  linAlg.axpy(N, -1.f, o_Ax, 1.f, o_r);

  // u = M*r [ Sanan notation ]
  precon.Operator(o_r, o_u);

  // p = u
  o_p.copyFrom(o_u);

  // w = A*p
  solver.Operator(o_p, o_w);

  // gamma = u.r
  // delta = u.w
  Update0NBFPCG(o_r);

  precon.Operator(o_w, o_m);

  solver.Operator(o_m, o_n);
  o_s.copyFrom(o_w);
  o_q.copyFrom(o_m);
  o_z.copyFrom(o_n);

  MPI_Wait(&request, &status);
  gamma0 = globaldots[0]; // udotr
  delta0 = globaldots[1]; // udotw
  rdotr0 = globaldots[2]; // rdotr
  eta0   = delta0;
  alpha0 = gamma0/eta0;

  dfloat TOL = mymax(tol*tol*rdotr0,tol*tol);

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
    Update1NBFPCG(alpha0, o_x, o_r);

    // n <= (w-r)
    linAlg.zaxpy(N, 1.0, o_w, -1.0, o_r, o_n);

    // m <= M*(w-r)
    precon.Operator(o_n, o_m);

    // m <= u + M*(w-r)
    linAlg.axpy(N, 1.0, o_u, 1.0, o_m);

    // n = A*m
    solver.Operator(o_m, o_n);

    // block for delta
    MPI_Wait(&request, &status);
    gamma0 = globaldots[0];       //  u.r
    beta0  = -globaldots[1]/eta0; // -u.s/eta
    delta0 = globaldots[2];       //  u.w
    rdotr0 = globaldots[3];       // r.r

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

void nbfpcg::Update0NBFPCG(occa::memory &o_r){

  // (u.r)
  // (u.w)
  // (r.r)
  int Nblocks = (N+NBFPCG_BLOCKSIZE-1)/NBFPCG_BLOCKSIZE;
  Nblocks = (Nblocks>NBFPCG_BLOCKSIZE) ? NBFPCG_BLOCKSIZE : Nblocks; //limit to NBFPCG_BLOCKSIZE entries

  update0NBFPCGKernel(N, Nblocks, o_u, o_r, o_w, o_tmpdots);

  o_tmpdots.copyTo(tmpdots, 3*Nblocks*sizeof(dfloat));

  localdots[0] = 0;
  localdots[1] = 0;
  localdots[2] = 0;
  for(int n=0;n<Nblocks;++n) {
    localdots[0] += tmpdots[0+3*n];
    localdots[1] += tmpdots[1+3*n];
    localdots[2] += tmpdots[2+3*n];
  }

  globaldots[0] = 0;
  globaldots[1] = 0;
  globaldots[2] = 0;
  MPI_Iallreduce(localdots, globaldots, 3, MPI_DFLOAT, MPI_SUM, comm, &request);
}

void nbfpcg::Update1NBFPCG(const dfloat alpha, occa::memory &o_x, occa::memory &o_r){

  // p <= z + beta*p
  // s <= Z + beta*s
  // dot(p,s)
  int Nblocks = (N+NBFPCG_BLOCKSIZE-1)/NBFPCG_BLOCKSIZE;
  Nblocks = (Nblocks>NBFPCG_BLOCKSIZE) ? NBFPCG_BLOCKSIZE : Nblocks; //limit to NBFPCG_BLOCKSIZE entries

  update1NBFPCGKernel(N, Nblocks, o_p, o_s, o_q, o_z, alpha, o_x, o_r, o_u, o_w, o_tmpdots);

  o_tmpdots.copyTo(tmpdots, 4*Nblocks*sizeof(dfloat));

  localdots[0] = 0;
  localdots[1] = 0;
  localdots[2] = 0;
  localdots[3] = 0;
  for(int n=0;n<Nblocks;++n) {
    localdots[0] += tmpdots[0+4*n];
    localdots[1] += tmpdots[1+4*n];
    localdots[2] += tmpdots[2+4*n];
    localdots[3] += tmpdots[3+4*n];
  }

  globaldots[0] = 0;
  globaldots[1] = 0;
  globaldots[2] = 0;
  globaldots[3] = 0;
  MPI_Iallreduce(localdots, globaldots, 4, MPI_DFLOAT, MPI_SUM, comm, &request);
}

nbfpcg::~nbfpcg() {
  update0NBFPCGKernel.free();
  update1NBFPCGKernel.free();
}

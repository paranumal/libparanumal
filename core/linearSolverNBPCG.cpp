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

#define NBPCG_BLOCKSIZE 512

nbpcg::nbpcg(solver_t& _solver):
  linearSolver_t(_solver) {};

nbpcg::~nbpcg() {
  update1NBPCGKernel.free();
  update2NBPCGKernel.free();
}

void nbpcg::Init(int _weighted, occa::memory& o_weight_) {

  N = mesh.Np*mesh.Nelements;
  dlong Nhalo  = mesh.Np*mesh.totalHaloPairs;
  dlong Ntotal = N + Nhalo;

  /*aux variables */
  o_p  = device.malloc(Ntotal*sizeof(dfloat));
  o_s  = device.malloc(Ntotal*sizeof(dfloat));
  o_S  = device.malloc(Ntotal*sizeof(dfloat));
  o_z  = device.malloc(Ntotal*sizeof(dfloat));
  o_Z  = device.malloc(Ntotal*sizeof(dfloat));
  o_Ax = device.malloc(Ntotal*sizeof(dfloat));

  localdots  = (dfloat*) calloc(3, sizeof(dfloat));
  globaldots = (dfloat*) calloc(3, sizeof(dfloat));

  weighted = _weighted;
  o_weight = o_weight_;

  //pinned tmp buffer for reductions
  occa::properties mprops;
  mprops["mapped"] = true;
  h_tmpdots = device.malloc(3*NBPCG_BLOCKSIZE*sizeof(dfloat), mprops);
  tmpdots = (dfloat*) h_tmpdots.ptr(mprops);
  o_tmpdots = device.malloc(3*NBPCG_BLOCKSIZE*sizeof(dfloat));

  /* build kernels */
  occa::properties kernelInfo = props; //copy base properties

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)NBPCG_BLOCKSIZE;

  // combined NBPCG update kernels
  update1NBPCGKernel = buildKernel(device,
                                LIBP_DIR "/core/okl/linearSolverUpdateNBPCG.okl",
                                "update1NBPCG", kernelInfo, comm);
  update2NBPCGKernel = buildKernel(device,
                                LIBP_DIR "/core/okl/linearSolverUpdateNBPCG.okl",
                                "update2NBPCG", kernelInfo, comm);
}

int nbpcg::Solve(solver_t& solver, precon_t& precon,
                 occa::memory &o_x, occa::memory &o_r,
                 const dfloat tol, const int MAXIT, const int verbose) {

  int rank = mesh.rank;

  // register scalars
  dfloat zdotz0 = 0;
  dfloat rdotr0 = 0;

  dfloat alpha0 = 0;
  dfloat beta0  = 0;
  dfloat gamma0 = 0;
  dfloat delta0 = 0;

  dfloat gamma1 = 0; // history gamma

  // compute A*x
  solver.Operator(o_x, o_Ax);

  // subtract r = r - A*x
  linAlg.axpy(N, -1.f, o_Ax, 1.f, o_r);

   // z = M*r [ Gropp notation ]
  precon.Operator(o_r, o_z);

  // set alpha = 0 to get
  // r.z and z.z
  alpha0 = 0;
  Update2NBPCG(alpha0, o_r);

  solver.Operator(o_z, o_Z);

  MPI_Wait(&request, &status);
  gamma0 = globaldots[0]; // rdotz
  zdotz0 = globaldots[1];
  rdotr0 = globaldots[2];

  dfloat TOL = mymax(tol*tol*rdotr0,tol*tol);

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
    MPI_Wait(&request, &status);
    delta0 = globaldots[0];

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
    solver.Operator(o_z, o_Z);

    // block for delta
    MPI_Wait(&request, &status);
    gamma1 = gamma0;
    gamma0 = globaldots[0]; // gamma = r.z
    zdotz0 = globaldots[1]; //
    rdotr0 = globaldots[2]; //

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
  Nblocks = (Nblocks>NBPCG_BLOCKSIZE) ? NBPCG_BLOCKSIZE : Nblocks; //limit to NBPCG_BLOCKSIZE entries

  update1NBPCGKernel(N, Nblocks, weighted, o_weight, o_z, o_Z, beta, o_p, o_s, o_tmpdots);

  o_tmpdots.copyTo(tmpdots, Nblocks*sizeof(dfloat));

  localdots[0] = 0;
  for(int n=0;n<Nblocks;++n)
    localdots[0] += tmpdots[n];

  globaldots[0] = 0;
  MPI_Iallreduce(localdots, globaldots, 1, MPI_DFLOAT, MPI_SUM, comm, &request);
}

void nbpcg::Update2NBPCG(const dfloat alpha, occa::memory &o_r){

  // r <= r - alpha*s
  // z <= z - alpha*S
  // dot(r,z)
  // dot(z,z)
  // dot(r,r)
  int Nblocks = (N+NBPCG_BLOCKSIZE-1)/NBPCG_BLOCKSIZE;
  Nblocks = (Nblocks>NBPCG_BLOCKSIZE) ? NBPCG_BLOCKSIZE : Nblocks; //limit to NBPCG_BLOCKSIZE entries

  update2NBPCGKernel(N, Nblocks, weighted, o_weight, o_s, o_S, alpha, o_r, o_z, o_tmpdots);

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
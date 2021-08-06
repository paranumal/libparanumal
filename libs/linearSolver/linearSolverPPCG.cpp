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

#define PPCG_BLOCKSIZE 512
#define PPCG_NREDUCTIONS 7

ppcg::ppcg(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, MPI_Comm _comm):
  linearSolver_t(_N, _Nhalo, _platform, _settings,  _comm) {

  dlong Ntotal = N+Nhalo;

  printf("ppcg::ppcg N=%d\n", Ntotal);
  
  // Make sure LinAlg has the necessary kernels
  platform.linAlg.InitKernels({"axpy", "amxpy",
	"innerProd", "norm2", "set"});

  Nblocks = (Ntotal+PPCG_BLOCKSIZE-1)/PPCG_BLOCKSIZE;

  /*aux variables */
  dfloat *dummy = (dfloat *) calloc(Ntotal,sizeof(dfloat)); //need this to avoid uninitialized memory warnings

  o_p  = platform.malloc(N*sizeof(dfloat),dummy);
  o_v  = platform.malloc(N*sizeof(dfloat),dummy);
  o_z  = platform.malloc(N*sizeof(dfloat),dummy);
  o_invM  = platform.malloc(N*sizeof(dfloat),dummy);
  
  free(dummy);

  //pinned tmp buffer for reductions
  reductionTmps = (dfloat*) platform.hostMalloc(PPCG_NREDUCTIONS*Nblocks*sizeof(dfloat), NULL, h_reductionTmps);
  o_reductionTmps = platform.malloc(PPCG_NREDUCTIONS*Nblocks*sizeof(dfloat));
  
  globalTmps = (dfloat*) calloc(PPCG_NREDUCTIONS, sizeof(dfloat));
  localTmps = (dfloat*) calloc(PPCG_NREDUCTIONS, sizeof(dfloat));
  
  /* build kernels */
  occa::properties kernelInfo = platform.props; //copy base properties

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)PPCG_BLOCKSIZE;

  // combined PPCG reductions kernel
  reductionsPPCGKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/linearSolverPPCG.okl",
					      "reductionsPPCG", kernelInfo);
  // combined update kernel
  updatePPCGKernel = platform.buildKernel(LINEARSOLVER_DIR "/okl/linearSolverPPCG.okl",
					  "updatePPCG", kernelInfo);
}

void ppcg::SetupPreconditioner(occa::memory &_o_invM) {
  o_invM.copyFrom(_o_invM);
}

int ppcg::Solve(solver_t& solver, precon_t& precon,
               occa::memory &o_x, occa::memory &o_r,
               const dfloat tol, const int MAXIT, const int verbose) {

  int rank;
  MPI_Comm_rank(comm, &rank);
  linAlg_t &linAlg = platform.linAlg;

  // register scalars
  dfloat a, b, c, d, e, f, g;
  dfloat alpha, alpha_old;
  dfloat beta, beta_old;
  dfloat TOL = 0.0;

  // Comput norm of RHS (for stopping tolerance).
  if (settings.compareSetting("LINEAR SOLVER STOPPING CRITERION", "ABS/REL-RHS-2NORM")) {
    dfloat normb = linAlg.norm2(N, o_r, comm);
    TOL = mymax(tol*tol*normb*normb, tol*tol);
  }

  // compute A*x
  solver.Operator(o_x, o_Ax);

  // subtract r = r - A*x
  linAlg.axpy(N, -1.f, o_Ax, 1.f, o_r);

  g = linAlg.norm2(N, o_r, comm);
  g = g*g;

  if (settings.compareSetting("LINEAR SOLVER STOPPING CRITERION", "ABS/REL-INITRESID")) {
    TOL = mymax(tol*tol*g,tol*tol);
  }

  if (verbose&&(rank==0))
    printf("PPCG: initial res norm %12.12f \n", sqrt(g));

  // TW TO DO -
  // 1. alpha = (r.(M\r))/(r.r);
  // 2. beta = ?
  // 3. check  if(fabs(g-2.*alpha*b+alpha*alpha*c)<TOL){ .. } 
  // 4. need to initialize p, v
  
  int iter = 0;

  // do a warm start PCG iterations =====>
  linAlg.amxpy(N, 1.0, o_invM, o_r, 0.0, o_z); // z = M\r ( a bit lazy )

  o_z.copyTo(o_p, N*sizeof(dfloat)); // p = z

  e = linAlg.innerProd(N, o_r, o_z, comm); // e = r.z

  for(iter=0;iter<2;++iter){
  
    solver.Operator(o_p, o_v); // v = A*p
    
    a = linAlg.innerProd(N, o_p, o_v, comm); // e = r.z

    alpha_old = alpha;
    alpha = e/a; // alpha = r.z/p.v
    
    linAlg.axpy(N,  alpha, o_p, 1.0, o_x); // x <= x + alpha*p
    linAlg.axpy(N, -alpha, o_v, 1.0, o_r); // r <= r - alpha*v
    
    g = linAlg.norm2(N, o_r, comm);
    g = g*g;
    
    linAlg.amxpy(N, 1.0, o_invM, o_r, 0.0, o_z); // z <= M\r ( a bit lazy )
    
    dfloat enew = linAlg.innerProd(N, o_r, o_z, comm); // enw = r.z
    beta_old = beta;
    beta  = enew/e; // beta = r1.z1/r0.z0
    e = enew; // e = r1.z1
    
    linAlg.axpy(N, 1.0, o_z, beta, o_p); // p <= z + beta*p
  }
  // <==== end PCG iteration here
  
  // do PPCG iterations (start with iter 2)
  for(;iter<MAXIT;++iter){

    if (verbose&&(rank==0)) {
      // some diagnostics
      printf("CG: it %d, r norm %12.12le, alpha = %le \n", iter, sqrt(g), alpha);
    }
    
    /* 
       M. Kronbichler  - "efficient matrix-free methods with deal.ii" CEED Annual Meeting 5, 2021

       update kernel:
       
       if(k even){
         x <= x + alpha*p + (alpha + alpha_old/beta_old)*p - (alpha_old/beta_old)*M\r;
       }
       
       r <= r - alpha*v;
       p <= M\r + beta*p;

       [ needs alpha and beta at first iteration (assuming k is odd at start) ]
    */

    int updatex = (iter%2);
    updatePPCGKernel(N, (int)updatex, alpha, beta, (dfloat)(alpha_old/beta_old), o_invM, o_p, o_r, o_v, o_x);
    
    /* 
       mat-vec 
       v = A*p
    */

    solver.Operator(o_p, o_v);
    
    /* 
       M. Kronbichler  - "efficient matrix-free methods with deal.ii" CEED Annual Meeting 5, 2021
       
       reduction kernel:
       a = p.v;
       b = r.v;
       c = v.v;
       d = r.(M\r);
       e = r.(M\v);
       f = v.(M\v);
       g = r.r;

    */

    ReductionsPPCG(o_r, &a, &b, &c, &d, &e, &f, &g);

    // TW: check this?
    alpha_old = alpha;
    if(a!=0){
      alpha = d/a;
    }else{
      printf("ppcg::Solve a=zero\n");
      exit(-1);
    }

    // TW check this
    beta_old = beta;
    beta = (d-2.*alpha*e+alpha*alpha*f)/d;

    if(fabs(g-2.*alpha*b+alpha*alpha*c)<TOL){ // notice comparing square
      printf("incrementing x in ppcg solve\n");
      linAlg.axpy(N, alpha, o_p, 1., o_x); // x += alpha*p
      break;
    }

    if (verbose&&(rank==0)) {
      if(g<0)
        printf("WARNING CG: rdotr = %17.15lf\n", g);
    }
  }

  return iter;
}

/*
 Kronbichler merged reduction kernel calling:

       a = p.v;
       b = r.v;
       c = v.v;
       d = r.(M\r);
       e = r.(M\v);
       f = v.(M\v);
       g = r.r;

*/

void ppcg::ReductionsPPCG(occa::memory &o_r,
			  dfloat *a, dfloat *b, dfloat *c,
			  dfloat *d, dfloat *e, dfloat *f,
			  dfloat *g){

  // invoke reductions kernel
  reductionsPPCGKernel(N, Nblocks, o_invM, o_r, o_p, o_v, o_reductionTmps);

  // copy partial sums to HOST
  o_reductionTmps.copyTo(reductionTmps);

  // initialize accumulators
  for(int fld=0;fld<PPCG_NREDUCTIONS;++fld){
    localTmps[fld] = 0;
  }

  // finalize local accumulations
  dlong id = 0;
  for(int n=0;n<Nblocks;++n){
    for(int fld=0;fld<PPCG_NREDUCTIONS;++fld){  // matches PPCG_NREDUCTIONS
      localTmps[fld] += reductionTmps[id++];
    }
  }

  // globalize accumulations
  MPI_Allreduce(localTmps, globalTmps, PPCG_NREDUCTIONS, MPI_DFLOAT, MPI_SUM, comm);
  
  *a = globalTmps[0];
  *b = globalTmps[1];
  *c = globalTmps[2];
  *d = globalTmps[3];
  *e = globalTmps[4];
  *f = globalTmps[5];
  *g = globalTmps[6];

}

ppcg::~ppcg() {
  //  updatePPCGKernel.free();
  //  reductionsPPCGKernel.free();
  o_p.free();
  o_v.free();
  o_z.free();
  o_Ax.free();
  o_reductionTmps.free();
  h_reductionTmps.free();

  free(localTmps);
  free(globalTmps);
}

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

#define PPCG_BLOCKSIZE 256
#define PPCG_NREDUCTIONS 7

ppcg::ppcg(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, terminator_t &_terminator, MPI_Comm _comm):
  linearSolver_t(_N, _Nhalo, _platform, _settings, _terminator, _comm) {

  dlong Ntotal = N;

  Nblocks = (_N+PPCG_BLOCKSIZE-1)/PPCG_BLOCKSIZE;

  /*aux variables */
  dfloat *dummy = (dfloat *) calloc(Ntotal,sizeof(dfloat)); //need this to avoid uninitialized memory warnings

  o_p  = platform.malloc(Ntotal*sizeof(dfloat),dummy);
  o_r  = platform.malloc(Ntotal*sizeof(dfloat),dummy);
  o_v  = platform.malloc(Ntotal*sizeof(dfloat),dummy);
  
  free(dummy);

  //pinned tmp buffer for reductions
  reductionTmps = (dfloat*) platform.hostMalloc(PPCG_NREDUCTIONS*PPCG_BLOCKSIZE*sizeof(dfloat),
					       NULL, h_reductionTmps);
  o_reductionTmps = platform.malloc(PPCG_NREDUCTIONS*PPCG_BLOCKSIZE*sizeof(dfloat));

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

int ppcg::Solve(solver_t& solver, precon_t& precon,
               occa::memory &o_x, occa::memory &o_r,
               const dfloat tol, const int MAXIT, const int verbose) {

  int rank;
  MPI_Comm_rank(comm, &rank);
  linAlg_t<dfloat> &linAlg = platform.linAlg;

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

  rdotr = linAlg.norm2(N, o_r, comm);
  rdotr = rdotr*rdotr;

  if (settings.compareSetting("LINEAR SOLVER STOPPING CRITERION", "ABS/REL-INITRESID")) {
    TOL = mymax(tol*tol*rdotr,tol*tol);
  }

  if (verbose&&(rank==0))
    printf("PPCG: initial res norm %12.12f \n", sqrt(rdotr));

  // reset terminator state for new solve
  terminator.reset();
  
  int iter;
  for(iter=0;iter<MAXIT;++iter){

    // Exit if tolerance is reached, taking at least one step.
    if ((iter == 0) && (g == 0.0)) break;

    if(iter>0)
      if(terminator.stopTest(o_x, g, TOL)){
	// need to check if iter is odd (or is that even)
	break;
      }

    int updatex = !(iter%2);

    /* 
       M. Kronbichler  - "efficient matrix-free methods with deal.ii" CEED Annual Meeting 5, 2021

       update kernel:
       
       if(k even){
         x <= x + alpha*p + (alpha + alpha_old/beta_old)*p - (alpha_old/beta_old)*M\r;
       }
       
       r <= r - alpha*v;
       p <= M\r + beta*p;
       
    */
    
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

    reductionsPPCG(N, o_invM, o_r, o_p, o_v, &a, &b, &c, &d, &e, &f, &g);

    alpha_old = alpha;
    alpha = d/a;


    if(fabs(g-2.*alpha*b+alpha*alpha*c)<TOL){ // notice comparing square
      linAlg.axpy(N, alpha, o_p, 1., o_x);
    }
    beta_old = beta;
    beta = (d-2.*alpha*e+alpha*apha*f)/d;

    if (verbose&&(rank==0)) {
      if(rdotr<0)
        printf("WARNING CG: rdotr = %17.15lf\n", rdotr);
      
      printf("CG: it %d, r norm %12.12le, alpha = %le \n", iter+1, sqrt(g), alpha);
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

void ppcg::ReductionsPPCG(const dlong N,
			  occa::memory &o_invM,
			  occa::memory &o_r,
			  occa::memory &o_p,
			  occa::memory &o_v,
			  dfloat *a, dfloat *b, dfloat *c,
			  dfloat *d, dfloat *e, dfloat *f,
			  dfloat *g){
  
  reductionsPPCGKernel(N, Nblocks, o_invM, o_r, o_p, o_v, o_reductionTmps);
  
  o_reductionTmps.copyTo(reductionTmps, PPCG_NREDUCTIONS*Nblocks*sizeof(dfloat));

  dlong id = 0;
  for(int fld=0;fld<PPCG_NREDUCTIONS;++fld){
    localTmps[fld] = 0;
  }
  
  for(int n=0;n<Nblocks;++n){
    for(int fld=0;fld<PPCG_NREDUCTIONS;++fld){  // matches PPCG_NREDUCTIONS
      localTmps[fld] += reductionTmps[id++];
    }
  }
  
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
  updatePPCGKernel.free();
}

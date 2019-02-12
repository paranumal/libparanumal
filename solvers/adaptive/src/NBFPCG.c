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

#include "adaptive.h"

void serialTic();
double serialElapsed();

// Non-blocking flexible preconditioned conjugate gradient
// "Pipelined, Flexible Krylov Subspace Methods", P. Sanan, S.M. Schnepp, D.A. May
// https://arxiv.org/pdf/1511.07226.pdf

int nbfpcg(adaptive_t* adaptive, dfloat lambda, 
	  occa::memory &o_r, occa::memory &o_x, 
	  const dfloat tol, const int MAXIT){
  
  mesh_t *mesh = adaptive->mesh;

  setupAide &options = adaptive->options;
  int verbose = options.compareArgs("VERBOSE", "TRUE");

  int fixedIterationCountFlag = 0;
  if(options.compareArgs("FIXED ITERATION COUNT", "TRUE")){
    fixedIterationCountFlag = 1;
  }

  // register scalars
  dfloat alpha0 = 0;
  dfloat beta0  = 0;
  dfloat gamma0 = 0;
  dfloat delta0 = 0;
  dfloat eta0   = 0;
  dfloat rdotr0 = 0;

  dfloat one = 1, zero = 0;
  
  /*aux variables */
  occa::memory &o_u  = adaptive->o_pcgWork[0];
  occa::memory &o_p  = adaptive->o_pcgWork[1];
  occa::memory &o_w  = adaptive->o_pcgWork[2];
  occa::memory &o_n  = adaptive->o_pcgWork[3];
  occa::memory &o_m  = adaptive->o_pcgWork[4];
  occa::memory &o_s  = adaptive->o_pcgWork[5];
  occa::memory &o_z  = adaptive->o_pcgWork[6];
  occa::memory &o_q  = adaptive->o_pcgWork[7];
  occa::memory &o_Ax = adaptive->o_pcgWork[8];

  MPI_Request request;
  MPI_Status  status;

  dfloat *localdots = (dfloat*) calloc(4, sizeof(dfloat));
  dfloat *globaldots = (dfloat*) calloc(4, sizeof(dfloat));

#if 0
  // blocks
  dfloat normB = adaptiveWeightedNorm2(adaptive, adaptive->o_invDegree, o_r);
  dfloat TOL = mymax(normB*tol*tol, tol*tol);
#endif
  
  // Ax = A*x
  adaptiveOperator(adaptive, lambda, o_x, o_Ax, dfloatString); 
  
  // subtract r = b - A*x
  adaptiveScaledAdd(adaptive, -one, o_Ax, one, o_r);

  // u = M*r [ Sanan notation ]
  adaptivePreconditioner(adaptive, lambda, o_r, o_u);
  
  // p = u
  o_p.copyFrom(o_u);

  // w = A*p
  adaptiveOperator(adaptive, lambda, o_p, o_w, dfloatString); 
  
  // gamma = u.r
  // delta = u.w
  // [
  adaptiveNonBlockingUpdate0NBFPCG(adaptive, o_u, o_r, o_w, localdots, globaldots, &request);

  adaptivePreconditioner(adaptive, lambda, o_w, o_m);
  
  adaptiveOperator(adaptive, lambda, o_m, o_n, dfloatString);

  o_s.copyFrom(o_w);
  o_q.copyFrom(o_m);
  o_z.copyFrom(o_n);

  MPI_Wait(&request, &status);
  gamma0 = globaldots[0]; // udotr
  delta0 = globaldots[1]; // udotw
  rdotr0 = globaldots[2]; // rdotr
  eta0   = delta0;
  alpha0 = gamma0/eta0;
  // ]

  dfloat TOL = mymax(rdotr0*tol*tol, tol*tol);
  
  int iter;

  for(iter=1;iter<=MAXIT;++iter){
    
    // x <= x + alpha*p
    // r <= r - alpha*s
    // u <= u - alpha*q
    // w <= w - alpha*z
    // gamma <= u.r
    // beta  <= -u.s/eta
    // delta <= u.w
    adaptiveNonBlockingUpdate1NBFPCG(adaptive, o_p, o_s, o_q, o_z, alpha0,
				     o_x, o_r, o_u, o_w, localdots, globaldots, &request);

    // n <= (w-r)
    //    adaptive->subKernel(mesh->Nelements*mesh->Np, o_w, o_r, o_n);
    adaptiveScaledAdd(adaptive,  one, o_w, zero, o_n);
    adaptiveScaledAdd(adaptive, -one, o_r,  one, o_n);

    // m <= M*(w-r)
    adaptivePreconditioner(adaptive, lambda, o_n, o_m);

    // m <= u + M*(w-r)
    adaptiveScaledAdd(adaptive, one, o_u, one, o_m);

    // n = A*m
    adaptiveOperator(adaptive, lambda, o_m, o_n, dfloatString);    

    // block for delta
    MPI_Wait(&request, &status);
    gamma0 = globaldots[0];       //  u.r
    beta0  = -globaldots[1]/eta0; // -u.s/eta
    delta0 = globaldots[2];       //  u.w
    rdotr0 = globaldots[3];       // r.r
    
    //  p <= u + beta*p
    adaptiveScaledAdd(adaptive, one, o_u, beta0, o_p);

    //  s <= w + beta*s
    adaptiveScaledAdd(adaptive, one, o_w, beta0, o_s);

    //  q <= m + beta*q
    adaptiveScaledAdd(adaptive, one, o_m, beta0, o_q);

    //  z <= n + beta*z
    adaptiveScaledAdd(adaptive, one, o_n, beta0, o_z);

    // eta = delta - beta^2*eta
    eta0 = delta0 - beta0*beta0*eta0;

    // alpha = gamma/eta
    alpha0 = gamma0/eta0;

    if (verbose&&(mesh->rank==0)) {

      if(gamma0<0)
	printf("WARNING CG: gamma = %17.15lf\n", gamma0);
      
      printf("CG: it %d z alpha = %12.12le beta = %le "
	     "rdotr = %le gamma = %le delta = %le, eta = %le\n",
	     iter, alpha0, beta0, rdotr0, gamma0, delta0, eta0);
    }

    if(rdotr0<=TOL && !fixedIterationCountFlag) break;
    
  }

  free(localdots);
  free(globaldots);
  
  return iter;
}



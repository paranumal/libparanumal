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

int pcg(adaptive_t* adaptive,
	dfloat lambda, 
        occa::memory &o_b,
	occa::memory &o_x, 
        const dfloat tol, const int MAXIT){

  level_t *level = adaptive->lvl;
  setupAide options = adaptive->options;

  int fixedIterationCountFlag = 0;
  int flexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");
  int verbose = options.compareArgs("VERBOSE", "TRUE");
  
  if(options.compareArgs("FIXED ITERATION COUNT", "TRUE"))
    fixedIterationCountFlag = 1;
  
  // register scalars
  dfloat rdotz1 = 0;
  dfloat rdotz2 = 0;

  // now initialized
  dfloat alpha = 0, beta = 0, pAp = 0;
  
  /*aux variables */
  occa::memory &o_p  = level->o_pcgWork[0];
  occa::memory &o_z  = level->o_pcgWork[1];
  occa::memory &o_Ap = level->o_pcgWork[2];
  occa::memory &o_Ax = level->o_pcgWork[3];
  occa::memory &o_r  = level->o_pcgWork[4];

  pAp = 0;
  rdotz1 = 1;

  dfloat rdotr0;

  // compute A*x
  adaptiveOperator(adaptive, level, lambda, o_x, o_Ax);

  // r <= b
  o_b.copyTo(o_r, level->Klocal*level->Np*sizeof(dfloat));
  
  // subtract r <= b - A*x
  adaptiveScaledAdd(adaptive, level, -1.f, o_Ax, 1.f, o_r);

  rdotr0 = adaptiveWeightedNorm2(adaptive, level, level->ogs->o_invDegree, o_r);

  dfloat TOL =  ASD_MAX(tol*tol*rdotr0,tol*tol);
  
  int iter;
  for(iter=1;iter<=MAXIT;++iter){

    // z = Precon^{-1} r
    adaptivePreconditioner(adaptive, lambda, o_r, o_z);
    
    rdotz2 = rdotz1;

    // r.z
    rdotz1 = adaptiveWeightedInnerProduct(adaptive, level, level->ogs->o_invDegree, o_r, o_z); 

    if(flexible){
      dfloat zdotAp;
      zdotAp = adaptiveWeightedInnerProduct(adaptive, level, level->ogs->o_invDegree, o_z, o_Ap);  

      beta = -alpha*zdotAp/rdotz2;
    }
    else{
      beta = (iter==1) ? 0:rdotz1/rdotz2;
    }

    // p = z + beta*p
    adaptiveScaledAdd(adaptive, level, 1.f, o_z, beta, o_p);
    
    // A*p
    adaptiveOperator(adaptive, level, lambda, o_p, o_Ap); 

    // dot(p,A*p)
    pAp =  adaptiveWeightedInnerProduct(adaptive, level, level->ogs->o_invDegree, o_p, o_Ap);

    alpha = rdotz1/pAp;

    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)
    
    dfloat rdotr = adaptiveUpdatePCG(adaptive, level, o_p, o_Ap, alpha, o_x, o_r);
	
    if (verbose&&(adaptive->rank==0)) {

      if(rdotr<0)
	printf("WARNING CG: rdotr = %17.15lf\n", rdotr);
      
      printf("CG: it %d r norm %12.12le alpha = %le \n", iter, sqrt(rdotr), alpha);    
    }
    
    if(rdotr<=TOL && !fixedIterationCountFlag) break;
    
  }

  return iter;
}



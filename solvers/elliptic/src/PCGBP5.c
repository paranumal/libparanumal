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

#define TIMER 0

#include "elliptic.h"

#if 1
void pcgBP5(elliptic_t* elliptic, dfloat lambda, 
        occa::memory &o_r, occa::memory &o_x, const int MAXIT) {

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  /*aux variables */
  occa::memory &o_p  = elliptic->o_p;
  occa::memory &o_z  = elliptic->o_z;
  occa::memory &o_Ap = elliptic->o_Ap;
  occa::memory &o_Ax = elliptic->o_Ax;

  // Precon^{-1} (b-A*x)
  ellipticPreconditioner(elliptic, lambda, o_r, o_z);

  // p = z
  o_p.copyFrom(o_z); // PCG

  // dot(r,z)
  dfloat rdotz0 = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  dfloat rdotr0 = 0; 
  dfloat alpha, beta, pAp = 0;

  for(int iter=0; iter<MAXIT; iter++) {

    // [
    // A*p
    ellipticOperator(elliptic, lambda, o_p, o_Ap, dfloatString);

    // dot(p,A*p)
    pAp =  ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_p, o_Ap);
    // ]
    
    // alpha = dot(r,z)/dot(p,A*p)
    alpha = rdotz0/pAp;

    // x <= x + alpha*p
    ellipticScaledAdd(elliptic,  alpha, o_p,  1.f, o_x);

#if TIMER
    occaTimerTic(mesh->device,"Residual update");
#endif
    // [
    // r <= r - alpha*A*p
    ellipticScaledAdd(elliptic, -alpha, o_Ap, 1.f, o_r);

    // dot(r,r)
#if 1
    rdotr1 = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_r);
#else
    rdotr1 = ellipticWeightedNorm2(elliptic, elliptic->o_invDegree, o_r);
#endif
    // ]

#if TIMER
  occaTimerToc(mesh->device,"Residual update");
#endif



#if TIMER
    occaTimerTic(mesh->device,"Preconditioner");
#endif

    // [
    // z = Precon^{-1} r
    ellipticPreconditioner(elliptic, lambda, o_r, o_z);

    // dot(r,z)
    rdotz1 = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
    // ]
    
    // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
    // if(options.compareArgs("KRYLOV SOLVER", "PCG+FLEXIBLE") ||
       // options.compareArgs("KRYLOV SOLVER", "PCG,FLEXIBLE")) {
      // dfloat zdotAp = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_z, o_Ap);
      // beta = -alpha*zdotAp/rdotz0;
    // } else {
      beta = rdotz1/rdotz0;
    // }

    // p = z + beta*p
    ellipticScaledAdd(elliptic, 1.f, o_z, beta, o_p);

    // switch rdotz0 <= rdotz1
    rdotz0 = rdotz1;

#if TIMER
    occaTimerToc(mesh->device,"Preconditioner");
#endif

    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;
  }
}



#else

void pcgBP5(elliptic_t* elliptic, dfloat lambda, 
        occa::memory &o_r, occa::memory &o_x, const int MAXIT) {

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int DEBUG_ENABLE_REDUCTIONS = 1;
  options.getArgs("DEBUG ENABLE REDUCTIONS", DEBUG_ENABLE_REDUCTIONS);
  
  // register scalars
  dfloat rdotz0 = 0;
  dfloat rdotr0 = 0;
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  dfloat zdotAp = 0;
  
  dfloat alpha, beta, pAp = 0;
  
  /*aux variables */
  occa::memory &o_p  = elliptic->o_p;
  occa::memory &o_z  = elliptic->o_z;
  occa::memory &o_Ap = elliptic->o_Ap;
  occa::memory &o_Ax = elliptic->o_Ax;

  // Precon^{-1} (b-A*x)
  ellipticPreconditioner(elliptic, lambda, o_r, o_z);

  // p = z
  o_p.copyFrom(o_z); // PCG

  // dot(r,z)
#if 0
  rdotz0 = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
#else
  rdotz0 = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
#endif

  for(int iter=0; iter<MAXIT; iter++) {

    // [
    // A*p
    ellipticOperator(elliptic, lambda, o_p, o_Ap, dfloatString);
    
    // dot(p,A*p)
    if(DEBUG_ENABLE_REDUCTIONS==1){
#if 0
      pAp =  ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_p, o_Ap);
#else
      pAp =  ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_p, o_Ap);
#endif
    }
    else
      pAp = 1;
    // ]
    
    // alpha = dot(r,z)/dot(p,A*p)
    alpha = rdotz0/pAp;

    // TO DO:
    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)
    //
    rdotr1 = ellipticUpdatePCG(elliptic, o_p, o_Ap, alpha, o_x, o_r);
    
    //    if (options.compareArgs("VERBOSE", "TRUE")&&(mesh->rank==0)) 
    //  printf("CG: it %d r norm %12.12f alpha = %f \n",Niter, sqrt(rdotr1), alpha);

    // if(rdotr1 < TOL) {
    //   rdotr0 = rdotr1;
    //  break;
    // }

    // [
    // z = Precon^{-1} r
    ellipticPreconditioner(elliptic, lambda, o_r, o_z);

    // dot(r,z)
    if(DEBUG_ENABLE_REDUCTIONS==1){
#if 0
      rdotz1 = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
#else
      rdotz1 = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
#endif
    }
    else
      rdotz1 = 1;
    
    // ]
    
    // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
    //if(options.compareArgs("KRYLOV SOLVER", "PCG+FLEXIBLE") ||
    //  options.compareArgs("KRYLOV SOLVER", "PCG,FLEXIBLE")) {
    
    //      if(DEBUG_ENABLE_REDUCTIONS==1){
    // #if 0
    //	zdotAp = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_z, o_Ap);
    // #else
    //  zdotAp = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_z, o_Ap);
    //#endif
    //  }
    //  else
    //	zdotAp = 1;
      
    //  beta = -alpha*zdotAp/rdotz0;
    // } else {
      beta = rdotz1/rdotz0;
      // }

    // p = z + beta*p
    ellipticScaledAdd(elliptic, 1.f, o_z, beta, o_p);

    // switch rdotz0 <= rdotz1
    rdotz0 = rdotz1;

    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;
  }
}

#if 0
dfloat ellipticUpdatePCG(elliptic_t *elliptic,
			 occa::memory &o_p, occa::memory &o_Ap, dfloat alpha,
			 occa::memory &o_x, occa::memory &o_r){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int DEBUG_ENABLE_REDUCTIONS = 1;
  options.getArgs("DEBUG ENABLE REDUCTIONS", DEBUG_ENABLE_REDUCTIONS);
  
  dfloat rdotr1 = 0;
  
  if(!options.compareArgs("DISCRETIZATION", "CONTINUOUS")){
    
    // x <= x + alpha*p
    ellipticScaledAdd(elliptic,  alpha, o_p,  1.f, o_x);
    
    // [
    // r <= r - alpha*A*p
    ellipticScaledAdd(elliptic, -alpha, o_Ap, 1.f, o_r);
    
    // dot(r,r)
    if(DEBUG_ENABLE_REDUCTIONS==1){
#if 0
      rdotr1 = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_r);
#else
      rdotr1 = ellipticWeightedNorm2(elliptic, elliptic->o_invDegree, o_r);
#endif
    }
    else
      rdotr1 = 1;
  }else{
    
    // x <= x + alpha*p
    // r <= r - alpha*A*p
    // dot(r,r)
    elliptic->updatePCGKernel(mesh->Nelements*mesh->Np, elliptic->NblocksUpdatePCG,
			      elliptic->o_invDegree, o_p, o_Ap, alpha, o_x, o_r, elliptic->o_tmpNormr);

    elliptic->o_tmpNormr.copyTo(elliptic->tmpNormr);

    rdotr1 = 0;
    for(int n=0;n<elliptic->NblocksUpdatePCG;++n){
      rdotr1 += elliptic->tmpNormr[n];
    }

    
    dfloat globalrdotr1 = 0;
    MPI_Allreduce(&rdotr1, &globalrdotr1, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);


    rdotr1 = globalrdotr1;
    
  }

  return rdotr1;
}

#endif

#endif

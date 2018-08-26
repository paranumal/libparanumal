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

#include "elliptic.h"

int pcg(elliptic_t* elliptic, dfloat lambda, 
        occa::memory &o_r, occa::memory &o_x, 
        const dfloat tol, const int MAXIT) {

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  /*aux variables */
  occa::memory &o_p  = elliptic->o_p;
  occa::memory &o_z  = elliptic->o_z;
  occa::memory &o_Ap = elliptic->o_Ap;
  occa::memory &o_Ax = elliptic->o_Ax;

  /*compute norm b, set the tolerance */
#if 1
  dfloat normB = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_r);
#else
  dfloat normB = ellipticWeightedNorm2(elliptic, elliptic->o_invDegree, o_r);
#endif

  dfloat TOL =  mymax(tol*tol*normB,tol*tol);
  // compute A*x
  ellipticOperator(elliptic, lambda, o_x, elliptic->o_Ax, dfloatString);

  // subtract r = b - A*x
  ellipticScaledAdd(elliptic, -1.f, o_Ax, 1.f, o_r);

#if 1
  dfloat rdotr0 = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_r);
#else
  dfloat rdotr0 = ellipticWeightedNorm2(elliptic, elliptic->o_invDegree, o_r);
#endif

  dfloat rdotz0 = 0;
  int Niter = 0;

  //sanity check
  if (rdotr0<1E-20) {
    if (options.compareArgs("VERBOSE", "TRUE")&&(mesh->rank==0)){
      printf("converged in ZERO iterations. Stopping.\n");}
    return 0;
  } 

  if (options.compareArgs("VERBOSE", "TRUE")&&(mesh->rank==0)) 
    printf("CG: initial res norm %12.12f WE NEED TO GET TO %12.12f \n", sqrt(rdotr0), sqrt(TOL));

  // Precon^{-1} (b-A*x)
  ellipticPreconditioner(elliptic, lambda, o_r, o_z);

  // p = z
  o_p.copyFrom(o_z); // PCG

  // dot(r,z)
  rdotz0 = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;

  dfloat alpha, beta, pAp = 0;

  while((Niter <MAXIT)) {

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

    occaTimerTic(mesh->device,"Residual update");
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
    occaTimerToc(mesh->device,"Residual update");
    
    if (options.compareArgs("VERBOSE", "TRUE")&&(mesh->rank==0)) 
      printf("CG: it %d r norm %12.12f alpha = %f \n",Niter, sqrt(rdotr1), alpha);

    if(rdotr1 < TOL) {
      rdotr0 = rdotr1;
      break;
    }

    occaTimerTic(mesh->device,"Preconditioner");

    // [
    // z = Precon^{-1} r
    ellipticPreconditioner(elliptic, lambda, o_r, o_z);

    // dot(r,z)
    rdotz1 = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
    // ]
    
    // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
    if(options.compareArgs("KRYLOV SOLVER", "PCG+FLEXIBLE") ||
       options.compareArgs("KRYLOV SOLVER", "PCG,FLEXIBLE")) {
      dfloat zdotAp = ellipticCascadingWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_z, o_Ap);
      beta = -alpha*zdotAp/rdotz0;
    } else {
      beta = rdotz1/rdotz0;
    }

    // p = z + beta*p
    ellipticScaledAdd(elliptic, 1.f, o_z, beta, o_p);

    // switch rdotz0 <= rdotz1
    rdotz0 = rdotz1;

    occaTimerToc(mesh->device,"Preconditioner");

    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;

    ++Niter;
  }
  return Niter;
}

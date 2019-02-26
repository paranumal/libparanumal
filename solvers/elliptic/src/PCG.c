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

#if 1
// FROM NEKBONE: not appropriate since it assumes zero initial data
int pcg(elliptic_t* elliptic, dfloat lambda,
        occa::memory &o_r, occa::memory &o_x,
        const dfloat tol, const int MAXIT){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int fixedIterationCountFlag = 0;
  int enableGatherScatters = 1;
  int enableReductions = 1;
  int flexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");
  int verbose = options.compareArgs("VERBOSE", "TRUE");

  // options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);
  // options.getArgs("DEBUG ENABLE OGS", enableGatherScatters);
  // if(options.compareArgs("FIXED ITERATION COUNT", "TRUE"))
  //   fixedIterationCountFlag = 1;

  // register scalars
  dfloat rdotz1 = 0;
  dfloat rdotz2 = 0;


  // now initialized
  dfloat alpha = 0, beta = 0, pAp = 0;

  /*aux variables */
  occa::memory &o_p  = elliptic->o_p;
  occa::memory &o_z  = elliptic->o_z;
  occa::memory &o_Ap = elliptic->o_Ap;
  occa::memory &o_Ax = elliptic->o_Ax;

#if 0
  dfloat normB = ellipticWeightedNorm2(elliptic, elliptic->o_invDegree, o_r);
#endif

  pAp = 0;
  rdotz1 = 1;

  dfloat rdotr0;

  // compute A*x
  ellipticOperator(elliptic, lambda, o_x, elliptic->o_Ax, dfloatString);

  // subtract r = b - A*x
  ellipticScaledAdd(elliptic, -1.f, o_Ax, 1.f, o_r);

  if(enableReductions)
    rdotr0 = ellipticWeightedNorm2(elliptic, elliptic->o_invDegree, o_r);
  else
    rdotr0 = 1;

  dfloat TOL =  mymax(tol*tol*rdotr0,tol*tol);

  int iter;
  for(iter=1;iter<=MAXIT;++iter){

    // z = Precon^{-1} r
    ellipticPreconditioner(elliptic, lambda, o_r, o_z);

    rdotz2 = rdotz1;

    // r.z
    if(enableReductions){
      rdotz1 = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
    }
    else
      rdotz1 = 1;

    if(flexible){
      dfloat zdotAp;
      if(enableReductions){
	zdotAp = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_z, o_Ap);
      }
      else
	zdotAp = 1;

      beta = -alpha*zdotAp/rdotz2;
    }
    else{
      beta = (iter==1) ? 0:rdotz1/rdotz2;
    }

    // p = z + beta*p
    ellipticScaledAdd(elliptic, 1.f, o_z, beta, o_p);

    // A*p
    ellipticOperator(elliptic, lambda, o_p, o_Ap, dfloatString);

    // dot(p,A*p)
    if(enableReductions){
      pAp =  ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_p, o_Ap);
    }
    else
      pAp = 1;

    alpha = rdotz1/pAp;

    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)

    dfloat rdotr = ellipticUpdatePCG(elliptic, o_p, o_Ap, alpha, o_x, o_r);

    if (verbose&&(mesh->rank==0)) {

      if(rdotr<0)
	printf("WARNING CG: rdotr = %17.15lf\n", rdotr);

      printf("CG: it %d r norm %12.12le alpha = %le \n", iter, sqrt(rdotr), alpha);
    }

    if(rdotr<=TOL && !fixedIterationCountFlag) break;

  }

  return iter;
}


#else

// original pcg
int pcg(elliptic_t* elliptic, dfloat lambda,
        occa::memory &o_r, occa::memory &o_x,
        const dfloat tol, const int MAXIT) {

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int fixedIterationCountFlag = 0;
  int enableGatherScatters = 1;
  int enableReductions = 1;
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int flexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");
  int verbose = options.compareArgs("VERBOSE", "TRUE");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");

  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);
  options.getArgs("DEBUG ENABLE OGS", enableGatherScatters);
  if(options.compareArgs("FIXED ITERATION COUNT", "TRUE")){
    fixedIterationCountFlag = 1;
  }


  // register scalars
  dfloat rdotz0 = 0;
  dfloat rdotr0 = 0;
  int Niter = 0;
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  dfloat zdotAp = 0;

  dfloat alpha, beta, pAp = 0;
  dfloat TOL, normB;

  /*aux variables */
  occa::memory &o_p  = elliptic->o_p;
  occa::memory &o_z  = elliptic->o_z;
  occa::memory &o_Ap = elliptic->o_Ap;
  occa::memory &o_Ax = elliptic->o_Ax;

  /*compute norm b, set the tolerance */
  normB = ellipticWeightedNorm2(elliptic, elliptic->o_invDegree, o_r);

  TOL =  mymax(tol*tol*normB,tol*tol);

  // compute A*x
  ellipticOperator(elliptic, lambda, o_x, elliptic->o_Ax, dfloatString);

  // subtract r = b - A*x
  ellipticScaledAdd(elliptic, -1.f, o_Ax, 1.f, o_r);

  if(enableReductions)
    rdotr0 = ellipticWeightedNorm2(elliptic, elliptic->o_invDegree, o_r);
  else
    rdotr0 = 1;

  //sanity check
  if (rdotr0<1E-20) {
    if (verbose&&(mesh->rank==0)){
      printf("converged in ZERO iterations. Stopping.\n");}
    return 0;
  }

  if (verbose&&(mesh->rank==0))
    printf("CG: initial res norm %12.12f WE NEED TO GET TO %12.12f \n", sqrt(rdotr0), sqrt(TOL));

  // Precon^{-1} (b-A*x)
  ellipticPreconditioner(elliptic, lambda, o_r, o_z);

  // p = z
  o_p.copyFrom(o_z); // PCG

  // dot(r,z)
  if(enableReductions)
    rdotz0 = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
  else
    rdotz0 = 1;

  while((Niter <MAXIT)) {

    // A*p
    ellipticOperator(elliptic, lambda, o_p, o_Ap, dfloatString);

    // dot(p,A*p)
    if(enableReductions)
      pAp =  ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_p, o_Ap);
    else
      pAp = 1;

    // alpha = dot(r,z)/dot(p,A*p)
    alpha = rdotz0/pAp;

    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)
    rdotr1 = ellipticUpdatePCG(elliptic, o_p, o_Ap, alpha, o_x, o_r);

    if (verbose&&(mesh->rank==0))
      printf("CG: it %d r norm %12.12le alpha = %le \n",Niter, sqrt(rdotr1), alpha);

    if(rdotr1 < TOL) {
      rdotr0 = rdotr1;
      break;
    }

    // z = Precon^{-1} r
    ellipticPreconditioner(elliptic, lambda, o_r, o_z);

    // dot(r,z)
    if(enableReductions)
      rdotz1 = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
    else
      rdotz1 = 1;

    // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
    if(flexible){

      if(enableReductions)
	zdotAp = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_z, o_Ap);
      else
	zdotAp = 1;

      beta = -alpha*zdotAp/rdotz0;
    } else {
      beta = rdotz1/rdotz0;
    }

    // p = z + beta*p
    ellipticScaledAdd(elliptic, 1.f, o_z, beta, o_p);

    // switch rdotz0 <= rdotz1
    rdotz0 = rdotz1;

    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;

    ++Niter;
  }

  return Niter;
}

#endif

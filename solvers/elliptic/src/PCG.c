#include "elliptic.h"

int pcg(elliptic_t* elliptic, dfloat lambda, 
        occa::memory &o_r, occa::memory &o_x, 
        const dfloat tol, const int MAXIT) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  /*aux variables */
  occa::memory &o_p  = elliptic->o_p;
  occa::memory &o_z  = elliptic->o_z;
  occa::memory &o_Ap = elliptic->o_Ap;
  occa::memory &o_Ax = elliptic->o_Ax;

  /*compute norm b, set the tolerance */
  dfloat normB = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_r);
  dfloat TOL =  mymax(tol*tol*normB,tol*tol);
  // compute A*x
  ellipticOperator(elliptic, lambda, o_x, elliptic->o_Ax);

  // subtract r = b - A*x
  ellipticScaledAdd(elliptic, -1.f, o_Ax, 1.f, o_r);

  dfloat rdotr0 = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_r);

  dfloat rdotz0 = 0;
  int Niter = 0;

  //sanity check
  if (rdotr0<1E-20) {
    if (options.compareArgs("VERBOSE", "TRUE")&&(rank==0)){
      printf("converged in ZERO iterations. Stopping.\n");}
    return 0;
  } 

  if (options.compareArgs("VERBOSE", "TRUE")&&(rank==0)) 
    printf("CG: initial res norm %12.12f WE NEED TO GET TO %12.12f \n", sqrt(rdotr0), sqrt(TOL));

  // Precon^{-1} (b-A*x)
  ellipticPreconditioner(elliptic, lambda, o_r, o_z);

  // p = z
  o_p.copyFrom(o_z); // PCG

  // dot(r,z)
  rdotz0 = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;

  dfloat alpha, beta, pAp = 0;

  while((Niter <MAXIT)) {

    // [
    // A*p
    ellipticOperator(elliptic, lambda, o_p, o_Ap);

    // dot(p,A*p)
    pAp =  ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree,o_p, o_Ap);
    // ]
    
    // alpha = dot(r,z)/dot(p,A*p)
    alpha = rdotz0/pAp;

    // x <= x + alpha*p
    ellipticScaledAdd(elliptic,  alpha, o_p,  1.f, o_x);

    // [
    // r <= r - alpha*A*p
    ellipticScaledAdd(elliptic, -alpha, o_Ap, 1.f, o_r);

    // dot(r,r)
    rdotr1 = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_r);
    // ]
    
    if (options.compareArgs("VERBOSE", "TRUE")&&(rank==0)) 
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
    rdotz1 = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_r, o_z);
    // ]
    
    // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
    if(options.compareArgs("KRYLOV SOLVER", "PCG+FLEXIBLE") ||
       options.compareArgs("KRYLOV SOLVER", "PCG,FLEXIBLE")) {
      dfloat zdotAp = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_z, o_Ap);
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

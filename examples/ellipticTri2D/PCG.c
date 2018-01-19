#include "ellipticTri2D.h"

int pcg(solver_t* solver, const char* options, dfloat lambda, 
    occa::memory &o_r, occa::memory &o_x, 
    const dfloat tol, const int MAXIT) {

  iint rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = solver->mesh;

  /*aux variables */
  occa::memory &o_p  = solver->o_p;
  occa::memory &o_z  = solver->o_z;
  occa::memory &o_Ap = solver->o_Ap;
  occa::memory &o_Ax = solver->o_Ax;

  /*compute norm b, set the tolerance */
  dfloat normB = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
  dfloat TOL =  tol*tol*normB;
  // compute A*x
  ellipticOperator2D(solver, lambda, o_x, solver->o_Ax, options);

  // subtract r = b - A*x
  ellipticScaledAdd(solver, -1.f, o_Ax, 1.f, o_r);

  dfloat rdotr0 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);

  dfloat rdotz0 = 0;
  int Niter = 0;

  //sanity check
  if (rdotr0<=(TOL)) {
    if (strstr(options, "VERBOSE")&&(rank==0)){
      printf("converged in ZERO iterations. Stopping.\n");}
    return 0;
  } 

  if (strstr(options,"VERBOSE")&&(rank==0)) 
    printf("CG: initial res norm^2 %12.12f WE NEED TO GET TO %12.12f \n", rdotr0, TOL);

  // Precon^{-1} (b-A*x)
  ellipticPreconditioner2D(solver, lambda, o_r, o_z, options);

  // p = z
  o_p.copyFrom(o_z); // PCG

  // dot(r,z)
  rdotz0 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_z, options);
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;

  dfloat alpha, beta, pAp = 0;

  while((rdotr0>TOL) && (Niter <MAXIT)) {

    // A*p
    ellipticOperator2D(solver, lambda, o_p, o_Ap, options);

    // dot(p,A*p)
    pAp =  ellipticWeightedInnerProduct(solver, solver->o_invDegree,o_p, o_Ap, options);

    // alpha = dot(r,z)/dot(p,A*p)
    alpha = rdotz0/pAp;

    // x <= x + alpha*p
    ellipticScaledAdd(solver,  alpha, o_p,  1.f, o_x);

    // r <= r - alpha*A*p
    ellipticScaledAdd(solver, -alpha, o_Ap, 1.f, o_r);

    // dot(r,r)
    rdotr1 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);

    if (strstr(options,"VERBOSE")&&(rank==0)) 
      printf("CG: it %d r norm^2 %12.12f alpha = %f \n",Niter, rdotr1, alpha);


    if(rdotr1 < TOL) {
      rdotr0 = rdotr1;
      break;
    }

    occaTimerTic(mesh->device,"Preconditioner");

    // z = Precon^{-1} r
    ellipticPreconditioner2D(solver, lambda, o_r, o_z, options);

    // dot(r,z)
    rdotz1 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_z, options);

    // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
    if(strstr(options,"FLEXIBLE")) {
      dfloat zdotAp = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_z, o_Ap, options);
      beta = -alpha*zdotAp/rdotz0;
    } else {
      beta = rdotz1/rdotz0;
    }

    // p = z + beta*p
    ellipticScaledAdd(solver, 1.f, o_z, beta, o_p);

    // switch rdotz0 <= rdotz1
    rdotz0 = rdotz1;

    occaTimerToc(mesh->device,"Preconditioner");

    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;

    ++Niter;
  }
  return Niter;
}
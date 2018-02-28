#include "ellipticTri2D.h"

//written by KS
// This code is based on Henk A. van Der Worst book "Iterative Krylov methods for large linear systems"
int pbicgstab(solver_t* solver, const char* options, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const dfloat tol, const int MAXIT) {

  // aux variables
  dfloat bnorm, snorm, rnorm;
  dfloat omega, alpha, beta=0;
  dfloat res;
  dfloat rho0, rho00;
  dfloat aux, aux2;

  occa::memory &o_Ax = solver->o_Ax;
  occa::memory &o_p  = solver->o_p;
  occa::memory &o_s  = solver->o_z;
  occa::memory &o_t = solver->o_rtmp;
  occa::memory &o_rtilde = solver->o_Ap;
  occa::memory &o_phat = solver->o_res;
  occa::memory &o_shat = solver->o_Sres;

  if (strstr(options,"LEFT")){
    ellipticPreconditioner2D(solver, lambda, o_r, o_s , options);
    o_r.copyFrom(o_s);
  }


  /* set the tolerance */
  bnorm = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
  dfloat tolrel =  tol*tol*bnorm;

  /* compute the initial residual */

  // compute A*x
  ellipticOperator2D(solver, lambda, o_x, o_Ax, options);

  // subtract r = b - A*x
  ellipticScaledAdd(solver, -1.f, o_Ax, 1.f, o_r);

  rnorm = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
  printf("initial res norm %16.18f tol rel = %15.18f\n", bnorm, tolrel); 
  if (rnorm<tolrel) {
    if (strstr(options, "VERBOSE")){    
      printf("BiCGstab: converged in 0 iterations \n");}
    return 0;
  }
  o_rtilde.copyFrom(o_r);

  int conv = 0;
  int error = 0;
  o_rtilde.copyFrom(o_r);

  if (strstr(options,"VERBOSE")) 
    printf("BiCGStab: initial res norm^2 %12.12f WE NEED TO GET TO %12.12f \n", rnorm, tolrel);

  int it = 0;
  while ((!conv)&&(!error)&&(it<MAXIT)) {
    rho0 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_rtilde, o_r, options);
    if (rho0 == 0) {
      conv = 0;
      error = -1;
    }

    if (it == 0) {
      o_p.copyFrom(o_r);
    } else {
      beta = (rho0/rho00) * (alpha/omega);

      /* p = r + beta * (p - omega * o_Ax);*/
      // ellipticScaledAdd(solver, -1.f*omega, o_Ax, 1.0, o_p);
      // ellipticScaledAdd(solver, 1.f, o_r, beta, o_p);
      ellipticScaledAdd(solver, -1.f*omega, o_Ax, beta, o_p);
      ellipticScaledAdd(solver, 1.f, o_r, 1.0f, o_p);


    }

    if (strstr(options, "LEFT")){
      ellipticOperator2D(solver, lambda, o_p, o_phat, options);
      ellipticPreconditioner2D(solver, lambda, o_phat, o_Ax, options);

    }
    else{
      ellipticPreconditioner2D(solver, lambda, o_p, o_phat, options);
      ellipticOperator2D(solver, lambda, o_phat, o_Ax, options);

      // o_phat.copyFrom(o_p);
    }    

    /*alpha = rho0/(rtilde^T*Ax)*/
    aux = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_rtilde, o_Ax, options);
    alpha = rho0/aux;
    /* s = r - alpha*Ax */
    o_s.copyFrom(o_r);
    ellipticScaledAdd(solver, -1.0*alpha, o_Ax, 1.0, o_s);
    //check for early convergence
    snorm = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_s, o_s, options);

    if (strstr(options,"VERBOSE")) 
      printf("BiCGStab: it %d sres norm^2 %12.12f alpha = %f beta = %f omega = %f  \n",it, snorm, alpha, beta, omega);

    if (snorm<tolrel) {
      /* o_x = o_x + alpha*o_p*/
      ellipticScaledAdd(solver, alpha, o_phat, 1.0f, o_x);
      //resid = ||s|| / bnorm;
      res = snorm/bnorm;
      res = sqrt(res);
      conv = 1;
      error = 1; // early convergence
    } else {
      /* t = A*s */ 
      if(strstr(options, "LEFT")){
        ellipticOperator2D(solver, lambda, o_s, o_shat, options);
        ellipticPreconditioner2D(solver, lambda, o_shat, o_t, options);

      }
      else{
        ellipticPreconditioner2D(solver, lambda, o_s, o_shat, options);
        ellipticOperator2D(solver, lambda, o_shat, o_t, options);

        //  o_shat.copyFrom(o_s);
      }

      /*omega = s^Tt / t^Tt*/
      aux = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_s, o_t, options);
      aux2 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_t, o_t, options);
      omega = aux/aux2;

      /* update x */
      /*x = x+alpha p*/
      ellipticScaledAdd(solver, alpha, o_phat, 1.0f, o_x);
      /* x = x + omega s*/
      ellipticScaledAdd(solver, omega, o_shat, 1.0f, o_x);
      /*update r; r = s - omega t */
      o_r.copyFrom(o_s);
      ellipticScaledAdd(solver, -1.0*omega, o_t, 1.0f, o_r);

      //check for convergence
      rnorm = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
      if (strstr(options,"VERBOSE")) 
        printf("BiCGStab: it %d rres norm^2 %12.12f \n",it, rnorm);

      if (rnorm<tolrel) {
        rnorm = sqrt(rnorm);
        conv = 1;
      } else {
        if (omega == 0) {
          error = 2;
        }
      }         
    }

    it++;
    rho00= rho0;
  }

  //why did we quit
  if (conv) {
    if (strstr(options, "RIGHT")){
      // ellipticPreconditioner2D(solver, lambda, o_x, o_Ax, options);
      // o_x.copyFrom(o_Ax);
    }

    if (strstr(options,"VERBOSE"))  
      printf("BiCGstab: converged in %d iterations \n", it);

    if (error == 1) {
      if (strstr(options,"VERBOSE"))  {
        printf("actually, converged in %d .5 steps \n", it-1);
      }      
    }

    if (error==2) {
      if (strstr(options,"VERBOSE"))  {
        printf("BiCGstab: quited because of omega == 0\n");
      }
    }
    if (error==-1) {
      if (strstr(options,"VERBOSE"))  {
        printf("BiCGstab: failed to converge because r perp to rtilde\n");
      }
      it = -1;
      return -1;
    }
  }

  return it;
}

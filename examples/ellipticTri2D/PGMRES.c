//gmres(m)
//set m in options
//if no m, default is 10
//m is the restart frequency
// based on Henk van der Vorst book, p.67
// and cuda_itsol package by Li and Saad
// this is a space-saving version by KS
// in a standard version by Saad we need to keep an additional space of
// vectors w_{k+1} = M*w_k
#include "ellipticTri2D.h"

#define IDX2C(i,j,ld) (((i)*(ld))+(j))
#define ZERO 0.0
#define EPSILON   1.0e-18
#define EPSMAC    1.0e-16

int pgmresm(solver_t* solver, const char* options, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const dfloat tol, const int MAXIT) {
  
  mesh2D *mesh = solver->mesh;

  occa::memory &o_w = solver->o_Ax;
  occa::memory &o_Ax  = solver->o_Ap;
  occa::memory &o_b = solver->o_rtmp;

  occa::memory o_invDegree = solver->o_invDegree;
  occa::memory o_z = solver->o_z;
  occa::memory *o_V = solver->o_V;

  dfloat *HH = solver->HH;
  int m = solver->GMRESrestartFreq;

  dfloat * c = (dfloat*) calloc(m,   sizeof(dfloat));
  dfloat * s = (dfloat*) calloc(m,   sizeof(dfloat));
  dfloat * rs = (dfloat*) calloc(m+1,   sizeof(dfloat));

  dfloat beta, t, gam, ro;
  int k1;


  /*multiply by M if left preconditioning */
  if (strstr(options,"LEFT")){
    ellipticPreconditioner2D(solver, lambda, o_r, o_z , options);
    o_r.copyFrom(o_z);
  }
  
  /* set the tolerance */
  dfloat bnorm = ellipticWeightedInnerProduct(solver, o_invDegree, o_r, o_r, options);
  bnorm = sqrt(bnorm);

  dfloat tolrel =  tol*tol*bnorm;
  
  /*start iterating*/
  int conv = 0;
  int Niter = -1;
  o_b.copyFrom(o_r);
  while ((!conv)&&(Niter<MAXIT)) {

    /* compute the initial residual */

    // compute A*(M*x) NOTE: DIFFERENT than in standard implementation
    //since intial guess is 0, and you cannot call the preconditioner on zero vector, just DONT DO IT
    if (Niter!=-1) {
      if (strstr(options, "RIGHT")){
        ellipticPreconditioner2D(solver, lambda, o_x, o_w, options);
        ellipticOperator2D(solver, lambda, o_w, o_Ax, options);
      } else{
        ellipticOperator2D(solver, lambda, o_x, o_w, options);
        ellipticPreconditioner2D(solver, lambda, o_w, o_Ax, options);
      }
    }
  
    // subtract r = b - A*(M*x);
    o_r.copyFrom(o_b);
    ellipticScaledAdd(solver, -1.f, o_Ax, 1.f, o_r);

    /*beta  = ||r|| */
    beta = ellipticWeightedInnerProduct(solver, o_invDegree, o_r, o_r, options);
    beta = sqrt(beta);

    /*v(:,0) = r/beta */
    dfloat zero = 0.0f;
    ellipticScaledAdd(solver, 1./beta, o_r, zero , o_V[0]);

    /*hessenberg matrix */
    rs[0] = beta;
    /*z = M*v, w = A*z */

    /*inner cycle */
    int NinnerIt =-1;
    if (strstr(options,"VERBOSE")) printf("GMRES Restarted\n");
    for (int k=0; k< m; ++k) {
      NinnerIt++;
      Niter++;

      if (strstr(options, "RIGHT")){
        ellipticPreconditioner2D(solver, lambda, o_V[k], o_w, options);
        ellipticOperator2D(solver, lambda, o_w, o_V[k+1], options);
      }else{
        ellipticOperator2D(solver, lambda, o_V[k], o_w, options);
        ellipticPreconditioner2D(solver, lambda, o_w, o_V[k+1], options);
      }

      /* Gram-Schmidt */
      for (int j=0; j<= k; ++j) {
        //i->k, j ->j
        HH[IDX2C(k,j,m+1)] = ellipticWeightedInnerProduct(solver, o_invDegree, o_V[(k+1)], o_V[(j)], options);
        ellipticScaledAdd(solver, -1.f*HH[IDX2C(k,j,m+1)], o_V[(j)], 1, o_V[(k+1)]);
      }
      t = ellipticWeightedInnerProduct(solver,  o_invDegree,  o_V[(k+1)],  o_V[(k+1)], options);
      t = sqrt(t);

      HH[IDX2C(k,k+1,m+1)] = t;

      /* protect */
      if (fabs(t-ZERO) > EPSILON) {
        t = 1.0 / t;
        /*----------------*
          |    V(k+1)*t
         *----------------*/
        ellipticScaledAdd(solver, t, o_V[(k+1)], 0.0,o_V[(k+1)]);

      }

      if (k !=0 )
        for (int j=1; j<=k; j++) {
          //i->k,  k->j
          k1 = j-1;

          t  = HH[IDX2C(k,k1,m+1)];

          HH[IDX2C(k,k1,m+1)] =  c[k1]*t + s[k1]*HH[IDX2C(k,j,m+1)];
          HH[IDX2C(k,j, m+1)] = -s[k1]*t + c[k1]*HH[IDX2C(k,j,m+1)];
        }

      dfloat Hii  = HH[IDX2C(k,k,m+1)];
      dfloat Hii1 = HH[IDX2C(k,k+1,m+1)];

      gam = sqrt(Hii*Hii + Hii1*Hii1);

      if (fabs(gam-ZERO) <= EPSILON) gam = EPSMAC;

      /* next Given's rotation */
      c[k] = Hii  / gam;
      s[k] = Hii1 / gam;
      rs[k+1] = -s[k] * rs[k];
      rs[k]   =  c[k] * rs[k];

      /* residue norm */
      HH[IDX2C(k,k,m+1)] = c[k]*Hii + s[k]*Hii1;
      ro = fabs(rs[k+1]);

      if (strstr(options, "VERBOSE")) {
        printf("GMRES: iteration %d residual %15.15f\n", Niter, ro);
      }
      
      /* test convergence */
      if (ro < tolrel) {
        conv = 1;
        break;
      }
    }// inner cycle

    /*-------------------------------*
      | Solve upper triangular system  |
     *-------------------------------*/

    rs[NinnerIt] /= HH[IDX2C(NinnerIt,NinnerIt,m+1)];

    for (int ii=2; ii<=NinnerIt+1; ii++) {
      int k  = NinnerIt-ii+1;
      k1 = k+1;
      t  = rs[k];
      for (int j=k1; j<=NinnerIt; j++)
        t -= HH[IDX2C(j,k,m+1)]*rs[j];

      rs[k] = t / HH[IDX2C(k,k,m+1)];
    }

    /*---------------*
      |  Get solution  
     *---------------*/
    for (int j=0; j<=NinnerIt; j++)
      ellipticScaledAdd(solver, rs[j], o_V[j], 1, o_x);
  }

  if (strstr(options, "VERBOSE")) {
    if (conv) {
      printf("GMRES(%d) CONVERGED IN %d IT \n",m, Niter);
    } else {
      printf("GMRES(%d) stagnated residual %15.15f\n", m, ro);
    }
  }

  if (strstr(options, "RIGHT")){
    ellipticPreconditioner2D(solver, lambda, o_x, o_Ax, options);
    o_x.copyFrom(o_Ax);
  }
  
  return Niter;
}


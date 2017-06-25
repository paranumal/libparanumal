/* wall 1, inflow 2, outflow 3 */

// Weakly Impose Nonlinear term BCs
#define insAdvectionBoundaryConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
  {	dfloat nu = 0.025;\
    dfloat lambda = 1.f/(2.f * nu) - occaSqrt(1.f/(4.f*nu*nu) + 4.f*OCCA_PI*OCCA_PI);\
    if(bc==1){								\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
    } else if(bc==2){							\
      *(uB) = 1.f - occaExp(lambda*x)*occaCos(2.f*OCCA_PI*y);   \
      *(vB) =  lambda/(2.f*OCCA_PI)*occaExp(lambda*x)*occaSin(2.f*OCCA_PI*y);  	\
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;							\
    }									\
  }

#define insDivergenceBoundaryConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
  {		dfloat nu = 0.025;\
    dfloat lambda = 1.f/(2.f * nu) - occaSqrt(1.f/(4.f*nu*nu) + 4.f*OCCA_PI*OCCA_PI);\
    if(bc==1){								\
      *(uB)= 0.f;							\
      *(vB)= 0.f;							\
    } else if(bc==2){							\
      *(uB) = 1.f - occaExp(lambda*x)*occaCos(2.f*OCCA_PI*y);   		\
      *(vB) = lambda/(2.f*OCCA_PI)*occaExp(lambda*x)*occaSin(2.f*OCCA_PI*y);\
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;							\
    }									\
  }

// Gradient only applies to Pressure and Pressure Incremament
// Boundary Conditions are implemented in strong form
#define insGradientBoundaryConditions2D(bc,t,x,y,nx,ny,pM,pB)	\
  {	dfloat nu = 0.025;\
    dfloat lambda = 1.f/(2.f * nu) - occaSqrt(1.f/(4.f*nu*nu) + 4.f*OCCA_PI*OCCA_PI);\
    if(bc==1){							\
      *(pB) = pM;						\
    } else if(bc==2){						\
      *(pB) = pM;						\
    } else if(bc==3){						\
      *(pB) = 0.5f*(1.f- occaExp(2.f*lambda*x));\
    }								\
  }

#define insHelmholtzBoundaryConditionsIpdg2D(bc, t, x, y, nx, ny, uB, uxB, uyB, vB, vxB, vyB) \
  {	dfloat nu = 0.025;\
    dfloat lambda = 1.f/(2.f*nu)-occaSqrt(1.f/(4.f*nu*nu) + 4.f*OCCA_PI*OCCA_PI);\
    if((bc==1)||(bc==4)){						\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
									\
      *(uxB) = 0.f;							\
      *(uyB) = 0.f;							\
      *(vxB) = 0.f;							\
      *(vyB) = 0.f;							\
    } else if(bc==2){							\
									\
      *(uB) = 1.f - occaExp(lambda*x)*occaCos(2.f*OCCA_PI*y);   	\
      *(vB) = lambda/(2.f*OCCA_PI)*occaExp(lambda*x)*occaSin(2.f*OCCA_PI*y);	\
									\
      *(uxB) = 0.f;							\
      *(uyB) = 0.f;							\
      *(vxB) = 0.f;							\
      *(vyB) = 0.f;							\
    } else if(bc==3){							\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
      *(uxB) = -lambda*occaExp(lambda*x)*occaCos(2.f*OCCA_PI*y);\
      *(uyB) = 2.f*OCCA_PI*occaExp(lambda*x)*occaSin(2.f*OCCA_PI*y); \
      *(vxB) = lambda*lambda/(2.f*OCCA_PI)*occaExp(lambda*x)*occaSin(2.f*OCCA_PI*y);   \
      *(vyB) = lambda*occaExp(lambda*x)*occaCos(2.f*OCCA_PI*y);              \
    }									\
  }


// Compute bcs for P increment: if c0 = 0 give Pr BCs, zero if time independent
#define insPoissonBoundaryConditions2D(bc,t,dt,x,y,nx,ny,pB,pxB,pyB)	\
  {	dfloat nu = 0.025;\
    dfloat lambda = 1.f/(2.f * nu) - occaSqrt(1.f/(4.f*nu*nu) + 4.f*OCCA_PI*OCCA_PI);\
    if((bc==1)||(bc==4)){						\
      *(pB) = 0.f;							\
									\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    }									\
    if(bc==2){								\
      *(pB)  = 0.f;							\
									\
      *(pxB) =  0.f; \
      *(pyB) =  0.f; \
    }									\
    if(bc==3){								\
      *(pB) =  0.5f*(1.f- occaExp(2.f*lambda*x));\
									\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    }									\
  }





/* wall 1, inflow 2, outflow 3 */

// Weakly Impose Nonlinear term BCs
#define insAdvectionBoundaryConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
  {	dfloat nu = .01f ;								\
    if(bc==1){								\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
    } else if(bc==2){							\
      *(uB) = -occaSin(2.f*OCCA_PI*y)*occaExp(-nu*4.f*OCCA_PI*OCCA_PI*t);   \
      *(vB) =  occaSin(2.f*OCCA_PI*x)*occaExp(-nu*4.f*OCCA_PI*OCCA_PI*t);  	\
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;							\
    }									\
  }

#define insDivergenceBoundaryConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
  {		dfloat nu = .01f;							\
    if(bc==1){								\
      *(uB)= 0.f;							\
      *(vB)= 0.f;							\
    } else if(bc==2){							\
      *(uB) = -occaSin(2.f*OCCA_PI*y)*occaExp(-nu*4.f*OCCA_PI*OCCA_PI*t);   		\
      *(vB) =  occaSin(2.f*OCCA_PI*x)*occaExp(-nu*4.f*OCCA_PI*OCCA_PI*t);   		\
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;							\
    }									\
  }

// Gradient only applies to Pressure and Pressure Incremament
// Boundary Conditions are implemented in strong form
#define insGradientBoundaryConditions2D(bc,t,x,y,nx,ny,pM,pB)	\
  {								\
		dfloat nu = .01f;						\
    if(bc==1){							\
      *(pB) = pM;						\
    } else if(bc==2){						\
      *(pB) = pM;						\
    } else if(bc==3){						\
      *(pB) = -occaCos(2.f*OCCA_PI*y)*occaCos(2.f*OCCA_PI*x)*occaExp(-nu*8.f*OCCA_PI*OCCA_PI*t); 	\
    }								\
  }

#define insHelmholtzBoundaryConditionsIpdg2D(bc, t, x, y, nx, ny, uB, uxB, uyB, vB, vxB, vyB) \
  {		\
     dfloat nu = .01f; 					\
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
      *(uB) = -occaSin(2.f*OCCA_PI*y)*occaExp(-nu*4.f*OCCA_PI*OCCA_PI*t); 	\
      *(vB) =  occaSin(2.f*OCCA_PI*x)*occaExp(-nu*4.f*OCCA_PI*OCCA_PI*t);		\
									\
      *(uxB) = 0.f;							\
      *(uyB) = 0.f;							\
      *(vxB) = 0.f;							\
      *(vyB) = 0.f;							\
    } else if(bc==3){							\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
      *(uxB) = 0.f;							\
      *(uyB) =-2.f*OCCA_PI*occaCos(2.f*OCCA_PI*y)*occaExp(-nu*4.f*OCCA_PI*OCCA_PI*t);  \
      *(vxB) = 2.f*OCCA_PI*occaCos(2.f*OCCA_PI*x)*occaExp(-nu*4.f*OCCA_PI*OCCA_PI*t);   \
      *(vyB) = 0.f;							\
    }									\
  }


// Compute bcs for P increment
#define insPoissonBoundaryConditions2D(bc,t,dt,x,y,nx,ny,pB,pxB,pyB)	\
  {		dfloat nu = .01f ; 		\
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
      *(pB) = -occaCos(2.f*OCCA_PI*y)*occaCos(2.f*OCCA_PI*x)*(occaExp(-nu*8.f*OCCA_PI*OCCA_PI*t)-occaExp(-nu*8.f*OCCA_PI*OCCA_PI*(t-dt))); \
									\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    }									\
  }


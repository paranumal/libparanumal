/* wall 1, inflow 2, outflow 3 */

// Weakly Impose Nonlinear term BCs
#define insAdvectionBoundaryConditions3D(bc, t, x, y, z, nx, ny, nz, uM, vM, wM, uB, vB, wB) \
  {	\
    dfloat lambda = 1.f/(2.f * p_nu) - occaSqrt(1.f/(4.f*p_nu*p_nu) + 4.f*OCCA_PI*OCCA_PI);\
    if(bc==1){								\
      *(uB) = 0.f;							\
      *(vB) = 0.f;              \
      *(wB) = 0.f;              \
    } else if(bc==2){							\
      *(uB) = 1.f - occaExp(lambda*x)*occaCos(2.f*OCCA_PI*y);   \
      *(vB) =  lambda/(2.f*OCCA_PI)*occaExp(lambda*x)*occaSin(2.f*OCCA_PI*y);   \
      *(wB) =  0.f;   \
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;             \
      *(wB) = wM;             \
    }									\
  }

#define insDivergenceBoundaryConditions3D(bc, t, x, y, z, nx, ny, nz, uM, vM, wM, uB, vB, wB) \
  {	\
    dfloat lambda = 1.f/(2.f * p_nu) - occaSqrt(1.f/(4.f*p_nu*p_nu) + 4.f*OCCA_PI*OCCA_PI);\
    if(bc==1){                \
      *(uB) = 0.f;              \
      *(vB) = 0.f;              \
      *(wB) = 0.f;              \
    } else if(bc==2){             \
      *(uB) = 1.f - occaExp(lambda*x)*occaCos(2.f*OCCA_PI*y);   \
      *(vB) =  lambda/(2.f*OCCA_PI)*occaExp(lambda*x)*occaSin(2.f*OCCA_PI*y);   \
      *(wB) =  0.f;   \
    } else if(bc==3){             \
      *(uB) = uM;             \
      *(vB) = vM;             \
      *(wB) = wM;             \
    }                 \
  }

// Gradient only applies to Pressure and Pressure Incremament
// Boundary Conditions are implemented in strong form
#define insGradientBoundaryConditions3D(bc,t,x,y,z,nx,ny,nz,pM,pB)  \
  {	\
    dfloat lambda = 1.f/(2.f * p_nu) - occaSqrt(1.f/(4.f*p_nu*p_nu) + 4.f*OCCA_PI*OCCA_PI);\
    if(bc==1){							\
      *(pB) = pM;						\
    } else if(bc==2){						\
      *(pB) = pM;						\
    } else if(bc==3){						\
      *(pB) = 0.5f*(1.f- occaExp(2.f*lambda*x));\
    }								\
  }

#define insHelmholtzBoundaryConditionsIpdg3D(bc,t,x,y,z, nx,ny,nz, uB,uxB,uyB,uzB, vB,vxB,vyB,vzB, wB,wxB,wyB,wzB) \
  {	\
    dfloat lambda = 1.f/(2.f*p_nu)-occaSqrt(1.f/(4.f*p_nu*p_nu) + 4.f*OCCA_PI*OCCA_PI);\
    if((bc==1)||(bc==4)){						\
      *(uB) = 0.f;              \
      *(vB) = 0.f;              \
      *(wB) = 0.f;              \
                                \
      *(uxB) = 0.f;             \
      *(uyB) = 0.f;             \
      *(uzB) = 0.f;             \
      *(vxB) = 0.f;             \
      *(vyB) = 0.f;             \
      *(vzB) = 0.f;\
      *(wxB) = 0.f;             \
      *(wyB) = 0.f;             \
      *(wzB) = 0.f;             \
    } else if(bc==2){							\
									\
      *(uB) = 1.f - occaExp(lambda*x)*occaCos(2.f*OCCA_PI*y);   	\
      *(vB) = lambda/(2.f*OCCA_PI)*occaExp(lambda*x)*occaSin(2.f*OCCA_PI*y);	\
      *(wB) = 0.f;              \
      \
      *(uxB) = 0.f;             \
      *(uyB) = 0.f;             \
      *(uzB) = 0.f;             \
      *(vxB) = 0.f;             \
      *(vyB) = 0.f;             \
      *(vzB) = 0.f;             \
                                \
      *(wxB) = 0.f;             \
      *(wyB) = 0.f;             \
      *(wzB) = 0.f;             \
    } else if(bc==3){							\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
      *(uxB) = -lambda*occaExp(lambda*x)*occaCos(2.f*OCCA_PI*y);\
      *(uyB) = 2.f*OCCA_PI*occaExp(lambda*x)*occaSin(2.f*OCCA_PI*y); \
      *(uzB) = 0.f; \
      *(vxB) = lambda*lambda/(2.f*OCCA_PI)*occaExp(lambda*x)*occaSin(2.f*OCCA_PI*y);   \
      *(vyB) = lambda*occaExp(lambda*x)*occaCos(2.f*OCCA_PI*y);              \
      *(vzB) = 0.f;             \
      *(wxB) = 0.f; \
      *(wyB) = 0.f; \
      *(wzB) = 0.f; \
    }									\
  }


// Compute bcs for P increment: if c0 = 0 give Pr BCs, zero if time independent
#define insPoissonBoundaryConditions3D(bc,t,dt,x,y,z,nx,ny,nz,pB,pxB,pyB,pzB) \
  {	\
    dfloat lambda = 1.f/(2.f * p_nu) - occaSqrt(1.f/(4.f*p_nu*p_nu) + 4.f*OCCA_PI*OCCA_PI);\
    if((bc==1)||(bc==4)){						\
      *(pB) = 0.f;							\
									\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
      *(pzB) = 0.f;              \
    }									\
    if(bc==2){								\
      *(pB)  = 0.f;							\
									\
      *(pxB) =  0.f; \
      *(pyB) =  0.f; \
      *(pzB) = 0.f;              \
    }									\
    if(bc==3){								\
      *(pB) =  0.5f*(1.f- occaExp(2.f*lambda*x));\
									\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
      *(pzB) = 0.f;              \
    }									\
  }





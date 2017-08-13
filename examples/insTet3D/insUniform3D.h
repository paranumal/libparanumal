/* wall 1, inflow 2, outflow 3 */

// Weakly Impose Nonlinear term BCs
#define insAdvectionBoundaryConditions3D(bc, t, x, y, z, nx, ny, nz, uM, vM, wM, uB, vB, wB) \
  {									\
    if(bc==1){								\
      *(uB) = 0.f;							\
      *(vB) = 0.f;              \
      *(wB) = 0.f;              \
    } else if(bc==2){							\
      *(uB) = 1.5f*y*(6.0f-y)/(3.0f*3.0f); \
      *(vB) = 0.f;              \
      *(wB) = 0.f;              \
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;             \
      *(wB) = wM;             \
    }									\
  }

#define insDivergenceBoundaryConditions3D(bc, t, x, y, z, nx, ny, nz, uM, vM, wM, uB, vB, wB) \
  {									\
    if(bc==1){								\
      *(uB)= 0.f;							\
      *(vB)= 0.f;             \
      *(wB)= 0.f;             \
    } else if(bc==2){							\
      *(uB) = 1.5f*y*(6.0f-y)/(3.0f*3.0f); \
      *(vB) = 0.0f;             \
      *(wB) = 0.0f;             \
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;             \
      *(wB) = wM;             \
    }									\
  }

// Gradient only applies to Pressure and Pressure Incremament
// Boundary Conditions are implemented in strong form
#define insGradientBoundaryConditions3D(bc,t,x,y,z,nx,ny,nz,pM,pB)  \
  {								\
								\
    if(bc==1){							\
      *(pB) = pM;						\
    } else if(bc==2){						\
      *(pB) = pM;						\
    } else if(bc==3){						\
      *(pB) = 0.f;						\
    }								\
  }

#define insHelmholtzBoundaryConditionsIpdg3D(bc,t,x,y,z, nx,ny,nz, uB,uxB,uyB,uzB, vB,vxB,vyB,vzB, wB,wxB,wyB,wzB) \
  {									\
    if((bc==1)||(bc==4)){						\
      *(uB) = 0.f;							\
      *(vB) = 0.f;              \
      *(wB) = 0.f;              \
									\
      *(uxB) = 0.f;             \
      *(uyB) = 0.f;             \
      *(uzB) = 0.f;             \
                                \
      *(vxB) = 0.f;             \
      *(vyB) = 0.f;             \
      *(vzB) = 0.f;             \
                                \
      *(wxB) = 0.f;             \
      *(wyB) = 0.f;             \
      *(wzB) = 0.f;             \
    } else if(bc==2){							\
									\
      *(uB) = 1.5f*y*(6.0f-y)/(3.0f*3.0f); \
      *(vB) = 0.f;              \
      *(wB) = 0.f;              \
									\
      *(uxB) = 0.f;             \
      *(uyB) = 0.f;\
      *(uzB) = 0.f;             \
                                \
      *(vxB) = 0.f;             \
      *(vyB) = 0.f;             \
      *(vzB) = 0.f;             \
                                \
      *(wxB) = 0.f;             \
      *(wyB) = 0.f;             \
      *(wzB) = 0.f;             \
    } else if(bc==3){							\
      *(uB) = 0.f;							\
      *(vB) = 0.f;              \
      *(wB) = 0.f;              \
                                \
      *(uxB) = 0.f;             \
      *(uyB) = 0.f;             \
      *(uzB) = 0.f;             \
                                \
      *(vxB) = 0.f;             \
      *(vyB) = 0.f;             \
      *(vzB) = 0.f;             \
                                \
      *(wxB) = 0.f;             \
      *(wyB) = 0.f;             \
      *(wzB) = 0.f;             \
    }									\
  }


// Compute bcs for P increment
#define insPoissonBoundaryConditions3D(bc,t,dt,x,y,z,nx,ny,nz,pB,pxB,pyB,pzB)	\
  {									\
    if((bc==1)||(bc==4)){           \
      *(pB) = 0.f;              \
                  \
      *(pxB) = 0.f;             \
      *(pyB) = 0.f;             \
      *(pzB) = 0.f;              \
    }                 \
    if(bc==2){                \
      *(pB)  = 0.f;             \
                  \
      *(pxB) = 0.f; \
      *(pyB) = 0.f; \
      *(pzB) = 0.f; \
    }                 \
    if(bc==3){                \
      *(pB) = 0.f; \
                  \
      *(pxB) = 0.f;             \
      *(pyB) = 0.f;             \
      *(pzB) = 0.f;             \
    }                 \
  }


// Initial conditions 
#define insFlowField3D(t,x,y,z,u,v,w,p)		\
  {						\
    *(u) = 1.0f;				\
    *(v) = 0.0f;        \
    *(w) = 0.0f;        \
    *(p) = 0.0f;				\
  }						

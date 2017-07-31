/* wall 1, inflow 2, outflow 3 */

// Weakly Impose Nonlinear term BCs
#define insAdvectionBoundaryConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
  {									\
    if(bc==1){								\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
    } else if(bc==2){							\
      *(uB) = 1.0f;	\
      *(vB) = 0.f;							\
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;							\
    }									\
  }

#define insDivergenceBoundaryConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
  {									\
    if(bc==1){								\
      *(uB)= 0.f;							\
      *(vB)= 0.f;							\
    } else if(bc==2){							\
      *(uB) = 1.0f; \
      *(vB) = 0.0f;							\
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;							\
    }									\
  }

// Gradient only applies to Pressure and Pressure Incremament
// Boundary Conditions are implemented in strong form
#define insGradientBoundaryConditions2D(bc,t,x,y,nx,ny,pM,pB)	\
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

#define insHelmholtzBoundaryConditionsIpdg2D(bc, t, x, y, nx, ny, uB, uxB, uyB, vB, vxB, vyB) \
  {									\
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
      *(uB) = 1.0f; \
      *(vB) = 0.f;							\
									\
      *(uxB) = 0.f;							\
      *(uyB) = 0.f;							\
      *(vxB) = 0.f;							\
      *(vyB) = 0.f;							\
    } else if(bc==3){							\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
      *(uxB) = 0.f;							\
      *(uyB) = 0.f;							\
      *(vxB) = 0.f;							\
      *(vyB) = 0.f;							\
    }									\
  }


// Compute bcs for P increment
#define insPoissonBoundaryConditions2D(bc,t,dt,x,y,nx,ny,pB,pxB,pyB)	\
  {									\
    if((bc==1)||(bc==4)){						\
      *(pB) = 0.f;							\
									\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    }									\
    if(bc==2){								\
      *(pB)  = 0.f;							\
									\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    }									\
    if(bc==3){								\
      *(pB) = 0.f;							\
									\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    }									\
  }


// Initial conditions 
#define insFlowField2D(t,x,y,u,v,p)		\
  {						\
    *(u) = 1.0f;				\
    *(v) = 0.0f;				\
    *(p) = 0.0f;				\
  }						

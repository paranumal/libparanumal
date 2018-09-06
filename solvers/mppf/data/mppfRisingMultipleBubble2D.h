// Initial conditions : f and g are forcing for exact solution 
// messy manufactured solution
#define mppfFlowField2D(t,x,y,u,v,p) \
{                           \
  *(u) = p_ubar;            \
  *(v) = p_vbar;            \
  *(p) = p_pbar;            \
}


#define mppfPhaseField2D(t,h,x,y,phi) \
{                                                         \
  dfloat L     = 1.0f;                                     \
  dfloat xc1    = 0.f;                                      \
  dfloat yc1    = 1.25f*L;                                   \
  dfloat rad1  = 0.25f*L;                                  \
  dfloat xc2   = 0.0f;                                      \
  dfloat yc2   = 0.5f*L;                                   \
  dfloat rad2  = 0.25f*L;                                  \
  dfloat phiL1 = sqrt((x-xc1)*(x-xc1) + (y-yc1)*(y-yc1)) - rad1;\
  dfloat phiL2 = sqrt((x-xc2)*(x-xc2) + (y-yc2)*(y-yc2)) - rad2;\
  *(phi)      = 1.0f + tanh(-phiL1/(sqrt(2.f)*h)) +tanh(-phiL2/(sqrt(2.f)*h));               \
}
  // dfloat xc3   = -0.5f;                                      \
  // dfloat yc3   = 0.5f*L;                                   \
  // dfloat rad3  = 0.25f*L;                                  \
  dfloat phiL3 = sqrt((x-xc3)*(x-xc3) + (y-yc3)*(y-yc3)) - rad3;\
//  *(phi)      = 1.0f + tanh(phiL1/(sqrt(2.f)*h)) +tanh(phiL2/(sqrt(2.f)*h)) +tanh(phiL3/(sqrt(2.f)*h));               \

// This is only for manufactured solutions
#define mppfPhaseFieldSource2D(t,x,y,g) \
{                                       \
  *(g)           = 0.f;                 \
} 

// Only gravitational force 
#define mppfVelocitySource2D(t,x,y,fx, fy)\
{                                         \
  *(fx)          =  0.f;                  \
  *(fy)          = -9.8f;                 \
} 

#define mppfPhaseFieldDirichletConditions2D(bc, t, x, y, nx, ny, phiM, phiB) \
{                                   \
  if(bc==1){                        \
    *(phiB) = phiM;                 \
  } else if(bc==2){                 \
    *(phiB) = phiM;                 \
  } else if(bc==3){                 \
    *(phiB) = phiM;                 \
  } else if(bc==4||bc==5){          \
    *(phiB) = phiM;                 \
  }                                 \
}

#define mppfPhaseFieldNeumannConditions2D(bc, t, x, y, nx, ny, phixM, phiyM, phixB, phiyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(phixB) = phixM;                      \
    *(phiyB) = phiyM;                      \
  } else if(bc==3){                        \
    *(phixB) = 0.f;                        \
    *(phiyB) = 0.f;                        \
  } else if(bc==4||bc==5){                 \
    *(phixB) = nx*nx*phixM;                \
    *(phiyB) = ny*ny*phiyM;                \
  }                                        \
}


// Boundary conditions
#define mppfVelocityDirichletConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
{                                   \
  if(bc==1){                        \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
  } else if(bc==2){                 \
    *(uB) = p_ubar;                 \
    *(vB) = p_vbar;                 \
  } else if(bc==3){                 \
    *(uB) = uM;                     \
    *(vB) = vM;                     \
  }else if(bc==4||bc==5){           \
    *(uB) = uM - (nx*uM+ny*vM)*nx;  \
    *(vB) = vM - (nx*uM+ny*vM)*ny;  \
  }                                 \
}

#define mppfVelocityNeumannConditions2D(bc, t, x, y, nx, ny, uxM, uyM, vxM, vyM, uxB, uyB, vxB, vyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
    *(vxB) = vxM;                          \
    *(vyB) = vyM;                          \
  } else if(bc==3){                        \
    *(uxB) = 0.f;                          \
    *(uyB) = 0.f;                          \
    *(vxB) = 0.f;                          \
    *(vyB) = 0.f;                          \
  } else if(bc==4||bc==5){                 \
    *(uxB) = nx*nx*uxM;                    \
    *(uyB) = nx*nx*uyM;                    \
    *(vxB) = ny*ny*vxM;                    \
    *(vyB) = ny*ny*vyM;                    \
  }                                        \
}



#define mppfPressureDirichletConditions2D(bc, t, x, y, nx, ny, pM, pB) \
{                                   \
  if(bc==1 || bc==2){               \
    *(pB) = pM;                     \
  } else if(bc==3){                 \
    *(pB) = p_pbar;                 \
  } else if(bc==4||bc==5){          \
    *(pB) = pM;                     \
  }                                 \
}
// 

#define mppfPressureNeumannConditions2D(bc, t, x, y, nx, ny, pxM, pyM, pxB, pyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
  } else if(bc==3){                        \
    *(pxB) = pxM;                          \
    *(pyB) = pyM;                          \
  } else if(bc==4||bc==5){                 \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
  }                                        \
}

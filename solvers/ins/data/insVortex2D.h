// Initial conditions 
#define insFlowField2D(t,x,y,u,v,p)   \
  {                                   \
    *(u) = -sin(2.f*M_PI*y)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(v) =  sin(2.f*M_PI*x)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(p) = -cos(2.f*M_PI*y)*cos(2.f*M_PI*x)*exp(-p_nu*8.f*M_PI*M_PI*t);   \
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define insVelocityDirichletConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
{                                   \
  if(bc==1){                        \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
  } else if(bc==2){                 \
    *(uB) = -sin(2.f*M_PI*y)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(vB) =  sin(2.f*M_PI*x)*exp(-p_nu*4.f*M_PI*M_PI*t);\
  } else if(bc==3){                 \
    *(uB) = uM;                     \
    *(vB) = vM;                     \
  } else if(bc==4){                 \
    *(uB) = 0.f;                    \
    *(vB) = vM;                     \
  } else if(bc==5){                 \
    *(uB) = uM;                     \
    *(vB) = 0.f;                    \
  }                                 \
}

#define insVelocityNeumannConditions2D(bc, t, x, y, nx, ny, uxM, uyM, vxM, vyM, uxB, uyB, vxB, vyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
    *(vxB) = vxM;                          \
    *(vyB) = vyM;                          \
  } else if(bc==3){                        \
    *(uxB) = 0.f;                          \
    *(uyB) = -2.f*M_PI*cos(2.f*M_PI*y)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(vxB) =  2.f*M_PI*cos(2.f*M_PI*x)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(vyB) = 0.f;                          \
  } else if(bc==4){                        \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
    *(vxB) = 0.f;                          \
    *(vyB) = 0.f;                          \
  } else if(bc==5){                        \
    *(uxB) = 0.f;                          \
    *(uyB) = 0.f;                          \
    *(vxB) = vxM;                          \
    *(vyB) = vyM;                          \
  }                                        \
}


#define insPressureDirichletConditions2D(bc, t, x, y, nx, ny, pM, pB) \
{                                   \
  if(bc==1 || bc==2){               \
    *(pB) = pM;                     \
  } else if(bc==3){                 \
    *(pB) = -cos(2.f*M_PI*y)*cos(2.f*M_PI*x)*exp(-p_nu*8.f*M_PI*M_PI*t);\
  } else if(bc==4){                 \
    *(pB) = pM;                     \
  } else if(bc==5){                 \
    *(pB) = pM;                     \
  }                                 \
}

#define insPressureNeumannConditions2D(bc, t, x, y, nx, ny, pxM, pyM, pxB, pyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(pxB) = 2.f*M_PI*cos(2.f*M_PI*y)*sin(2.f*M_PI*x)*exp(-p_nu*8.f*M_PI*M_PI*t);\
    *(pyB) = 2.f*M_PI*sin(2.f*M_PI*y)*cos(2.f*M_PI*x)*exp(-p_nu*8.f*M_PI*M_PI*t);\
  } else if(bc==3){                        \
    *(pxB) = pxM;                          \
    *(pyB) = pyM;                          \
  } else if(bc==4){                        \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
  } else if(bc==5){                        \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
  }                                        \
}

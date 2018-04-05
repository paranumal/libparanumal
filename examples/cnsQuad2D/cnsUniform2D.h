
// Initial conditions 
#define cnsFlowField2D(t,x,y,u,v,p)   \
  {                                   \
    *(r) = p_rbar;                    \
    *(u) = p_ubar;                    \
    *(v) = p_vbar;                    \
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define cnsDirichletConditions2D(bc, t, x, y, nx, ny, rM, uM, vM, rB, uB, vB) \
{                                   \
  if(bc==1){                        \
    *(rB) = rM;                     \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
  } else if(bc==2){                 \
    *(rB) = p_rbar;                 \
    *(uB) = p_ubar;                 \
    *(vB) = p_vbar;                 \
  } else if(bc==3){                 \
    *(rB) = p_rbar;                 \
    *(uB) = uM;                     \
    *(vB) = vM;                     \
  } else if(bc==4||bc==5){          \
    *(rB) = rM;                     \
    *(uB) = uM - (nx*uM+ny*vM)*nx;  \
    *(vB) = vM - (nx*uM+ny*vM)*ny;  \
  }                                 \
}

#define cnsNeumannConditions2D(bc, t, x, y, nx, ny, uxM, uyM, vxM, vyM, uxB, uyB, vxB, vyB) \
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

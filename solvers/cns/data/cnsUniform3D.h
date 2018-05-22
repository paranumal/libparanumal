
// Initial conditions 
#define cnsFlowField3D(t,x,y,z,u,v,w,p)		\
  {                                   \
    *(r) = p_rbar;                    \
    *(u) = p_ubar;                    \
    *(v) = p_vbar;                    \
    *(w) = p_wbar;                    \
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define cnsDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, rM, uM, vM, wM, rB, uB, vB, wB) \
{                                   \
  if(bc==1){                        \
    *(rB) = rM;                     \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
    *(wB) = 0.f;                    \
  } else if(bc==2){                 \
    *(rB) = p_rbar;                 \
    *(uB) = p_ubar;                 \
    *(vB) = p_vbar;                 \
    *(wB) = p_wbar;                 \
  } else if(bc==3){                 \
    *(rB) = rM;                     \
    *(uB) = p_ubar;                 \
    *(vB) = p_vbar;                 \
    *(wB) = p_wbar;                 \
  } else if(bc==4||bc==5||bc==6){   \
    *(rB) = rM;                     \
    *(uB) = uM - (nx*uM+ny*vM+nz*wM)*nx;  \
    *(vB) = vM - (nx*uM+ny*vM+nz*wM)*ny;  \
    *(wB) = wM - (nx*uM+ny*vM+nz*wM)*nz;  \
  }                                 \
}

//needs updating for 3D
#define cnsNeumannConditions3D(bc, t, x, y, z, nx, ny, nz, uxM, uyM, uzM, vxM, vyM, vzM, uxB, uyB, uzB, vxB, vyB, vzB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
    *(uzB) = uzM;                          \
    *(vxB) = vxM;                          \
    *(vyB) = vyM;                          \
    *(vzB) = vzM;                          \
  } else if(bc==3){                        \
    *(uxB) = 0.f;                          \
    *(uyB) = 0.f;                          \
    *(uzB) = 0.f;                          \
    *(vxB) = 0.f;                          \
    *(vyB) = 0.f;                          \
    *(vzB) = 0.f;                          \
  } else if(bc==4){                        \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
    *(uzB) = uzM;                          \
    *(vxB) = 0.f;                          \
    *(vyB) = 0.f;                          \
    *(vzB) = 0.f;                          \
  } else if(bc==5){                        \
    *(uxB) = 0.f;                          \
    *(uyB) = 0.f;                          \
    *(uzB) = 0.f;                          \
    *(vxB) = vxM;                          \
    *(vyB) = vyM;                          \
    *(vzB) = vzM;                          \
  }                                        \
}

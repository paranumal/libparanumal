
// Initial conditions 
#define insFlowField3D(t,x,y,z, u,v,w,p) \
  {                                   \
    dfloat a = OCCA_PI/4.f; \
    dfloat d = OCCA_PI/2.f; \
    *(u) = -a*(occaExp(a*x)*occaSin(a*y+d*z)+occaExp(a*z)*occaCos(a*x+d*y))*occaExp(-d*d*t);\
    *(v) = -a*(occaExp(a*y)*occaSin(a*z+d*x)+occaExp(a*x)*occaCos(a*y+d*z))*occaExp(-d*d*t);\
    *(w) = -a*(occaExp(a*z)*occaSin(a*x+d*y)+occaExp(a*y)*occaCos(a*z+d*x))*occaExp(-d*d*t);\
    *(p) = -a*a*occaExp(-2.f*d*d*t)*(occaExp(2.f*a*x)+occaExp(2.f*a*y)+occaExp(2.f*a*z))*(occaSin(a*x+d*y)*occaCos(a*z+d*x)*occaExp(a*(y+z))+occaSin(a*y+d*z)*occaCos(a*x+d*y)*occaExp(a*(x+z))+occaSin(a*z+d*x)*occaCos(a*y+d*z)*occaExp(a*(x+y))); \
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
#define insVelocityDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, uM, vM, wM, uB, vB, wB) \
{                                   \
  dfloat a = OCCA_PI/4.f; \
  dfloat d = OCCA_PI/2.f; \
  if(bc==1){                        \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
    *(wB) = 0.f;                    \
  } else if(bc==2){                 \
    *(uB) = -a*(occaExp(a*x)*occaSin(a*y+d*z)+occaExp(a*z)*occaCos(a*x+d*y))*occaExp(-d*d*t);\
    *(vB) = -a*(occaExp(a*y)*occaSin(a*z+d*x)+occaExp(a*x)*occaCos(a*y+d*z))*occaExp(-d*d*t);\
    *(wB) = -a*(occaExp(a*z)*occaSin(a*x+d*y)+occaExp(a*y)*occaCos(a*z+d*x))*occaExp(-d*d*t);\
  } else if(bc==3){                 \
    *(uB) = uM;                     \
    *(vB) = vM;                     \
    *(wB) = wM;                     \
  } else if(bc==4||bc==5||bc==6){   \
    *(uB) = uM - (nx*uM+ny*vM+nz*wM)*nx;\
    *(vB) = vM - (nx*uM+ny*vM+nz*wM)*ny;\
    *(wB) = wM - (nx*uM+ny*vM+nz*wM)*nz;\
  }                                 \
}

#define insVelocityNeumannConditions3D(bc, t, x, y, z, nx, ny, nz, uxM, uyM, uzM, vxM, vyM, vzM, wxM, wyM, wzM, uxB, uyB, uzB, vxB, vyB, vzB, wxB, wyB, wzB) \
{                                          \
  dfloat a = OCCA_PI/4.f; \
  dfloat d = OCCA_PI/2.f; \
  if(bc==1 || bc==2){                      \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
    *(uzB) = uzM;                          \
    *(vxB) = vxM;                          \
    *(vyB) = vyM;                          \
    *(vzB) = vzM;                          \
    *(wxB) = wxM;                          \
    *(wyB) = wyM;                          \
    *(wzB) = wzM;                          \
  } else if(bc==3){                        \
    *(uxB) = -a*(a*occaExp(a*x)*occaSin(a*y+d*z)-a*occaExp(a*z)*occaSin(a*x+d*y))*occaExp(-d*d*t); \
    *(uyB) = -a*(a*occaExp(a*x)*occaCos(a*y+d*z)-d*occaExp(a*z)*occaSin(a*x+d*y))*occaExp(-d*d*t); \
    *(uzB) = -a*(d*occaExp(a*x)*occaCos(a*y+d*z)+a*occaExp(a*z)*occaCos(a*x+d*y))*occaExp(-d*d*t); \
    *(vxB) = -a*(d*occaExp(a*y)*occaCos(a*z+d*x)+a*occaExp(a*x)*occaCos(a*y+d*z))*occaExp(-d*d*t); \
    *(vyB) = -a*(a*occaExp(a*y)*occaSin(a*z+d*x)-a*occaExp(a*x)*occaSin(a*y+d*z))*occaExp(-d*d*t); \
    *(vzB) = -a*(a*occaExp(a*y)*occaCos(a*z+d*x)-d*occaExp(a*x)*occaSin(a*y+d*z))*occaExp(-d*d*t); \
    *(wxB) =  a*(a*occaExp(a*z)*occaCos(a*x+d*y)-d*occaExp(a*y)*occaSin(a*z+d*x))*occaExp(-d*d*t); \
    *(wyB) =  a*(d*occaExp(a*z)*occaCos(a*x+d*y)+a*occaExp(a*y)*occaCos(a*z+d*x))*occaExp(-d*d*t); \
    *(wzB) =  a*(a*occaExp(a*z)*occaSin(a*x+d*y)-a*occaExp(a*y)*occaSin(a*z+d*x))*occaExp(-d*d*t); \
  } else if(bc==4||bc==5||bc==6){          \
    *(uxB) = nx*nx*uxM;                    \
    *(uyB) = nx*nx*uyM;                    \
    *(uzB) = nx*nx*uzM;                    \
    *(vxB) = ny*ny*vxM;                    \
    *(vyB) = ny*ny*vyM;                    \
    *(vzB) = ny*ny*vzM;                    \
    *(wxB) = nz*nz*wxM;                    \
    *(wyB) = nz*nz*wyM;                    \
    *(wzB) = nz*nz*wzM;                    \
  }                                        \
}


#define insPressureDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, pM, pB) \
{                                   \
  dfloat a = OCCA_PI/4.f; \
  dfloat d = OCCA_PI/2.f; \
  if(bc==1 || bc==2){               \
    *(pB) = pM;                     \
  } else if(bc==3){                 \
    *(pB) = -a*a*occaExp(-2.f*d*d*t)*( occaExp(2.f*a*x)+occaExp(2.f*a*y)+occaExp(2.f*a*z))*(occaSin(a*x+d*y)*occaCos(a*z+d*x)*occaExp(a*(y+z))+occaSin(a*y+d*z)*occaCos(a*x+d*y)*occaExp(a*(x+z))+occaSin(a*z+d*x)*occaCos(a*y+d*z)*occaExp(a*(x+y))); \
  } else if(bc==4||bc==5||bc==6){   \
    *(pB) = pM;                     \
  }                                 \
}

#define insPressureNeumannConditions3D(bc, t, x, y, z, nx, ny, nz, pxM, pyM, pzM, pxB, pyB, pzB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
    *(pzB) = 0.f;                          \
  } else if(bc==3){                        \
    *(pxB) = pxM;                          \
    *(pyB) = pyM;                          \
    *(pzB) = pzM;                          \
  } else if(bc==4||bc==5||bc==6){          \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
    *(pzB) = 0.f;                          \
  }                                        \
}

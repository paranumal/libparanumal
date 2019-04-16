/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/


// Initial conditions 
#define insFlowField3D(t,x,y,z, u,v,w,p) \
  {                                   \
    dfloat a = M_PI/4.f; \
    dfloat d = M_PI/2.f; \
    *(u) = -a*exp(-p_nu*d*d*t)*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y));\
    *(v) = -a*exp(-p_nu*d*d*t)*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z));\
    *(w) = -a*exp(-p_nu*d*d*t)*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x));\
    *(p) = -0.5f*a*a*exp(-2.f*p_nu*d*d*t)*(2.f*exp(a*(z+y))*cos(a*z+d*x)*sin(a*x+d*y)+2.f*exp(a*(z+x))*cos(a*x+d*y)*sin(a*y+d*z)+2.f*exp(a*(y+x))*cos(a*y+d*z)*sin(a*z+d*x)+exp(2.f*a*z)+exp(2.f*a*y)+exp(2.f*a*x));\
  }   
  
//*(p) = -a*a*exp(-2.f*d*d*t)*(exp(2.f*a*x)+exp(2.f*a*y)+exp(2.f*a*z))*(sin(a*x+d*y)*cos(a*z+d*x)*exp(a*(y+z))+sin(a*y+d*z)*cos(a*x+d*y)*exp(a*(x+z))+sin(a*z+d*x)*cos(a*y+d*z)*exp(a*(x+y))); \

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
#define insVelocityDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, uM, vM, wM, uB, vB, wB) \
{                                   \
  dfloat a = M_PI/4.f; \
  dfloat d = M_PI/2.f; \
  if(bc==1){                        \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
    *(wB) = 0.f;                    \
  } else if(bc==2){                 \
    *(uB) = -a*exp(-p_nu*d*d*t)*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y));\
    *(vB) = -a*exp(-p_nu*d*d*t)*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z));\
    *(wB) = -a*exp(-p_nu*d*d*t)*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x));\
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
  dfloat a = M_PI/4.f; \
  dfloat d = M_PI/2.f; \
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
    *(uxB) = -a*(a*exp(a*x)*sin(a*y+d*z)-a*exp(a*z)*sin(a*x+d*y))*exp(-p_nu*d*d*t); \
    *(uyB) = -a*(a*exp(a*x)*cos(a*y+d*z)-d*exp(a*z)*sin(a*x+d*y))*exp(-p_nu*d*d*t); \
    *(uzB) = -a*(d*exp(a*x)*cos(a*y+d*z)+a*exp(a*z)*cos(a*x+d*y))*exp(-p_nu*d*d*t); \
    *(vxB) = -a*(d*exp(a*y)*cos(a*z+d*x)+a*exp(a*x)*cos(a*y+d*z))*exp(-p_nu*d*d*t); \
    *(vyB) = -a*(a*exp(a*y)*sin(a*z+d*x)-a*exp(a*x)*sin(a*y+d*z))*exp(-p_nu*d*d*t); \
    *(vzB) = -a*(a*exp(a*y)*cos(a*z+d*x)-d*exp(a*x)*sin(a*y+d*z))*exp(-p_nu*d*d*t); \
    *(wxB) =  a*(a*exp(a*z)*cos(a*x+d*y)-d*exp(a*y)*sin(a*z+d*x))*exp(-p_nu*d*d*t); \
    *(wyB) =  a*(d*exp(a*z)*cos(a*x+d*y)+a*exp(a*y)*cos(a*z+d*x))*exp(-p_nu*d*d*t); \
    *(wzB) =  a*(a*exp(a*z)*sin(a*x+d*y)-a*exp(a*y)*sin(a*z+d*x))*exp(-p_nu*d*d*t); \
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
  dfloat a = M_PI/4.f; \
  dfloat d = M_PI/2.f; \
  if(bc==1 || bc==2){               \
    *(pB) = pM;                     \
  } else if(bc==3){                 \
    *(pB) = -0.5f*a*a*exp(-2.f*p_nu*d*d*t)*(2.f*exp(a*(z+y))*cos(a*z+d*x)*sin(a*x+d*y)+2.f*exp(a*(z+x))*cos(a*x+d*y)*sin(a*y+d*z)+2.f*exp(a*(y+x))*cos(a*y+d*z)*sin(a*z+d*x)+exp(2.f*a*z)+exp(2.f*a*y)+exp(2.f*a*x));\
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

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
    *(u) = -sin(2.f*M_PI*y)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(v) =  sin(2.f*M_PI*x)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(w) =  0.f;\
    *(p) = -cos(2.f*M_PI*y)*cos(2.f*M_PI*x)*exp(-p_nu*8.f*M_PI*M_PI*t);   \
  }   
  
   // *(p) = -0.5f*a*a*exp(-2.f*p_nu*d*d*t)*(2.f*exp(a*(z+y))*cos(a*z+d*x)*sin(a*x+d*y)+2.f*exp(a*(z+x))*cos(a*x+d*y)*sin(a*y+d*z)+2.f*exp(a*(y+x))*cos(a*y+d*z)*sin(a*z+d*x)+exp(2.f*a*z)+exp(2.f*a*y)+exp(2.f*a*x));\
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
    *(uB) = -sin(2.f*M_PI*y)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(vB) =  sin(2.f*M_PI*x)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(wB) =  0.f; \
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
    *(uxB) = 0.f;                           \
    *(uyB) = -2.f*M_PI*cos(2.f*M_PI*y)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(uzB) = 0.f; \
    *(vxB) = 2.f*M_PI*cos(2.f*M_PI*x)*exp(-p_nu*4.f*M_PI*M_PI*t);\
    *(vyB) = 0.f; \
    *(vzB) = 0.f; \
    *(wxB) = 0.f; \
    *(wyB) = 0.f; \
    *(wzB) = 0.f; \
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
  if(bc==1 || bc==2){               \
    *(pB) = pM;                     \
  } else if(bc==3){                 \
    *(pB) = -cos(2.f*M_PI*y)*cos(2.f*M_PI*x)*exp(-p_nu*8.f*M_PI*M_PI*t);\
  } else if(bc==4||bc==5||bc==6){   \
    *(pB) = pM;                     \
  }                                 \
}

#define insPressureNeumannConditions3D(bc, t, x, y, z, nx, ny, nz, pxM, pyM, pzM, pxB, pyB, pzB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(pxB) =2.f*M_PI*cos(2.f*M_PI*y)*sin(2.f*M_PI*x)*exp(-p_nu*8.f*M_PI*M_PI*t);\
    *(pyB) =2.f*M_PI*sin(2.f*M_PI*y)*cos(2.f*M_PI*x)*exp(-p_nu*8.f*M_PI*M_PI*t);\
    *(pzB) =0.f; \
  } else if(bc==3){                        \
    *(pxB) = pxM; \
    *(pyB) = pyM; \
    *(pzB) = pzM; \
  } else if(bc==4||bc==5||bc==6){          \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
    *(pzB) = 0.f;                          \
  }                                        \
}





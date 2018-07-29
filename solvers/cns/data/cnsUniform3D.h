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
#define cnsFlowField3D(t,x,y,z,u,v,w,p)		\
  {                                   \
    *(r) = p_rbar;                    \
    *(u) = p_ubar;                    \
    *(v) = p_vbar;                    \
    *(w) = p_wbar;                    \
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define cnsDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, intfx, intfy, intfz, rM, uM, vM, wM, rB, uB, vB, wB) \
{                                   \
  if(bc==1){                        \
    *(rB) = rM;                     \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
    *(wB) = 0.f;                    \
  } else if(bc==2){                 \
    *(rB) = p_rbar;                 \
    *(uB) = intfx;                 \
    *(vB) = intfy;                 \
    *(wB) = intfz;                 \
  } else if(bc==3){                 \
    *(rB) = rM;                     \
    *(uB) = intfx;                 \
    *(vB) = intfy;                 \
    *(wB) = intfz;                 \
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

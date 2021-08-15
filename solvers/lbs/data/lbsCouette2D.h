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

//mean flow
#define CS 1.0/sqrt(3.0)
#define MA 0.1
#define RBAR 1.0
#define VBAR 0.0

// const dfloat L = 1.f; \
//   const dfloat U = MA*1.0/sqrt(3.0); \
//   dfloat un = U*y/L; \
//   for(int n=1; n<1000; ++n){ \
//     const dfloat lambda = n*M_PI/L; \
//     const dfloat sone = (n%2)==0 ? 1.f:-1.f; \
//     un += 2.f*U*sone/(lambda*L)*exp(-nu*lambda*lambda*t)*sin(lambda*y);\
//   }\

// Initial conditions
#define lbsInitialConditions2D(nu, t, x, y, r, u, v) \
{                                         \
  *(r) = RBAR;                            \
  *(u) = 0.f*MA*1.f/sqrt(3.f)*y;              \
  *(v) = VBAR;                            \
}

// Body force
#define lbsBodyForce2D(nu, t, x, y, r, u, v, fx, fy) \
{                                                   \
  *(fx) = 0.0;                                      \
  *(fy) = 0.0;                                      \
}

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define lbsBoundaryConditions2D(bc, nu, \
                                t, x, y, nx, ny, \
                                rM, uM, vM, \
                                rB, uB, vB) \
{                                      \
  if(bc==1){                           \
    *(rB) = RBAR;                      \
    *(uB) = 0.0;                     \
    *(vB) = 0.0;                     \
  } else if(bc==2){                    \
    *(rB) = RBAR;                      \
    *(uB) = MA*1.f/sqrt(3.f)*y;        \
    *(vB) = VBAR;                      \
  } else if(bc==3){                    \
    *(rB) = RBAR;                      \
    *(uB) = uM;                        \
    *(vB) = vM;                        \
  } else if(bc==4||bc==5){             \
    *(rB) = rM;                        \
    *(uB) = uM - (nx*uM+ny*vM)*nx;     \
    *(vB) = vM - (nx*uM+ny*vM)*ny;     \
  }                                    \
}

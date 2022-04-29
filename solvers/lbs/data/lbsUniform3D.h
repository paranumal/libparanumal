/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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
#define RBAR 1.0
#define UBAR 1.0
#define VBAR 0.0
#define WBAR 0.0


// Initial conditions
#define lbsInitialConditions3D(nu, t, x, y, z, \
                               r, u, v, w)        \
{                                         \
  *(r) = RBAR;                            \
  *(u) = UBAR;                            \
  *(v) = VBAR;                            \
  *(w) = WBAR;                            \
}

// Body force
#define lbsBodyForce3D(nu, t, x, y, z, r, u, v, w, fx, fy, fz) \
{                                                   \
  *(fx) = 0.0;                                      \
  *(fy) = 0.0;                                      \
  *(fz) = 0.0;                                      \
}

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define lbsBoundaryConditions3D(bc, c, nu, \
                                t, x, y, z, nx, ny, nz, \
                                rM, uM, vM, wM, \
                                rB, uB, vB, wB) \
{                                      \
  if(bc==1){                           \
    *(rB) = rM;                        \
    *(uB) = 0.0;                       \
    *(vB) = 0.0;                       \
    *(wB) = 0.0;                       \
  } else if(bc==2){                    \
    *(rB) = RBAR;                      \
    *(uB) = UBAR;                      \
    *(vB) = VBAR;                      \
    *(wB) = WBAR;                      \
  } else if(bc==3){                    \
    *(rB) = RBAR;                      \
    *(uB) = uM;                        \
    *(vB) = vM;                        \
    *(wB) = wM;                        \
  } else if(bc==4||bc==5||bc==6){      \
    *(rB) = rM;                        \
    *(uB) = uM - (nx*uM+ny*vM+nz*wM)*nx; \
    *(vB) = vM - (nx*uM+ny*vM+nz*wM)*ny; \
    *(wB) = wM - (nx*uM+ny*vM+nz*wM)*nz; \
  }                                    \
}

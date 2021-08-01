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
#define RBAR 1.0
#define UBAR 1.0
#define VBAR 0.0
#define WBAR 0.0


// Initial conditions
#define lbsInitialConditions3D(c, nu, t, x, y, z, \
                               r, u, v, w,        \
                               s11, s12, s13,     \
                               s22, s23, s33)     \
{                                         \
  *(r) = RBAR;                            \
  *(u) = UBAR;                            \
  *(v) = VBAR;                            \
  *(w) = WBAR;                            \
  *(s11) = 0.0;                           \
  *(s12) = 0.0;                           \
  *(s13) = 0.0;                           \
  *(s22) = 0.0;                           \
  *(s23) = 0.0;                           \
  *(s33) = 0.0;                           \
}

// Body force
#define lbsBodyForce3D(c, nu, t, x, y, z, r, u, v, w, fx, fy, fz) \
{                                                   \
  *(fx) = 0.0;                                      \
  *(fy) = 0.0;                                      \
  *(fz) = 0.0;                                      \
}

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define lbsBoundaryConditions3D(bc, c, nu, \
                                t, x, y, z, nx, ny, nz, \
                                rM, uM, vM, wM, s11M, s12M, s13M, s22M, s23M, s33M, \
                                rB, uB, vB, wB, s11B, s12B, s13B, s22B, s23B, s33B) \
{                                      \
  if(bc==1){                           \
    *(rB) = rM;                        \
    *(uB) = 0.0;                       \
    *(vB) = 0.0;                       \
    *(wB) = 0.0;                       \
    *(s11B) = 0.0;                     \
    *(s12B) = 0.0;                     \
    *(s13B) = 0.0;                     \
    *(s22B) = 0.0;                     \
    *(s23B) = 0.0;                     \
    *(s33B) = 0.0;                     \
  } else if(bc==2){                    \
    *(rB) = RBAR;                      \
    *(uB) = UBAR;                      \
    *(vB) = VBAR;                      \
    *(wB) = WBAR;                      \
    *(s11B) = 0.0;                     \
    *(s12B) = 0.0;                     \
    *(s13B) = 0.0;                     \
    *(s22B) = 0.0;                     \
    *(s23B) = 0.0;                     \
    *(s33B) = 0.0;                     \
  } else if(bc==3){                    \
    *(rB) = RBAR;                      \
    *(uB) = uM;                        \
    *(vB) = vM;                        \
    *(wB) = wM;                        \
    *(s11B) = s11M;                    \
    *(s12B) = s12M;                    \
    *(s13B) = s13M;                    \
    *(s22B) = s22M;                    \
    *(s23B) = s23M;                    \
    *(s33B) = s33M;                    \
  } else if(bc==4||bc==5||bc==6){      \
    *(rB) = rM;                        \
    *(uB) = uM - (nx*uM+ny*vM+nz*wM)*nx; \
    *(vB) = vM - (nx*uM+ny*vM+nz*wM)*ny; \
    *(wB) = wM - (nx*uM+ny*vM+nz*wM)*nz; \
    *(s11B) = s11M;                    \
    *(s12B) = s12M;                    \
    *(s13B) = s13M;                    \
    *(s22B) = s22M;                    \
    *(s23B) = s23M;                    \
    *(s33B) = s33M;                    \
  }                                    \
}

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
#define lbsInitialConditions2D(nu, t, x, y, r, u, v) \
{                                         \
  *(r) = 1 + exp(-3*(x*x+y*y));           \
  *(u) = exp(-3*(x*x+y*y));               \
  *(v) = exp(-3*(x*x+y*y));               \
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
    *(rB) = rM;                        \
    *(uB) = 0.0;                       \
    *(vB) = 0.0;                       \
  } else if(bc==2){                    \
    *(rB) = 1.0;                       \
    *(uB) = 0.0;                       \
    *(vB) = 0.0;                       \
  } else if(bc==3){                    \
    *(rB) = 1.0;                       \
    *(uB) = uM;                        \
    *(vB) = vM;                        \
  } else if(bc==4||bc==5){             \
    *(rB) = rM;                        \
    *(uB) = uM - (nx*uM+ny*vM)*nx;     \
    *(vB) = vM - (nx*uM+ny*vM)*ny;     \
  }                                    \
}

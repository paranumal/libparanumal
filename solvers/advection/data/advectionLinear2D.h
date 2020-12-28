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

#define ADVECTION_SPEED_X 1.0
#define ADVECTION_SPEED_Y 0.5

// Flux function
#define advectionFlux2D(t, x, y, q, cx, cy) \
{                                       \
  *(cx) = ADVECTION_SPEED_X*q;          \
  *(cy) = ADVECTION_SPEED_Y*q;          \
}

// max wavespeed (should be max eigen of Jacobian of flux function)
#define advectionMaxWaveSpeed2D(t, x, y, q, u, v) \
{                                                 \
  *(u) = ADVECTION_SPEED_X;                       \
  *(v) = ADVECTION_SPEED_Y;                       \
}

// Boundary conditions
/* wall 1, outflow 2 */
#define advectionDirichletConditions2D(bc, t, x, y, nx, ny, qM, qB) \
{                                       \
  if(bc==1){                            \
    *(qB) = 0.0;                        \
  } else if(bc==2){                     \
    *(qB) = qM;                         \
  }                                     \
}

// Initial conditions
#define advectionInitialConditions2D(t, x, y, q) \
{                                       \
  *(q) = exp(-3*(x*x+y*y));             \
}

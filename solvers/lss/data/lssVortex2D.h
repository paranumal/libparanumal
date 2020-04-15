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

// Level-Set function
#define lssInitialConditions2D(t, x, y, q) \
{                                       \
  const dfloat xc = 0.50; 			    \
  const dfloat yc = 0.75; 				\
  dfloat rc = 0.15; 					\
  *(q) = sqrt((x-xc+0.5)*(x-xc+0.5) + (y-yc+0.5)*(y-yc+0.5)) - rc; ; \
}

// LS Advective field
#define lssAdvectionField2D(t, x, y, q, u, v) \
{                                       \
  const dfloat PERIOD = 8.f; 			\
  const dfloat xs = x + 0.5; 			\
  const dfloat ys = y + 0.5; 			\
 *(u) = sin(M_PI*xs)*sin(M_PI*xs)*sin(2.f*M_PI*ys)*cos(M_PI*t/PERIOD); \
 *(v) =-sin(M_PI*ys)*sin(M_PI*ys)*sin(2.f*M_PI*xs)*cos(M_PI*t/PERIOD); \
}



// Boundary conditions
/* wall 1, outflow 2 */
#define lssDirichletConditions2D(bc, t, x, y, nx, ny, qM, qB) \
{                                       \
  if(bc==1){                            \
    *(qB) = qM;                        \
  } else if(bc==2){                     \
    *(qB) = qM;                         \
  }                                     \
}
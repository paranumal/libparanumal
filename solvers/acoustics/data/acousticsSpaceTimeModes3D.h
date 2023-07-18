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



// Boundary conditions
/* wall 1, inflow 2 */
#define acousticsDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, rM, uM, vM, wM, rB, uB, vB, wB) \
{                                                \
  if(bc==1){                                     \
    *(rB) = -rM;                                 \
    *(uB) = uM;                                  \
    *(vB) = vM;                                  \
    *(wB) = wM;                                  \
  } else if(bc==2){                              \
    *(rB) = rM;                                  \
    *(uB) = uM - 2.0*(nx*uM+ny*vM+nz*wM)*nx;     \
    *(vB) = vM - 2.0*(nx*uM+ny*vM+nz*wM)*ny;     \
    *(wB) = wM - 2.0*(nx*uM+ny*vM+nz*wM)*nz;     \
  }                                              \
}

#if 0
// Initial conditions
#define acousticsInitialConditions3D(t, x, y, z, r, u, v, w) \
{                                       \
  dfloat xc = 14, yc = 22, zc = 2;                               \
  *(r) = exp(-5.*((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)));    \
  *(u) = 0.0;                                                     \
  *(v) = 0.0;                                                     \
  *(w) = 0.0;                                                     \
}
#else

// Initial conditions
#define acousticsInitialConditions3D(t, x, y, z, r, u, v, w)            \
  {                                                                     \
    *(r) = 0; *(u) = 0; *(v) = 0; *(w)  = 0;                            \
    for(int mode=1;mode<4;++mode){                                      \
      dfloat mPI = mode*M_PI;                                           \
      dfloat sc = exp(-(dfloat)mode);                                   \
      *(r) += sc*sin(mPI*x)*sin(mPI*y)*sin(mPI*z)*sin(mPI*sqrt(3.)*t);  \
      *(u) += sc*cos(mPI*x)*sin(mPI*y)*sin(mPI*z)*cos(mPI*sqrt(3.)*t)/sqrt(3.); \
      *(v) += sc*sin(mPI*x)*cos(mPI*y)*sin(mPI*z)*cos(mPI*sqrt(3.)*t)/sqrt(3.); \
      *(w) += sc*sin(mPI*x)*sin(mPI*y)*cos(mPI*z)*cos(mPI*sqrt(3.)*t)/sqrt(3.); \
    }                                                                   \
  }

#endif

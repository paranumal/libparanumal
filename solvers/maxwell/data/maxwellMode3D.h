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
/* PEC=1 */
#define maxwellDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, HxM, HyM, HzM, ExM, EyM, EzM, HxB, HyB, HzB, ExB, EyB, EzB) \
  {									\
    if(bc==1){								\
      *(HxB) = HxM;							\
      *(HyB) = HyM;							\
      *(HzB) = HzM;							\
      *(ExB) = -ExM;							\
      *(EyB) = -EyM;							\
      *(EzB) = -EzM;							\
    }									\
  }

// Initial conditions
#define maxwellInitialConditions3D(t, x, y, z, Hx, Hy, Hz, Ex, Ey, Ez)	\
  {									\
    dfloat n = 1., m = 2.;						\
    const dfloat omega = M_PI*sqrt(n*n+m*m);				\
    *(Hx) = -M_PI*n*sin(m*M_PI*x)*cos(n*M_PI*y)*sin(omega*t)/omega;	\
    *(Hy) =  M_PI*m*cos(m*M_PI*x)*sin(n*M_PI*y)*sin(omega*t)/omega;	\
    *(Hz) = 0;								\
    *(Ex) = 0;								\
    *(Ey) = 0;								\
    *(Ez) = sin(m*M_PI*x)*sin(n*M_PI*y)*cos(omega*t);			\
  }

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

#define PI 3.14159265358979323846

#define hS 10

/* forcing function   */
#define ellipticForcing2D(x, y, lambda, f)				\
  {									\
    dfloat  spix = sin(PI*x);						\
    dfloat  spix2 = spix*spix;						\
    dfloat c2pix = cos(2*PI*x);						\
    dfloat u = spix*exp(-hS*spix2);					\
    dfloat  spiy = sin(PI*y);						\
    dfloat  spiy2 = spiy*spiy;						\
    dfloat c2piy = cos(2*PI*y);						\
    dfloat v = spiy*exp(-hS*spiy2);					\
    dfloat lapu =  -PI*PI*exp(-hS*spix2)*spix*(2.*hS + hS*hS*c2pix*c2pix + 4.*hS*c2pix - hS*hS + 1.); \
    dfloat lapv =  -PI*PI*exp(-hS*spiy2)*spiy*(2.*hS + hS*hS*c2piy*c2piy + 4.*hS*c2piy - hS*hS + 1.); \
    f = -(lapu*v + u*lapv) + lambda*u*v; \
  }

/* Dirichlet boundary condition   */
#define ellipticDirichletCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)	\
  {									\
    dfloat  spix = sin(PI*x);						\
    dfloat  spix2 = spix*spix;						\
    dfloat c2pix = cos(2*PI*x);						\
    dfloat u = spix*exp(-hS*spix2);					\
    dfloat  spiy = sin(PI*y);						\
    dfloat  spiy2 = spiy*spiy;						\
    dfloat c2piy = cos(2*PI*y);						\
    dfloat v = spiy*exp(-hS*spiy2);					\
    uB = u*v;								\
    uxB = uxM;								\
    uyB = uyM;								\
  }

/* Neumann boundary condition   */
#define ellipticNeumannCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)	\
  {									\
  dfloat  spix = sin(PI*x);						\
  dfloat  spix2 = spix*spix;						\
  dfloat  cpix = cos(PI*x);						\
  dfloat u = spix*exp(-hS*spix2);					\
  dfloat  spiy = sin(PI*y);						\
  dfloat  spiy2 = spiy*spiy;						\
  dfloat  cpiy = cos(PI*y);						\
  dfloat v = spiy*exp(-hS*spiy2);					\
  uB  = uM;								\
  uxB = PI*(cpix -hS*2.*spix2*cpix)*exp(-hS*spix2)*v;			\
  uyB = PI*(cpiy -hS*2.*spiy2*cpiy)*exp(-hS*spiy2)*u;			\
  }

#define ellipticExactSolution2D(x,y,lambda,exu,exux,exuy)		\
  {									\
    dfloat  spix = sin(PI*x);						\
    dfloat  spix2 = spix*spix;						\
    dfloat  cpix = cos(PI*x);						\
    dfloat c2pix = cos(2*PI*x);						\
    dfloat  spiy = sin(PI*y);						\
    dfloat  spiy2 = spiy*spiy;						\
    dfloat  cpiy = cos(PI*y);						\
    dfloat c2piy = cos(2*PI*y);						\
    dfloat u = spix*exp(-hS*spix2);					\
    dfloat v = spiy*exp(-hS*spiy2);					\
    dfloat dudx = PI*(cpix - hS*2.*spix2*cpix)*exp(-hS*spix2);		\
    dfloat dvdy = PI*(cpiy - hS*2.*spiy2*cpiy)*exp(-hS*spiy2);		\
    *exu = u*v;								\
    *exux = dudx*v;							\
    *exuy = dvdy*u;							\
  }

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

#define hornS 10

/* forcing function   */
#define ellipticForcing1D(x, lambda, f)					\
  {									\
    dfloat  spi = sin(PI*x);						\
    dfloat  spi2 = spi*spi;						\
    dfloat c2pi = cos(2*PI*x);						\
    dfloat u = spi*exp(-hornS*spi2);					\
    f = -( -PI*PI*exp(-hornS*spi2)*spi*(2.*hornS + hornS*hornS*c2pi*c2pi + 4.*hornS*c2pi - hornS*hornS + 1.) ) + lambda*u; \
  }

/* Dirichlet boundary condition   */
#define ellipticDirichletCondition1D(x,nx,uM,uxM,uB,uxB)		\
  {									\
    dfloat  spi = sin(PI*x);						\
    dfloat  spi2 = spi*spi;						\
    dfloat c2pi = cos(2*PI*x);						\
    uB = spi*exp(-hornS*spi2);						\
    uxB = uxM;								\
  }

/* Neumann boundary condition   */
#define ellipticNeumannCondition1D(x,nx,uM,uxM,uB,uxB)			\
  {									\
  dfloat  spi = sin(PI*x);						\
  dfloat  spi2 = spi*spi;						\
  dfloat  cpi = cos(PI*x);						\
  uB  = uM;								\
  uxB = PI*(cpi -hornS*2.*spi2*cpi)*exp(-hornS*spi2);			\
  }

#define ellipticExactSolution1D(x,lambda,exu,exux)			\
  {									\
    dfloat  spi = sin(PI*x);						\
    dfloat  spi2 = spi*spi;						\
    dfloat  cpi = cos(PI*x);						\
    dfloat c2pi = cos(2*PI*x);						\
    *exu = spi*exp(-hornS*spi2);					\
    *exux = PI*(cpi -hornS*2.*spi2*cpi)*exp(-hornS*spi2);		\
  }

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

#define ellipticForcing2D(x, y, lambda, f) \
  { \
  f = 0.;                                     \
  }


/* Homogeneous Dirichlet boundary condition   */
#define ellipticDirichletCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = 0.f;   \
    uxB = uxM;   \
    uyB = uyM;   \
  }

/* Homogeneous Neumann boundary condition   */
#define ellipticNeumannCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = uM;    \
    uxB = 0.f;   \
    uyB = 0.f;   \
  }

#define waveForcingFunction2D(x, y, sigma, omega, f)                   \
  {                                                                   \
    /* taken from https://arxiv.org/pdf/1910.10148.pdf (4.2.1) */     \
    dfloat tmp = -sigma*((xn-0.01)*(xn-0.01)+(yn-0.015)*(yn-0.015));  \
    f = -omega*omega*exp(tmp);                                       \
  }

#define waveInitialConditionsFunction2D(t, x, y, d, p)  \
  {                                                     \
    d = 0.;                                             \
    p = 0.;                                             \
  }


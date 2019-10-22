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

#define PI 3.14159265358979323846

/* coefficient function */
#define ellipticCoefficient2D(x, y, lambda_0, lambda_1)    \
  {                                         \
    lambda_0 = sin(PI*x)*sin(PI*y) + 1.0; \
    lambda_1 = sin(PI*x)*sin(PI*y) + 1.0; \
  }

/* forcing function   */
#define ellipticForcing2D(x, y, lambda, f)  \
  {                                         \
    dfloat sxx  = sin(PI*x)*sin(PI*x); \
    dfloat syy  = sin(PI*y)*sin(PI*y); \
    dfloat sxy  = sin(PI*x)*sin(PI*y); \
    f = -PI*PI*(sxx + syy) + sxy*(1.f + 2.f*PI*PI + sxy + 4.f*PI*PI*sxy);  \
  }

/* Dirichlet boundary condition   */
#define ellipticDirichletCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = sin(PI*x)*sin(PI*y);   \
    uxB = uxM;   \
    uyB = uyM;   \
  }

/* Neumann boundary condition   */
#define ellipticNeumannCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    dfloat lambda_0 = sin(PI*x)*sin(PI*y) + 1.0; \
    uB  = uM;    \
    uxB = -lambda_0*PI*cos(PI*x)*sin(PI*y); \
    uyB = -lambda_0*PI*sin(PI*x)*cos(PI*y); \
  }
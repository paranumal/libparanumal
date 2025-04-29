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

/* forcing function   */
#define stressForcing2D(x, y, lambda, f)  \
  {                                         \
    f[0]  = (2*PI*PI+lambda)*sin(PI*x)*sin(PI*y);   \
    f[1]  = 0;  \
  }

/* Dirichlet boundary condition   */
#define stressDirichletCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)	\
  {									\
    uB[0]  = sin(PI*x)*sin(PI*y);					\
    uB[1]  = 0;								\
    uxB[0] = uxM[0];							\
    uxB[1] = uxM[1];							\
    uyB[0] = uyM[0];							\
    uyB[1] = uyM[1];							\
  }

/* Neumann boundary condition   */
#define stressNeumannCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {								   \
    uB[0]  = uM[0];						   \
    uB[1]  = uM[1];						   \
    uxB[0] = -PI*cos(PI*x)*sin(PI*y);				   \
    uxB[1] = -PI*cos(PI*x)*sin(PI*y);				   \
    uyB[0] = -PI*sin(PI*x)*cos(PI*y);				   \
    uyB[1] = -PI*sin(PI*x)*cos(PI*y);				   \
  }

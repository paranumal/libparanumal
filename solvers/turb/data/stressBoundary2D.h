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

/* Dirichlet 1, Neumann 2, Robin 3 (defaulted to Neumann for now) */
#define stressBoundaryConditions2D(bc,x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
{									\
  if (bc==1) stressDirichletCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB) \
  else if(bc==2) stressNeumannCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB) \
  else           stressNeumannCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB) \
}


/*-----------------------------------------------------------------------------------------------*/
/* Homogeneuous Boundary conditions used in stressAx.
/*-----------------------------------------------------------------------------------------------*/

/* Homogeneous Dirichlet boundary condition   */
#define stressHomogeneousDirichlet2D(uM,uxM,uyM,uB,uxB,uyB)  \
  {							     \
    uB[0]  = 0.f;					     \
    uB[1]  = 0.f;					     \
    uxB[0] = uxM[0];					     \
    uyB[0] = uyM[0];					     \
    uxB[1] = uxM[1];					     \
    uyB[1] = uyM[1];					     \
  }

/* Homogeneous Neumann boundary condition   */
#define stressHomogeneousNeumann2D(uM,uxM,uyM,uB,uxB,uyB)	\
  {								\
  uB[0] = uM[0];						\
  uB[1] = uM[1];						\
  uxB[0] = 0.f;							\
  uxB[1] = 0.f;							\
  uyB[0] = 0.f;							\
  uyB[1] = 0.f;							\
  }

/* Dirichlet 1, Neumann 2, Robin 3 (defaulted to Neumann for now) */
#define stressHomogeneousBC2D(bc,uM,uxM,uyM,uB,uxB,uyB)			\
  {									\
    if     (bc==1) stressHomogeneousDirichlet2D(uM,uxM,uyM,uB,uxB,uyB)	\
    else if(bc==2) stressHomogeneousNeumann2D(uM,uxM,uyM,uB,uxB,uyB)	\
    else           stressHomogeneousNeumann2D(uM,uxM,uyM,uB,uxB,uyB)	\
  }


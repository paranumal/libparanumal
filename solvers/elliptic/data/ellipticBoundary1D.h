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
#define ellipticBoundaryConditions1D(bc,x,nx,uM,uxM,uB,uxB)  \
  {                 \
    if     (bc==1) ellipticDirichletCondition1D(x,nx,uM,uxM,uB,uxB) \
    else if(bc==2) ellipticNeumannCondition1D(x,nx,uM,uxM,uB,uxB)  \
    else           ellipticNeumannCondition1D(x,nx,uM,uxM,uB,uxB)  \
  }


/*-----------------------------------------------------------------------------------------------*/
/* Homogeneuous Boundary conditions used in ellipticAx.
/*-----------------------------------------------------------------------------------------------*/

/* Homogeneous Dirichlet boundary condition   */
#define ellipticHomogeneousDirichlet1D(uM,uxM,uB,uxB)  \
  {              \
    uB  = 0.f;   \
    uxB = uxM;   \
  }

/* Homogeneous Neumann boundary condition   */
#define ellipticHomogeneousNeumann1D(uM,uxM,uB,uxB)  \
  {              \
    uB = uM;     \
    uxB = 0.f;   \
  }

/* Dirichlet 1, Neumann 2, Robin 3 (defaulted to Neumann for now) */
#define ellipticHomogeneousBC1D(bc,uM,uxM,uB,uxB)  \
  {                 \
    if     (bc==1) ellipticHomogeneousDirichlet1D(uM,uxM,uB,uxB) \
    else if(bc==2) ellipticHomogeneousNeumann1D(uM,uxM,uB,uxB)  \
    else           ellipticHomogeneousNeumann1D(uM,uxM,uB,uxB)  \
  }


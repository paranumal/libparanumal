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

/* wall 1, inflow 2, outflow 3 */

// Weakly Impose Nonlinear term BCs
#define insAdvectionBoundaryConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
  {	\
    dfloat lambda = 1.f/(2.f * p_nu) - occaSqrt(1.f/(4.f*p_nu*p_nu) + 4.f*M_PI*M_PI);\
    if(bc==1){								\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
    } else if(bc==2){							\
      *(uB) = 1.f - exp(lambda*x)*cos(2.f*M_PI*y);   \
      *(vB) =  lambda/(2.f*M_PI)*exp(lambda*x)*sin(2.f*M_PI*y);  	\
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;							\
    }									\
  }

#define insDivergenceBoundaryConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
  {	\
    dfloat lambda = 1.f/(2.f * p_nu) - occaSqrt(1.f/(4.f*p_nu*p_nu) + 4.f*M_PI*M_PI);\
    if(bc==1){								\
      *(uB)= 0.f;							\
      *(vB)= 0.f;							\
    } else if(bc==2){							\
      *(uB) = 1.f - exp(lambda*x)*cos(2.f*M_PI*y);   		\
      *(vB) = lambda/(2.f*M_PI)*exp(lambda*x)*sin(2.f*M_PI*y);\
    } else if(bc==3){							\
      *(uB) = uM;							\
      *(vB) = vM;							\
    }									\
  }

// Gradient only applies to Pressure and Pressure Incremament
// Boundary Conditions are implemented in strong form
#define insGradientBoundaryConditions2D(bc,t,x,y,nx,ny,pM,pB)	\
  {	\
    dfloat lambda = 1.f/(2.f * p_nu) - occaSqrt(1.f/(4.f*p_nu*p_nu) + 4.f*M_PI*M_PI);\
    if(bc==1){							\
      *(pB) = pM;						\
    } else if(bc==2){						\
      *(pB) = pM;						\
    } else if(bc==3){						\
      *(pB) = 0.5f*(1.f- exp(2.f*lambda*x));\
    }								\
  }

#define insHelmholtzBoundaryConditionsIpdg2D(bc, t, x, y, nx, ny, uB, uxB, uyB, vB, vxB, vyB) \
  {	\
    dfloat lambda = 1.f/(2.f*p_nu)-occaSqrt(1.f/(4.f*p_nu*p_nu) + 4.f*M_PI*M_PI);\
    if((bc==1)||(bc==4)){						\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
									\
      *(uxB) = 0.f;							\
      *(uyB) = 0.f;							\
      *(vxB) = 0.f;							\
      *(vyB) = 0.f;							\
    } else if(bc==2){							\
									\
      *(uB) = 1.f - exp(lambda*x)*cos(2.f*M_PI*y);   	\
      *(vB) = lambda/(2.f*M_PI)*exp(lambda*x)*sin(2.f*M_PI*y);	\
									\
      *(uxB) = 0.f;							\
      *(uyB) = 0.f;							\
      *(vxB) = 0.f;							\
      *(vyB) = 0.f;							\
    } else if(bc==3){							\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
      *(uxB) = -lambda*exp(lambda*x)*cos(2.f*M_PI*y);\
      *(uyB) = 2.f*M_PI*exp(lambda*x)*sin(2.f*M_PI*y); \
      *(vxB) = lambda*lambda/(2.f*M_PI)*exp(lambda*x)*sin(2.f*M_PI*y);   \
      *(vyB) = lambda*exp(lambda*x)*cos(2.f*M_PI*y);              \
    }									\
  }


// Compute bcs for P increment: if c0 = 0 give Pr BCs, zero if time independent
#define insPoissonBoundaryConditions2D(bc,t,dt,x,y,nx,ny,pB,pxB,pyB)	\
  {	\
    dfloat lambda = 1.f/(2.f * p_nu) - occaSqrt(1.f/(4.f*p_nu*p_nu) + 4.f*M_PI*M_PI);\
    if((bc==1)||(bc==4)){						\
      *(pB) = 0.f;							\
									\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    }									\
    if(bc==2){								\
      *(pB)  = 0.f;							\
									\
      *(pxB) =  0.f; \
      *(pyB) =  0.f; \
    }									\
    if(bc==3){								\
      *(pB) =  0.5f*(1.f- exp(2.f*lambda*x));\
									\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    }									\
  }



// Compute bcs for P increment
#define insPoissonNeumannTimeDerivative2D(bc,t,x,y,dpdt)  \
  { \
    if((bc==1)||(bc==4)||(bc==2) ){           \
      *(dpdt) =0.f; \
    }                 \
    if(bc==3){                \
      *(dpdt) = 0.0; \
    }                 \
  }

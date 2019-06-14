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

// Initial conditions 
#define insFlowField2D(t,x,y,u,v,p)   \
  {                                   \
    dfloat lambda = 0.5*1.0/p_nu - sqrt(1.0/(4.0*p_nu*p_nu) + 4.f*M_PI*M_PI); \
    *(u) = 1.0 - exp(lambda*x)*cos(2.0*M_PI*y);\
    *(v) = 0.5*lambda/M_PI*exp(lambda*x)*sin(2.0*M_PI*y);\
    *(p) = 0.5*(1.0 - exp(2.0*lambda*x)); \
  }   

// Boundary conditions
// BC == 1 is already handled in the solver...  
// Default u+ = u-, modify if it is different  
#define insVelocityDirichletConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
{                                                       \
    dfloat lambda = 0.5*1.0/p_nu - sqrt(1.0/(4.0*p_nu*p_nu) + 4.0*M_PI*M_PI); \
  if(bc==2){                                            \
    *(uB) = 1.0 - exp(lambda*x)*cos(2.0*M_PI*y);\
    *(vB) = 0.5*lambda/M_PI*exp(lambda*x)*sin(2.0*M_PI*y);\
  }else if(bc==4){                  \
    *(uB) = 0.0;                    \
  } else if(bc==5){                 \
    *(vB) = 0.0;                    \
  }                                 \
}

// default dudx = 0.0; modify only if you have a specific outflow bc 
#define insVelocityNeumannConditions2D(bc, t, x, y, nx, ny, uxM, uyM, vxM, vyM, uxB, uyB, vxB, vyB) \
{                                          \
    dfloat lambda = 0.5*1.0/p_nu - sqrt(1.0/(4.0*p_nu*p_nu) + 4.0*M_PI*M_PI); \
  if(bc==3){                               \
    *(uxB) =-lambda*exp(lambda*x)*cos(2.0*M_PI*y);\
    *(uyB) = 2.0*M_PI*exp(lambda*x)*sin(2.0*M_PI*y);\
    *(vxB) = 0.5*lambda*lambda/M_PI*exp(lambda*x)*sin(2.0*M_PI*y);\
    *(vyB) = lambda*exp(lambda*x)*cos(2.0*M_PI*y); \
  }                                        \
}

// default is pB = pM; modify only if you have a specific outflow bc
#define insPressureDirichletConditions2D(bc, t, x, y, nx, ny, pM, pB) \
{                                   \
  dfloat lambda = 0.5*1.0/p_nu - sqrt(1.0/(4.0*p_nu*p_nu) + 4.0*M_PI*M_PI); \
  if(bc==3){                        \
    *(pB) = 0.5f*(1.f - exp(2.f*lambda*x)); \
   } \
}


// default is dPdx = 0.0; I think we dont need that....  
#define insPressureNeumannConditions2D(bc, t, x, y, nx, ny, pxM, pyM, pxB, pyB) \
{                                          \
  if(bc==3){                               \
  }                                        \
}




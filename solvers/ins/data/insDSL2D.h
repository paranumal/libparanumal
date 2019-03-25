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

// Contributed by Stefen Kerkemmeier

// Initial conditions 
#define insFlowField2D(t,x,y,u,v,p)		\
  {						\
    if(y<=0.5){					\
      *(u) = tanh(30.f*(y-0.25f));		\
    }else{					\
      *(u) = tanh(30.f*(0.75f-y));		\
    }						\
    *(v) = 0.05*sin(2*M_PI*x);			\
    *(p) = 0.f;					\
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
#define insVelocityDirichletConditions2D(bc, t, x, y,nx, ny,uM, vM, uB, vB) \
  {									\
  }

#define insVelocityNeumannConditions2D(bc, t, x, y, nx, ny, uxM, uyM, vxM, vyM,uxB, uyB,vxB, vyB) \
  {									\
  }


#define insPressureDirichletConditions2D(bc, t, x, y, nx, ny, pM, pB)	\
  {									\
  }

#define insPressureNeumannConditions2D(bc, t, x, y, nx, ny, pxM, pyM, pxB, pyB) \
  {									\
    *(pxB) = 0.f;							\
    *(pyB) = 0.f;							\
  }

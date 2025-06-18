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

#define p_ubar 1.0
#define p_vbar 0.0
#define p_pbar 1.0

#define p_kbar 1.0
#define p_taubar 1.0

#define copysign(a) ( ((a)<0) ? -1.f : 1.f )

// Initial conditions
#define insInitialConditions2D(nu,t,x,y,u,v,k,tau,p)	\
  {						\
    *(u) = p_ubar;				\
    *(v) = p_vbar;				\
    *(k) = p_kbar;				\
    *(tau) = p_taubar;				\
    *(p) = p_pbar;				\
  }

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define insVelocityDirichletConditions2D(bc, ndotUe, nu, t, x, y, nx, ny, uM, uB) \
  {									\
  if(bc==1){ /* wall */							\
    uB[0] = 0.f;							\
    uB[1] = 0.f;							\
    uB[2] = 0.f;							\
    uB[3] = 0.f;							\
  } else if(bc==2){ /* inflow */					\
    uB[0] = p_ubar;							\
    uB[1] = p_vbar;							\
    uB[2] = p_kbar;							\
    uB[3] = p_taubar;							\
  } else if(bc==3){							\
    /* penalize tangential part */					\
    uB[0] = (ndotUe>0) ? uM[0] : 0;					\
    uB[1] = (ndotUe>0) ? uM[1] : 0;					\
    uB[2] = (ndotUe>0) ? uM[2] : 0;					\
    uB[3] = (ndotUe>0) ? uM[3] : 0;					\
  }									\
  }

// no slip condition
#define insVelocityNeumannConditions2D(bc, nu, t, x, y, nx, ny, QxM, QyM, QxB, QyB) \
  {									\
    if(bc==1 || bc==2){							\
      for(int fld=0;fld<p_NVfields;++fld){				\
	QxB[fld] = QxM[fld];						\
	QyB[fld] = QyM[fld];						\
      }									\
    } else if(bc==3){							\
      for(int fld=0;fld<p_NVfields;++fld){				\
	QxB[fld] = 0.f;							\
	QyB[fld] = 0.f;							\
      }									\
    }									\
    }


#define insPressureDirichletConditions2D(bc, nu, t, x, y, nx, ny, pM, pB) \
  {									\
    if(bc==1 || bc==2){							\
      *(pB) = pM;							\
    } else if(bc==3){							\
      *(pB) = 0.0;							\
    } else if(bc==4){							\
      *(pB) = pM;							\
    } else if(bc==5){							\
      *(pB) = pM;							\
    }									\
  }

#define insPressureNeumannConditions2D(bc, nu, t, x, y, nx, ny, pxM, pyM, pxB, pyB) \
  {									\
    if(bc==1 || bc==2){							\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    } else if(bc==3){							\
      *(pxB) = pxM;							\
      *(pyB) = pyM;							\
    } else if(bc==4){							\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    } else if(bc==5){							\
      *(pxB) = 0.f;							\
      *(pyB) = 0.f;							\
    }									\
  }

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
#define cdsFlowField2D(t,x,y,u,v) \
  {				       \
    *(u) = -y;		       \
    *(v) =  x;		       \
  }   
#define cdsScalarField2D(t,x,y,s)	\
  {					\
    dfloat mtime = t + M_PI; \
    dfloat cond  = 0.001; \
    dfloat xc = 0.00f;			\
    dfloat yc = 0.50f;     \
    dfloat xt = xc*cos(mtime)  - yc*sin(mtime);      \
    dfloat yt = -xc*sin(mtime) + yc*cos(mtime);     \
    dfloat r2 = (x-xt)*(x-xt) + (y-yt)*(y-yt);      \
    dfloat sExact =  1.f / (4.f*M_PI*cond*mtime) * exp(-r2/ (4.f*cond*mtime));	\
    *(s) = sExact;				\
  }   


// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
#define cdsDirichletConditions2D(bc, t, x, y, nx, ny, sM, sB) \
{                                   \
    dfloat mtime = t + M_PI/2; \
    dfloat cond  = 0.001; \
    dfloat xc = 0.00f;      \
    dfloat yc = 0.50f;     \
    dfloat xt = xc*cos(mtime) - yc*sin(mtime);      \
    dfloat yt = -xc*sin(mtime) + yc*cos(mtime);     \
    dfloat r2 = (x-xt)*(x-xt) + (y-yt)*(y-yt);      \
    dfloat sExact =  1.f / (4.f*M_PI*cond*mtime) * exp(-r2/ (4.f*cond*mtime));  \
  if(bc==1){                        \
    *(sB) = 0.f;                    \
  } else if(bc==2){                 \
    *(sB) = sExact;		    \
  } else if(bc==3){                 \
    *(sB) = sM;                     \
  } else if(bc==4||bc==5||bc==6){   \
    *(sB) = sM; \
  }                                 \
}

#define cdsNeumannConditions2D(bc, t, x, y, nx, ny, sxM, syM, sxB, syB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(sxB) = sxM;                          \
    *(syB) = syM;                          \
  } else if(bc==3){                        \
    *(sxB) = 0.f;                          \
    *(syB) = 0.f;                          \
  } else if(bc==4||bc==5||bc==6){          \
    *(sxB) = nx*nx*sxM;                    \
    *(syB) = ny*ny*syM;                    \
  }                                        \
}

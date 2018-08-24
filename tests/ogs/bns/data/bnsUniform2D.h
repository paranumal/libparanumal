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


// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define boundaryConditionsPML2D(bc, t, x, y, nx, ny, intfx, intfy, q1M, q2M, q3M, q4M, q5M, q6M, q1B, q2B, q3B, q4B, q5B, q6B) \
  {									\
    if(bc==1){								\
      *(q1B) =  q1M;							\
      *(q2B) = -q2M;							\
      *(q3B) = -q3M;							\
      *(q4B) =  q4M;							\
      *(q5B) =  q5M;							\
      *(q6B) =  q6M;							\
    } else if(bc==4||bc==5){						\
      *(q1B) = q1M;							\
      *(q2B) = q2M-2.f*(nx*q2M+ny*q3M)*nx;				\
      *(q3B) = q3M-2.f*(nx*q2M+ny*q3M)*ny;				\
      *(q4B) =  q4M;							\
      *(q5B) =  q5M;							\
      *(q6B) =  q6M;							\
    }									\
    else {								\
      *(q1B) = 2.f*p_q1bar - q1M;					\
      *(q2B) = 2.f*p_q1bar*intfx*p_isqrtRT - q2M;			\
      *(q3B) = 2.f*p_q1bar*intfy*p_isqrtRT - q3M;			\
      *(q4B) = 2.f*p_q1bar*intfx*intfy*p_isqrtRT*p_isqrtRT-q4M;		\
      *(q5B) = 2.f*p_q1bar*intfx*intfx*p_isqrtRT*p_isqrtRT*p_invsqrt2-q5M; \
      *(q6B) = 2.f*p_q1bar*intfy*intfy*p_isqrtRT*p_isqrtRT*p_invsqrt2-q6M; \
    }									\
  }

// Boundary conditions: Did not check yet
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define boundaryConditions2D(bc, t, x, y, nx, ny, intfx, intfy, q1M, q2M, q3M, q4M, q5M, q6M, q1B, q2B, q3B, q4B, q5B, q6B) \
  {									\
    if(bc==1){								\
      *(q1B) =  q1M;							\
      *(q2B) = -q2M;							\
      *(q3B) = -q3M;							\
      *(q4B) =  q4M;							\
      *(q5B) =  q5M;							\
      *(q6B) =  q6M;							\
    } else if(bc==2 || bc==3){						\
      *(q1B) = 2.f*p_q1bar - q1M;					\
      *(q2B) = 2.f*p_q1bar*intfx*p_isqrtRT - q2M;			\
      *(q3B) = 2.f*p_q1bar*intfy*p_isqrtRT - q3M;			\
      *(q4B) = 2.f*p_q1bar*intfx*intfy*p_isqrtRT*p_isqrtRT-q4M;		\
      *(q5B) = 2.f*p_q1bar*intfx*intfx*p_isqrtRT*p_isqrtRT*p_invsqrt2-q5M; \
      *(q6B) = 2.f*p_q1bar*intfy*intfy*p_isqrtRT*p_isqrtRT*p_invsqrt2-q6M; \
    } else if(bc==4||bc==5){						\
      *(q1B) = q1M;							\
      *(q2B) = q2M-2.f*(nx*q2M+ny*q3M)*nx;				\
      *(q3B) = q3M-2.f*(nx*q2M+ny*q3M)*ny;				\
      *(q4B) = q4M;							\
      *(q5B) = q5M;							\
      *(q6B) = q6M;							\
    }									\
  }

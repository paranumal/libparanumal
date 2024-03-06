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

#define p_gamma (1.4)

//  (0.82)
// C2 was 1.5
// N=4,cubN=4
// C2=1.21 and tau=2 worked
// C2=0.82 and tau=2 worked
// C2=0.60 and tau=2 ok
// C2=0.45 and tau=3 ok
// N5
// C2 = 1.21 with tau=2 ok
// C2 = 0.9 with tau=2 ok
// C2 = 0.45 with tau=6 ok (switching to 2 after t=0.1 ok)
// C2 = 0.25 not stab
// art.visc N=5 with C2=0.16, 0.11 ok

//#define C2    (0.6)
//#define C2      (4)
#if 0
#define ub    1.0
#define vb    0
#define pb    (C2/p_gamma)
#define rb    1.0
#endif
#define esdgInitialConditions2D(gamma, t, cx, cy, x, y, r, u, v, p)	\
  {									\
    *(r) = p_rbar;								\
    *(u) = p_ubar;								\
    *(v) = p_vbar;								\
    *(p) = p_pbar;								\
  }

// Body force
#define esdgBodyForce2D(gamma, t, x, y, r, u, v, p, fx, fy)	\
  {								\
    *(fx) = 0.0;						\
    *(fy) = 0.0;						\
  }

// Boundary conditions (actually run with periodic) 
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define p_WALL 1
#define p_INFLOW 2
#define p_OUTFLOW 3

#define esdgBoundaryConditions2D(bc, gamma,			\
                                  t, x, y, nx, ny,		\
                                  rM, uM, vM, pM,		\
                                  rB, uB, vB, pB)		\
  {								\
    switch(bc){							\
    case p_WALL:{ /* wall */					\
      *(rB) = rM;						\
      *(uB) = uM-2.*(uM*nx+vM*ny)*nx;				\
      *(vB) = vM-2.*(uM*nx+vM*ny)*ny;				\
      *(pB) = pM;						\
      break;							\
    }								\
    case p_INFLOW:{ /* inflow */				\
      *(rB) = p_rbar;						\
      *(uB) = p_ubar;						\
      *(vB) = p_vbar;						\
      *(pB) = p_pbar;						\
      break;							\
    }								\
    case p_OUTFLOW:{ /* outflow */				\
      *(rB) = p_rbar;						\
      *(uB) = p_ubar;						\
      *(vB) = p_vbar;						\
      *(pB) = p_pbar;						\
      break;							\
    }								\
    default:{							\
      printf("ARGGH - UNHANDLED BC \n");			\
    }								\
    }								\
  }

/*
      *(rB) = rM;						\
      *(uB) = uM;						\
      *(vB) = vM;						\
      *(pB) = p_pbar;						\
      */

#define esdgGradientBoundaryConditions2D(bc, gamma,		\
					  t, x, y, nx, ny,	\
					  v1M, v2M, v3M, v4M,	\
					  v1B, v2B, v3B, v4B)	\
  { /* only apply bcs for fields 1:3 */				\
    switch(bc){							\
    case p_WALL:{ /* wall */					\
      *(v2B) =  -v2M;						\
      *(v3B) =  -v3M;						\
      *(v4B) =  v4M;						\
      break;    						\
    }								\
    case p_INFLOW:{ /* inflow */				\
      *(v2B) =  v2M;						\
      *(v3B) =  v3M;						\
      *(v4B) =  v4M;						\
      break;							\
    }								\
    case p_OUTFLOW:{ /* outflow */				\
      *(v2B) =  v2M;						\
      *(v3B) =  v3M;						\
      *(v4B) =  v4M;						\
      break;							\
    }								\
    default:{							\
      printf("ARGGH - UNHANDLED BC \n");			\
    }								\
    }								\
  }


#define esdgDivergenceBoundaryConditions2D(bc, gamma,			\
					    t, x, y, nx, ny,		\
					    dv1dnM, dv2dnM, dv3dnM, dv4dnM, \
					    dv1dnB, dv2dnB, dv3dnB, dv4dnB) \
  {									\
    switch(bc){								\
    case p_WALL:{ /* wall */						\
      *(dv2dnB) =  dv2dnM;						\
      *(dv3dnB) =  dv3dnM;						\
      *(dv4dnB) = -dv4dnM;						\
      break;								\
    }									\
    case p_INFLOW:{ /* inflow */					\
      *(dv2dnB) = -dv2dnM;						\
      *(dv3dnB) = -dv3dnM;						\
      *(dv4dnB) = -dv4dnM;						\
      break;								\
    }									\
    case p_OUTFLOW:{ /* outflow */					\
      *(dv2dnB) = dv2dnM;						\
      *(dv3dnB) = dv3dnM;						\
      *(dv4dnB) = dv4dnM;						\
      break;								\
    }									\
    default:{								\
      printf("ARGGH - UNHANDLED BC \n");				\
    }									\
    }									\
  }

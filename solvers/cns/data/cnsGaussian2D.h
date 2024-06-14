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

#define p_RBAR  1.0000000
#define p_UBAR  1.0000000
#define p_VBAR  0.0
#define p_PBAR  71.428571428571431
#define p_TBAR  1.0000000


// Initial conditions (p is ignored for isothermal)
#define cnsInitialConditions2D(gamma, mu, t, x, y, r, u, v, p)                   \
{                                                                                \
  *(r) = p_RBAR*(1 + exp(-3*(x*x+y*y)));                                         \
  *(u) = p_UBAR*exp(-3*(x*x+y*y));                                               \
  *(v) = p_UBAR*exp(-3*(x*x+y*y));                                               \
  *(p) = p_PBAR*(1 + exp(-3*(x*x+y*y)));                                         \
}    



// Body force
#define cnsBodyForce2D(gamma, mu, t, x, y, r, u, v, p, fx, fy)                    \
{                                                                                 \
  *(fx) = 0.0;                                                                    \
  *(fy) = 0.0;                                                                    \
}


// Viscous Riemann Solve BC's
// ************************************************************************
#define cnsViscousBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny,   \
                                      rM, ruM, rvM, reM, rB, ruB, rvB, reB,  \
                                      drrdxM, drrdyM, drudxM, drudyM,        \
                                      drvdxM, drvdyM, dredxM, dredyM,        \
                                      drrdxB, drrdyB, drudxB, drudyB,        \
                                      drvdxB, drvdyB, dredxB, dredyB)        \
{                                                                            \
const dfloat uM  = ruM/rM;                                                   \
const dfloat vM =  rvM/rM;                                                   \
const dfloat pM  = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));               \
const dfloat unM = uM*nx + vM*ny;                                            \
const dfloat cM  = sqrt(gamma*pM/rM);                                        \
const dfloat mM  = fabs(unM/cM);                                             \
const dfloat keREF = 0.5*p_RBAR*(p_UBAR*p_UBAR + p_VBAR*p_VBAR);             \
  if(bc==11){                                                                \
    *( rB) = rM;                                                             \
    *(ruB) = 0.0;                                                            \
    *(rvB) = 0.0;                                                            \
    *(reB) = rM*CP*p_TBAR/gamma;                                             \
    *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
    *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
    *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
    *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
   }else if(bc==12){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = 0.0;                                                            \
    *(rvB) = 0.0;                                                            \
    *(reB) = reM - 0.5*(ruM*ruM + rvM*rvM)/rM;                               \
    const dfloat uxM = drudxM - uM*drrdxM;                                   \
    const dfloat uyM = drudyM - uM*drrdyM;                                   \
    const dfloat vxM = drvdxM - vM*drrdxM;                                   \
    const dfloat vyM = drvdyM - vM*drrdyM;                                   \
    const dfloat TxM = dredxM - (drrdxM*reM/rM + uM*uxM + vM*vxM);          \
    const dfloat TyM = dredyM - (drrdyM*reM/rM + uM*uyM + vM*vyM);          \
    *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
    *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
    *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
    *(dredxB) = dredxM - (nx*nx*TxM + nx*ny*TyM);                            \                                                     \
    *(dredyB) = dredyM - (nx*ny*TxM + ny*ny*TyM);                            \
  }else if(bc==13){                                                          \
    *( rB) =  rM;                                                            \
    *(ruB) =  ruM - (nx*ruM+ny*rvM)*nx;                                      \
    *(rvB) =  rvM - (nx*ruM+ny*rvM)*ny;                                      \
    *(reB) =  reM;                                                           \
    *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
    *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
    *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
    *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
  }else if(bc==20){                                                          \
    if(unM<=0){                                                              \
      if(mM > 1.0){                                                          \
        *( rB) = p_RBAR;                                                     \
        *(ruB) = p_RBAR*p_UBAR;                                              \
        *(rvB) = p_RBAR*p_VBAR;                                              \
        *(reB) = p_PBAR/(gamma -1.0) + keREF;                                \
        *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                     \
        *(drudxB) = 0.0;*(drudyB) = 0.0;                                     \
        *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                     \
        *(dredxB) = 0.0;*(dredyB) = 0.0;                                     \
      }else{                                                                 \
        *( rB) = p_RBAR;                                                     \
        *(ruB) = p_RBAR*p_UBAR;                                              \
        *(rvB) = p_RBAR*p_VBAR;                                              \
        *(reB) = reM-0.5*rM*(uM*uM+vM*vM)+keREF;                             \
        *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                    \
        *(drudxB) = 0.0; *(drudyB) = 0.0;                                    \
        *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                    \
        *(dredxB) = 0.0; *(dredyB) = 0.0;                                    \
      }                                                                      \
    }else{                                                                   \
      if(mM > 1.0){                                                          \
        *( rB) = rM;                                                         \
        *(ruB) = ruM;                                                        \
        *(rvB) = rvM;                                                        \
        *(reB) = reM;                                                        \
        *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                               \
        *(drudxB) = drudxM;*(drudyB) = drudyM;                               \
        *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                               \
        *(dredxB) = dredxM;*(dredyB) = dredyM;                               \
        }else{                                                               \
        *( rB) = rM;                                                         \
        *(ruB) = ruM;                                                        \
        *(rvB) = rvM;                                                        \
        *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);               \
        *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                    \
        *(drudxB) = 0.0; *(drudyB) = 0.0;                                    \
        *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                    \
        *(dredxB) = 0.0; *(dredyB) = 0.0;                                    \
      }                                                                      \
    }                                                                        \
  }else if(bc==21){                                                          \
    *( rB) = p_RBAR;                                                         \
    *(ruB) = p_RBAR*p_UBAR;                                                  \
    *(rvB) = p_RBAR*p_VBAR;                                                  \
    *(reB) = reM-0.5*rM*(uM*uM+vM*vM)+keREF;                                 \
    *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
    *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
    *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
    *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
  }else if(bc==22){                                                          \
    *( rB) = p_RBAR;                                                         \
    *(ruB) = p_RBAR*p_UBAR;                                                  \
    *(rvB) = p_RBAR*p_VBAR;                                                  \
    *(reB) = p_PBAR/(gamma -1.0) + keREF;                                    \
    *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                         \
    *(drudxB) = 0.0;*(drudyB) = 0.0;                                         \
    *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                         \
    *(dredxB) = 0.0;*(dredyB) = 0.0;                                         \
  }else if(bc==31){                                                          \
    *( rB) = rM;                                                             \
    *(ruB) = ruM;                                                            \
    *(rvB) = rvM;                                                            \
    *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                   \
    *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
    *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
    *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
    *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
  }else if(bc==32){                                                          \
    *( rB) = rM;                                                             \
    *(ruB) = ruM;                                                            \
    *(rvB) = rvM;                                                            \
    *(reB) = reM;                                                            \
    *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
    *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
    *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
    *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
  }                                                                          \
}



// // Initial conditions (p is ignored for isothermal)
// #define cnsInitialConditions2D(gamma, mu, t, x, y, r, u, v, p)                   \
// {                                                                                \
//   *(r) = 1.0;                                         \
//   *(u) = 2*x;                                               \
//   *(v) = 3*y;                                               \
//   *(p) = x*x*x + y*y*y;                                         \
// }    

// // Body force
// #define cnsBodyForce2D(gamma, mu, t, x, y, r, u, v, p, fx, fy)                    \
// {                                                                                 \
//   *(fx) = 0.0;                                                                    \
//   *(fy) = 0.0;                                                                    \
// }


// //Viscous Riemann Solve BC's
// // ************************************************************************
// #define cnsViscousBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny,   \
//                                       rM, ruM, rvM, reM, rB, ruB, rvB, reB,  \
//                                       drrdxM, drrdyM, drudxM, drudyM,        \
//                                       drvdxM, drvdyM, dredxM, dredyM,        \
//                                       drrdxB, drrdyB, drudxB, drudyB,        \
//                                       drvdxB, drvdyB, dredxB, dredyB)        \
// {                                                                            \
//   if(bc==11){                                                                \
//     *( rB) = 1.0;                                                            \
//     *(ruB) = x*x*x;                                                            \
//     *(rvB) = y*y;                                                              \
//     *(reB) = (x*x*x + y*y*y)/(gamma-1.0) + 0.5*1.0*(4*x*x + 9*y*y) ; \
//     *(drrdxB) = 0.0;     *(drrdyB) = 0.0;                                    \
//     *(drudxB) = 2.0; *(drudyB) = 0.0;                                    \
//     *(drvdxB) = 0.0;     *(drvdyB) = 3.0;                                    \
//     *(dredxB) = dredxM;  *(dredyB) = dredyM;                                   \
//    }                                                                         \
// }


// // Initial conditions (p is ignored for isothermal)
// #define cnsInitialConditions2D(gamma, mu, t, x, y, r, u, v, p)                   \
// {                                                                                \
//   *(r) = x+y >0 ? 1:0.1;                                         \
//   *(u) = 2*x;                                               \
//   *(v) = 3*y;                                               \
//   *(p) = x*x*x + y*y*y;                                         \
// }    

// // Body force
// #define cnsBodyForce2D(gamma, mu, t, x, y, r, u, v, p, fx, fy)                    \
// {                                                                                 \
//   *(fx) = 0.0;                                                                    \
//   *(fy) = 0.0;                                                                    \
// }


// //Viscous Riemann Solve BC's
// // ************************************************************************
// #define cnsViscousBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny,   \
//                                       rM, ruM, rvM, reM, rB, ruB, rvB, reB,  \
//                                       drrdxM, drrdyM, drudxM, drudyM,        \
//                                       drvdxM, drvdyM, dredxM, dredyM,        \
//                                       drrdxB, drrdyB, drudxB, drudyB,        \
//                                       drvdxB, drvdyB, dredxB, dredyB)        \
// {                                                                            \
//   if(bc==11){                                                                \
//     *( rB) = 1.0;                                                            \
//     *(ruB) = x*x*x;                                                            \
//     *(rvB) = y*y;                                                              \
//     *(reB) = (x*x*x + y*y*y)/(gamma-1.0) + 0.5*1.0*(4*x*x + 9*y*y) ; \
//     *(drrdxB) = 0.0;     *(drrdyB) = 0.0;                                    \
//     *(drudxB) = 2.0; *(drudyB) = 0.0;                                    \
//     *(drvdxB) = 0.0;     *(drvdyB) = 3.0;                                    \
//     *(dredxB) = dredxM;  *(dredyB) = dredyM;                                   \
//    }                                                                         \
// }
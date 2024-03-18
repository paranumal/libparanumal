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

/* ************************************************************************ */
#define MACH 0.4
/* ************************************************************************ */
#if MACH==0.1
  #define p_RBAR 1.4
  #define p_UBAR 1.0
  #define p_VBAR 0.0
  #define p_PBAR 100
  #define p_TBAR 1.00000
#elif MACH==0.15
  #define p_RBAR 1.4
  #define p_UBAR 0.984807753012208
  #define p_VBAR 0.173648177666930
  #define p_PBAR 44.44444444444444
  #define p_TBAR 1.00000
#elif MACH==0.2
  #define p_RBAR 1.0
  #define p_UBAR 0.965925826289068
  #define p_VBAR 0.258819045102521
  #define p_PBAR 17.857142857142858
  #define p_TBAR 1.00000
#elif MACH==0.4
  #define p_RBAR 1.4
  #define p_UBAR 1.0
  #define p_VBAR 0.0
  #define p_PBAR 6.25
  #define p_TBAR 1.00000
#elif MACH==0.8
  #define p_RBAR 1.79200
  #define p_UBAR 1.00000
  #define p_VBAR 0.00000
  #define p_PBAR 2.00000
  #define p_TBAR 1.00000
#elif MACH==1.2
  #define p_RBAR 2.01600
  #define p_UBAR 0.984807753012208
  #define p_VBAR 0.173648177666930
  #define p_PBAR 1.00000
  #define p_TBAR 1.00000
#elif MACH==2.0
  #define p_RBAR 5.60000
  #define p_UBAR 0.984807753012208
  #define p_VBAR 0.173648177666930
  #define p_PBAR 1.00000
  #define p_TBAR 1.00000
#elif MACH==3.0
  #define p_RBAR 12.6000
  #define p_UBAR 0.984807753012208
  #define p_VBAR 0.173648177666930
  #define p_PBAR 1.00000
  #define p_TBAR 1.00000
#elif MACH==4.0
  #define p_RBAR 22.4000
  #define p_UBAR 0.984807753012208
  #define p_VBAR 0.173648177666930
  #define p_PBAR 1.00000
  #define p_TBAR 1.00000
#elif MACH==5.0
  #define p_RBAR 35.0
  #define p_UBAR 1.0
  #define p_VBAR 0.0
  #define p_PBAR 1.0
  #define p_TBAR 1.00000
#else
  #error "Reference State is not defined in the BC data file."
#endif
/***************************************************************************/
// BCs Summary
// 1-1 : No-Slip Wall : Isothermal: Navier-Stokes Solver
// 1-2 : No-Slip Wall : Adiabatic : Navier-Stokes Solver
// 1-3 : Slip Wall    : Adiabatic : Euler Wall 
// 1-3 : Symmetry     : Navier-Stokes Symmetry 
// 2-0 : FarField     : Supersonic/Subsonic Inlet/Outlet
// 2-1 : Inflow       : Subsonic Inlet
// 2-2 : Inflow       : Supersonic Inlet
// 3-1 : Outflow      : Subsonic Outlet
// 3-2 : Outflow      : Supersonic Outlet
// 4-1 : WeakDrichlet : Define
/***************************************************************************/

/***************************************************************************/
// ************************************************************************
// Pressure Riemann Solver on BC based on NASA Report
// AK: We dont use in the current form!!!!
/* ***********************************************************************/
#define PressureRiemann2D(gamma, R, CP, CV, UN, CN, rM, uM, vM, pM, PR){     \
const dfloat pfunc= 2.0*gamma/(gamma -1.0);                                  \
const dfloat AR   = 2.0/(rM*(gamma+1.0));                                    \
const dfloat BR   = pM*(gamma -1.0)/(gamma+1.0);                             \
const dfloat PR1  = pM*pow(max(0.0001, 1.0+0.5*(gamma -1.0)*UN/CN), pfunc);  \
const dfloat PR2  = pM+0.5*UN/AR*(UN+sqrt(UN*UN+4.0*AR*(pM+BR)));            \
*(PR)   = (UN<=0) ? PR1 : PR2;                                               \
}
/* ************************************************************************ */


// Initial conditions (p is ignored for isothermal)
// ************************************************************************
#define cnsInitialConditions2D(gamma, mu, t, x, y, r, u, v, p) \
{                                         \
  *(r) = p_RBAR;           \
  *(u) = p_UBAR;           \
  *(v) = p_VBAR;           \
  *(p) = p_PBAR;           \
}

// Body force
// ************************************************************************
#define cnsBodyForce2D(gamma, mu, t, x, y, r, u, v, p, fx, fy) \
{                                                   \
  *(fx) = 0.0;                                      \
  *(fy) = 0.0;                                      \
}


    // *(reB) = rM*CP*p_TBAR/gamma;                                             


// Viscous Riemann Solve BC's
// ************************************************************************
#define cnsViscousBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny, \
                                      rM, ruM, rvM, reM, rB, ruB, rvB, reB,\
                                      drrdxM, drrdyM, drudxM, drudyM,        \
                                      drvdxM, drvdyM, dredxM, dredyM,        \
                                      drrdxB, drrdyB, drudxB, drudyB,        \
                                      drvdxB, drvdyB, dredxB, dredyB)        \
{                                                                            \
const dfloat uM  = ruM/rM;                                                   \
const dfloat vM  =  rvM/rM;                                                   \
const dfloat pM  = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));               \
const dfloat unM = uM*nx + vM*ny;                                            \
const dfloat cM  = sqrt(gamma*pM/rM);                                        \
const dfloat mM  = fabs(unM/cM);                                             \
const dfloat keREF = 0.5*p_RBAR*(p_UBAR*p_UBAR + p_VBAR*p_VBAR);             \
  if(bc==11){                                                                \
    *( rB) = rM;                                                             \
    *(ruB) = -ruM;                                                           \
    *(rvB) = -rvM;                                                           \
    *(reB) = rM*CV*p_TBAR + 0.5*(ruM*ruM + rvM*rvM)/rM;                      \
    *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
    *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
    *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
    *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
   }else if(bc==12){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = -ruM;                                                           \
    *(rvB) = -rvM;                                                           \
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
  } else if(bc==41){                                                         \
    *( rB) = p_RBAR;                                                         \
    *(ruB) = p_RBAR*p_UBAR;                                                  \
    *(rvB) = p_RBAR*p_VBAR;                                                  \
    *(reB) = p_PBAR/(gamma -1.0) + keREF;                                    \
    *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
    *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
    *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
    *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
  }                                                                          \
}



// // Inviscid Riemann Solve BC's
// // ************************************************************************
// #define cnsInviscidBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny, \
//                                         rM, ruM, rvM, reM, rB, ruB, rvB, reB) \
// {                                                                            \
// const dfloat uM    = ruM/rM;                                                 \
// const dfloat vM    = rvM/rM;                                                 \
// const dfloat pM    = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));             \
// const dfloat unM   = uM*nx + vM*ny;                                          \
// const dfloat cM    = sqrt(gamma*pM/rM);                                      \
// const dfloat mM    = fabs(unM/cM);                                           \
// const dfloat keREF = 0.5*p_RBAR*(p_UBAR*p_UBAR + p_VBAR*p_VBAR);             \
//   if(bc==11){                                                                \
//     *( rB) = rM;                                                             \
//     *(ruB) = -ruM;                                                           \
//     *(rvB) = -rvM;                                                           \
//     *(reB) = rM*CP*p_TBAR/gamma + 0.5*(ruM*ruM + rvM*rvM)/rM;                \
//    }else if(bc==12){                                                         \
//     *( rB) = rM;                                                             \
//     *(ruB) = -ruM;                                                           \
//     *(rvB) = -rvM;                                                           \
//     *(reB) = reM;                                                            \
//   }else if(bc==13){                                                          \
//     *( rB) =  rM;                                                            \
//     *(ruB) =  ruM - 2.0*(nx*ruM+ny*rvM)*nx;                                  \
//     *(rvB) =  rvM - 2.0*(nx*ruM+ny*rvM)*ny;                                  \
//     *(reB) =  reM;                                                           \
//    }else if(bc==20){                                                         \
//     if(unM<=0){                                                              \
//       if(mM > 1.0){                                                          \
//         *( rB) = p_RBAR;                                                     \
//         *(ruB) = p_RBAR*p_UBAR;                                              \
//         *(rvB) = p_RBAR*p_VBAR;                                              \
//         *(reB) = p_PBAR/(gamma -1.0) + keREF;                                \
//       }else{                                                                 \
//         *( rB) = p_RBAR;                                                     \
//         *(ruB) = p_RBAR*p_UBAR;                                              \
//         *(rvB) = p_RBAR*p_VBAR;                                              \
//         *(reB) = reM-0.5*rM*(uM*uM+vM*vM)+keREF ;                            \
//       }                                                                      \
//     }else{                                                                   \
//       if(mM > 1.0){                                                          \
//       *( rB) = rM;                                                           \
//       *(ruB) = ruM;                                                          \
//       *(rvB) = rvM;                                                          \
//       *(reB) = reM;                                                          \
//       }else{                                                                 \
//       *( rB) = rM;                                                           \
//       *(ruB) = ruM;                                                          \
//       *(rvB) = rvM;                                                          \
//       *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                 \
//       }                                                                      \
//     }                                                                        \
//    }else if(bc==21){                                                         \
//     *( rB) = p_RBAR;                                                         \
//     *(ruB) = p_RBAR*p_UBAR;                                                  \
//     *(rvB) = p_RBAR*p_VBAR;                                                  \
//     *(reB) = reM-0.5*rM*(uM*uM+vM*vM)+keREF ;                                \
//    }else if(bc==22){                                                         \
//     *( rB) = p_RBAR;                                                         \
//     *(ruB) = p_RBAR*p_UBAR;                                                  \
//     *(rvB) = p_RBAR*p_VBAR;                                                  \
//     *(reB) = p_PBAR/(gamma -1.0) + keREF;                                    \
//    }else if(bc==31){                                                         \
//     *( rB) = rM;                                                             \
//     *(ruB) = ruM;                                                            \
//     *(rvB) = rvM;                                                            \
//     *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                   \
//    }else if(bc==32){                                                         \
//     *( rB) = rM;                                                             \
//     *(ruB) = ruM;                                                            \
//     *(rvB) = rvM;                                                            \
//     *(reB) = reM;                                                            \
//   }else if(bc==41){                                                          \
//     *( rB) = p_RBAR;                                                         \
//     *(ruB) = p_RBAR*p_UBAR;                                                  \
//     *(rvB) = p_RBAR*p_VBAR;                                                  \
//     *(reB) = p_PBAR/(gamma -1.0) + keREF;                                    \
//  }                                                                           \
// }


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
/***************************************************************************/
// BCs Summary
// 1-1 : No-Slip Isothermall Wall: Isothermal: Navier-Stokes Solver
// 1-2 : No-Slip Wall    : Adiabatic : Navier-Stokes Solver
// 1-3 : Slip Wall       : Adiabatic : Euler Wall 
// 1-3 : Symmetry        : Navier-Stokes Symmetry 
// 2-0 : FarField        : Riemann Farfield / Euler-Navier-Stokes
// 2-1 : Inflow          : Subsonic Inlet
// 2-2 : Inflow          : Supersonic Inlet
// 3-1 : Outflow         : Subsonic Outlet
// 3-2 : Outflow         : Supersonic Outlet
// 3-3 : Pressure-Outlet : Pressure Outlet:Navier-Stokes
// 4-1 : WeakDrichlet    : Define
/***************************************************************************/

/***************************************************************************/
// ************************************************************************
// Pressure Riemann Solver on BC based on NASA Report
// AK: We dont use in the current form!!!!
// /* ***********************************************************************/
// #define PressureRiemann2D(gamma, R, CP, CV, UN, CN, rM, uM, vM, pM, PR){     \
// const dfloat pfunc= 2.0*gamma/(gamma -1.0);                                  \
// const dfloat AR   = 2.0/(rM*(gamma+1.0));                                    \
// const dfloat BR   = pM*(gamma -1.0)/(gamma+1.0);                             \
// const dfloat PR1  = pM*pow(max(0.0001, 1.0+0.5*(gamma -1.0)*UN/CN), pfunc);  \
// const dfloat PR2  = pM+0.5*UN/AR*(UN+sqrt(UN*UN+4.0*AR*(pM+BR)));            \
// *(PR)   = (UN<=0) ? PR1 : PR2;                                               \
// }
/* ************************************************************************ */


// Initial conditions (p is ignored for isothermal)
// ************************************************************************
#define cnsInitialConditions2D(gamma, mu, t, R, RBAR, UBAR, VBAR, PBAR, x, y, r, u, v, p) \
{                                         \
  *(r) = RBAR;           \
  *(u) = UBAR;           \
  *(v) = VBAR;           \
  *(p) = PBAR;           \
}

// Body force
// ************************************************************************
#define cnsBodyForce2D(gamma, mu, t, x, y, r, u, v, p, fx, fy) \
{                                                   \
  *(fx) = 0.0;                                      \
  *(fy) = 0.0;                                      \
}


// Inviscid Riemann Solve BC's
// ************************************************************************
#define cnsInviscidBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t,         \
                                        RBAR, UBAR, VBAR, PBAR,              \
                                        x, y, nx, ny,                        \
                                        rM, ruM, rvM, reM,                   \
                                        rB, ruB, rvB, reB)                   \
{                                                                            \
const dfloat uM    = ruM/rM;                                                 \
const dfloat vM    = rvM/rM;                                                 \
const dfloat pM    = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));             \
const dfloat unM   = uM*nx + vM*ny;                                          \
const dfloat cM    = sqrt(gamma*pM/rM);                                      \
const dfloat mM    = fabs(unM/cM);                                           \
const dfloat keREF = 0.5*RBAR*(UBAR*UBAR + VBAR*VBAR);                       \
  if(bc==11){                                                                \
    *( rB) = rM;                                                             \
    *(ruB) = -ruM;                                                           \
    *(rvB) = -rvM;                                                           \
    *(reB) = reM;                                                            \
  }else if(bc==12){                                                          \
    *( rB) = rM;                                                             \
    *(ruB) = -ruM;                                                           \
    *(rvB) = -rvM;                                                           \
    *(reB) = reM;                                                            \
   }else if(bc==13){                                                         \
    *( rB) =  rM;                                                            \
    *(ruB) =  rM*(uM - 2.0*(nx*uM+ny*vM)*nx);                                \
    *(rvB) =  rM*(vM - 2.0*(nx*uM+ny*vM)*ny);                                \
    *(reB) =  reM;                                                           \
   }else if(bc==20){                                                         \                                                              \
    *( rB) = RBAR;                                                           \
    *(ruB) = RBAR*UBAR;                                                      \
    *(rvB) = RBAR*VBAR;                                                      \
    *(reB) = PBAR/(gamma-1)+keREF ;                                          \
   }else if(bc==21){                                                         \
    *( rB) = RBAR;                                                           \
    *(ruB) = RBAR*UBAR;                                                      \
    *(rvB) = RBAR*VBAR;                                                      \
    *(reB) = pM/(gamma-1.0)+keREF ;                                          \
   }else if(bc==22){                                                         \
    *( rB) = RBAR;                                                           \
    *(ruB) = RBAR*UBAR;                                                      \
    *(rvB) = RBAR*VBAR;                                                      \
    *(reB) = PBAR/(gamma-1)+keREF ;                                          \
   }else if(bc==31){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = ruM;                                                            \
    *(rvB) = rvM;                                                            \
    *(reB) = (2.0*PBAR-pM)/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);             \
   }else if(bc==32){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = ruM;                                                            \
    *(rvB) = rvM;                                                            \
    *(reB) = reM;                                                            \
  }else if(bc==40){                                                          \
   if(unM>0.0){                                                              \
    const dfloat pB = mM > 1.0 ? pM :  PBAR;                                 \
    *( rB) = rM;                                                             \
    *(ruB) = ruM;                                                            \
    *(rvB) = rvM;                                                            \
    *(reB) = pB/(gamma-1)+0.5*rM*(uM*uM + vM*vM) ;                           \
   }else{                                                                    \
    *( rB) = RBAR;                                                           \
    *(ruB) = RBAR*UBAR;                                                      \
    *(rvB) = RBAR*VBAR;                                                      \
    *(reB) = PBAR/(gamma-1)+keREF ;                                          \
   }                                                                         \
 }                                                                           \
}


// Viscous Riemann Solve BC's
// ************************************************************************
#define cnsViscousBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t,          \
                                      RBAR, UBAR, VBAR, PBAR,                \
                                      x, y, nx, ny,                          \
                                      rM, ruM, rvM, reM, rB, ruB, rvB, reB,  \
                                      drrdxM, drrdyM, drudxM, drudyM,        \
                                      drvdxM, drvdyM, dredxM, dredyM,        \
                                      drrdxB, drrdyB, drudxB, drudyB,        \
                                      drvdxB, drvdyB, dredxB, dredyB)        \
{                                                                            \
const dfloat uM  =  ruM/rM;                                                   \
const dfloat vM  =  rvM/rM;                                                  \
const dfloat pM  = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));               \
const dfloat tM  = pM/(rM*R);                                                \
const dfloat unM = uM*nx + vM*ny;                                            \
const dfloat cM  = sqrt(gamma*pM/rM);                                        \
const dfloat mM  = fabs(unM/cM);                                             \
const dfloat keREF = 0.5*RBAR*(UBAR*UBAR + VBAR*VBAR);                       \
  if(bc==11){                                                                \
    const dfloat TBAR  = PBAR/(R*RBAR);                                      \
    *( rB) = rM;                                                             \
    *(ruB) = 0.0;                                                            \
    *(rvB) = 0.0;                                                            \
    *(reB) = (rM*R*TBAR)/(gamma-1.0);                                        \
    *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
    *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
    *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
    *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
   }else if(bc==12){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = 0.0;                                                            \
    *(rvB) = 0.0;                                                            \
    *(reB) = reM - 0.5*rM*(uM*uM+vM*vM);                                     \
    const dfloat uxM = drudxM - uM*drrdxM;                                   \
    const dfloat uyM = drudyM - uM*drrdyM;                                   \
    const dfloat vxM = drvdxM - vM*drrdxM;                                   \
    const dfloat vyM = drvdyM - vM*drrdyM;                                   \
    const dfloat TxM = dredxM - (drrdxM*reM/rM + uM*uxM + vM*vxM);           \
    const dfloat TyM = dredyM - (drrdyM*reM/rM + uM*uyM + vM*vyM);           \
    *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
    *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
    *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
    *(dredxB) = dredxM - nx*(nx*TxM + ny*TyM);                               \                                                     \
    *(dredyB) = dredyM - ny*(nx*TxM + ny*TyM);                               \
  }else if(bc==13){                                                          \
    *( rB) =  rM;                                                            \
    *(ruB) =  rM*(uM - 2.0*(nx*uM+ny*vM)*nx);                                \
    *(rvB) =  rM*(vM - 2.0*(nx*uM+ny*vM)*ny);                                \
    *(reB) =  reM;                                                           \
    *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                         \
    *(drudxB) = 0.0;*(drudyB) = 0.0;                                         \
    *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                         \
    *(dredxB) = 0.0;*(dredyB) = 0.0;                                         \
  }else if(bc==20||bc==40){                                                  \
    *( rB) = rM;                                                             \
    *(ruB) = rM*uM;                                                          \
    *(rvB) = rM*vM;                                                          \
    *(reB) = reM;                                                            \
    *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
    *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
    *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
    *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
  }else if(bc==21){                                                          \
    *( rB) = RBAR;                                                           \
    *(ruB) = RBAR*UBAR;                                                      \
    *(rvB) = RBAR*VBAR;                                                      \
    *(reB) = reM-0.5*rM*(uM*uM+vM*vM)+keREF ;                                \
    *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                         \
    *(drudxB) = 0.0;*(drudyB) = 0.0;                                         \
    *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                         \
    *(dredxB) = 0.0;*(dredyB) = 0.0;                                         \
  }else if(bc==22){                                                          \
    *( rB) = RBAR;                                                           \
    *(ruB) = RBAR*UBAR;                                                      \
    *(rvB) = RBAR*VBAR;                                                      \
    *(reB) = PBAR/(gamma -1.0) + keREF;                                      \
    *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                         \
    *(drudxB) = 0.0;*(drudyB) = 0.0;                                         \
    *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                         \
    *(dredxB) = 0.0;*(dredyB) = 0.0;                                         \
  }else if(bc==31){                                                          \
    *( rB) = rM;                                                             \
    *(ruB) = ruM;                                                            \
    *(rvB) = rvM;                                                            \
    *(reB) = PBAR/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                     \
    *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                         \
    *(drudxB) = 0.0;*(drudyB) = 0.0;                                         \
    *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                         \
    *(dredxB) = 0.0;*(dredyB) = 0.0;                                         \
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









// // Viscous Riemann Solve BC's
// // ************************************************************************
// #define cnsViscousBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t,          \
//                                       RBAR, UBAR, VBAR, PBAR, TBAR,           \
//                                       x, y, nx, ny,                          \
//                                       rM, ruM, rvM, reM, rB, ruB, rvB, reB,  \
//                                       drrdxM, drrdyM, drudxM, drudyM,        \
//                                       drvdxM, drvdyM, dredxM, dredyM,        \
//                                       drrdxB, drrdyB, drudxB, drudyB,        \
//                                       drvdxB, drvdyB, dredxB, dredyB)        \
// {                                                                            \
// const dfloat uM  = ruM/rM;                                                   \
// const dfloat vM  = rvM/rM;                                                   \
// const dfloat pM  = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));               \
// const dfloat unM = uM*nx + vM*ny;                                            \
// const dfloat cM  = sqrt(gamma*pM/rM);                                        \
// const dfloat mM  = fabs(unM/cM);                                             \
// const dfloat keREF = 0.5*RBAR*(UBAR*UBAR + VBAR*VBAR);                       \
// dfloat PR;                                                                   \
// PressureRiemann2D(gamma, R, CP, CV, unM, cM, rM, uM, vM, pM, &PR);           \
// const dfloat TM = (reM-0.5*rM*(uM*uM+vM*vM))/(rM*CV);                       \
//   if(bc==11){                                                                \
//     *( rB) = PR/(TBAR*R);                                                    \
//     *(ruB) = 0.0;                                                            \
//     *(rvB) = 0.0;                                                            \
//     *(reB) = PR/(gamma-1.0);                                                 \
//     *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
//     *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
//     *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
//     *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
//    }else if(bc==12){                                                         \
//     *( rB) = PR/(TM*R);                                                      \
//     *(ruB) = 0.0;                                                            \
//     *(rvB) = 0.0;                                                            \
//     *(reB) = PR/(gamma-1.0);                                                 \
//     const dfloat uxM = drudxM - uM*drrdxM;                                   \
//     const dfloat uyM = drudyM - uM*drrdyM;                                   \
//     const dfloat vxM = drvdxM - vM*drrdxM;                                   \
//     const dfloat vyM = drvdyM - vM*drrdyM;                                   \
//     const dfloat TxM = dredxM - (drrdxM*reM/rM + uM*uxM + vM*vxM);           \
//     const dfloat TyM = dredyM - (drrdyM*reM/rM + uM*uyM + vM*vyM);           \
//     *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
//     *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
//     *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
//     *(dredxB) = dredxM - nx*(nx*TxM + ny*TyM);                               \                                                     \
//     *(dredyB) = dredyM - ny*(nx*TxM + ny*TyM);                               \
//   }else if(bc==13){                                                          \
//     *( rB) =  rM;                                                            \
//     *(ruB) =  ruM - (nx*ruM+ny*rvM)*nx;                                      \
//     *(rvB) =  rvM - (nx*ruM+ny*rvM)*ny;                                      \
//     *(reB) =  reM;                                                           \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }else if(bc==20){                                                          \
//     *( rB) = RBAR;                                                           \
//     *(ruB) = RBAR*UBAR;                                                      \
//     *(rvB) = RBAR*VBAR;                                                      \
//     *(reB) = PBAR/(gamma-1) + keREF;                                           \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \                                                                  \
//   }else if(bc==21){                                                          \
//     *( rB) = RBAR;                                                           \
//     *(ruB) = RBAR*UBAR;                                                      \
//     *(rvB) = RBAR*VBAR;                                                      \
//     *(reB) = pM/(gamma-1) + keREF;                                           \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }else if(bc==22){                                                          \
//     *( rB) = RBAR;                                                           \
//     *(ruB) = RBAR*UBAR;                                                      \
//     *(rvB) = RBAR*VBAR;                                                      \
//     *(reB) = PBAR/(gamma -1.0) + keREF;                                      \
//     *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                         \
//     *(drudxB) = 0.0;*(drudyB) = 0.0;                                         \
//     *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                         \
//     *(dredxB) = 0.0;*(dredyB) = 0.0;                                         \
//   }else if(bc==31){                                                          \
//     *( rB) = rM;                                                             \
//     *(ruB) = ruM;                                                            \
//     *(rvB) = rvM;                                                            \
//     *(reB) = 4.4642857/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                     \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }else if(bc==32){                                                          \
//     *( rB) = rM;                                                             \
//     *(ruB) = ruM;                                                            \
//     *(rvB) = rvM;                                                            \
//     *(reB) = reM;                                                            \
//     *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
//     *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
//     *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
//     *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
//   } else if(bc==41){                                                         \
//     *( rB) = RBAR;                                                           \
//     *(ruB) = RBAR*UBAR;                                                      \
//     *(rvB) = RBAR*VBAR;                                                      \
//     *(reB) = PBAR/(gamma -1.0) + keREF;                                      \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }                                                                          \
// }

// // Viscous Riemann Solve BC's
// // ************************************************************************
// #define cnsViscousBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t,          \
//                                       RBAR, UBAR, VBAR, PBAR, TBAR,           \
//                                       x, y, nx, ny,                          \
//                                       rM, ruM, rvM, reM, rB, ruB, rvB, reB,  \
//                                       drrdxM, drrdyM, drudxM, drudyM,        \
//                                       drvdxM, drvdyM, dredxM, dredyM,        \
//                                       drrdxB, drrdyB, drudxB, drudyB,        \
//                                       drvdxB, drvdyB, dredxB, dredyB)        \
// {                                                                            \
// const dfloat uM  = ruM/rM;                                                   \
// const dfloat vM  = rvM/rM;                                                   \
// const dfloat pM  = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));               \
// const dfloat unM = uM*nx + vM*ny;                                            \
// const dfloat cM  = sqrt(gamma*pM/rM);                                        \
// const dfloat mM  = fabs(unM/cM);                                             \
// const dfloat keREF = 0.5*RBAR*(UBAR*UBAR + VBAR*VBAR);                       \
// dfloat PR;                                                                   \
// PressureRiemann2D(gamma, R, CP, CV, unM, cM, rM, uM, vM, pM, &PR);           \
// const dfloat TM = (reM-0.5*rM*(uM*uM+vM*vM))/(rM*CV);                       \
//   if(bc==11){                                                                \
//     *( rB) = PR/(TBAR*R);                                                    \
//     *(ruB) = 0.0;                                                            \
//     *(rvB) = 0.0;                                                            \
//     *(reB) = PR/(gamma-1.0);                                                 \
//     *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
//     *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
//     *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
//     *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
//    }else if(bc==12){                                                         \
//     *( rB) = PR/(TM*R);                                                      \
//     *(ruB) = 0.0;                                                            \
//     *(rvB) = 0.0;                                                            \
//     *(reB) = PR/(gamma-1.0);                                                 \
//     const dfloat uxM = drudxM - uM*drrdxM;                                   \
//     const dfloat uyM = drudyM - uM*drrdyM;                                   \
//     const dfloat vxM = drvdxM - vM*drrdxM;                                   \
//     const dfloat vyM = drvdyM - vM*drrdyM;                                   \
//     const dfloat TxM = dredxM - (drrdxM*reM/rM + uM*uxM + vM*vxM);           \
//     const dfloat TyM = dredyM - (drrdyM*reM/rM + uM*uyM + vM*vyM);           \
//     *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
//     *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
//     *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
//     *(dredxB) = dredxM - nx*(nx*TxM + ny*TyM);                               \                                                     \
//     *(dredyB) = dredyM - ny*(nx*TxM + ny*TyM);                               \
//   }else if(bc==13){                                                          \
//     *( rB) =  rM;                                                            \
//     *(ruB) =  ruM - (nx*ruM+ny*rvM)*nx;                                      \
//     *(rvB) =  rvM - (nx*ruM+ny*rvM)*ny;                                      \
//     *(reB) =  reM;                                                           \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }else if(bc==20){                                                          \
//     if(unM<=0){                                                              \
//       if(mM > 1.0){                                                          \
//         *( rB) = RBAR;                                                       \
//         *(ruB) = RBAR*UBAR;                                                  \
//         *(rvB) = RBAR*VBAR;                                                  \
//         *(reB) = PBAR/(gamma -1.0) + keREF;                                  \
//         *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                     \
//         *(drudxB) = 0.0;*(drudyB) = 0.0;                                     \
//         *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                     \
//         *(dredxB) = 0.0;*(dredyB) = 0.0;                                     \
//       }else{                                                                 \
//         *( rB) = RBAR;                                                       \
//         *(ruB) = RBAR*UBAR;                                                  \
//         *(rvB) = RBAR*VBAR;                                                  \
//         *(reB) = pM/(gamma-1) + keREF;                                       \
//         *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                    \
//         *(drudxB) = 0.0; *(drudyB) = 0.0;                                    \
//         *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                    \
//         *(dredxB) = 0.0; *(dredyB) = 0.0;                                    \
//       }                                                                      \
//     }else{                                                                   \
//       if(mM > 1.0){                                                          \
//         *( rB) = rM;                                                         \
//         *(ruB) = ruM;                                                        \
//         *(rvB) = rvM;                                                        \
//         *(reB) = reM;                                                        \
//         *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                               \
//         *(drudxB) = drudxM;*(drudyB) = drudyM;                               \
//         *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                               \
//         *(dredxB) = dredxM;*(dredyB) = dredyM;                               \
//         }else{                                                               \
//         *( rB) = rM;                                                         \
//         *(ruB) = ruM;                                                        \
//         *(rvB) = rvM;                                                        \
//         *(reB) = 4.4642857/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                 \
//         *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                    \
//         *(drudxB) = 0.0; *(drudyB) = 0.0;                                    \
//         *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                    \
//         *(dredxB) = 0.0; *(dredyB) = 0.0;                                    \
//       }                                                                      \
//     }                                                                        \
//   }else if(bc==21){                                                          \
//     *( rB) = RBAR;                                                           \
//     *(ruB) = RBAR*UBAR;                                                      \
//     *(rvB) = RBAR*VBAR;                                                      \
//     *(reB) = pM/(gamma-1) + keREF;                                           \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }else if(bc==22){                                                          \
//     *( rB) = RBAR;                                                           \
//     *(ruB) = RBAR*UBAR;                                                      \
//     *(rvB) = RBAR*VBAR;                                                      \
//     *(reB) = PBAR/(gamma -1.0) + keREF;                                      \
//     *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                         \
//     *(drudxB) = 0.0;*(drudyB) = 0.0;                                         \
//     *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                         \
//     *(dredxB) = 0.0;*(dredyB) = 0.0;                                         \
//   }else if(bc==31){                                                          \
//     *( rB) = rM;                                                             \
//     *(ruB) = ruM;                                                            \
//     *(rvB) = rvM;                                                            \
//     *(reB) = 4.4642857/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                     \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }else if(bc==32){                                                          \
//     *( rB) = rM;                                                             \
//     *(ruB) = ruM;                                                            \
//     *(rvB) = rvM;                                                            \
//     *(reB) = reM;                                                            \
//     *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
//     *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
//     *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
//     *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
//   } else if(bc==41){                                                         \
//     *( rB) = RBAR;                                                           \
//     *(ruB) = RBAR*UBAR;                                                      \
//     *(rvB) = RBAR*VBAR;                                                      \
//     *(reB) = PBAR/(gamma -1.0) + keREF;                                      \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }                                                                          \
// }


// // Viscous Riemann Solve BC's
// // ************************************************************************
// #define cnsViscousBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t,          \
//                                       RBAR, UBAR, VBAR, PBAR, TBAR,          \
//                                       x, y, nx, ny,                          \
//                                       rM, ruM, rvM, reM, rB, ruB, rvB, reB,  \
//                                       drrdxM, drrdyM, drudxM, drudyM,        \
//                                       drvdxM, drvdyM, dredxM, dredyM,        \
//                                       drrdxB, drrdyB, drudxB, drudyB,        \
//                                       drvdxB, drvdyB, dredxB, dredyB)        \
// {                                                                            \
// const dfloat uM  = ruM/rM;                                                   \
// const dfloat vM  = rvM/rM;                                                   \
// const dfloat pM  = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));               \
// const dfloat tM  = pM/(rM*R);                                                \
// const dfloat unM = uM*nx + vM*ny;                                            \
// const dfloat cM  = sqrt(gamma*pM/rM);                                        \
// const dfloat mM  = fabs(unM/cM);                                             \
// const dfloat keREF = 0.5*RBAR*(UBAR*UBAR + VBAR*VBAR);                       \
// dfloat PR; PressureRiemann2D(gamma, R, CP, CV, unM, cM, rM, uM, vM, pM, &PR);\
//   if(bc==11){                                                                \
//     *( rB) = PR/(TBAR*R);                                                    \
//     *(ruB) = 0.0;                                                            \
//     *(rvB) = 0.0;                                                            \
//     *(reB) = PR/(gamma-1.0);                                                 \
//     *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
//     *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
//     *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
//     *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
//    }else if(bc==12){                                                         \
//     *( rB) = PR/(tM*R);                                                      \
//     *(ruB) = 0.0;                                                            \
//     *(rvB) = 0.0;                                                            \
//     *(reB) = PR/(gamma-1.0);                                                 \
//     *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                         \
//     *(drudxB) = 0.0;*(drudyB) = 0.0;                                         \
//     *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                         \
//     *(dredxB) = 0.0;*(dredyB) = 0.0;                                         \
//   }else if(bc==13){                                                          \
//     const dfloat ruS = ruM - (nx*ruM+ny*rvM)*nx;                                 \
//     const dfloat rvS = rvM - (nx*ruM+ny*rvM)*ny;                                 \
//     *( rB) =  rM;                                                            \
//     *(ruB) =  ruM - (nx*ruM+ny*rvM)*nx;                                      \
//     *(rvB) =  rvM - (nx*ruM+ny*rvM)*ny;                                      \
//     *(reB) =  PR/(gamma -1.0)+0.5*(ruS*ruS + rvS*rvS)/rM;                    \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }else if(bc==20){                                                          \
//     *( rB) = RBAR;                                                           \
//     *(ruB) = RBAR*UBAR;                                                      \
//     *(rvB) = RBAR*VBAR;                                                      \
//     *(reB) = PBAR/(gamma -1.0) + keREF;                                      \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }else if(bc==23){                                                          \
//     const dfloat vRmag = sqrt(UBAR*UBAR+VBAR*VBAR);                          \                                                           \
//     const dfloat cR    = sqrt(gamma*PBAR/RBAR);                              \
//     const dfloat MaOut = vRmag/cR;                                           \
//     const dfloat vMmag = sqrt(uM*uM+vM*vM);                                  \                                                           \
//     const dfloat cM    = sqrt(gamma*pM/rM);                                  \
//     const dfloat MaM   = vMmag/cM;                                           \
//     const dfloat cB    = vMmag/MaOut;                                        \
//     dfloat pB  = 1.0;                                                        \
//     if(MaM<1.0){                                                                 \
//       dfloat pT = pM*pow((1.0+0.5*(gamma-1.0)*MaM*MaM),     gamma/(gamma-1.0));  \
//       pB        = pT*pow((1.0+0.5*(gamma-1.0)*MaOut*MaOut),-gamma/(gamma-1.0));  \
//     }else{                                                                   \
//       pB = pM + 0.5*rM*(uM*uM+vM*vM);                                      \
//     }                                                                        \
//     const dfloat rb = gamma*pB/(cB*cB);                                      \
//     *( rB) = rb;                                                             \
//     *(ruB) = rb*uM;                                                          \
//     *(rvB) = rb*vM;                                                          \
//     *(reB) = pB/(gamma-1) + 0.5*rb*(uM*uM+vM*vM);                            \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }else if(bc==22){                                                          \
//     *( rB) = RBAR;                                                           \
//     *(ruB) = RBAR*UBAR;                                                      \
//     *(rvB) = RBAR*VBAR;                                                      \
//     *(reB) = PBAR/(gamma -1.0) + keREF;                                      \
//     *(drrdxB) = 0.0;*(drrdyB) = 0.0;                                         \
//     *(drudxB) = 0.0;*(drudyB) = 0.0;                                         \
//     *(drvdxB) = 0.0;*(drvdyB) = 0.0;                                         \
//     *(dredxB) = 0.0;*(dredyB) = 0.0;                                         \
//   }else if(bc==31){                                                          \
//     *( rB) = rM;                                                             \
//     *(ruB) = ruM;                                                            \
//     *(rvB) = rvM;                                                            \
//     *(reB) = 4.4642857/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }else if(bc==32){                                                          \
//     *( rB) = rM;                                                             \
//     *(ruB) = ruM;                                                            \
//     *(rvB) = rvM;                                                            \
//     *(reB) = reM;                                                            \
//     *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;                                   \
//     *(drudxB) = drudxM;*(drudyB) = drudyM;                                   \
//     *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;                                   \
//     *(dredxB) = dredxM;*(dredyB) = dredyM;                                   \
//   } else if(bc==41){                                                         \
//     *( rB) = RBAR;                                                           \
//     *(ruB) = RBAR*UBAR;                                                      \
//     *(rvB) = RBAR*VBAR;                                                      \
//     *(reB) = PBAR/(gamma -1.0) + keREF;                                      \
//     *(drrdxB) = 0.0; *(drrdyB) = 0.0;                                        \
//     *(drudxB) = 0.0; *(drudyB) = 0.0;                                        \
//     *(drvdxB) = 0.0; *(drvdyB) = 0.0;                                        \
//     *(dredxB) = 0.0; *(dredyB) = 0.0;                                        \
//   }                                                                          \
// }

// 
// 
// 

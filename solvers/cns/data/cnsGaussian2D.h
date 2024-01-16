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
#define p_VBAR  0.0000000
#define p_PBAR  71.428571428571431
#define p_TBAR  1.0000000

// Pressure Riemann Solver on BC based on NASA Report
/* ************************************************************************ */
#define PressureRiemann2D(gamma, R, CP, CV, UN, CN, rM, uM, vM, pM, PR){          \
const dfloat pfunc= 2.0*gamma/(gamma -1.0);                                       \
const dfloat AR   = 2.0/(rM*(gamma+1.0));                                         \
const dfloat BR   = pM*(gamma -1.0)/(gamma+1.0);                                  \
const dfloat PR1  = pM*pow(max(0.0001, 1.0+0.5*(gamma -1.0)*UN/CN), pfunc);       \
const dfloat PR2  = pM+0.5*UN/AR*(UN+sqrt(UN*UN+4.0*AR*(pM+BR)));                 \
*(PR)   = (UN<=0) ? PR1 : PR2;                                                    \
}

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

#define cnsInviscidBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny, \
                                        rM, ruM, rvM, reM, rB, ruB, rvB, reB) \
{                                                                            \
const dfloat uM  = ruM/rM;                                                   \
const dfloat vM =  rvM/rM;                                                   \
const dfloat pM  = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));              \
const dfloat unM = uM*nx + vM*ny;                                            \
const dfloat cM  = sqrt(gamma*pM/rM);                                        \
const dfloat mM  = fabs(unM/cM);                                             \
  if(bc==11){                                                                \
    *( rB) = rM;                                                             \
    *(ruB) = -ruM;                                                           \
    *(rvB) = -rvM;                                                           \
    *(reB) = rM*CP*p_TBAR/gamma + 0.5*(ruM*ruM + rvM*rvM)/rM;                \
   }else if(bc==12){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = -ruM;                                                           \
    *(rvB) = -rvM;                                                           \
    *(reB) = rM*CP*p_TBAR/gamma + 0.5*(ruM*ruM + rvM*rvM)/rM;                \
   }else if(bc==20){                                                         \
    if(unM<=0){                                                              \
      if(mM > 1.0){                                                          \
        *( rB) = p_RBAR;                                                     \
        *(ruB) = p_RBAR*p_UBAR;                                              \
        *(rvB) = p_RBAR*p_VBAR;                                              \
        *(reB) = p_PBAR/(gamma -1.0) + 0.5*p_RBAR*(p_UBAR*p_UBAR + p_VBAR*p_VBAR);\
      }else{                                                                 \
        *( rB) = p_RBAR;                                                     \
        *(ruB) = p_RBAR*p_UBAR;                                              \
        *(rvB) = p_RBAR*p_VBAR;                                              \
        *(reB) = pM/(gamma -1.0) + 0.5*p_RBAR*(p_UBAR*p_UBAR + p_VBAR*p_VBAR);\
      }                                                                       \
    }else{                                                                    \
      if(mM > 1.0){                                                           \
      *( rB) = rM;                                                            \
      *(ruB) = uM;                                                            \
      *(rvB) = vM;                                                            \
      *(reB) = reM;                                                           \
      }else{                                                                  \
      *( rB) = rM;                                                            \
      *(ruB) = ruM;                                                           \
      *(rvB) = rvM;                                                           \
      *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                  \
      }                                                                       \
    }                                                                         \
  }                                                                           \
}


#define cnsViscousBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny, \
                                      rM, ruM, rvM, reM, rB, ruB, rvB, reB) \
{                                                                            \
const dfloat uM  = ruM/rM;                                                   \
const dfloat vM =  rvM/rM;                                                   \
const dfloat pM  = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));              \
const dfloat unM = uM*nx + vM*ny;                                            \
const dfloat cM  = sqrt(gamma*pM/rM);                                        \
const dfloat mM  = fabs(unM/cM);                                             \
  if(bc==11){                                                                \
    *( rB) = rM;                                                             \
    *(ruB) = 0.0;                                                            \
    *(rvB) = 0.0;                                                            \
    *(reB) = rM*CP*p_TBAR/gamma;                                             \
   }else if(bc==12){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = -ruM;                                                           \
    *(rvB) = -rvM;                                                           \
    *(reB) = rM*CP*p_TBAR/gamma + 0.5*(ruM*ruM + rvM*rvM)/rM;                \
  }else if(bc==20){                                                          \
    if(unM<=0){                                                              \
      if(mM > 1.0){                                                          \
        *( rB) = p_RBAR;                                                     \
        *(ruB) = p_RBAR*p_UBAR;                                              \
        *(rvB) = p_RBAR*p_VBAR;                                              \
        *(reB) = p_PBAR/(gamma -1.0) + 0.5*p_RBAR*(p_UBAR*p_UBAR + p_VBAR*p_VBAR);\
      }else{                                                                 \
        *( rB) = p_RBAR;                                                     \
        *(ruB) = p_RBAR*p_UBAR;                                              \
        *(rvB) = p_RBAR*p_VBAR;                                              \
        *(reB) = pM/(gamma -1.0) + 0.5*p_RBAR*(p_UBAR*p_UBAR + p_VBAR*p_VBAR);\
      }                                                                       \
    }else{                                                                    \
      if(mM > 1.0){                                                           \
      *( rB) = rM;                                                            \
      *(ruB) = uM;                                                            \
      *(rvB) = vM;                                                            \
      *(reB) = reM;                                                           \
      }else{                                                                  \
      *( rB) = rM;                                                            \
      *(ruB) = ruM;                                                           \
      *(rvB) = rvM;                                                           \
      *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                  \
      }                                                                       \
    }                                                                         \
  }                                                                           \
}


//   Walls Start with 1
// 11: No-Slip Isothermal Wall
// 12: No-Slip Adiabatic Wall 
// 13: Slip Wall or Euler Wall 
// 20: Far Field with subsonic Inlet-Outlet
    // *(reB) = rM*CP*p_TBAR/gamma;                                             

// #define cnsBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny, \
//                                 rM, ruM, rvM, reM, rB, ruB, rvB, reB) \
// {                                                                            \
// const dfloat uM  = ruM/rM;                                                   \
// const dfloat vM =  rvM/rM;                                                   \
// const dfloat pM  = (gamma-1.0)*(reM - 0.5*rM*(uM*uM + vM*vM));              \
// const dfloat unM = uM*nx + vM*ny;                                            \
// const dfloat cM  = sqrt(gamma*pM/rM);                                        \
// const dfloat mM  = fabs(unM/cM);                                             \
//   if(bc==11){                                                                \
//     *( rB) = rM;                                                             \
//     *(ruB) = 0.0;                                                            \
//     *(rvB) = 0.0;                                                            \
//     *(reB) = reM;                                      \
//    }else if(bc==18){                                                         \
//     *( rB) = rM;                                                             \
//     *(ruB) = rM*p_WALL_VELOCITY;                                             \
//     *(rvB) = rM*p_WALL_VELOCITY;                                             \
//     *(reB) = rM*(CP*p_TBAR/gamma + 0.5*(p_WALL_VELOCITY*p_WALL_VELOCITY));   \
//   } else if(bc==20){                                                         \
//     if(unM<=0){                                                              \
//       if(mM > 1.0){                                                          \
//         *( rB) = p_RBAR;                                                     \
//         *(ruB) = p_RBAR*p_UBAR;                                              \
//         *(rvB) = p_RBAR*p_VBAR;                                              \
//         *(reB) = p_PBAR/(gamma -1.0) + 0.5*p_RBAR*(p_UBAR*p_UBAR + p_VBAR*p_VBAR);\
//       }else{                                                                 \
//         *( rB) = p_RBAR;                                                     \
//         *(ruB) = p_RBAR*p_UBAR;                                              \
//         *(rvB) = p_RBAR*p_VBAR;                                              \
//         *(reB) = pM/(gamma -1.0) + 0.5*p_RBAR*(p_UBAR*p_UBAR + p_VBAR*p_VBAR);\
//       }                                                                       \
//     }else{                                                                    \
//       if(mM > 1.0){                                                           \
//       *( rB) = rM;                                                            \
//       *(ruB) = uM;                                                            \
//       *(rvB) = vM;                                                            \
//       *(reB) = reM;                                                           \
//       }else{                                                                  \
//       *( rB) = rM;                                                            \
//       *(ruB) = ruM;                                                           \
//       *(rvB) = rvM;                                                           \
//       *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM);                  \
//       }                                                                       \
//     }                                                                         \
//   }                                                                           \
// }




// // Boundary conditions
// // 11: Adiabatic Wall: Diffusion-> density->inside, U=0, rhoE = Inside
// /* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
// // *(rB) = PR/(p_TBAR*R);                                                   

// #define cnsBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny, \
//                                   rM, uM, vM, pM, uxM, uyM, vxM, vyM, \
//                                   rB, uB, vB, pB, uxB, uyB, vxB, vyB) \
// {                                    \
// const dfloat uin  = uM*nx + vM*ny;                                           \
// const dfloat cin  = sqrt(gamma*pM/rM);                                       \
// const dfloat min  = fabs(uin/cin);                                           \
//   if(bc==11){                                                                \
//     *(rB) =  rM;                                                             \
//     *(uB) = -uM;                                                             \
//     *(vB) = -vM;                                                             \
//     *(pB) = PR;                                                              \
//     *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
//   } else if(bc==12){                                                         \
//     dfloat PR = 0;                                                           \
//     PressureRiemann2D(gamma, R, CP, CV, uin, cin, rM, uM, vM, pM, &PR);      \
//     const dfloat TB = pM/(rM*R);                                             \
//     *(rB) = PR/(TB*R);                                                       \
//     *(uB) = 0.0;                                                             \
//     *(vB) = 0.0;                                                             \
//     *(pB) = PR;                                                              \
//     *(uxB) = 0.0; *(uyB) = 0.0;*(vxB) = 0.0; *(vyB) = 0.0;                   \
//   }else if(bc==13){                                                          \
//     dfloat PR = 0;                                                           \
//     PressureRiemann2D(gamma, R, CP, CV, uin, cin, rM, uM, vM, pM, &PR);      \
//     *(rB) = rM;                                                              \
//     *(uB) = uM - (nx*uM+ny*vM)*nx;                                           \
//     *(vB) = vM - (nx*uM+ny*vM)*ny;                                           \
//     *(pB) = PR;                                                              \
//     *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
//   } else if(bc==20){                                                         \
//     if(uin<=0){                                                              \
//       if(min <=1.0){                                                         \
//         *(rB) = p_RBAR; *(uB) = p_UBAR;*(vB) = p_VBAR;                       \
//         *(pB) = pM;                                                          \
//       }else{                                                                 \
//         *(rB) = p_RBAR; *(uB) = p_UBAR; *(vB) = p_VBAR;*(pB) = p_PBAR;       \
//       }                                                                      \
//      *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                  \
//     }else{                                                                   \
//       if(min <=1.0){                                                         \
//       *(rB) = rM; *(uB) = uM; *(vB) = vM; *(pB) = p_PBAR;                    \
//       }else{                                                                 \
//        *(rB) = rM; *(uB) = uM; *(vB) = vM; *(pB) = pM;                       \
//       }                                                                      \
//     *(uxB) = 0.0; *(uyB) = 0.0;*(vxB) = 0.0; *(vyB) = 0.0;                   \
//     }                                                                        \
//   }else if(bc==41){                                                          \
//     *(rB) = rM;                                                              \
//     *(uB) = uM - (nx*uM+ny*vM)*nx;                                           \
//     *(vB) = vM - (nx*uM+ny*vM)*ny;                                           \
//     *(pB) = pM;                                                              \
//     *(uyB) = uyM;                                                            \
//     *(vxB) = vxM;                                                            \
//     *(vyB) = vyM;                                                            \
//     *(uxB) = uxM;                                                            \
//   }  else if(bc==42){                                                        \
//     *(rB) = rM;                                                              \
//     *(uB) = uM;                                                              \
//     *(vB) = vM;                                                              \
//     *(pB) = pM;                                                              \
//     *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
//   }                                                                          \
// }






// // Boundary conditions
// // 11: Adiabatic Wall: Diffusion-> density->inside, U=0, rhoE = Inside
// /* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
// // *(rB) = PR/(p_TBAR*R);                                                   

// #define cnsBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny, \
//                                   rM, uM, vM, pM, uxM, uyM, vxM, vyM, \
//                                   rB, uB, vB, pB, uxB, uyB, vxB, vyB) \
// {                                    \
// const dfloat uin  = uM*nx + vM*ny;                                           \
// const dfloat cin  = sqrt(gamma*pM/rM);                                       \
// const dfloat min  = fabs(uin/cin);                                           \
//   if(bc==11){                                                                \
//     dfloat PR = 0;                                                           \
//     PressureRiemann2D(gamma, R, CP, CV, uin, cin, rM, uM, vM, pM, &PR);      \
//     *(rB) = PR/(p_TBAR*R);                                                   \
//     *(uB) = 0.0;                                                             \
//     *(vB) = 0.0;                                                             \
//     *(pB) = PR;                                                              \
//     *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
//   } else if(bc==12){                                                         \
//     dfloat PR = 0;                                                           \
//     PressureRiemann2D(gamma, R, CP, CV, uin, cin, rM, uM, vM, pM, &PR);      \
//     const dfloat TB = pM/(rM*R);                                             \
//     *(rB) = PR/(TB*R);                                                       \
//     *(uB) = 0.0;                                                             \
//     *(vB) = 0.0;                                                             \
//     *(pB) = PR;                                                              \
//     *(uxB) = 0.0; *(uyB) = 0.0;*(vxB) = 0.0; *(vyB) = 0.0;                   \
//   }else if(bc==13){                                                          \
//     dfloat PR = 0;                                                           \
//     PressureRiemann2D(gamma, R, CP, CV, uin, cin, rM, uM, vM, pM, &PR);      \
//     *(rB) = rM;                                                              \
//     *(uB) = uM - (nx*uM+ny*vM)*nx;                                           \
//     *(vB) = vM - (nx*uM+ny*vM)*ny;                                           \
//     *(pB) = PR;                                                              \
//     *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
//   } else if(bc==20){                                                         \
//     if(uin<=0){                                                              \
//       if(min <=1.0){                                                         \
//         *(rB) = p_RBAR; *(uB) = p_UBAR;*(vB) = p_VBAR;                       \
//         *(pB) = pM;                                                          \
//       }else{                                                                 \
//         *(rB) = p_RBAR; *(uB) = p_UBAR; *(vB) = p_VBAR;*(pB) = p_PBAR;       \
//       }                                                                      \
//      *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                  \
//     }else{                                                                   \
//       if(min <=1.0){                                                         \
//       *(rB) = rM; *(uB) = uM; *(vB) = vM; *(pB) = p_PBAR;                    \
//       }else{                                                                 \
//        *(rB) = rM; *(uB) = uM; *(vB) = vM; *(pB) = pM;                       \
//       }                                                                      \
//     *(uxB) = 0.0; *(uyB) = 0.0;*(vxB) = 0.0; *(vyB) = 0.0;                   \
//     }                                                                        \
//   }else if(bc==41){                                                          \
//     *(rB) = rM;                                                              \
//     *(uB) = uM - (nx*uM+ny*vM)*nx;                                           \
//     *(vB) = vM - (nx*uM+ny*vM)*ny;                                           \
//     *(pB) = pM;                                                              \
//     *(uyB) = uyM;                                                            \
//     *(vxB) = vxM;                                                            \
//     *(vyB) = vyM;                                                            \
//     *(uxB) = uxM;                                                            \
//   }  else if(bc==42){                                                        \
//     *(rB) = rM;                                                              \
//     *(uB) = uM;                                                              \
//     *(vB) = vM;                                                              \
//     *(pB) = pM;                                                              \
//     *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
//   }                                                                          \
// }


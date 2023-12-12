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

//Mean flow
// #define p_RBAR 5.6
// #define p_UBAR 1.0
// #define p_VBAR 0.0
// #define p_PBAR 6.25

// MAch 2
// #define p_RBAR 11.2
// #define p_UBAR 0.984807753012208
// #define p_VBAR 0.173648177666930
// #define p_PBAR 2.0

// Mach 3
#define p_RBAR 44.8
#define p_UBAR 0.984807753012208
#define p_VBAR 0.173648177666930
#define p_PBAR 3.555555555555555


// // Ma = 0.8
// #define p_RBAR 1.4
// #define p_UBAR 0.984807753012208
// #define p_VBAR 0.173648177666930
// #define p_PBAR 1.5625
// Define Temperature for Isothermal Wall T = P /(rho * R)
#define p_TBAR 1.00000

/* ************************************************************************ */
/* ************************************************************************ */


// Pressure Riemann Solver on BC based on NASA Report
/* ************************************************************************ */
#define PressureRiemann2D(gamma, R, CP, CV, UN, CN, rM, uM, vM, pM, PR){          \
const dfloat pfunc= 2.0*gamma/(gamma -1.0);                                     \
const dfloat AR   = 2.0/(rM*(gamma+1.0));                                       \
const dfloat BR   = pM*(gamma -1.0)/(gamma+1.0);                                \
const dfloat PR1  = pM*pow(max(0.0001, 1.0+0.5*(gamma -1.0)*UN/CN), pfunc);     \
const dfloat PR2  = pM+0.5*UN/AR*(UN+sqrt(UN*UN+4.0*AR*(pM+BR)));                \
*(PR)   = (UN<=0) ? PR1 : PR2;                                                   \
}

/* ************************************************************************ */
/* ************************************************************************ */



// ************************************************************************
// Initial conditions (p is ignored for isothermal)
#define cnsInitialConditions2D(gamma, mu, t, x, y, r, u, v, p) \
{                                         \
  *(r) = p_RBAR;           \
  *(u) = p_UBAR;           \
  *(v) = p_VBAR;           \
  *(p) = p_PBAR;           \
}

// Body force
#define cnsBodyForce2D(gamma, mu, t, x, y, r, u, v, p, fx, fy) \
{                                                   \
  *(fx) = 0.0;                                      \
  *(fy) = 0.0;                                      \
}


/*
// Isothermall Wall 1-1, all gradients from interior, temperature from reference state
// Adibatic Wall 1-2, all gradients are zero, temperature from internal state
// Euler Slip Wall 1-3, normal velocity zero, gradient from inside !!!!, T = P/(rho R)
// FarFiled (Dirichlet) 20, free-stream based on characteristics, 
//   subsonic   - inlet -> pressure from inside outlet->pressure is freestream (or density)
//   supersonic - inlet -> all freestream       outlet-> all from inside
// Subsonic inflow 21, if supersonic use total pressure to compute density
// Supersonic inflow 22, ceck horses2D
// Pressure Inlet 23, ceck horses2D
// Subsonic outflow 31, 
// Supersonic outflow 32,
// Symmetry 41  normal velocity zero, normal of grad of other fields are zero, i.e. pB = pM and rB = rM 
// Periodic 42  copy state qB = qM and nabla(qB) = nabla(qM) 
// Pressure Outflow 33,
*/
#define cnsBoundaryConditions2D(bc, gamma, R, CP, CV, mu, t, x, y, nx, ny, \
                                  rM, uM, vM, pM, uxM, uyM, vxM, vyM, \
                                  rB, uB, vB, pB, uxB, uyB, vxB, vyB) \
{                                                                            \
const dfloat uin  = uM*nx + vM*ny;                                           \
const dfloat cin  = sqrt(gamma*pM/rM);                                       \
const dfloat min  = fabs(uin/cin);                                           \
  if(bc==11){                                                                \
    dfloat PR = 0;                                                           \
    PressureRiemann2D(gamma, R, CP, CV, uin, cin, rM, uM, vM, pM, &PR);      \
    *(rB) = PR/(p_TBAR*R);                                                   \
    *(uB) = 0.0;                                                             \
    *(vB) = 0.0;                                                             \
    *(pB) = PR;                                                              \
    *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
  } else if(bc==12){                                                         \
    dfloat PR = 0;                                                           \
    PressureRiemann2D(gamma, R, CP, CV, uin, cin, rM, uM, vM, pM, &PR);      \
    const dfloat TB = pM/(rM*R);                                             \
    *(rB) = PR/(TB*R);                                                       \
    *(uB) = 0.0;                                                             \
    *(vB) = 0.0;                                                             \
    *(pB) = PR;                                                              \
    *(uxB) = 0.0; *(uyB) = 0.0;*(vxB) = 0.0; *(vyB) = 0.0;                   \
  }else if(bc==13){                                                          \
    dfloat PR = 0;                                                           \
    PressureRiemann2D(gamma, R, CP, CV, uin, cin, rM, uM, vM, pM, &PR);      \
    *(rB) = rM;                                                              \
    *(uB) = uM - (nx*uM+ny*vM)*nx;                                           \
    *(vB) = vM - (nx*uM+ny*vM)*ny;                                           \
    *(pB) = PR;                                                              \
    *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
  } else if(bc==20){                                                         \
    if(uin<=0){                                                              \
      if(min <=1.0){                                                         \
        *(rB) = p_RBAR; *(uB) = p_UBAR;*(vB) = p_VBAR;                       \
        *(pB) = pM;                                                          \
      }else{                                                                 \
        *(rB) = p_RBAR; *(uB) = p_UBAR; *(vB) = p_VBAR;*(pB) = p_PBAR;       \
      }                                                                      \
     *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                  \
    }else{                                                                   \
      if(min <=1.0){                                                         \
      *(rB) = rM; *(uB) = uM; *(vB) = vM; *(pB) = p_PBAR;                    \
      }else{                                                                 \
       *(rB) = rM; *(uB) = uM; *(vB) = vM; *(pB) = pM;                       \
      }                                                                      \
    *(uxB) = 0.0; *(uyB) = 0.0;*(vxB) = 0.0; *(vyB) = 0.0;                   \
    }                                                                        \
  }else if(bc==41){                                                          \
    *(rB) = rM;                                                              \
    *(uB) = uM - (nx*uM+ny*vM)*nx;                                           \
    *(vB) = vM - (nx*uM+ny*vM)*ny;                                           \
    *(pB) = pM;                                                              \
    *(uyB) = uyM;                                                            \
    *(vxB) = vxM;                                                            \
    *(vyB) = vyM;                                                            \
    *(uxB) = uxM;                                                            \
  }  else if(bc==42){                                                        \
    *(rB) = rM;                                                              \
    *(uB) = uM;                                                              \
    *(vB) = vM;                                                              \
    *(pB) = pM;                                                              \
    *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
  }                                                                          \
}


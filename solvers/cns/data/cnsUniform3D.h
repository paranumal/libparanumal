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
// ffmpeg -r 20 -f image2 -s 1538x784 -i mach%*.png -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -vframes 400 -crf 15 -pix_fmt yuv420p test.mp4
/**************************************************************************/
//Mean flow n*1.4/(1.4*Ma^2) = P
// #define p_RBAR 5.6
// #define p_UBAR 1.0
// #define p_VBAR 0.0
// #define p_WBAR 0.0
// #define p_PBAR 6.25

// // Mach = 0.6
// #define p_RBAR 5.6000000
// #define p_UBAR 1.0
// #define p_VBAR 0.0
// #define p_WBAR 0.0
// #define p_PBAR 11.111111

// // Mach = 2.0
// #define p_RBAR 11.2
// #define p_UBAR 1.0
// #define p_VBAR 0.0
// #define p_WBAR 0.0
// #define p_PBAR 2.0

// // Mach 3
// #define p_RBAR 44.8
// #define p_UBAR 1.0
// #define p_VBAR 0.0
// #define p_WBAR 0.0
// #define p_PBAR 3.555555555555555

// Ma = 0.8
#define p_RBAR 1.4
#define p_UBAR 1.0
#define p_VBAR 0.0
#define p_WBAR 0.0
#define p_PBAR 1.5625


#define p_TBAR 1.00000
/**************************************************************************/


// Pressure Riemann Solver on BC based on NASA Report
/**************************************************************************/
#define PressureRiemann3D(gamma, R, CP, CV, UN, CN, rM, uM, vM, wM, pM, PR){    \
const dfloat pfunc= 2.0*gamma/(gamma -1.0);                                     \
const dfloat AR   = 2.0/(rM*(gamma+1.0));                                       \
const dfloat BR   = pM*(gamma -1.0)/(gamma+1.0);                                \
const dfloat PR1  = pM*pow(max(0.0001, 1.0+0.5*(gamma -1.0)*UN/CN), pfunc);     \
const dfloat PR2  = pM+0.5*UN/AR*(UN+sqrt(UN*UN+4.0*AR*(pM+BR)));               \
*(PR)   = (UN<=0) ? PR1 : PR2;                                                  \
}


// Initial conditions (p is ignored for isothermal)
/*************************************************************************/
#define cnsInitialConditions3D(gamma, mu, t, x, y, z, r, u, v, w, p) \
{                                         \
  *(r) = p_RBAR;           \
  *(u) = p_UBAR;           \
  *(v) = p_VBAR;           \
  *(w) = p_WBAR;           \
  *(p) = p_PBAR;           \
}

// Body force
/**************************************************************************/
#define cnsBodyForce3D(gamma, mu, t, x, y, z, r, u, v, w, p, fx, fy, fz) \
{                                                   \
  *(fx) = 0.0;                                      \
  *(fy) = 0.0;                                      \
  *(fz) = 0.0;                                      \
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
#define cnsBoundaryConditions3D(bc, gamma, R, CP, CV, mu, t, x, y, z, nx, ny, nz, \
                                  rM, uM, vM, wM, pM, uxM, uyM, uzM, vxM, vyM, vzM, wxM, wyM, wzM, \
                                  rB, uB, vB, wB, pB, uxB, uyB, uzB, vxB, vyB, vzB, wxB, wyB, wzB) \
{                                                                            \
const dfloat uin  = uM*nx + vM*ny + wM*nz;                                   \
const dfloat cin  = sqrt(gamma*pM/rM);                                       \
const dfloat min  = fabs(uin/cin);                                           \
  if(bc==11){                                                                \
    dfloat PR = 0;                                                           \
    PressureRiemann3D(gamma, R, CP, CV, uin, cin, rM, uM, vM, wM, pM, &PR);  \
    const dfloat TB = pM/(rM*R);                                             \
    *(rB) = PR/(TB*R);                                                       \
    *(uB) = 0.0;                                                             \
    *(vB) = 0.0;                                                             \
    *(wB) = 0.0;                                                             \
    *(pB) = PR;                                                              \
    *(uxB) = uxM;*(uyB) = uyM; *(uzB) = uzM;                                 \
    *(vxB) = vxM;*(vyB) = vyM; *(vzB) = vzM;                                 \
    *(wxB) = wxM;*(wyB) = wyM; *(wzB) = wzM;                                 \
  } else if(bc==12){                                                         \
    dfloat PR = 0;                                                           \
    PressureRiemann3D(gamma, R, CP, CV, uin, cin, rM, uM, vM, wM, pM, &PR);  \
    const dfloat TB = pM/(rM*R);                                             \
    *(rB) = PR/(TB*R);                                                       \
    *(uB) = 0.0;                                                             \
    *(vB) = 0.0;                                                             \
    *(wB) = 0.0;                                                             \
    *(pB) = PR;                                                              \
    *(uxB) = 0.f;*(uyB) = 0.f; *(uzB) = 0.f;                                 \
    *(vxB) = 0.f;*(vyB) = 0.f; *(vzB) = 0.f;                                 \
    *(wxB) = 0.f;*(wyB) = 0.f; *(wzB) = 0.f;                                 \
  }else if(bc==13){                                                          \
    dfloat PR = 0;                                                           \
    PressureRiemann3D(gamma, R, CP, CV, uin, cin, rM, uM, vM, wM, pM, &PR);  \
    *(rB) = rM;                                                              \
    *(uB) = uM - (nx*uM+ny*vM+nz*wM)*nx;                                     \
    *(vB) = vM - (nx*uM+ny*vM+nz*wM)*ny;                                     \
    *(wB) = wM - (nx*uM+ny*vM+nz*wM)*nz;                                     \
    *(pB) = PR;                                                              \
    *(uxB) = uxM;*(uyB) = uyM; *(uzB) = uzM;                                 \
    *(vxB) = vxM;*(vyB) = vyM; *(vzB) = vzM;                                 \
    *(wxB) = wxM;*(wyB) = wyM; *(wzB) = wzM;                                 \
  } else if(bc==20){                                                         \
    if(uin<=0){                                                              \
      if(min <=1.0){                                                         \
        *(rB) = p_RBAR; *(uB) = p_UBAR;*(vB) = p_VBAR; *(wB) = p_WBAR;       \
        *(pB) = pM;                                                          \
      }else{                                                                 \
        *(rB) = p_RBAR; *(uB) = p_UBAR; *(vB) = p_VBAR;*(wB) = p_WBAR;        \
        *(pB) = p_PBAR;                                                      \
      }                                                                      \
    *(uxB) = uxM;*(uyB) = uyM; *(uzB) = uzM;                                 \
    *(vxB) = vxM;*(vyB) = vyM; *(vzB) = vzM;                                 \
    *(wxB) = wxM;*(wyB) = wyM; *(wzB) = wzM;                                 \
    }else{                                                                   \
      if(min <=1.0){                                                         \
      *(rB) = rM; *(uB) = uM; *(vB) = vM; *(wB) = wM;                        \
      *(pB) = p_PBAR;                                                        \
      }else{                                                                 \
       *(rB) = rM; *(uB) = uM; *(vB) = vM; *(wB) = wM;*(pB) = pM;            \
      }                                                                      \
     *(uxB) = 0.0;*(uyB) = 0.0; *(uzB) = 0.0;                                \
     *(vxB) = 0.0;*(vyB) = 0.0; *(vzB) = 0.0;                                \
     *(wxB) = 0.0;*(wyB) = 0.0; *(wzB) = 0.0;                                \
    }                                                                        \
  }else if(bc==41){                                                          \
    *(rB) = rM;                                                              \
    *(uB) = uM - (nx*uM+ny*vM+nz*wM)*nx;                                     \
    *(vB) = vM - (nx*uM+ny*vM+nz*wM)*ny;                                     \
    *(wB) = wM - (nx*uM+ny*vM+nz*wM)*nz;                                     \
    *(pB) = pM;                                                              \
    *(uxB) = uxM;*(uyB) = uyM; *(uzB) = uzM;                                 \
    *(vxB) = vxM;*(vyB) = vyM; *(vzB) = vzM;                                 \
    *(wxB) = wxM;*(wyB) = wyM; *(wzB) = wzM;                                 \
  }  else if(bc==42){                                                        \
    *(rB) = rM;                                                              \
    *(uB) = uM;                                                              \
    *(vB) = vM;                                                              \
    *(wB) = wM;                                                              \
    *(pB) = pM;                                                              \
    *(uxB) = uxM;*(uyB) = uyM; *(uzB) = uzM;                                 \
    *(vxB) = vxM;*(vyB) = vyM; *(vzB) = vzM;                                 \
    *(wxB) = wxM;*(wyB) = wyM; *(wzB) = wzM;                                 \
  }                                                                          \
}



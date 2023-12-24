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

#define p_TBAR 1.00000

// Initial conditions (p is ignored for isothermal)
#define cnsInitialConditions3D(gamma, mu, t, x, y, z, r, u, v, w, p) \
{                                         \
  *(r) = 1 + exp(-3*(x*x+y*y+z*z));       \
  *(u) = exp(-3*(x*x+y*y+z*z));           \
  *(v) = exp(-3*(x*x+y*y+z*z));           \
  *(w) = exp(-3*(x*x+y*y+z*z));           \
  *(p) = 71.428571428571431*(1 + exp(-3*(x*x+y*y+z*z)));       \
}

// Body force
#define cnsBodyForce3D(gamma, mu, t, x, y, z, r, u, v, w, p, fx, fy, fz) \
{                                                   \
  *(fx) = 0.0;                                      \
  *(fy) = 0.0;                                      \
  *(fz) = 0.0;                                      \
}


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

/* ************************************************************************ */
#define cnsDiffusionBoundaryConditions3D(bc, gamma, R, CP, CV, mu, t, x, y, z, nx, ny, nz, \
                                  rM, uM, vM, wM, pM, uxM, uyM, uzM, vxM, vyM, vzM, wxM, wyM, wzM, \
                                  rB, uB, vB, wB, pB, uxB, uyB, uzB, vxB, vyB, vzB, wxB, wyB, wzB) \
{                                      \
const dfloat uin  = uM*nx + vM*ny;                                           \
const dfloat cin  = sqrt(gamma*pM/rM);                                       \
const dfloat min  = fabs(uin/cin);                                           \
  if(bc==11){                                                                \
    *(rB) = rM;                                                              \
    *(uB) = 0.0;                                                             \
    *(vB) = 0.0;                                                             \
    *(wB) = 0.0;                                                             \
    *(pB) = (gamma-1.0)*rM*CV*p_TBAR;                                        \
    *(uxB) = uxM;*(uyB) = uyM; *(uzB) = uzM;                                 \
    *(vxB) = vxM;*(vyB) = vyM; *(vzB) = vzM;                                 \
    *(wxB) = wxM;*(wyB) = wyM; *(wzB) = wzM;                                 \
  } else if(bc==12){                                                         \
    *(rB) = rM;                                                              \
    *(uB) = 0.0;                                                             \
    *(vB) = 0.0;                                                             \
    *(wB) = 0.0;                                                             \
    *(pB) = pM + 0.5*(gamma -1.0)*rM*(uM*uM + vM*vM + wM*wM);                \                                                              \
    *(uxB) = 0.0; *(uyB) = 0.0;*(vxB) = 0.0; *(vyB) = 0.0;                   \
  }else if(bc==13){                                                          \
    *(rB) = rM;                                                              \
    *(uB) = uM - (nx*uM+ny*vM)*nx;                                           \
    *(vB) = vM - (nx*uM+ny*vM)*ny;                                           \
    *(pB) = pM;                                                              \                                                              \
    *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
  }else{                                                                     \
    *(rB) = rM;                                                              \
    *(uB) = uM;                                                              \
    *(vB) = vM;                                                              \
    *(pB) = pM;                                                               \
    *(uxB) = uxM;*(uyB) = uyM; *(vxB) = vxM; *(vyB) = vyM;                   \
  }\
}

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
    // dfloat PR = 0;                                                           
    // PressureRiemann3D(gamma, R, CP, CV, uin, cin, rM, uM, vM, wM, pM, &PR);  
#define cnsBoundaryConditions3D(bc, gamma, R, CP, CV, mu, t, x, y, z, nx, ny, nz, \
                                  rM, uM, vM, wM, pM, uxM, uyM, uzM, vxM, vyM, vzM, wxM, wyM, wzM, \
                                  rB, uB, vB, wB, pB, uxB, uyB, uzB, vxB, vyB, vzB, wxB, wyB, wzB) \
{                                      \
const dfloat uin  = uM*nx + vM*ny + wM*nz;                                   \
const dfloat cin  = sqrt(gamma*pM/rM);                                       \
const dfloat min  = fabs(uin/cin);                                           \
  if(bc==11){                                                                \
    dfloat PR = 0;                                                           \
    PressureRiemann3D(gamma, R, CP, CV, uin, cin, rM, uM, vM, wM, pM, &PR);  \
    *(rB) = PR/(p_TBAR*R);                                                   \
    *(uB) = 0.0;                                                             \
    *(vB) = 0.0;                                                             \
    *(wB) = 0.0;                                                             \
    *(pB) = PR;                                                              \
    *(uxB) = uxM;*(uyB) = uyM; *(uzB) = uzM;                                 \
    *(vxB) = vxM;*(vyB) = vyM; *(vzB) = vzM;                                 \
    *(wxB) = wxM;*(wyB) = wyM; *(wzB) = wzM;                                 \
  } else if(bc==2){                    \
    dfloat PR = 0;                                                           \
    PressureRiemann3D(gamma, R, CP, CV, uin, cin, rM, uM, vM, wM, pM, &PR);  \
    const dfloat TB = pM/(rM*R);                                             \
    *(rB) = PR/(TB*R);                                                   \
    *(uB) = 0.0;                                                             \
    *(vB) = 0.0;                                                             \
    *(wB) = 0.0;                                                             \
    *(pB) = PR;                                                              \
    *(uxB) = 0,0;                      \
    *(uyB) = 0,0;                      \
    *(uzB) = 0,0;                      \
    *(vxB) = 0,0;                      \
    *(vyB) = 0,0;                      \
    *(vzB) = 0,0;                      \
    *(wxB) = 0,0;                      \
    *(wyB) = 0,0;                      \
    *(wzB) = 0,0;                      \
  } else if(bc==3){                    \
    *(rB) = rM;                        \
    *(uB) = uM;                        \
    *(vB) = vM;                        \
    *(wB) = wM;                        \
    *(pB) = pM;                        \
    *(uxB) = 0.0;                      \
    *(uyB) = 0.0;                      \
    *(uzB) = 0.0;                      \
    *(vxB) = 0.0;                      \
    *(vyB) = 0.0;                      \
    *(vzB) = 0.0;                      \
    *(wxB) = 0.0;                      \
    *(wyB) = 0.0;                      \
    *(wzB) = 0.0;                      \
  } else if(bc==4||bc==5){             \
    *(rB) = rM;                        \
    *(uB) = uM - (nx*uM+ny*vM+nz*wM)*nx; \
    *(vB) = vM - (nx*uM+ny*vM+nz*wM)*ny; \
    *(wB) = wM - (nx*uM+ny*vM+nz*wM)*nz; \
    *(pB) = pM;                        \
    *(uxB) = uxM;                      \
    *(uyB) = uyM;                      \
    *(uzB) = uzM;                      \
    *(vxB) = vxM;                      \
    *(vyB) = vyM;                      \
    *(vzB) = vzM;                      \
    *(wxB) = wxM;                      \
    *(wyB) = wyM;                      \
    *(wzB) = wzM;                      \
  }                                    \
}

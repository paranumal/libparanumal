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

#define p_RBAR 1.4
#define p_UBAR 0.984807753012208
#define p_VBAR 0.173648177666930
#define p_WBAR 0.173648177666930
#define p_PBAR 100
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

// Inviscid Riemann Solve BC's
// ************************************************************************
#define cnsInviscidBoundaryConditions3D(bc, gamma, R, CP, CV, mu, t, x, y, z,\
                                        nx, ny, nz,                          \
                                        rM, ruM, rvM, rwM, reM,              \
                                        rB, ruB, rvB, rwB, reB)              \
{                                                                            \
const dfloat uM    = ruM/rM;                                                 \
const dfloat vM    = rvM/rM;                                                 \
const dfloat wM    = rwM/rM;                                                 \
const dfloat pM    = (gamma-1.0)*(reM - 0.5*rM*(uM*uM+vM*vM+wM*wM));         \
const dfloat unM   = uM*nx + vM*ny + wM*nz;                                  \
const dfloat cM    = sqrt(gamma*pM/rM);                                      \
const dfloat mM    = fabs(unM/cM);                                           \
const dfloat keREF = 0.5*p_RBAR*(p_UBAR*p_UBAR+p_VBAR*p_VBAR+p_WBAR*p_WBAR); \
  if(bc==11){                                                                \
    *( rB) = rM;                                                             \
    *(ruB) = -ruM;                                                           \
    *(rvB) = -rvM;                                                           \
    *(rwB) = -rwM;                                                           \
    *(reB) = rM*CP*p_TBAR/gamma + 0.5*(ruM*ruM + rvM*rvM + rwM*rwM)/rM;      \
   }else if(bc==12){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = -ruM;                                                           \
    *(rvB) = -rvM;                                                           \
    *(rwB) = -rwM;                                                           \
    *(reB) = reM;                                                            \
  }else if(bc==13){                                                          \
    *( rB) =  rM;                                                            \
    *(ruB) =  ruM - 2.0*(nx*ruM+ny*rvM+nz*rwM)*nx;                           \
    *(rvB) =  rvM - 2.0*(nx*ruM+ny*rvM+nz*rwM)*ny;                           \
    *(rwB) =  rwM - 2.0*(nx*ruM+ny*rvM+nz*rwM)*nz;                           \
    *(reB) =  reM;                                                           \
   }else if(bc==20){                                                         \
    if(unM<=0){                                                              \
      if(mM > 1.0){                                                          \
        *( rB) = p_RBAR;                                                     \
        *(ruB) = p_RBAR*p_UBAR;                                              \
        *(rvB) = p_RBAR*p_VBAR;                                              \
        *(rvB) = p_RBAR*p_WBAR;                                              \
        *(reB) = p_PBAR/(gamma -1.0) + keREF;                                \
      }else{                                                                 \
        *( rB) = p_RBAR;                                                     \
        *(ruB) = p_RBAR*p_UBAR;                                              \
        *(rvB) = p_RBAR*p_VBAR;                                              \
        *(rwB) = p_RBAR*p_WBAR;                                              \
        *(reB) = reM-0.5*rM*(uM*uM+vM*vM+wM*wM)+keREF;                       \
      }                                                                      \
    }else{                                                                   \
      if(mM > 1.0){                                                          \
      *( rB) = rM;                                                           \
      *(ruB) = ruM;                                                          \
      *(rvB) = rvM;                                                          \
      *(rwB) = rwM;                                                          \
      *(reB) = reM;                                                          \
      }else{                                                                 \
      *( rB) = rM;                                                           \
      *(ruB) = ruM;                                                          \
      *(rvB) = rvM;                                                          \
      *(rwB) = rwM;                                                          \
      *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM+vM*vM+wM*wM);             \
      }                                                                      \
    }                                                                        \
   }else if(bc==21){                                                         \
    *( rB) = p_RBAR;                                                         \
    *(ruB) = p_RBAR*p_UBAR;                                                  \
    *(rvB) = p_RBAR*p_VBAR;                                                  \
    *(rwB) = p_RBAR*p_WBAR;                                                  \
    *(reB) = reM-0.5*rM*(uM*uM+vM*vM+wM*wM)+keREF ;                          \
   }else if(bc==22){                                                         \
    *( rB) = p_RBAR;                                                         \
    *(ruB) = p_RBAR*p_UBAR;                                                  \
    *(rvB) = p_RBAR*p_VBAR;                                                  \
    *(rwB) = p_RBAR*p_WBAR;                                                  \
    *(reB) = p_PBAR/(gamma -1.0) + keREF;                                    \
   }else if(bc==31){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = ruM;                                                            \
    *(rvB) = rvM;                                                            \
    *(rwB) = rwM;                                                            \
    *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM+wM*wM);             \
   }else if(bc==32){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = ruM;                                                            \
    *(rvB) = rvM;                                                            \
    *(rwB) = rwM;                                                            \
    *(reB) = reM;                                                            \
 }                                                                           \
}

// Viscous Riemann Solve BC's
// ************************************************************************
#define cnsViscousBoundaryConditions3D(bc, gamma, R, CP, CV, mu, t, x, y, z, \
                                      nx, ny, nz,                            \
                                      rM, ruM, rvM, rwM, reM,                \
                                      rB, ruB, rvB, rwB, reB,                \
                                      drrdxM, drrdyM, drrdzM,                \
                                      drudxM, drudyM, drudzM,                \
                                      drvdxM, drvdyM, drvdzM,                \
                                      drwdxM, drwdyM, drwdzM,                \
                                      dredxM, dredyM, dredzM,                \
                                      drrdxB, drrdyB, drrdzB,                \
                                      drudxB, drudyB, drudzB,                \
                                      drvdxB, drvdyB, drvdzB,                \
                                      drwdxB, drwdyB, drwdzB,                \
                                      dredxB, dredyB, dredzB)                \
{                                                                            \
const dfloat uM  = ruM/rM;                                                   \
const dfloat vM =  rvM/rM;                                                   \
const dfloat wM =  rwM/rM;                                                   \
const dfloat pM  = (gamma-1.0)*(reM-0.5*rM*(uM*uM+vM*vM+wM*wM));             \
const dfloat unM = uM*nx+vM*ny+wM*nz;                                        \
const dfloat cM  = sqrt(gamma*pM/rM);                                        \
const dfloat mM  = fabs(unM/cM);                                             \
const dfloat keREF = 0.5*p_RBAR*(p_UBAR*p_UBAR+p_VBAR*p_VBAR+p_WBAR*p_WBAR); \
  if(bc==11){                                                                \
    *( rB) = rM;                                                             \
    *(ruB) = 0.0;                                                            \
    *(rvB) = 0.0;                                                            \
    *(rwB) = 0.0;                                                            \
    *(reB) = rM*CP*p_TBAR/gamma;                                             \
    *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;*(drrdzB) = drrdzM;                \
    *(drudxB) = drudxM;*(drudyB) = drudyM;*(drudzB) = drudzM;                \
    *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;*(drvdzB) = drvdzM;                \
    *(drwdxB) = drwdxM;*(drwdyB) = drwdyM;*(drwdzB) = drwdzM;                \
    *(dredxB) = dredxM;*(dredyB) = dredyM;*(dredzB) = dredzM;                \
   }else if(bc==12){                                                         \
    *( rB) = rM;                                                             \
    *(ruB) = 0.0;                                                            \
    *(rvB) = 0.0;                                                            \
    *(rwB) = 0.0;                                                            \
    *(reB) = reM-0.5*(ruM*ruM+rvM*rvM+rwM*rwM)/rM;                           \
    const dfloat uxM = drudxM - uM*drrdxM;                                   \
    const dfloat uyM = drudyM - uM*drrdyM;                                   \
    const dfloat uzM = drudzM - uM*drrdzM;                                   \
    const dfloat vxM = drvdxM - vM*drrdxM;                                   \
    const dfloat vyM = drvdyM - vM*drrdyM;                                   \
    const dfloat vzM = drvdzM - vM*drrdzM;                                   \
    const dfloat wxM = drwdxM - wM*drrdxM;                                   \
    const dfloat wyM = drwdyM - wM*drrdyM;                                   \
    const dfloat wzM = drwdzM - wM*drrdzM;                                   \
    const dfloat TxM = dredxM - (drrdxM*reM/rrM + uM*uxM + vM*vxM+ wM*wxM);  \
    const dfloat TyM = dredyM - (drrdyM*reM/rrM + uM*uyM + vM*vyM+ wM*wyM);  \
    const dfloat TzM = dredzM - (drrdzM*reM/rrM + uM*uzM + vM*vzM+ wM*wzM);  \
    *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;*(drrdzB) = drrdzM;                \
    *(drudxB) = drudxM;*(drudyB) = drudyM;*(drudzB) = drudzM;                \
    *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;*(drvdzB) = drvdzM;                \
    *(drwdxB) = drwdxM;*(drwdyB) = drwdyM;*(drwdzB) = drwdzM;                \
    *(dredxB) = dredxM - (nx*nx*TxM + nx*ny*TyM+ nx*nz*TzM);                 \                                                     \
    *(dredyB) = dredyM - (ny*nx*TxM + ny*ny*TyM+ ny*nz*TzM);                 \
    *(dredzB) = dredzM - (nz*nx*TxM + nz*ny*TyM+ nz*nz*TzM);                 \
  }else if(bc==13){                                                          \
    *( rB) =  rM;                                                            \
    *(ruB) =  ruM - (nx*ruM+ny*rvM+nz*rwM)*nx;                               \
    *(rvB) =  rvM - (nx*ruM+ny*rvM+nz*rwM)*ny;                               \
    *(rwB) =  rwM - (nx*ruM+ny*rvM+nz*rwM)*nz;                               \
    *(reB) =  reM;                                                           \
    *(drrdxB) = 0.0;*(drrdyB) = 0.0;*(drrdzB) = 0.0;                         \
    *(drudxB) = 0.0;*(drudyB) = 0.0;*(drudzB) = 0.0;                         \
    *(drvdxB) = 0.0;*(drvdyB) = 0.0;*(drvdzB) = 0.0;                         \
    *(drwdxB) = 0.0;*(drwdyB) = 0.0;*(drwdzB) = 0.0;                         \
    *(dredxB) = 0.0;*(dredyB) = 0.0;*(dredzB) = 0.0;                         \
  }else if(bc==20){                                                          \
    if(unM<=0){                                                              \
      if(mM > 1.0){                                                          \
        *( rB) = p_RBAR;                                                     \
        *(ruB) = p_RBAR*p_UBAR;                                              \
        *(rvB) = p_RBAR*p_VBAR;                                              \
        *(rwB) = p_RBAR*p_WBAR;                                              \
        *(reB) = p_PBAR/(gamma -1.0) + keREF;                                \
        *(drrdxB) = 0.0;*(drrdyB) = 0.0;*(drrdzB) = 0.0;                     \
        *(drudxB) = 0.0;*(drudyB) = 0.0;*(drudzB) = 0.0;                     \
        *(drvdxB) = 0.0;*(drvdyB) = 0.0;*(drvdzB) = 0.0;                     \
        *(drwdxB) = 0.0;*(drwdyB) = 0.0;*(drwdzB) = 0.0;                     \
        *(dredxB) = 0.0;*(dredyB) = 0.0;*(dredzB) = 0.0;                     \
      }else{                                                                 \
        *( rB) = p_RBAR;                                                     \
        *(ruB) = p_RBAR*p_UBAR;                                              \
        *(rvB) = p_RBAR*p_VBAR;                                              \
        *(rwB) = p_RBAR*p_WBAR;                                              \
        *(reB) = reM-0.5*rM*(uM*uM+vM*vM+wM*wM)+keREF;                       \
        *(drrdxB) = 0.0;*(drrdyB) = 0.0;*(drrdzB) = 0.0;                     \
        *(drudxB) = 0.0;*(drudyB) = 0.0;*(drudzB) = 0.0;                     \
        *(drvdxB) = 0.0;*(drvdyB) = 0.0;*(drvdzB) = 0.0;                     \
        *(drwdxB) = 0.0;*(drwdyB) = 0.0;*(drwdzB) = 0.0;                     \
        *(dredxB) = 0.0;*(dredyB) = 0.0;*(dredzB) = 0.0;                     \
      }                                                                      \
    }else{                                                                   \
      if(mM > 1.0){                                                          \
        *( rB) = rM;                                                         \
        *(ruB) = ruM;                                                        \
        *(rvB) = rvM;                                                        \
        *(rwB) = rwM;                                                        \
        *(reB) = reM;                                                        \
        *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;*(drrdzB) = drrdzM;            \
        *(drudxB) = drudxM;*(drudyB) = drudyM;*(drudzB) = drudzM;            \
        *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;*(drvdzB) = drvdzM;            \
        *(drwdxB) = drwdxM;*(drwdyB) = drwdyM;*(drwdzB) = drwdzM;            \
        *(dredxB) = dredxM;*(dredyB) = dredyM;*(dredzB) = dredzM;            \
        }else{                                                               \
        *( rB) = rM;                                                         \
        *(ruB) = ruM;                                                        \
        *(rvB) = rvM;                                                        \
        *(rwB) = rwM;                                                        \
        *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM+vM*vM+wM*wM);           \
        *(drrdxB) = 0.0;*(drrdyB) = 0.0;*(drrdzB) = 0.0;                     \
        *(drudxB) = 0.0;*(drudyB) = 0.0;*(drudzB) = 0.0;                     \
        *(drvdxB) = 0.0;*(drvdyB) = 0.0;*(drvdzB) = 0.0;                     \
        *(drwdxB) = 0.0;*(drwdyB) = 0.0;*(drwdzB) = 0.0;                     \
        *(dredxB) = 0.0;*(dredyB) = 0.0;*(dredzB) = 0.0;                     \
      }                                                                      \
    }                                                                        \
  }else if(bc==21){                                                          \
    *( rB) = p_RBAR;                                                         \
    *(ruB) = p_RBAR*p_UBAR;                                                  \
    *(rvB) = p_RBAR*p_VBAR;                                                  \
    *(rwB) = p_RBAR*p_WBAR;                                                  \
    *(reB) = reM-0.5*rM*(uM*uM+vM*vM+wM*wM)+keREF;                           \
    *(drrdxB) = 0.0;*(drrdyB) = 0.0;*(drrdzB) = 0.0;                         \
    *(drudxB) = 0.0;*(drudyB) = 0.0;*(drudzB) = 0.0;                         \
    *(drvdxB) = 0.0;*(drvdyB) = 0.0;*(drvdzB) = 0.0;                         \
    *(drwdxB) = 0.0;*(drwdyB) = 0.0;*(drwdzB) = 0.0;                         \
    *(dredxB) = 0.0;*(dredyB) = 0.0;*(dredzB) = 0.0;                         \
  }else if(bc==22){                                                          \
    *( rB) = p_RBAR;                                                         \
    *(ruB) = p_RBAR*p_UBAR;                                                  \
    *(rvB) = p_RBAR*p_VBAR;                                                  \
    *(rwB) = p_RBAR*p_WBAR;                                                  \
    *(reB) = p_PBAR/(gamma -1.0) + keREF;                                    \
    *(drrdxB) = 0.0;*(drrdyB) = 0.0;*(drrdzB) = 0.0;                         \
    *(drudxB) = 0.0;*(drudyB) = 0.0;*(drudzB) = 0.0;                         \
    *(drvdxB) = 0.0;*(drvdyB) = 0.0;*(drvdzB) = 0.0;                         \
    *(drwdxB) = 0.0;*(drwdyB) = 0.0;*(drwdzB) = 0.0;                         \
    *(dredxB) = 0.0;*(dredyB) = 0.0;*(dredzB) = 0.0;                         \
  }else if(bc==31){                                                          \
    *( rB) = rM;                                                             \
    *(ruB) = ruM;                                                            \
    *(rvB) = rvM;                                                            \
    *(rwB) = rwM;                                                            \
    *(reB) = p_PBAR/(gamma -1.0) + 0.5*rM*(uM*uM + vM*vM + wM*wM);           \
    *(drrdxB) = 0.0;*(drrdyB) = 0.0;*(drrdzB) = 0.0;                         \
    *(drudxB) = 0.0;*(drudyB) = 0.0;*(drudzB) = 0.0;                         \
    *(drvdxB) = 0.0;*(drvdyB) = 0.0;*(drvdzB) = 0.0;                         \
    *(drwdxB) = 0.0;*(drwdyB) = 0.0;*(drwdzB) = 0.0;                         \
    *(dredxB) = 0.0;*(dredyB) = 0.0;*(dredzB) = 0.0;                         \
  }else if(bc==32){                                                          \
    *( rB) = rM;                                                             \
    *(ruB) = ruM;                                                            \
    *(rvB) = rvM;                                                            \
    *(rwB) = rwM;                                                            \
    *(reB) = reM;                                                            \
    *(drrdxB) = drrdxM;*(drrdyB) = drrdyM;*(drrdzB) = drrdzM;                \
    *(drudxB) = drudxM;*(drudyB) = drudyM;*(drudzB) = drudzM;                \
    *(drvdxB) = drvdxM;*(drvdyB) = drvdyM;*(drvdzB) = drvdzM;                \
    *(drwdxB) = drwdxM;*(drwdyB) = drwdyM;*(drwdzB) = drwdzM;                \
    *(dredxB) = dredxM;*(dredyB) = dredyM;*(dredzB) = dredzM;                \
  }                                                                          \
}
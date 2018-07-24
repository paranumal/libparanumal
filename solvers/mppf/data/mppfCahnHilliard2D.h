// Initial conditions : f and g are forcing for exact solution 
// messy manufactured solution
#define mppfFlowField2D(t,x,y,u,v,p,phi) \
  {                         \
    dfloat epsilon = 0.1;   \
    dfloat lambda  = 0.001; \
    dfloat M       = 0.001; \
    *(u)    =  occaCos(OCCy)*occaSin(OCCA_PI*x)*sin(t);\
    *(v)    = -occaSin(OCCy)*occaCos(OCCA_PI*x)*sin(t);\
    *(p)    =  occaSin(OCCy)*occaSin(OCCA_PI*x)*cos(t);\
    *(phi)  =  occaCos(OCCx)*occaCos(OCCA_PI*y)*sin(t);\
  } 

// #define mppfSourceField2D(t,x,y,f,g) \
//   {                                    \
//     dfloat epsilon = 0.1; \
//     dfloat lambda  = 0.001; \
//     dfloat M       = 0.001; \
//     dfloat u       =  cos(OCCA_PI*y)*sin(OCCA_PI*x)*sin(t);\
//     dfloat v       = -sin(OCCA_PI*y)*cos(OCCA_PI*x)*sin(t);\
//     dfloat p       =  sin(OCCA_PI*y)*sin(OCCA_PI*x)*cos(t);\
//     dfloat phi     =  cos(OCCA_PI*x)*cos(OCCA_PI*y)*sin(t);\
//     dfloat phit    = cos(OCCA_PI*x)*cos(OCCA_PI*y)*cos(t); \
//     dfloat phix    = -OCCA_PI*cos(OCCA_PI*y)*sin(OCCA_PI*x)*sin(t); \
//     dfloat phiy    = -OCCA_PI*cos(OCCA_PI*x)*sin(OCCA_PI*y)*sin(t); \
//     dfloat h       = phi*(phi*phi-1.f)/(epsilon*epsilon); \
//     dfloat phi2xx = 2.f*OCCA_PI*OCCA_PI*OCCA_PI*OCCA_PI*cos(OCCA_PI*x)*cos(OCCA_PI*y)*sin(t); \
//     dfloat phi2yy = 2.f*OCCA_PI*OCCA_PI*OCCA_PI*OCCA_PI*cos(OCCA_PI*x)*cos(OCCA_PI*y)*sin(t); \
//     dfloat hxx1    = cos(OCCA_PI*x)*cos(OCCA_PI*y)*sin(t); \
//     dfloat hxx2    = cos(OCCA_PI*y)*sin(OCCA_PI*x)*sin(t); \
//     dfloat hxx     = (OCCA_PI*OCCA_PI*hxx1*(6.f*hxx2*hxx2 - 3.f*hxx1*hxx1 + 1.f))/(epsilon*epsilon); \
//     dfloat hyy1    = cos(OCCA_PI*x)*cos(OCCA_PI*y)*sin(t); \
//     dfloat hyy2    = cos(OCCA_PI*x)*sin(OCCA_PI*y)*sin(t) ; \
//     dfloat hyy     = (OCCA_PI*OCCA_PI*hyy1*(6.f*hyy2*hyy2 - 3.f*hyy1*hyy1 + 1.f))/(epsilon*epsilon); \
//     *(g)           = phit + u*phix + v*phiy + lambda*M*(phi2xx + phi2yy - hxx -hyy);\
//     *(f)           = 0.f;\
//   } 

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
/*  Inflow condition should be phiM to create zero normal derivative ? */
#define mppfPhaseFieldDirichletConditions2D(bc, t, x, y, nx, ny, phiM, phiB) \
{                                   \
  if(bc==1){                        \
    *(phiB) = phiM;                  \
  } else if(bc==2){                 \
    *(phiB) = cos(OCCA_PI*x)*cos(OCCA_PI*y)*sin(t);\
  } else if(bc==3){                 \
    *(phiB) = phiM;                 \
  } else if(bc==4){                 \
    *(phiB) = phiM;                 \
  } else if(bc==5){                 \
    *(phiB) = phiM;                 \
  }                                 \
}

#define mppfPhaseFieldNeumannConditions2D(bc, t, x, y, nx, ny, phixM, phiyM, phixB, phiyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(phixB) = 0.f;                        \
    *(phiyB) = 0.f;                        \
  } else if(bc==3){                        \
    *(phixB) =0.f;                         \
    *(phiyB) =0.f;                         \
  } else if(bc==4){                        \
    *(phixB) = 0.f;                        \
    *(phiyB) = 0.f;                        \
  } else if(bc==5){                        \
    *(phixB) = 0.f;                        \
    *(phiyB) = 0.f;                        \
  }                                        \
}


// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define mppfVelocityDirichletConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
{                                   \
  if(bc==1){                        \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
  } else if(bc==2){                 \
    *(u)    =  cos(OCCA_PI*y)*sin(OCCA_PI*x)*sin(t);\
    *(v)    = -sin(OCCA_PI*y)*cos(OCCA_PI*x)*sin(t);\
  } else if(bc==3){                 \
    *(uB) = uM;                     \
    *(vB) = vM;                     \
  } else if(bc==4){                 \
    *(uB) = 0.f;                    \
    *(vB) = vM;                     \
  } else if(bc==5){                 \
    *(uB) = uM;                     \
    *(vB) = 0.f;                    \
  }                                 \
}

#define mppfVelocityNeumannConditions2D(bc, t, x, y, nx, ny, uxM, uyM, vxM, vyM, uxB, uyB, vxB, vyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
    *(vxB) = vxM;                          \
    *(vyB) = vyM;                          \
  } else if(bc==3){                        \
    *(uxB) = OCCA_PI*cos(OCCA_PI*x)*cos(OCCA_PI*y)*sin(t);\
    *(uyB) =-OCCA_PI*sin(OCCA_PI*x)*sin(OCCA_PI*y)*sin(t);\
    *(vxB) = OCCA_PI*sin(OCCA_PI*x)*sin(OCCA_PI*y)*sin(t);\
    *(vyB) =-OCCA_PI*cos(OCCA_PI*x)*cos(OCCA_PI*y)*sin(t);\
  } else if(bc==4){                        \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
    *(vxB) = 0.f;                          \
    *(vyB) = 0.f;                          \
  } else if(bc==5){                        \
    *(uxB) = 0.f;                          \
    *(uyB) = 0.f;                          \
    *(vxB) = vxM;                          \
    *(vyB) = vyM;                          \
  }                                        \
}


#define mppfPressureDirichletConditions2D(bc, t, x, y, nx, ny, pM, pB) \
{                                   \
  if(bc==1 || bc==2){               \
    *(pB) = pM;                     \
  } else if(bc==3){                 \
    *(pB) = sin(OCCA_PI*y)*sin(OCCA_PI*x)*cos(t);\
  } else if(bc==4){                 \
    *(pB) = pM;                     \
  } else if(bc==5){                 \
    *(pB) = pM;                     \
  }                                 \
}

#define mppfPressureNeumannConditions2D(bc, t, x, y, nx, ny, pxM, pyM, pxB, pyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(pxB) = OCCA_PI*cos(OCCA_PI*x)*sin(OCCA_PI*y)*cos(t);\
    *(pyB) = OCCA_PI*cos(OCCA_PI*y)*sin(OCCA_PI*x)*cos(t);\
  } else if(bc==3){                        \
    *(pxB) = pxM;                          \
    *(pyB) = pyM;                          \
  } else if(bc==4){                        \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
  } else if(bc==5){                        \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
  }                                        \
}

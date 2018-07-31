// Initial conditions : f and g are forcing for exact solution 
// messy manufactured solution
#define mppfFlowField2D(t,x,y,u,v,p) \
  {                         \
    *(u)    =  cos(M_PI*y)*sin(M_PI*x)*sin(t);\
    *(v)    = -sin(M_PI*y)*cos(M_PI*x)*sin(t);\
    *(p)    =  sin(M_PI*y)*sin(M_PI*x)*cos(t);\
  }


#define mppfPhaseField2D(t,h,x,y,phi) \
  {                         \
    *(phi)  =  cos(M_PI*x)*cos(M_PI*y)*sin(t);\
  }

#define mppfPhaseFieldSource2D(t,x,y,g)\
  {                                    \
    dfloat eta     = 0.1; \
    dfloat lambda  = 0.001; \
    dfloat M       = 0.001; \
    dfloat eta2    = eta*eta;\
    dfloat u       =  cos(M_PI*y)*sin(M_PI*x)*sin(t);\
    dfloat v       = -sin(M_PI*y)*cos(M_PI*x)*sin(t);\
    dfloat p       =  sin(M_PI*y)*sin(M_PI*x)*cos(t);\
    dfloat phi     =  cos(M_PI*x)*cos(M_PI*y)*sin(t);\
    dfloat phit    = cos(M_PI*x)*cos(M_PI*y)*cos(t);\
    dfloat phix    = -M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(t);\
    dfloat phiy    = -M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(t);\
    dfloat dt1     = sin(M_PI*x)*sin(t);\
    dfloat dt2     = sin(M_PI*y)*sin(t); \
    dfloat dt3     = sin(M_PI*x)*sin(M_PI*y)*sin(t);\
    dfloat dif     = (2.f*M*M_PI*M_PI*lambda*phi*(3.f*sin(t)*sin(t) - 6.f*dt1*dt1 - 6.f*dt2*dt2 + 2.f*M_PI*M_PI*eta2 + 9.f*dt3*dt3 - 1.f))/eta2; \
    *(g)           = phit + u*phix + v*phiy + dif;\
  } 




  #define mppfPhaseFieldFixTermNu2D(t,x,y,nu)\
  {                                    \
    dfloat eta     = 0.1; \
    dfloat lambda  = 0.001; \
    dfloat M       = 0.001; \
    dfloat eta2    = eta*eta;\
    dfloat u       =  cos(M_PI*y)*sin(M_PI*x)*sin(t);\
    dfloat v       = -sin(M_PI*y)*cos(M_PI*x)*sin(t);\
    dfloat p       =  sin(M_PI*y)*sin(M_PI*x)*cos(t);\
    dfloat phi     =  cos(M_PI*x)*cos(M_PI*y)*sin(t);\
    dfloat phix    = -M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(t);\
    dfloat phiy    = -M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(t);\
    *(nu)          =  u*phix + v*phiy;\
  } 



   #define mppfPhaseFieldFixTermLap2D(t,x,y,inveta2, chSeta2, lap)\
  {                                    \
    dfloat eta     = 0.1; \
    dfloat lambda  = 0.001; \
    dfloat M       = 0.001; \
    dfloat eta2    = eta*eta;\
    dfloat phi     =  cos(M_PI*x)*cos(M_PI*y)*sin(t);\
    dfloat dt1     = cos(M_PI*y)*sin(M_PI*x)*sin(t);\
    *(lap)         =  (2.f*M_PI*M_PI*phi*(chSeta2/inveta2 - 3.f*phi*phi + 6.f*dt1*dt1 + 1.f))*inveta2;\
  } 


// (2*M_PI^2*cos(M_PI*x)*cos(M_PI*y)*sin(t)*(chSeta2*eta2 - 3*cos(M_PI*x)^2*cos(M_PI*y)^2*sin(t)^2 + 6*cos(M_PI*y)^2*sin(M_PI*x)^2*sin(t)^2 + 1))/eta2
 

// dfloat term1   = 2.f*M*M_PI*M_PI*M_PI*M_PI*lambda*cos(M_PI*x)*cos(M_PI*y)*sin(t);\
    dfloat h1      = sin(M_PI*x)*sin(M_PI*y)*sin(t);\
    dfloat term2   = -(2.f*M*lambda*M_PI*M_PI*phi*(3.f*sin(t)*sin(t) + 9.f*h1*h1 + 1.f))/eta2;\
// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
/*  Inflow condition should be phiM to create zero normal derivative ? */



#define mppfPhaseFieldDirichletConditions2D(bc, t, x, y, nx, ny, phiM, phiB) \
{                                   \
  if(bc==1){                        \
    *(phiB) = phiM;                 \
  } else if(bc==2){                 \
    *(phiB) = phiM;                 \
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
    *(u)    =  cos(M_PI*y)*sin(M_PI*x)*sin(t);\
    *(v)    = -sin(M_PI*y)*cos(M_PI*x)*sin(t);\
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
    *(uxB) = M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(t);\
    *(uyB) =-M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(t);\
    *(vxB) = M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(t);\
    *(vyB) =-M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(t);\
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
    *(pB) = sin(M_PI*y)*sin(M_PI*x)*cos(t);\
  } else if(bc==4){                 \
    *(pB) = pM;                     \
  } else if(bc==5){                 \
    *(pB) = pM;                     \
  }                                 \
}

#define mppfPressureNeumannConditions2D(bc, t, x, y, nx, ny, pxM, pyM, pxB, pyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(pxB) = M_PI*cos(M_PI*x)*sin(M_PI*y)*cos(t);\
    *(pyB) = M_PI*cos(M_PI*y)*sin(M_PI*x)*cos(t);\
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


// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define boundaryConditionsPML3D(bc, t, x, y, z, nx, ny, nz, intfx, intfy, intfz, q1M, q2M, q3M, q4M, q5M, q6M, q7M, q8M, q9M, q10M, q1B, q2B, q3B, q4B, q5B, q6B, q7B, q8B, q9B, q10B) \
{                                           \
  if(bc==1){                                \
    *(q1B)  =  q1M;                         \
    *(q2B)  = -q2M;                         \
    *(q3B)  = -q3M;                         \
    *(q4B)  = -q4M;                         \
    *(q5B)  =  q5M;                         \
    *(q6B)  =  q6M;                         \
    *(q7B)  =  q7M;                         \
    *(q8B)  =  q8M;                         \
    *(q9B)  =  q9M;                         \
    *(q10B) =  q10M;                        \
  } else if(bc==4||bc==5){                  \
    *(q1B)  = q1M;                           \
    *(q2B)  = q2M-2.f*(nx*q2M+ny*q3M*nz*q4M)*nx;\
    *(q3B)  = q3M-2.f*(nx*q2M+ny*q3M*nz*q4M)*ny;\
    *(q4B)  = q4M-2.f*(nx*q2M+ny*q3M*nz*q4M)*nz;\
    *(q5B)  = q5M;                           \
    *(q6B)  = q6M;                           \
    *(q7B)  = q7M;                           \
    *(q8B)  = q8M;                           \
    *(q9B)  = q9M;                           \
    *(q10B) = q10M;                          \
  }else{								\
    *(q1B) = 2.f*p_q1bar - q1M;						\
    *(q2B) = 2.f*p_q1bar*intfx*p_isqrtRT - q2M;				\
    *(q3B) = 2.f*p_q1bar*intfy*p_isqrtRT - q3M;				\
    *(q4B) = 2.f*p_q1bar*intfz*p_isqrtRT - q4M;				\
    *(q5B) = 2.f*p_q1bar*intfx*intfy*p_isqrtRT*p_isqrtRT-q5M;		\
    *(q6B) = 2.f*p_q1bar*intfx*intfz*p_isqrtRT*p_isqrtRT-q6M;		\
    *(q7B) = 2.f*p_q1bar*intfy*intfz*p_isqrtRT*p_isqrtRT-q7M;		\
    *(q8B) = 2.f*p_q1bar*intfx*intfx*p_isqrtRT*p_isqrtRT*p_invsqrt2-q8M; \
    *(q9B) = 2.f*p_q1bar*intfy*intfy*p_isqrtRT*p_isqrtRT*p_invsqrt2-q9M; \
    *(q10B) = 2.f*p_q1bar*intfz*intfz*p_isqrtRT*p_isqrtRT*p_invsqrt2-q10M; \
  }									\
}

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define boundaryConditions3D(bc, t, x, y, z, nx, ny, nz, intfx, intfy, intfz, q1M, q2M, q3M, q4M, q5M, q6M, q7M, q8M, q9M, q10M, q1B, q2B, q3B, q4B, q5B, q6B, q7B, q8B, q9B, q10B) \
{                                           \
  if(bc==1){                                \
    *(q1B)  =  q1M;                         \
    *(q2B)  = -q2M;                         \
    *(q3B)  = -q3M;                         \
    *(q4B)  = -q4M;                         \
    *(q5B)  =  q5M;                         \
    *(q6B)  =  q6M;                         \
    *(q7B)  =  q7M;                         \
    *(q8B)  =  q8M;                         \
    *(q9B)  =  q9M;                         \
    *(q10B) =  q10M;                        \
  } else if(bc==2 || bc==3){		    \
    *(q1B) = 2.f*p_q1bar - q1M;						\
    *(q2B) = 2.f*p_q1bar*intfx*p_isqrtRT - q2M;				\
    *(q3B) = 2.f*p_q1bar*intfy*p_isqrtRT - q3M;				\
    *(q4B) = 2.f*p_q1bar*intfz*p_isqrtRT - q4M;				\
    *(q5B) = 2.f*p_q1bar*intfx*intfy*p_isqrtRT*p_isqrtRT-q5M;		\
    *(q6B) = 2.f*p_q1bar*intfx*intfz*p_isqrtRT*p_isqrtRT-q6M;		\
    *(q7B) = 2.f*p_q1bar*intfy*intfz*p_isqrtRT*p_isqrtRT-q7M;		\
    *(q8B) = 2.f*p_q1bar*intfx*intfx*p_isqrtRT*p_isqrtRT*p_invsqrt2-q8M; \
    *(q9B) = 2.f*p_q1bar*intfy*intfy*p_isqrtRT*p_isqrtRT*p_invsqrt2-q9M; \
    *(q10B) = 2.f*p_q1bar*intfz*intfz*p_isqrtRT*p_isqrtRT*p_invsqrt2-q10M; \
  } else if(bc==4||bc==5){                  \
    *(q1B)  = q1M;                           \
    *(q2B)  = q2M-2.f*(nx*q2M+ny*q3M*nz*q4M)*nx;\
    *(q3B)  = q3M-2.f*(nx*q2M+ny*q3M*nz*q4M)*ny;\
    *(q4B)  = q4M-2.f*(nx*q2M+ny*q3M*nz*q4M)*nz;\
    *(q5B)  = q5M;                           \
    *(q6B)  = q6M;                           \
    *(q7B)  = q7M;                           \
    *(q8B)  = q8M;                           \
    *(q9B)  = q9M;                           \
    *(q10B) = q10M;                          \
  }                                         \
}

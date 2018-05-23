
// Initial conditions: not in use cirrently 
#define bnsFlowField3D(t,x,y,ramp,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10)  \
  {                                                   \
    *(q1)  = q1bar;                                    \
    *(q2)  = ramp*q2bar;                               \
    *(q3)  = ramp*q3bar;                               \
    *(q4)  = ramp*q4bar;                               \
    *(q5)  = ramp*ramp*q5bar;                          \
    *(q6)  = ramp*ramp*q6bar;                          \
    *(q7)  = ramp*ramp*q7bar;                          \
    *(q8)  = ramp*ramp*q8bar;                          \
    *(q9)  = ramp*ramp*q9bar;                          \
    *(q10) = ramp*ramp*q10bar;                         \
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define boundaryConditionsPML3D(bc, t, x, y, z, nx, ny, nz, ramp, q1M, q2M, q3M, q4M, q5M, q6M, q7M, q8M, q9M, q10M, q1B, q2B, q3B, q4B, q5B, q6B, q7B, q8B, q9B, q10B) \
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
    }else{                                   \
    *(q1B)  = 2.f*p_q1bar - q1M;            \
    *(q2B)  = 2.f*ramp*p_q2bar - q2M;       \
    *(q3B)  = 2.f*ramp*p_q3bar - q3M;       \
    *(q4B)  = 2.f*ramp*p_q4bar - q4M;       \
    *(q5B)  = 2.f*ramp*ramp*p_q5bar  - q5M; \
    *(q6B)  = 2.f*ramp*ramp*p_q6bar  - q6M; \
    *(q7B)  = 2.f*ramp*ramp*p_q7bar  - q7M; \
    *(q8B)  = 2.f*ramp*ramp*p_q8bar  - q8M; \
    *(q9B)  = 2.f*ramp*ramp*p_q9bar  - q9M; \
    *(q10B) = 2.f*ramp*ramp*p_q10bar - q10M;\
  }                                         \
}

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define boundaryConditions3D(bc, t, x, y, z, nx, ny, nz, ramp, q1M, q2M, q3M, q4M, q5M, q6M, q7M, q8M, q9M, q10M, q1B, q2B, q3B, q4B, q5B, q6B, q7B, q8B, q9B, q10B) \
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
  } else if(bc==2){                         \
    *(q1B)  = 2.f*p_q1bar - q1M;            \
    *(q2B)  = 2.f*ramp*p_q2bar - q2M;       \
    *(q3B)  = 2.f*ramp*p_q3bar - q3M;       \
    *(q4B)  = 2.f*ramp*p_q4bar - q4M;       \
    *(q5B)  = 2.f*ramp*ramp*p_q5bar  - q5M; \
    *(q6B)  = 2.f*ramp*ramp*p_q6bar  - q6M; \
    *(q7B)  = 2.f*ramp*ramp*p_q7bar  - q7M; \
    *(q8B)  = 2.f*ramp*ramp*p_q8bar  - q8M; \
    *(q9B)  = 2.f*ramp*ramp*p_q9bar  - q9M; \
    *(q10B) = 2.f*ramp*ramp*p_q10bar - q10M;\
  } else if(bc==3){                         \
    *(q1B)  = 2.f*q1M - q1M;                \
    *(q2B)  = 2.f*ramp*p_q2bar - q2M;       \
    *(q3B)  = 2.f*ramp*p_q3bar - q3M;       \
    *(q4B)  = 2.f*ramp*p_q4bar - q4M;       \
    *(q5B)  = 2.f*ramp*ramp*p_q5bar  - q5M; \
    *(q6B)  = 2.f*ramp*ramp*p_q6bar  - q6M; \
    *(q7B)  = 2.f*ramp*ramp*p_q7bar  - q7M; \
    *(q8B)  = 2.f*ramp*ramp*p_q8bar  - q8M; \
    *(q9B)  = 2.f*ramp*ramp*p_q9bar  - q9M; \
    *(q10B) = 2.f*ramp*ramp*p_q10bar - q10M;\
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
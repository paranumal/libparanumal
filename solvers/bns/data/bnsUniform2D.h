
// Initial conditions: not in use cirrently 
#define bnsFlowField2D(t,x,y,ramp,q1,q2,q3,q4,q5,q6)  \
  {                                                   \
    *(q1) = q1bar;                                    \
    *(q2) = ramp*q2bar;                               \
    *(q3) = ramp*q3bar;                               \
    *(q4) = ramp*ramp*q4bar;                          \
    *(q5) = ramp*ramp*q5bar;                          \
    *(q6) = ramp*ramp*q6bar;                          \
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define boundaryConditionsPML2D(bc, t, x, y, nx, ny, ramp, q1M, q2M, q3M, q4M, q5M, q6M, q1B, q2B, q3B, q4B, q5B, q6B) \
{                                           \
  if(bc==1){                                \
    *(q1B) =  q1M;                          \
    *(q2B) = -q2M;                          \
    *(q3B) = -q3M;                          \
    *(q4B) =  q4M;                          \
    *(q5B) =  q5M;                          \
    *(q6B) =  q6M;                          \
  } else if(bc==4||bc==5){                  \
    *(q1B) = q1M;                           \
    *(q2B) = q2M-2.f*(nx*q2M+ny*q3M)*nx;    \
    *(q3B) = q3M-2.f*(nx*q2M+ny*q3M)*ny;    \
    *(q4B) =  q4M;                          \
    *(q5B) =  q5M;                          \
    *(q6B) =  q6M;                          \
  }                                         \
  else {                                    \
    *(q1B) = 2.f*p_q1bar - q1M;             \
    *(q2B) = 2.f*ramp*p_q2bar - q2M;        \
    *(q3B) = 2.f*ramp*p_q3bar - q3M;        \
    *(q4B) = 2.f*ramp*ramp*p_q4bar - q4M;   \
    *(q5B) = 2.f*ramp*ramp*p_q5bar - q5M;   \
    *(q6B) = 2.f*ramp*ramp*p_q6bar - q6M;   \
  }                                         \
}

// Boundary conditions: Did not check yet
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define boundaryConditions2D(bc, t, x, y, nx, ny, ramp, q1M, q2M, q3M, q4M, q5M, q6M, q1B, q2B, q3B, q4B, q5B, q6B) \
{                                           \
  if(bc==1){                                \
    *(q1B) =  q1M;                          \
    *(q2B) = -q2M;                          \
    *(q3B) = -q3M;                          \
    *(q4B) =  q4M;                          \
    *(q5B) =  q5M;                          \
    *(q6B) =  q6M;                          \
  } else if(bc==2){                         \
    *(q1B) = 2.f*p_q1bar - q1M;             \
    *(q2B) = 2.f*ramp*p_q2bar - q2M;        \
    *(q3B) = 2.f*ramp*p_q3bar - q3M;        \
    *(q4B) = 2.f*ramp*ramp*p_q4bar - q4M;   \
    *(q5B) = 2.f*ramp*ramp*p_q5bar - q5M;   \
    *(q6B) = 2.f*ramp*ramp*p_q6bar - q6M;   \
  } else if(bc==3){                         \
    *(q1B) = 2.f*q1M - q1M;                 \
    *(q2B) = 2.f*ramp*p_q2bar - q2M;        \
    *(q3B) = 2.f*ramp*p_q3bar - q3M;        \
    *(q4B) = 2.f*ramp*ramp*p_q4bar - q4M;   \
    *(q5B) = 2.f*ramp*ramp*p_q5bar - q5M;   \
    *(q6B) = 2.f*ramp*ramp*p_q6bar - q6M;   \
  } else if(bc==4||bc==5){                  \
    *(q1B) = q1M;                           \
    *(q2B) = q2M-2.f*(nx*q2M+ny*q3M)*nx;    \
    *(q3B) = q3M-2.f*(nx*q2M+ny*q3M)*ny;    \
    *(q4B) = q4M;                           \
    *(q5B) = q5M;                           \
    *(q6B) = q6M;                           \
  }                                         \
}


// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define boundaryConditionsPML2D(bc, t, x, y, nx, ny, intfx, intfy, q1M, q2M, q3M, q4M, q5M, q6M, q1B, q2B, q3B, q4B, q5B, q6B) \
  {									\
    if(bc==1){								\
      *(q1B) =  q1M;							\
      *(q2B) = -q2M;							\
      *(q3B) = -q3M;							\
      *(q4B) =  q4M;							\
      *(q5B) =  q5M;							\
      *(q6B) =  q6M;							\
    } else if(bc==4||bc==5){						\
      *(q1B) = q1M;							\
      *(q2B) = q2M-2.f*(nx*q2M+ny*q3M)*nx;				\
      *(q3B) = q3M-2.f*(nx*q2M+ny*q3M)*ny;				\
      *(q4B) =  q4M;							\
      *(q5B) =  q5M;							\
      *(q6B) =  q6M;							\
    }									\
    else {								\
      *(q1B) = 2.f*p_q1bar - q1M;					\
      *(q2B) = 2.f*p_q1bar*intfx*p_isqrtRT - q2M;			\
      *(q3B) = 2.f*p_q1bar*intfy*p_isqrtRT - q3M;			\
      *(q4B) = 2.f*p_q1bar*intfx*intfy*p_isqrtRT*p_isqrtRT-q4M;		\
      *(q5B) = 2.f*p_q1bar*intfx*intfx*p_isqrtRT*p_isqrtRT*p_invsqrt2-q5M; \
      *(q6B) = 2.f*p_q1bar*intfy*intfy*p_isqrtRT*p_isqrtRT*p_invsqrt2-q6M; \
    }									\
  }

// Boundary conditions: Did not check yet
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define boundaryConditions2D(bc, t, x, y, nx, ny, intfx, intfy, q1M, q2M, q3M, q4M, q5M, q6M, q1B, q2B, q3B, q4B, q5B, q6B) \
  {									\
    if(bc==1){								\
      *(q1B) =  q1M;							\
      *(q2B) = -q2M;							\
      *(q3B) = -q3M;							\
      *(q4B) =  q4M;							\
      *(q5B) =  q5M;							\
      *(q6B) =  q6M;							\
    } else if(bc==2 || bc==3){						\
      *(q1B) = 2.f*p_q1bar - q1M;					\
      *(q2B) = 2.f*p_q1bar*intfx*p_isqrtRT - q2M;			\
      *(q3B) = 2.f*p_q1bar*intfy*p_isqrtRT - q3M;			\
      *(q4B) = 2.f*p_q1bar*intfx*intfy*p_isqrtRT*p_isqrtRT-q4M;		\
      *(q5B) = 2.f*p_q1bar*intfx*intfx*p_isqrtRT*p_isqrtRT*p_invsqrt2-q5M; \
      *(q6B) = 2.f*p_q1bar*intfy*intfy*p_isqrtRT*p_isqrtRT*p_invsqrt2-q6M; \
    } else if(bc==4||bc==5){						\
      *(q1B) = q1M;							\
      *(q2B) = q2M-2.f*(nx*q2M+ny*q3M)*nx;				\
      *(q3B) = q3M-2.f*(nx*q2M+ny*q3M)*ny;				\
      *(q4B) = q4M;							\
      *(q5B) = q5M;							\
      *(q6B) = q6M;							\
    }									\
  }

#define p_RDIMX 32
#define p_RDIMY 8

// used a macro since I am not sure what happens with exclusive variables in OpenMP mode
#define twoPhaseReduction(r_ip, s_ip, s_res, g_ip)		\
								\
  barrier(localMemFence);					\
								\
  for(int ty=0;ty<p_RDIMY;++ty;inner1){				\
    for(int tx=0;tx<p_RDIMX;++tx;inner0){			\
      s_ip[ty][tx] = r_ip;					\
      if(tx>  1*p_RDIMX/2  -1) s_ip[ty][tx] += s_ip[ty][tx-p_RDIMX/2];	\
      if(tx>  3*p_RDIMX/4  -1) s_ip[ty][tx] += s_ip[ty][tx-p_RDIMX/4];	\
      if(tx>  7*p_RDIMX/8  -1) s_ip[ty][tx] += s_ip[ty][tx-p_RDIMX/8];	\
      if(tx> 15*p_RDIMX/16 -1) s_ip[ty][tx] += s_ip[ty][tx-p_RDIMX/16];	\
      if(tx> 31*p_RDIMX/32 -1) s_ip[ty][tx] += s_ip[ty][tx-p_RDIMX/32];	\
      if(tx==p_RDIMX-1) s_res[ty] = s_ip[ty][tx];			\
    }								\
  }								\
								\
  barrier(localMemFence);					\
							\
  for(int ty=0;ty<p_RDIMY;++ty;inner1){				\
    for(int tx=0;tx<p_RDIMX;++tx;inner0){			\
      if(tx==0){						\
      	if(ty>  1*p_RDIMY/2  -1) s_res[ty] += s_res[ty-p_RDIMY/2];	\
      	if(ty>  3*p_RDIMY/4  -1) s_res[ty] += s_res[ty-p_RDIMY/4];	\
      	if(ty>  7*p_RDIMY/8  -1) s_res[ty] += s_res[ty-p_RDIMY/8];	\
      	if(ty==p_RDIMY-1) { \
          dfloat val = s_res[p_RDIMY-1]; \
          atomicAdd(g_ip, val);     \
        } \
      }								\
    }								\
  }


#define workgroup_scan_upsweep(sPrefix)			\
							\
  for(int offset= 1, d = 2*bdim>>1;  d>0; d>>=1){	\
							\
    barrier(localMemFence);				\
							\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t < d){					\
	const int ai = offset*(2*t+1)-1;		\
	const int bi = offset*(2*t+2)-1;		\
							\
	sPrefix[bi] += sPrefix[ai];			\
      }							\
    }							\
							\
    offset *= 2;					\
  }						

#define workgroup_scan_downsweep(sPrefix)	\
						\
  barrier(localMemFence);			\
						\
  for(int t=0;t<innerDim0;++t;inner0){	\
    if(t == bdim-1){				\
      sPrefix[2*bdim-1] = 0;			\
    }						\
  }						\
{    						\
  int offset = 2*bdim;				\
						\
  for(int d=1; d<2*bdim; d *= 2){		\
						\
    offset >>= 1;				\
						\
    barrier(localMemFence);			\
						\
    for(int t=0;t<innerDim0;++t;inner0){	\
      if(t < d){				\
	const int ai = offset*(2*t+1)-1;	\
	const int bi = offset*(2*t+2)-1;	\
						\
	int t = sPrefix[ai];			\
	sPrefix[ai] = sPrefix[bi];		\
	sPrefix[bi] += t;			\
      }						\
    }						\
  }						\
  }




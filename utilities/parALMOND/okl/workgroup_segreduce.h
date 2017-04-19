#define workgroup_segreduce(left, s_key, s_val)		\
  							\
  barrier(localMemFence);			\
						\
  if(bdim >= 1){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if(t >= 1 && s_key[t] == s_key[t-1]){		\
	left = s_val[t-1];				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 1){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t >= 1 && s_key[t] == s_key[t-1]){		\
	s_val[t] += left;				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 2){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if( t >= 2 && s_key[t] == s_key[t-2]){		\
	left = s_val[t-2];				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 2){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if( t >= 2 && s_key[t] == s_key[t-2]){		\
	s_val[t] += left;				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 4){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if(t >= 4 && s_key[t] == s_key[t-4]){		\
	left = s_val[t-4];				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 4){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t >= 4 && s_key[t] == s_key[t-4]){		\
	s_val[t] += left;				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 8){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if(t >= 8 && s_key[t] == s_key[t-8]){		\
	left = s_val[t-8];				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 8){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t >= 8 && s_key[t] == s_key[t-8]){		\
	s_val[t] += left;				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 16){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if(t >= 16 && s_key[t] == s_key[t-16]){	\
	left = s_val[t-16];				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 16){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t >= 16 && s_key[t] == s_key[t-16]){	\
	s_val[t] += left;				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 32){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if(t >= 32 && s_key[t] == s_key[t-32]){	\
	left = s_val[t-32];				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 32){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t >= 32 && s_key[t] == s_key[t-32]){	\
	s_val[t] += left;				\
      }							\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 64){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if(t >= 64 && s_key[t] == s_key[t-64])		\
	left = s_val[t-64];				\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 64){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t >= 64 && s_key[t] == s_key[t-64])		\
	s_val[t] += left;				\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 128){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if(t >= 128 && s_key[t] == s_key[t-128])	\
	left = s_val[t-128];				\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 128){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t >= 128 && s_key[t] == s_key[t-128])	\
	s_val[t] += left;				\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 256){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if(t >= 256 && s_key[t] == s_key[t-256])	\
	left = s_val[t-256];				\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 256){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t >= 256 && s_key[t] == s_key[t-256])	\
	s_val[t] += left;				\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 512){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if(t >= 512 && s_key[t] == s_key[t-512])	\
	left = s_val[t-512];				\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 512){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t >= 512 && s_key[t] == s_key[t-512])	\
	s_val[t] += left;				\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 1024){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      left = 0.;					\
      if(t >= 1024 && s_key[t] == s_key[t-1024])	\
	left = s_val[t-1024];			\
    }							\
  }							\
							\
  barrier(localMemFence);				\
							\
  if(bdim >= 1024){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t >= 1024 && s_key[t] == s_key[t-1024])	\
	s_val[t] += left;				\
    }							\
  }							\
  barrier(localMemFence);				

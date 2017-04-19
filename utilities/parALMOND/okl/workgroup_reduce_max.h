#define mymax(a,b) ((a>b)?a:b)

// assumes power of two - and lots of compiler defines
#define  workgroup_reduce_max(s_a)			\
							\
  barrier(localMemFence);				\
							\
  if(bdim>=1024){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t<512) s_a[t] = mymax(s_a[t],s_a[t+512]);	\
    }							\
  }							\
  barrier(localMemFence);				\
  if(bdim>=512){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t<256) s_a[t] = mymax(s_a[t],s_a[t+256]);	\
    }							\
  }							\
  barrier(localMemFence);				\
							\
  if(bdim>=256){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t<128) s_a[t] = mymax(s_a[t],s_a[t+128]);	\
    }							\
  }							\
							\
  barrier(localMemFence);				\
  if(bdim>=128){					\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t<64)  s_a[t] = mymax(s_a[t],s_a[t+64]);	\
    }							\
  }							\
							\
  barrier(localMemFence);				\
  if(bdim>=64){						\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t<32) s_a[t] = mymax(s_a[t],s_a[t+32]);	\
    }							\
  }							\
  barrier(localMemFence);				\
							\
  if(bdim>=32){						\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t<16) s_a[t] = mymax(s_a[t],s_a[t+16]);	\
    }							\
  }							\
  barrier(localMemFence);				\
							\
  if(bdim>=16){						\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t<8)  s_a[t] = mymax(s_a[t],s_a[t+8]);	\
    }							\
  }							\
  barrier(localMemFence);				\
							\
  if(bdim>=8){						\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t<4)  s_a[t] = mymax(s_a[t],s_a[t+4]);	\
    }							\
  }							\
  barrier(localMemFence);				\
							\
  if(bdim>=4){						\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t<2)  s_a[t] = mymax(s_a[t],s_a[t+2]);	\
    }							\
  }							\
  barrier(localMemFence);				\
							\
  if(bdim>=2){						\
    for(int t=0;t<innerDim0;++t;inner0){		\
      if(t==0) s_a[t] = mymax(s_a[t],s_a[1]);	\
    }							\
  }							\
  barrier(localMemFence);				\
  									

#define  workgroup_reduce(s_a)			\
  						\
  barrier(localMemFence);			\
  for(int it=1024;it>=2;it/=2){			\
    for(int t=0;t<innerDim0;++t;inner0){	\
      if(bdim>=it){				\
	if(t<(it/2)) s_a[t] += s_a[t+(it/2)];	\
      }						\
    }						\
    barrier(localMemFence);			\
  }						

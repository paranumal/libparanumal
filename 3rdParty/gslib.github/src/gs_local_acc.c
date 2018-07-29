/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/


/*------------------------------------------------------------------------------
  The basic gather kernel
------------------------------------------------------------------------------*/
#define DEFINE_GATHER_ACC(T,OP) \
static void gather_##T##_##OP##_acc( \
  T* out, const T* in, const unsigned in_stride,           \
  const uint *restrict map)                                                  \
{                                                                            \
  uint i,j;                                                                  \
  while((i=*map++)!=-(uint)1) {                                              \
    T t=out[i];                                                              \
    j=*map++; do GS_DO_##OP(t,in[j*in_stride]); while((j=*map++)!=-(uint)1); \
    out[i]=t;                                                                \
  }                                                                          \
}

/*------------------------------------------------------------------------------
  The basic scatter kernel
------------------------------------------------------------------------------*/
#define DEFINE_SCATTER_ACC(T) \
static void scatter_##T##_acc( \
  T* out, const unsigned out_stride,                      \
  const T* in, const unsigned in_stride,                  \
  const uint *restrict map)                                        \
{                                                                  \
  uint i,j;                                                        \
  while((i=*map++)!=-(uint)1) {                                    \
    T t=in[i*in_stride];                                           \
    j=*map++; do out[j*out_stride]=t; while((j=*map++)!=-(uint)1); \
  }                                                                \
}

/*------------------------------------------------------------------------------
  The basic initialization kernel
------------------------------------------------------------------------------*/
#define DEFINE_INIT_ACC(T) \
static void init_##T##_acc(T* out, const uint* restrict map, gs_op op) \
{                                                       \
  uint i; const T e = gs_identity_##T[op];              \
  while((i=*map++)!=-(uint)1) out[i]=e;                 \
}


#define DEFINE_PROCS(T) \
  GS_FOR_EACH_OP(T,DEFINE_GATHER_ACC) \
  DEFINE_SCATTER_ACC(T) \
  DEFINE_INIT_ACC(T)

GS_FOR_EACH_DOMAIN(DEFINE_PROCS)


/*------------------------------------------------------------------------------
  Multiple array kernels
------------------------------------------------------------------------------*/
void gs_gather_many_acc(void** out, const void** in, const unsigned vn,
			const uint* restrict map, gs_dom dom, gs_op op)
{
  uint k;
#define WITH_OP(T,OP) for(k=0;k<vn;++k) gather_##T##_##OP##_acc(out[k],in[k],1,map)
#define WITH_DOMAIN(T) SWITCH_OP(T,op)
  SWITCH_DOMAIN(dom);
#undef  WITH_DOMAIN
#undef  WITH_OP
}

void gs_scatter_many_acc(void** out, const void** in, const unsigned vn,
			 const uint* restrict map, gs_dom dom)
{
  uint k;
#define WITH_DOMAIN(T) for(k=0;k<vn;++k) scatter_##T##_acc(out[k],1,in[k],1,map)
  SWITCH_DOMAIN(dom);
#undef  WITH_DOMAIN
}

void gs_init_many_acc(void** out, const unsigned vn, const uint* restrict map,
		      gs_dom dom, gs_op op)
{
  uint k;
#define WITH_DOMAIN(T) for(k=0;k<vn;++k) init_##T##_acc(out[k],map,op)
  SWITCH_DOMAIN(dom);
#undef  WITH_DOMAIN
}


/*------------------------------------------------------------------------------
  Gather from strided array -> multiple arrays
  Scatter from multiple arrays -> strided array,
  Scatter from strided array -> multiple arrays,
------------------------------------------------------------------------------*/
void gs_gather_vec_to_many_acc(void *out, const void *in, const unsigned vn,
                           const uint *map, gs_dom dom, gs_op op)
{
  unsigned i; const unsigned unit_size = gs_dom_size[dom];
  typedef void *ptr_to_void;
  const ptr_to_void *p = out; const char *q = in;
#define WITH_OP(T,OP) \
  for(i=vn;i;--i) gather_##T##_##OP##_acc(*p++,(const T*)q,vn,map), q+=unit_size
#define WITH_DOMAIN(T) SWITCH_OP(T,op)
  SWITCH_DOMAIN(dom);
#undef  WITH_DOMAIN
#undef  WITH_OP
}

void gs_scatter_many_to_vec_acc(void *out, const void *in, const unsigned vn,
                            const uint *map, gs_dom dom)
{
  unsigned i; const unsigned unit_size = gs_dom_size[dom];
  typedef const void *ptr_to_const_void;
  char *p = out; const ptr_to_const_void *q = in;
#define WITH_DOMAIN(T) \
  for(i=vn;i;--i) scatter_##T##_acc((T*)p,vn,*q++,1,map), p+=unit_size
  SWITCH_DOMAIN(dom);
#undef  WITH_DOMAIN
}


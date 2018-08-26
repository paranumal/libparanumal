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


// used a macro since I am not sure what happens with @exclusive variables in OpenMP mode
#define twoPhaseReduction(r_ip, s_ip, s_res, g_ip)			\
									\
  @barrier("local");						\
									\
  for(int ty=0;ty<p_RDIMY;++ty;@inner(1)){					\
    for(int tx=0;tx<p_RDIMX;++tx;@inner(0)){				\
      s_ip[ty][tx] = r_ip;						\
      if(tx>=  1*p_RDIMX/2) s_ip[ty][tx] += s_ip[ty][tx-p_RDIMX/2];	\
      if(tx>=  3*p_RDIMX/4) s_ip[ty][tx] += s_ip[ty][tx-p_RDIMX/4];	\
      if(tx>=  7*p_RDIMX/8) s_ip[ty][tx] += s_ip[ty][tx-p_RDIMX/8];	\
      if(tx>= 15*p_RDIMX/16) s_ip[ty][tx] += s_ip[ty][tx-p_RDIMX/16];	\
      if(tx>= 31*p_RDIMX/32) s_ip[ty][tx] += s_ip[ty][tx-p_RDIMX/32];	\
      if(tx==(p_RDIMX-1)) s_res[ty] = s_ip[ty][tx];			\
    }									\
  }									\
									\
  @barrier("local");						\
									\
  for(int ty=0;ty<p_RDIMY;++ty;@inner(1)){					\
    for(int tx=0;tx<p_RDIMX;++tx;@inner(0)){				\
      if(ty==0 && tx<p_RDIMY){						\
      	if(tx >= 1*p_RDIMY/2) s_res[tx] += s_res[tx-p_RDIMY/2];		\
      	if(tx >= 3*p_RDIMY/4) s_res[tx] += s_res[tx-p_RDIMY/4];		\
      	if(tx >= 7*p_RDIMY/8) s_res[tx] += s_res[tx-p_RDIMY/8];		\
      	if(tx==(p_RDIMY-1)) {						\
          g_ip = s_res[p_RDIMY-1];				\
        }								\
      }									\
    }									\
  }

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

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

#define USE_ODD_EVEN 1

static const int p_Nq = comp_Nq;
static const int p_cubNq = comp_cubNq;

static const int p_halfNq = ((comp_Nq+1)/2);
static const int p_halfCubNq = ((comp_cubNq+1)/2);

static const int p_padNq = (p_Nq%4) ? 0:1;
static const int p_padCubNq = (p_cubNq%4) ? 0:1;

static const int p_MAX_ITERATIONS= comp_MAX_ITERATIONS;

#define p_Nq2 (p_Nq*p_Nq)
#define p_Np  (p_Nq*p_Nq*p_Nq)


#define p_cubNq2 (p_cubNq*p_cubNq)
#define p_cubNp  (p_cubNq*p_cubNq*p_cubNq)


#define p_Nvgeo 1
#define p_JWID 0

#define p_Nwarps ((p_Nq2+32-1)/32)

#define dlong int
#define hlong dlong
#define dfloat double

#if USE_ODD_EVEN==1
__constant__ dfloat const_I[p_halfCubNq*p_Nq];
__constant__ dfloat const_IT[p_cubNq*p_halfNq];
#else
__constant__ dfloat const_I[p_cubNq*p_Nq];
__constant__ dfloat const_IT[p_cubNq*p_Nq];
#endif


// do C op first
__device__ void advectionMassMatrixMultiply(const dlong element,
					    dfloat * __restrict__ r_p,
					    dfloat s_WJ[p_cubNq][p_cubNq][p_cubNq],
					    dfloat s_Ap[p_cubNq][p_cubNq][p_cubNq+p_padCubNq],
					    dfloat * __restrict__ r_Ap){

  dfloat r_tmp[p_cubNq];
  
  const int t = threadIdx.x;

  __syncthreads();

  // transform in 'c'
  {
    const int a = t%p_Nq;
    const int b = t/p_Nq;
    
    if(b<p_Nq){

#if USE_2D_CONSTANT_MEMORY==0
      const dfloat * __restrict__ cI = const_I;
#endif

#pragma unroll p_cubNq
      for(int k=0;k<p_cubNq;++k){
	dfloat res = 0;

#pragma unroll p_Nq
	for(int c=0;c<p_Nq;++c){
#if USE_2D_CONSTANT_MEMORY==1
	  res += const_I[k][c]*r_p[c];
#else
	  res += *(cI++)*r_p[c];
#endif
	}
	s_Ap[k][b][a] = res;
      }
    }
  }
  
  __syncthreads();

  // transform in 'b'
  {
    const int a = t%p_Nq;
    const int k = t/p_Nq;
    
    if(k<p_cubNq){
      for(int b=0;b<p_Nq;++b){
	r_tmp[b] = s_Ap[k][b][a];
      }

#if USE_2D_CONSTANT_MEMORY==0
      const dfloat * __restrict__ cI = const_I;
#endif
      
#pragma unroll p_cubNq
      for(int j=0;j<p_cubNq;++j){
	dfloat res = 0;
	
#pragma unroll p_Nq
	for(int b=0;b<p_Nq;++b){
#if USE_2D_CONSTANT_MEMORY==1
	  res += const_I[j][b]*r_tmp[b];
#else
	  res += *(cI++)*r_tmp[b];
#endif
	}
	s_Ap[k][j][a] = res;
      }
    }
  }

  __syncthreads();

  // transform in 'a'
  {
    const int j = t%p_cubNq;
    const int k = t/p_cubNq;

    for(int a=0;a<p_Nq;++a){
      r_tmp[a] = s_Ap[k][j][a];
    }

#if USE_2D_CONSTANT_MEMORY==0
    const dfloat * __restrict__ cI = const_I;
#endif
	  
#pragma unroll p_cubNq
    for(int i=0;i<p_cubNq;++i){
      dfloat res = 0;

#pragma unroll p_Nq
      for(int a=0;a<p_Nq;++a){
#if USE_2D_CONSTANT_MEMORY==1
	res += const_I[i][a]*r_tmp[a];
#else
	//	res += const_I[cnt++]*r_tmp[a];
	res += *(cI++)*r_tmp[a];
#endif
      }
    
      //      const dlong gid = element*p_cubNp*p_Nvgeo + t + k*p_cubNq*p_cubNq + p_JWID*p_cubNp;
      //      const dfloat WJ = cubvgeo[gid];
      const dfloat WJ = s_WJ[k][j][i];
      
      s_Ap[k][j][i] = WJ*res;
    }
  }

  __syncthreads();

  // test in 'a'
  {
    const int j = t%p_cubNq;
    const int k = t/p_cubNq;
    
    for(int i=0;i<p_cubNq;++i){
      r_tmp[i] = s_Ap[k][j][i];
    }

    const dfloat * __restrict__ cIT = const_IT;

#pragma unroll p_Nq
    for(int a=0;a<p_Nq;++a){
      dfloat res = 0;

#pragma unroll p_cubNq
      for(int i=0;i<p_cubNq;++i){
	res += *(cIT++)*r_tmp[i];
      }

      s_Ap[k][j][a] = res;
    }
  }

  __syncthreads();

  // test in 'b'
  {
    const int a = t%p_Nq;
    const int k = t/p_Nq;

    const dfloat * __restrict__ cIT = const_IT;
    
    if(k<p_cubNq){
      for(int j=0;j<p_cubNq;++j){
	r_tmp[j] = s_Ap[k][j][a];
      }

#pragma unroll p_Nq
      for(int b=0;b<p_Nq;++b){
	dfloat res = 0;

#pragma unroll p_cubNq
	for(int j=0;j<p_cubNq;++j){
	  res += *(cIT++)*r_tmp[j];
	}

	s_Ap[k][b][a] = res;
      }
    }
  }
  
  __syncthreads();

  // test in 'c'
  {
    const int a = t%p_Nq;
    const int b = t/p_Nq;
    
    if(t<p_Nq2){
      for(int k=0;k<p_cubNq;++k){
	r_tmp[k] = s_Ap[k][b][a];
      }
      
      const dfloat * __restrict__ cIT = const_IT;
      
#pragma unroll p_Nq
      for(int c=0;c<p_Nq;++c){
	dfloat res = 0;
	
#pragma unroll p_cubNq
	for(int k=0;k<p_cubNq;++k){
	  res += *(cIT++)*r_tmp[k];
	}
	
	r_Ap[c] = res;
      }
    }
  }
}


__forceinline__ __device__ dfloat advectionMassMatrixMultiplyOddEven(const dlong element,
								     dfloat * __restrict__ r_p,
								     dfloat s_WJ[p_cubNq][p_cubNq][p_cubNq],
								     dfloat s_Ap[p_cubNq][p_cubNq][p_cubNq+p_padCubNq],
								     dfloat * __restrict__ r_Ap){
  
  dfloat r_tmpOdd[p_halfCubNq];
  dfloat r_tmpEven[p_halfCubNq];
  
  const int t = threadIdx.x;

  // assumes barrier before s_Ap was used last
  
  // transform in 'c'
  {
    const int a = t%p_Nq;
    const int b = t/p_Nq;
    
    if(b<p_Nq){
      
#if USE_2D_CONSTANT_MEMORY==0
      const dfloat * __restrict__ cI = const_I;
#endif
      
#pragma unroll p_halfNq
      for(int c=0;c<p_halfNq;++c){
	r_tmpOdd[c]  = r_p[c] + r_p[p_Nq-1-c];
	r_tmpEven[c] = r_p[c] - r_p[p_Nq-1-c];
      }
      
#pragma unroll p_halfCubNq
      for(int k=0;k<p_halfCubNq;++k){
	dfloat resOdd = 0, resEven = 0;
      
#pragma unroll p_halfNq
	for(int c=0;c<p_halfNq;++c){
	  resOdd += *(cI++)*r_tmpOdd[c];
	  resEven += *(cI++)*r_tmpEven[c];
	}
	
	s_Ap[k][b][a]           = resOdd+resEven;
	s_Ap[p_cubNq-1-k][b][a] = resOdd-resEven;
      }
    }
  }
  
  __syncthreads();

  // transform in 'b'
  {
    const int a = t%p_Nq;
    const int k = t/p_Nq;
    
    if(k<p_cubNq){

#pragma unroll p_halfNq
      for(int b=0;b<p_halfNq;++b){
	dfloat ApOdd  = s_Ap[k][b][a];
	dfloat ApEven = s_Ap[k][p_Nq-1-b][a];
	r_tmpOdd[b]  = ApOdd + ApEven;
	r_tmpEven[b] = ApOdd - ApEven;
      }      

      const dfloat * __restrict__ cI = const_I;
      
#pragma unroll p_halfCubNq
      for(int j=0;j<p_halfCubNq;++j){
	dfloat resOdd = 0, resEven = 0;
	
#pragma unroll p_halfNq
	for(int b=0;b<p_halfNq;++b){
	  resOdd += *(cI++)*r_tmpOdd[b];
	  resEven += *(cI++)*r_tmpEven[b];
	}
	
	s_Ap[k][j][a]               = resOdd+resEven;
	s_Ap[k][p_halfCubNq-1-j][a] = resOdd-resEven;
      }
    }
  }

  __syncthreads();

  // transform in 'a'
  {
    const int j = t%p_cubNq;
    const int k = t/p_cubNq;

#pragma unroll p_halfNq
    for(int a=0;a<p_halfNq;++a){
      dfloat ApOdd  = s_Ap[k][j][a];
      dfloat ApEven = s_Ap[k][j][p_Nq-1-a];
      r_tmpOdd[a]  = ApOdd + ApEven;
      r_tmpEven[a] = ApOdd - ApEven;
    }
    
    const dfloat * __restrict__ cI = const_I;
	  
#pragma unroll p_halfCubNq
    for(int i=0;i<p_halfCubNq;++i){
      dfloat resOdd = 0, resEven = 0;

#pragma unroll p_halfNq
      for(int a=0;a<p_halfNq;++a){
	resOdd  += *(cI++)*r_tmpOdd[a];
	resEven += *(cI++)*r_tmpEven[a];
      }
      
      r_Ap[i]           =           s_WJ[k][j][i]*(resOdd + resEven);
      r_Ap[p_cubNq-1-i] = s_WJ[k][j][p_cubNq-1-i]*(resOdd - resEven);
    }
  }

  __syncthreads();

  // test in 'a'
  {
    const int j = t%p_cubNq;
    const int k = t/p_cubNq;

    for(int i=0;i<p_halfCubNq;++i){
      dfloat ApOdd  = r_Ap[i];
      dfloat ApEven = r_Ap[p_cubNq-1-i];
      r_tmpOdd[i]  = ApOdd + ApEven;
      r_tmpEven[i] = ApOdd - ApEven;
    }
    
    const dfloat * __restrict__ cIT = const_IT;

#pragma unroll p_halfNq
    for(int a=0;a<p_halfNq;++a){
      dfloat resOdd = 0, resEven = 0;
      
#pragma unroll p_halfCubNq
      for(int i=0;i<p_halfCubNq;++i){
	resOdd  += *(cIT++)*r_tmpOdd[i];
	resEven += *(cIT++)*r_tmpEven[i];
      }
      
      s_Ap[k][j][a]        = resOdd+resEven;
      s_Ap[k][j][p_Nq-1-a] = resOdd-resEven;
    }
  }

  __syncthreads();

  // test in 'b'
  {
    const int a = t%p_Nq;
    const int k = t/p_Nq;

    const dfloat * __restrict__ cIT = const_IT;
    
    if(k<p_cubNq){

      for(int j=0;j<p_halfCubNq;++j){
	dfloat ApOdd  = s_Ap[k][j][a];
	dfloat ApEven = s_Ap[k][p_cubNq-1-j][a];
	r_tmpOdd[j]  = ApOdd + ApEven;
	r_tmpEven[j] = ApOdd - ApEven;
      }

#pragma unroll p_halfNq
      for(int b=0;b<p_halfNq;++b){
	dfloat resOdd = 0, resEven = 0;

#pragma unroll p_halfCubNq
	for(int j=0;j<p_halfCubNq;++j){
	  resOdd  += *(cIT++)*r_tmpOdd[j];
	  resEven += *(cIT++)*r_tmpEven[j];
	}
	
	s_Ap[k][b][a]        = resOdd+resEven;
	s_Ap[k][p_Nq-1-b][a] = resOdd-resEven;
      }
    }
  }
  
  __syncthreads();

  dfloat pAp_ab = 0;
  
  // test in 'c'
  {
    const int a = t%p_Nq;
    const int b = t/p_Nq;
    
    if(t<p_Nq2){

      for(int k=0;k<p_halfCubNq;++k){
	dfloat ApOdd  = s_Ap[k][b][a];
	dfloat ApEven = s_Ap[p_cubNq-1-k][b][a];
	r_tmpOdd[k]  = ApOdd + ApEven;
	r_tmpEven[k] = ApOdd - ApEven;
      }

      const dfloat * __restrict__ cIT = const_IT;
      
#pragma unroll p_halfNq
      for(int c=0;c<p_halfNq;++c){
	dfloat resOdd = 0, resEven = 0;

#pragma unroll p_halfCubNq
	for(int k=0;k<p_halfCubNq;++k){
	  resOdd  += *(cIT++)*r_tmpOdd[k];
	  resEven += *(cIT++)*r_tmpEven[k];
	}

	r_Ap[c]        = resOdd + resEven;
	r_Ap[p_Nq-1-c] = resOdd - resEven;

	pAp_ab += r_Ap[c]*r_p[c];
	pAp_ab += r_Ap[p_Nq-1-c]*r_p[p_Nq-1-c];
      }
    }
  }

  return pAp_ab;
}


__forceinline__ __device__ dfloat dotProduct(const dfloat a, volatile dfloat *s_a, volatile dfloat *s_warpa){
  
  const int t = threadIdx.x;
  const int w = t>>5; // divide by 32
  const int n = t%32;
  
  if(t<p_Nq2){
    s_a[t] = a;
    if(n<16 && t+16<p_Nq2) s_a[t] += s_a[t+16];
    if(n< 8 && t+ 8<p_Nq2) s_a[t] += s_a[t+ 8];
    if(n< 4 && t+ 4<p_Nq2) s_a[t] += s_a[t+ 4];
    if(n< 2 && t+ 2<p_Nq2) s_a[t] += s_a[t+ 2];
    if(n==0) s_warpa[w] = s_a[t] + s_a[t+1]; // dangerous (p_Nq2 == 33)
  }

  __syncthreads();

  if(w==0 && t<p_Nwarps){
#if (p_Nwarps)>16
    if(t<16 && t+16<p_Nwarps) s_warpa[t] += s_warpa[t+16];
#endif
#if (p_Nwarps)>8    
    if(t< 8 && t+ 8<p_Nwarps) s_warpa[t] += s_warpa[t+ 8];
#endif
#if (p_Nwarps)>4
    if(t< 4 && t+ 4<p_Nwarps) s_warpa[t] += s_warpa[t+ 4];
#endif
    if(t< 2 && t+ 2<p_Nwarps) s_warpa[t] += s_warpa[t+ 2];
    if(t==0 && t+1 <p_Nwarps) s_warpa[0] += s_warpa[1]; 
  }

  __syncthreads();

  dfloat res = s_warpa[0];

  return res;
}

// M*q = rhsq
// blocks: Nelements
// threads: cubNq x cubNq
__global__ void advectionInvertMassMatrixKernel(const dlong Nelements,
						const dlong  * __restrict__ elementIds,
						const dfloat dt,
						const dfloat rka,
						const dfloat rkb,
						const dfloat                tol,
						const dfloat * __restrict__ cubvgeo,
						const dfloat * __restrict__ cubI,
						const dfloat * __restrict__ precon,
					        dfloat * __restrict__ resq,
						const dfloat * __restrict__ rhsq,
						const dfloat * __restrict__ q,
					        dfloat * __restrict__ qnew){

  __shared__ dfloat s_tmp1[p_cubNq][p_cubNq][p_cubNq+p_padCubNq];

  volatile __shared__ dfloat s_tmp2[p_Nq2];
  volatile __shared__ dfloat s_tmpWarp[p_Nwarps];

  __shared__ dfloat s_WJ[p_cubNq][p_cubNq][p_cubNq];
  
  //  __shared__ dfloat s_precon[p_Nq][p_Nq][p_Nq];
  dfloat r_precon[p_Nq];
  
  dfloat r_r[p_Nq], r_z[p_Nq], r_x[p_Nq];
  dfloat r_p[p_Nq], r_Ap[p_cubNq];
  
  const dlong e = blockIdx.x;

  const dlong element = elementIds[e];

  const int t = threadIdx.x;
  
  const int a = t%p_Nq;
  const int b = t/p_Nq;

  dfloat rdotz_ab = 0;

  for(int i=0;i<p_cubNq;++i){
    const dlong gid = element*p_cubNp*p_Nvgeo + t + i*p_cubNq*p_cubNq + p_JWID*p_cubNp; //  (i slowest, k middle, j fastest)
    s_WJ[i][0][t] = cubvgeo[gid];
  }

  int n = t;
  while(n<p_Np){
    const int ta = n%p_Nq;
    const int tb = (n/p_Nq)%p_Nq;
    const int tc = (n/p_Nq2);

    const dlong id = n + element*p_Np;
    
    s_tmp1[tc][tb][ta]  = rhsq[id];
    n+= p_cubNq2;
  }

  __syncthreads();
  
  if(t<p_Nq2){    
    for(int c=0;c<p_Nq;++c){

      const dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;

      const dfloat prec = precon[id];
      r_precon[c] = prec;
      
      //      r_r[c] = rhsq[id]; // not great when cubNq>Nq
      r_r[c] = s_tmp1[c][b][a];
      r_z[c] = prec*r_r[c];
      r_p[c] = r_z[c];

      rdotz_ab += r_r[c]*r_z[c];
    }
  }

  dfloat rdotz = dotProduct(rdotz_ab, s_tmp2, s_tmpWarp);

  for(int it=0;it<p_MAX_ITERATIONS;++it){

    __syncthreads();
    
#if USE_ODD_EVEN==0
    advectionMassMatrixMultiply(element, r_p, s_WJ, s_tmp1, r_Ap);
#else
    dfloat pAp_ab = advectionMassMatrixMultiplyOddEven(element, r_p, s_WJ, s_tmp1, r_Ap);
#endif
    
    const dfloat pAp = dotProduct(pAp_ab, s_tmp2, s_tmpWarp);
    
    const dfloat alpha = rdotz/pAp;

    rdotz_ab = 0;
    
    if(t<p_Nq2){
#pragma unroll p_Nq
      for(int c=0;c<p_Nq;++c){
	r_x[c] += alpha*r_p[c];
	r_r[c] -= alpha*r_Ap[c];

	r_z[c] = r_precon[c]*r_r[c];
	
	rdotz_ab += r_r[c]*r_z[c];
      }
    }

    __syncthreads();
    
    // r.z
    const dfloat rdotz_new = dotProduct(rdotz_ab, s_tmp2, s_tmpWarp);
    
    const dfloat beta = rdotz_new/rdotz;

    rdotz = rdotz_new;

    if(t<p_Nq2){
#pragma unroll p_Nq
      for(int c=0;c<p_Nq;++c){
	r_p[c] = r_z[c] + beta*r_p[c];
      }
    }
  } // end iterations

#if 0
  if(t<p_Nq2){
#pragma unroll p_Nq
    for(int c=0;c<p_Nq;++c){
      dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;
#if 1
      dfloat r_qcba = q[id];
      
      dfloat r_resq = resq[id];
      r_resq = rka*r_resq + dt*r_x[c];
      r_qcba += rkb*r_resq;
      
      resq[id] = r_resq;
      
      qnew[id] = r_qcba;
#else
      qnew[id] = r_x[c];
#endif
      
    }
  }
#else
  
  if(t<p_Nq2){
#pragma unroll p_Nq
    for(int c=0;c<p_Nq;++c){      
      s_tmp1[c][b][a] = r_x[c];
    }
  }

  __syncthreads();

  n = t;
  while(n<p_Np){
    const int ta = n%p_Nq;
    const int tb = (n/p_Nq)%p_Nq;
    const int tc = (n/p_Nq2);
    
    const dlong id = n + element*p_Np;
    
    dfloat r_qcba = q[id];
    
    dfloat r_resq = resq[id];
    r_resq = rka*r_resq + dt*s_tmp1[tc][tb][ta];
    r_qcba += rkb*r_resq;
    
    resq[id] = r_resq;
    
    qnew[id] = r_qcba;

    n+= p_cubNq2;
  }

  
#endif

  
}
	

__global__ void advectionMassMatrixMultiplyKernel(const dlong Nelements,
						  const dlong  * __restrict__ elementIds,
						  const dfloat * __restrict__ cubvgeo,
						  const dfloat * __restrict__ cubI,
						  const dfloat * __restrict__ q,
						  dfloat * __restrict__ qnew){
  
  __shared__ dfloat s_tmp1[p_cubNq][p_cubNq][p_cubNq+p_padCubNq];

  volatile __shared__ dfloat s_tmp2[p_Nq2];
  volatile __shared__ dfloat s_tmpWarp[p_Nwarps];

  dfloat r_q[p_Nq], r_Aq[p_Nq];
  
  __shared__ dfloat s_WJ[p_cubNq][p_cubNq][p_cubNq];

  const dlong e = blockIdx.x;

  const dlong element = elementIds[e];

  const unsigned int t = threadIdx.x;

  const unsigned int a = t%p_Nq;
  const unsigned int b = t/p_Nq;

  for(int i=0;i<p_cubNq;++i){
    const dlong gid = element*p_cubNp*p_Nvgeo + t + i*p_cubNq*p_cubNq + p_JWID*p_cubNp; //  (i slowest, k middle, j fastest)
    s_WJ[i][0][t] = cubvgeo[gid];
  }
  
  if(t<p_Nq2){    
    for(int c=0;c<p_Nq;++c){

      dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;

      r_q[c] = q[id];
    }
  }

#if USE_ODD_EVEN==0
  advectionMassMatrixMultiply(element, r_q, s_WJ, s_tmp1, r_Aq);
#else
  advectionMassMatrixMultiplyOddEven(element, r_q, s_WJ, s_tmp1, r_Aq);
#endif

  if(t<p_Nq2){
#pragma unroll p_Nq
    for(int c=0;c<p_Nq;++c){
      dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;
      qnew[id] = r_Aq[c];
    }
  }
}
	


void dfloatRandAlloc(int N, dfloat **h_a, dfloat **c_a){

  *h_a = (dfloat*) calloc(N, sizeof(dfloat));

  for(int n=0;n<N;++n)
    h_a[0][n] = drand48();

  cudaMalloc(c_a, N*sizeof(dfloat));

  cudaMemcpy(c_a[0], h_a[0], N*sizeof(dfloat), cudaMemcpyHostToDevice);

}
				  
int main(int argc, char **argv){

  if(argc!=2){
    printf("Usage: ./advectionInvertMassMatrix Nelements\n");
    exit(-1);
  }

  dfloat tol = 1e-8;
  
  dfloat *h_cubvgeo, *c_cubvgeo;
  dfloat *h_precon,  *c_precon;
  dfloat *h_rhsq,    *c_rhsq;
  dfloat *h_resq,    *c_resq;
  dfloat *h_qnew,    *c_qnew;
  dfloat *h_q,       *c_q;
  dfloat *h_I,       *c_I;

  int    *h_elementIds, *c_elementIds;

  hlong Nelements = atoi(argv[argc-1]);
  
  int    Ntotal = Nelements*p_Np;
  int cubNtotal = Nelements*p_cubNp;

  // list of elements
  h_elementIds = (int*) calloc(Nelements, sizeof(int));
  for(int e=0;e<Nelements;++e)
    h_elementIds[e] = e;
  cudaMalloc(&c_elementIds, Nelements*sizeof(int));
  cudaMemcpy(c_elementIds, h_elementIds, Nelements*sizeof(int), cudaMemcpyHostToDevice);
  
  // float fields
  dfloatRandAlloc(cubNtotal*p_Nvgeo, &h_cubvgeo, &c_cubvgeo);
  dfloatRandAlloc(Ntotal,       &h_q, &c_q);
  dfloatRandAlloc(Ntotal,       &h_rhsq, &c_rhsq);
  dfloatRandAlloc(Ntotal,       &h_resq, &c_resq);
  dfloatRandAlloc(Ntotal,       &h_qnew, &c_qnew);
  dfloatRandAlloc(Ntotal,       &h_precon, &c_precon);
  dfloatRandAlloc(p_Nq*p_cubNq, &h_I, &c_I);

  // populate constant basis matrix
  int NconstantI  = (USE_ODD_EVEN==1) ? p_halfCubNq*p_Nq:p_cubNq*p_Nq;
  int NconstantIT = (USE_ODD_EVEN==1) ? p_cubNq*p_halfNq:p_cubNq*p_Nq;
  cudaMemcpyToSymbol(const_I,  h_I, NconstantI*sizeof(dfloat));
  cudaMemcpyToSymbol(const_IT, h_I, NconstantIT*sizeof(dfloat));

  // flush L2 ??
  dfloat *h_garbage, *c_garbage;
  int sz = 32*1024*1024; // 32MB
  dfloatRandAlloc(sz, &h_garbage, &c_garbage);
  
  // LSERK fake constants
  dfloat rka = 1.3, rkb = 3.2, dt = 0.01;

  cudaEvent_t start, end;

  cudaEventCreate(&start);
  cudaEventCreate(&end);	

  // call matrix inverse
  dim3 G(Nelements, 1, 1);
  dim3 B(p_cubNq*p_cubNq, 1, 1);

  {
    cudaEventRecord(start);
    
    int Ntests = 10;
    for(int test=0;test<Ntests;++test)
      advectionInvertMassMatrixKernel <<< G, B >>>
	(Nelements, c_elementIds, dt, rka, rkb, tol, c_cubvgeo, c_I, c_precon, c_resq, c_rhsq, c_q, c_qnew);
    
    cudaEventRecord(end);
    
    cudaEventSynchronize(end);
    
    float elapsed;
    cudaEventElapsedTime(&elapsed, start, end);
    elapsed /= 1000.;
    elapsed /= (double) Ntests;
    
    printf("%d %d %d %d %lg %lg %lg %%%% "
	   "[InvertMassMatrix: N, Nelements, Ndofs, MAX_ITERATIONS, elapsed, dofsPerSecond, nodeIterationsPerSecond]\n",
	   p_Nq-1, Nelements, p_Np*Nelements, p_MAX_ITERATIONS, elapsed, Nelements*(p_Np/elapsed),
	   Nelements*(p_Np*p_MAX_ITERATIONS/elapsed));
  }
  
  cudaDeviceSynchronize();

  {
    cudaEventRecord(start);
    
    int Ntests = 10;
    for(int test=0;test<Ntests;++test)
      advectionMassMatrixMultiplyKernel <<< G, B >>>
	(Nelements, c_elementIds, c_cubvgeo, c_I, c_q, c_qnew);
    
    cudaEventRecord(end);
    
    cudaEventSynchronize(end);
    
    float elapsed;
    cudaEventElapsedTime(&elapsed, start, end);
    elapsed /= 1000.;
    elapsed /= (double) Ntests;
    
    printf("%d %d %d %lg %lg %%%% [MassMatrixMultiply: N, Nelements, Ndofs, elapsed, dofsPerSecond]\n", p_Nq-1, Nelements, p_Np*Nelements, elapsed, Nelements*(p_Np/elapsed));
  }
  
  cudaEventDestroy(start);
  cudaEventDestroy(end);	
  
  return 0;

}

/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton

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
#include <cuda_runtime.h>

#define USE_GRAPH 1


static const int p_Nq = comp_Nq;
static const int p_cubNq = comp_cubNq;

static const int p_halfNq = ((comp_Nq+1)/2);
static const int p_halfCubNq = ((comp_cubNq+1)/2);

static const int p_padCubNq = (p_cubNq%4) ? 0:1;

#define p_Nq2 (p_Nq*p_Nq)
#define p_Np  (p_Nq*p_Nq*p_Nq)

#define p_cubNq2 (p_cubNq*p_cubNq)
#define p_cubNp  (p_cubNq*p_cubNq*p_cubNq)

#define p_Nvgeo 1
#define p_JWID 0

#define p_Nwarps ((p_Nq2+32-1)/32)

#if comp_Nq<=2
#define p_Nblock 8
#warning "using 8 elements per block"
#elif comp_Nq<=4
#define p_Nblock 2
#warning "using 2 elements per block"
#elif comp_Nq<10
#define p_Nblock 1
#warning "using 1 elements per block"
#else
#define p_Nblock 1
#warning "using 1 elements per block"
#endif

#define dlong int
#define hlong dlong
#define dfloat double

__constant__ dfloat const_oddI[p_halfCubNq][p_halfNq];
__constant__ dfloat const_evenI[p_halfCubNq][p_halfNq];


void dfloatRandAlloc(int N, dfloat **h_a, dfloat **c_a){

  *h_a = (dfloat*) calloc(N, sizeof(dfloat));

  for(int n=0;n<N;++n)
    h_a[0][n] = drand48();

  cudaMalloc(c_a, N*sizeof(dfloat));

  cudaMemcpy(c_a[0], h_a[0], N*sizeof(dfloat), cudaMemcpyHostToDevice);

}

__global__ void nothingKernel(){  }

__forceinline__ __device__
void massMatrixMultiplyDevice(const dlong Nelements,
			      const dlong element,
			      const dlong elementId,
			      const dfloat * __restrict__ cubvgeo,
			      const dfloat r_oddI[p_halfCubNq][p_halfNq],
			      const dfloat r_evenI[p_halfCubNq][p_halfNq],
			      dfloat s_Ap[p_Nblock][p_cubNq][p_cubNq][p_cubNq+p_padCubNq],
			      dfloat * __restrict__ r_Ap){
  
  dfloat r_tmpOdd[p_halfCubNq];
  dfloat r_tmpEven[p_halfCubNq];
  
  const int t = threadIdx.x;
  const int blk = threadIdx.y;
  
  // assumes barrier before s_Ap was used last
  
  // transform in 'c'
  {
    const int a = t%p_Nq;
    const int b = t/p_Nq;
    
#pragma unroll p_halfNq
    for(int c=0;c<p_halfNq;++c){
      r_tmpOdd[c]  = r_Ap[c] + r_Ap[p_Nq-1-c];
      r_tmpEven[c] = r_Ap[c] - r_Ap[p_Nq-1-c];
    }
    
#pragma unroll p_halfCubNq
    for(int k=0;k<p_halfCubNq;++k){
      dfloat resOdd = 0, resEven = 0;
      
#pragma unroll p_halfNq
      for(int c=0;c<p_halfNq;++c){
	
	resOdd += r_oddI[k][c]*r_tmpOdd[c];
	resEven += r_evenI[k][c]*r_tmpEven[c];
	
      }
      s_Ap[blk][k][b][a]           = resOdd + resEven;
      s_Ap[blk][p_cubNq-1-k][b][a] = resOdd - resEven;
    }
    
  }
  
  __syncthreads();

  // transform in 'b'
  {
    for(int n=t;n<p_Nq*p_cubNq;n+=p_Nq2){
      const int a = n%p_Nq;
      const int k = n/p_Nq;

#pragma unroll p_halfNq
      for(int b=0;b<p_halfNq;++b){
	dfloat ApOdd  = s_Ap[blk][k][b][a];
	dfloat ApEven = s_Ap[blk][k][p_Nq-1-b][a];
	r_tmpOdd[b]  = ApOdd + ApEven;
	r_tmpEven[b] = ApOdd - ApEven;
      }      
      
#pragma unroll p_halfCubNq
      for(int j=0;j<p_halfCubNq;++j){
	dfloat resOdd = 0, resEven = 0;
	
#pragma unroll p_halfNq
	for(int b=0;b<p_halfNq;++b){
	  resOdd += r_oddI[j][b]*r_tmpOdd[b];
	  resEven += r_evenI[j][b]*r_tmpEven[b];
	}
	
	s_Ap[blk][k][j][a]           = resOdd+resEven;
	s_Ap[blk][k][p_cubNq-1-j][a] = resOdd-resEven;
	
      }
    }

  }
  
  __syncthreads();

  // transform in 'a'
  {
    for(int n=t;n<p_cubNq2;n+=p_Nq2){
      const int j = n%p_cubNq;
      const int k = n/p_cubNq;
      
#pragma unroll p_halfNq
      for(int a=0;a<p_halfNq;++a){
	dfloat ApOdd  = s_Ap[blk][k][j][a];
	dfloat ApEven = s_Ap[blk][k][j][p_Nq-1-a];
	r_tmpOdd[a]  = ApOdd + ApEven;
	r_tmpEven[a] = ApOdd - ApEven;
      }
      
#pragma unroll p_halfCubNq
      for(int i=0;i<p_halfCubNq;++i){
	dfloat resOdd = 0, resEven = 0;
	
#pragma unroll p_halfNq
	for(int a=0;a<p_halfNq;++a){
	  resOdd  += r_oddI[i][a]*r_tmpOdd[a];
	  resEven += r_evenI[i][a]*r_tmpEven[a];
	}

	dlong gid1 = elementId*p_cubNp + k*p_cubNq*p_cubNq + j*p_cubNq + i;
	dlong gid2 = elementId*p_cubNp + k*p_cubNq*p_cubNq + j*p_cubNq + p_cubNq-1-i;
	
	dfloat WJ1 = (element<Nelements) ? cubvgeo[gid1]: 0;
	dfloat WJ2 = (element<Nelements) ? cubvgeo[gid2]: 0;
	
	dfloat ApOdd = WJ1*(resOdd + resEven);
	dfloat ApEven = WJ2*(resOdd - resEven);

	r_Ap[i] = ApOdd + ApEven;
	r_Ap[p_cubNq-1-i] = ApOdd - ApEven;
      }

#pragma unroll p_halfNq
      for(int a=0;a<p_halfNq;++a){
	dfloat resOdd = 0, resEven = 0;
	
#pragma unroll p_halfCubNq
	for(int i=0;i<p_halfCubNq;++i){
	  resOdd  += r_oddI[i][a]*r_Ap[i];
	  resEven += r_evenI[i][a]*r_Ap[p_cubNq-1-i];
	}
	
	s_Ap[blk][k][j][a]        = resOdd + resEven;
	s_Ap[blk][k][j][p_Nq-1-a] = resOdd - resEven;
      }
    }
  }
  
  __syncthreads();

  
  // test in 'b'
  {

    for(int n=t;n<p_Nq*p_cubNq;n+=p_Nq2){
      const int a = n%p_Nq;
      const int k = n/p_Nq;

      for(int j=0;j<p_halfCubNq;++j){
	dfloat ApOdd  = s_Ap[blk][k][j][a];
	dfloat ApEven = s_Ap[blk][k][p_cubNq-1-j][a];
	r_tmpOdd[j]  = ApOdd + ApEven;
	r_tmpEven[j] = ApOdd - ApEven;
      }

#pragma unroll p_halfNq
      for(int b=0;b<p_halfNq;++b){
	dfloat resOdd = 0, resEven = 0;
	
#pragma unroll p_halfCubNq
	for(int j=0;j<p_halfCubNq;++j){
	  resOdd  += r_oddI[j][b]*r_tmpOdd[j];
	  resEven += r_evenI[j][b]*r_tmpEven[j];
	}
	
	s_Ap[blk][k][b][a]        = resOdd + resEven;
	s_Ap[blk][k][p_Nq-1-b][a] = resOdd - resEven;
      }
    }
  }
  
  __syncthreads();

  // test in 'c'
  {
    const int a = t%p_Nq;
    const int b = t/p_Nq;

    for(int k=0;k<p_halfCubNq;++k){
      dfloat ApOdd  = s_Ap[blk][k][b][a];
      dfloat ApEven = s_Ap[blk][p_cubNq-1-k][b][a];
      r_tmpOdd[k]  = ApOdd + ApEven;
      r_tmpEven[k] = ApOdd - ApEven;
    }
    
#pragma unroll p_halfNq
    for(int c=0;c<p_halfNq;++c){
      dfloat resOdd = 0, resEven = 0;
      
#pragma unroll p_halfCubNq
      for(int k=0;k<p_halfCubNq;++k){
	resOdd  += r_oddI[k][c]*r_tmpOdd[k];
	resEven += r_evenI[k][c]*r_tmpEven[k];
      }
      
      r_Ap[c]        = resOdd + resEven;
      r_Ap[p_Nq-1-c] = resOdd - resEven;
      
    }
  }

}


__global__ void massMatrixMultiplyRegisterKernel(const dlong Nelements,
						 const dlong  * __restrict__ elementIds,
						 const dfloat * __restrict__ cubvgeo,
						 const dfloat * __restrict__ oddI,
						 const dfloat * __restrict__ evenI,
						 const dfloat * __restrict__ q,
						 dfloat * __restrict__ qnew){
  
  __shared__ dfloat s_tmp1[p_Nblock][p_cubNq][p_cubNq][p_cubNq+p_padCubNq];
  __shared__ dfloat s_oddI[p_halfNq*p_halfCubNq];
  __shared__ dfloat s_evenI[p_halfCubNq*p_halfNq];

  dfloat r_oddI[p_halfCubNq][p_halfNq];
  dfloat r_evenI[p_halfCubNq][p_halfNq];

  dfloat r_Aq[p_cubNq];

  const unsigned int t = threadIdx.x;
  const int blk = threadIdx.y;
  
  const dlong e = blockIdx.x*p_Nblock + blk;

  const dlong element = (e<Nelements) ? elementIds[e]: 0;
  
  const unsigned int a = t%p_Nq;
  const unsigned int b = t/p_Nq;

  for(int c=0;c<p_Nq;++c){
    
    dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;
    
    r_Aq[c] = q[id];
  }

  for(int n=t;n<p_halfNq*p_halfCubNq;n+=p_Nq*p_Nq){
    s_oddI[n] = oddI[n];
    s_evenI[n] = evenI[n];
    n+=p_Nq*p_Nq;
  }

  __syncthreads();
  
  for(int n=0;n<p_halfNq*p_halfCubNq;++n){
    r_oddI[0][n] = s_oddI[n];
    r_evenI[0][n] = s_evenI[n];
  }
  
  massMatrixMultiplyDevice(Nelements, e, element, cubvgeo, r_oddI, r_evenI, s_tmp1, r_Aq);
  
  if(e<Nelements){
#pragma unroll p_Nq
    for(int c=0;c<p_Nq;++c){
      dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;
      qnew[id] = r_Aq[c];
    }
  }
}



__global__ void massMatrixMultiplySharedKernel(const dlong Nelements,
						 const dlong  * __restrict__ elementIds,
						 const dfloat * __restrict__ cubvgeo,
						 const dfloat * __restrict__ oddI,
						 const dfloat * __restrict__ evenI,
						 const dfloat * __restrict__ q,
						 dfloat * __restrict__ qnew){
  
  __shared__ dfloat s_tmp1[p_Nblock][p_cubNq][p_cubNq][p_cubNq+p_padCubNq];
  __shared__ dfloat s_oddI[p_halfCubNq][p_halfNq];
  __shared__ dfloat s_evenI[p_halfCubNq][p_halfNq];

  dfloat r_Aq[p_cubNq];

  const unsigned int t = threadIdx.x;
  const int blk = threadIdx.y;
  
  const dlong e = blockIdx.x*p_Nblock + blk;

  const dlong element = (e<Nelements) ? elementIds[e]: 0;
  
  const unsigned int a = t%p_Nq;
  const unsigned int b = t/p_Nq;

  for(int c=0;c<p_Nq;++c){
    
    dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;
    
    r_Aq[c] = q[id];
  }

  for(int n=t;n<p_halfNq*p_halfCubNq;n+=p_Nq*p_Nq){
    s_oddI[0][n] = oddI[n];
    s_evenI[0][n] = evenI[n];
    n+=p_Nq*p_Nq;
  }

  __syncthreads();
  
  massMatrixMultiplyDevice(Nelements, e, element, cubvgeo, s_oddI, s_evenI, s_tmp1, r_Aq);
  
  if(e<Nelements){
#pragma unroll p_Nq
    for(int c=0;c<p_Nq;++c){
      dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;
      qnew[id] = r_Aq[c];
    }
  }
}


__global__ void massMatrixMultiplyConstantKernel(const dlong Nelements,
						 const dlong  * __restrict__ elementIds,
						 const dfloat * __restrict__ cubvgeo,
						 const dfloat * __restrict__ oddI,
						 const dfloat * __restrict__ evenI,
						 const dfloat * __restrict__ q,
						 dfloat * __restrict__ qnew){
  
  __shared__ dfloat s_tmp1[p_Nblock][p_cubNq][p_cubNq][p_cubNq+p_padCubNq];

  dfloat r_Aq[p_cubNq];

  const unsigned int t = threadIdx.x;
  const int blk = threadIdx.y;
  
  const dlong e = blockIdx.x*p_Nblock + blk;

  const dlong element = (e<Nelements) ? elementIds[e]: 0;
  
  const unsigned int a = t%p_Nq;
  const unsigned int b = t/p_Nq;

  for(int c=0;c<p_Nq;++c){
    
    dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;
    
    r_Aq[c] = q[id];
  }

  __syncthreads();
  
  massMatrixMultiplyDevice(Nelements, e, element, cubvgeo, const_oddI, const_evenI, s_tmp1, r_Aq);
  
  if(e<Nelements){
#pragma unroll p_Nq
    for(int c=0;c<p_Nq;++c){
      dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;
      qnew[id] = r_Aq[c];
    }
  }
}








void massMatrixMultiplyHost(const dlong Nelements,
			    const dlong  * __restrict__ elementIds,
			    const dfloat * __restrict__ cubvgeo,
			    const dfloat * __restrict__ cubI,
			    const dfloat * __restrict__ q,
			    dfloat * __restrict__ qnew){


  dfloat qXXX[p_Nq][p_Nq][p_Nq];
  dfloat qIXX[p_cubNq][p_Nq][p_Nq];
  dfloat qIIX[p_cubNq][p_cubNq][p_Nq];
  dfloat qIII[p_cubNq][p_cubNq][p_cubNq];
    
  for(dlong e=0;e<Nelements;++e){

    for(int c=0;c<p_Nq;++c){
      for(int b=0;b<p_Nq;++b){
	for(int a=0;a<p_Nq;++a){
	  int id = e*p_Np + c*p_Nq2 + b*p_Nq + a;
	  qXXX[c][b][a] = q[id];
	}
      }
    }
    
    for(int k=0;k<p_cubNq;++k){
      for(int b=0;b<p_Nq;++b){
	for(int a=0;a<p_Nq;++a){
	  
	  dfloat res = 0;
	  
	  for(int c=0;c<p_Nq;++c){
	    dfloat Ikc = cubI[k*p_Nq+c];
	    res += Ikc*qXXX[c][b][a];
	  }
	  
	  qIXX[k][b][a] = res;
	}
      }
    }
    
    // interpolate in b
    for(int k=0;k<p_cubNq;++k){
      for(int j=0;j<p_cubNq;++j){
	for(int a=0;a<p_Nq;++a){
	  
	  dfloat res = 0;
	  
	  for(int b=0;b<p_Nq;++b){
	    dfloat Ijb = cubI[j*p_Nq+b];
	    res += Ijb*qIXX[k][b][a];
	  }
	  
	  qIIX[k][j][a] = res;
	}
      }
    }

    // interpolate in a
    for(int k=0;k<p_cubNq;++k){
      for(int j=0;j<p_cubNq;++j){
	for(int i=0;i<p_cubNq;++i){

	  dfloat res = 0;
	  
	  for(int a=0;a<p_Nq;++a){
	    dfloat Iia = cubI[i*p_Nq+a];
	    res += Iia*qIIX[k][j][a];
	  }
	  
	  int gid = e*p_cubNp + k*p_cubNq2 + j*p_cubNq + i;
	  
	  dfloat JW = cubvgeo[gid];

	  qIII[k][j][i] = res*JW;
	}
      }
    }


    // project in a
    for(int k=0;k<p_cubNq;++k){
      for(int j=0;j<p_cubNq;++j){
	for(int a=0;a<p_Nq;++a){

	  dfloat res = 0;
	  
	  for(int i=0;i<p_cubNq;++i){
	    dfloat Iia = cubI[i*p_Nq+a];
	    res += Iia*qIII[k][j][i];
	  }

	  qIIX[k][j][a] = res;
	}
      }
    }


    // project in b
    for(int k=0;k<p_cubNq;++k){
      for(int b=0;b<p_Nq;++b){
	for(int a=0;a<p_Nq;++a){

	  dfloat res = 0;

	  for(int j=0;j<p_cubNq;++j){
	    dfloat Ijb = cubI[j*p_Nq+b];
	    res += Ijb*qIIX[k][j][a];
	  }
	  
	  qIXX[k][b][a] = res;

	}
      }
    }


    // project in c
    for(int c=0;c<p_Nq;++c){
      for(int b=0;b<p_Nq;++b){
	for(int a=0;a<p_Nq;++a){

	  dfloat res = 0;

	  for(int k=0;k<p_cubNq;++k){
	    dfloat Ikc = cubI[k*p_Nq+c];
	    res += Ikc*qIXX[k][b][a];
	  }

	  int id = e*p_Np + c*p_Nq2 + b*p_Nq + a;
	  qnew[id] = res;
	}
      }
    }
  }
    
  
}


void buildInterpMatrices(dfloat *h_I,  dfloat **c_oddI, dfloat **c_evenI){

#if 0
  // now overwrite h_I and copy to c_I
  printf("I = [\n");
  for(int i=0;i<p_cubNq;++i){
    for(int a=0;a<p_Nq;++a){
      printf("% .4e ", h_I[i*p_Nq+a]);
    }
    printf("\n");
  }
  printf("];\n");
#endif
  
  //  cudaMemcpy(*c_I, *h_I, p_Nq*p_cubNq*sizeof(dfloat), cudaMemcpyHostToDevice);

  dfloat *X = (dfloat*) calloc(p_Nq*p_Nq, sizeof(dfloat));
  dfloat *invX = (dfloat*) calloc(p_Nq*p_Nq, sizeof(dfloat));

  dfloat *cubX = (dfloat*) calloc(p_cubNq*p_cubNq, sizeof(dfloat));
  dfloat *cubInvX = (dfloat*) calloc(p_cubNq*p_cubNq, sizeof(dfloat));

  for(int n=0;n<p_cubNq;++n){
    cubX[n*p_cubNq + n] = 1;
    cubInvX[n*p_cubNq + n] = 0.5;

    if(n<p_cubNq/2){
      cubX[n*p_cubNq + p_cubNq-1-n] = -1;
      cubInvX[n*p_cubNq + p_cubNq-1-n] = +0.5;
    }
    
    if(n>=(p_cubNq/2)){
      cubX[n*p_cubNq + p_cubNq-1-n] = +1;
      cubInvX[n*p_cubNq + p_cubNq-1-n] = -0.5;
    }
  }

  for(int n=0;n<p_Nq;++n){
    X[n*p_Nq + n] = 1;
    invX[n*p_Nq + n] = 0.5;

    if(n<p_Nq/2){
      X[n*p_Nq + p_Nq-1-n] = 1;
      invX[n*p_Nq + p_Nq-1-n] = -0.5;
    }
    
    if(n>=p_Nq/2){
      X[n*p_Nq + p_Nq-1-n] = -1;
      invX[n*p_Nq + p_Nq-1-n] = 0.5;
    }
  }

  if(p_Nq%2) invX[(p_Nq)*(p_Nq)/2] = 1;
  if(p_cubNq%2) cubInvX[(p_cubNq+1)*(p_cubNq+1)/2] = 1;
  
  dfloat *IinvX = (dfloat*) calloc(p_Nq*p_cubNq, sizeof(dfloat));
  dfloat *cubInvXIinvX = (dfloat*) calloc(p_Nq*p_cubNq, sizeof(dfloat));

  // post multiply by invX
  for(int i=0;i<p_cubNq;++i){
    for(int a=0;a<p_Nq;++a){
      dfloat res = 0;
      for(int n=0;n<p_Nq;++n){
	res += h_I[i*p_Nq+n]*invX[n*p_Nq+a];
      }
      IinvX[i*p_Nq+a] = res;
    }
  }

  // pre multiply by invX
  for(int i=0;i<p_cubNq;++i){
    for(int a=0;a<p_Nq;++a){
      dfloat res = 0;
      for(int n=0;n<p_cubNq;++n){
	res += cubInvX[i*p_cubNq+n]*IinvX[n*p_Nq + a];
      }
      cubInvXIinvX[i*p_Nq+a] = res;
    }
  }

  // now interleave the two non-zero blocks
  // [ A 0 ]  => [ A[0][0] B[0][0] A[0][1] B[0][1] .. A[0][p_halfNq-1] B[0][p_halfNq-1] .. 
  // [ 0 B ] 

  dfloat *oddI  = (dfloat*) calloc(p_cubNq*p_halfCubNq, sizeof(dfloat));
  dfloat *evenI = (dfloat*) calloc(p_cubNq*p_halfCubNq, sizeof(dfloat));
  
  for(int i=0;i<p_halfCubNq;++i){
    for(int a=0;a<p_halfNq;++a){

      oddI[i*p_halfNq+a] = cubInvXIinvX[i*p_Nq+a];
      evenI[i*p_halfNq+a] = cubInvXIinvX[(p_cubNq-1-i)*p_Nq + p_Nq-1-a];
      
    }
  }
      
  int NoddI = p_halfCubNq*p_halfNq;
  int NevenI = p_halfCubNq*p_halfNq;
  
  cudaMalloc(c_oddI, NoddI*sizeof(dfloat));
  cudaMalloc(c_evenI, NevenI*sizeof(dfloat));
  
  cudaMemcpy(*c_oddI,  oddI,  NoddI*sizeof(dfloat),  cudaMemcpyHostToDevice);
  cudaMemcpy(*c_evenI, evenI, NoddI*sizeof(dfloat), cudaMemcpyHostToDevice);
  
  cudaMemcpyToSymbol(const_oddI,  oddI,  NoddI*sizeof(dfloat));
  cudaMemcpyToSymbol(const_evenI, evenI, NoddI*sizeof(dfloat));
}


int main(int argc, char **argv){

  cudaStream_t stream;
  cudaStreamCreate(&stream);
  
  if(argc!=2){
    printf("Usage: ./massMatrixMultiplyVT Nelements\n");
    exit(-1);
  }

  // read number of elements
  hlong Nelements = atoi(argv[argc-1]);
  
  int    Ntotal = Nelements*p_Np;
  int cubNtotal = Nelements*p_cubNp;

  dfloat *h_cubvgeo, *c_cubvgeo;
  dfloat *h_qnew,    *c_qnew;
  dfloat *h_q,       *c_q;
  dfloat *h_I,       *c_I;
  dfloat *c_oddI,    *c_evenI;
  dfloat *h_garbage, *c_garbage;
  int    *h_elementIds, *c_elementIds;

  // list of elements
  h_elementIds = (int*) calloc(Nelements, sizeof(int));
  for(int e=0;e<Nelements;++e)
    h_elementIds[e] = e;
  cudaMalloc(&c_elementIds, Nelements*sizeof(int));
  cudaMemcpy(c_elementIds, h_elementIds, Nelements*sizeof(int), cudaMemcpyHostToDevice);
  
  // float fields
  dfloatRandAlloc(cubNtotal*p_Nvgeo, &h_cubvgeo, &c_cubvgeo);

  for(int e=0;e<Nelements;++e){
    for(int n=0;n<p_cubNp;++n){
      h_cubvgeo[e*p_cubNp+n] = drand48();
    }
  }

  cudaMemcpy(c_cubvgeo, h_cubvgeo, p_Nvgeo*Nelements*p_cubNp*sizeof(dfloat), cudaMemcpyHostToDevice);
  
  dfloatRandAlloc(Ntotal,       &h_q, &c_q);
  dfloatRandAlloc(Ntotal,       &h_qnew, &c_qnew);

  dfloatRandAlloc(p_Nq*p_cubNq, &h_I, &c_I);
  
  // give I the correct symmetry
  for(int i=0;i<p_halfCubNq;++i){
    for(int a=0;a<p_Nq;++a){
      h_I[(p_cubNq-1-i)*p_Nq + p_Nq-1-a] = h_I[i*p_Nq+a];
    }
  }

  // create Odd-even packed storage for I and transpose(I) and push to constant memory
  buildInterpMatrices(h_I, &c_oddI, &c_evenI);

  // flush L2 ??
  int sz = 32*1024*1024; // 32MB
  dfloatRandAlloc(sz, &h_garbage, &c_garbage);
  
  cudaEvent_t start, end;
  cudaEventCreate(&start);
  cudaEventCreate(&end);	

  int Ntests = 100;
  // KERNEL GRID
  dim3 G((Nelements+p_Nblock-1)/p_Nblock, 1, 1);
  dim3 B(p_Nq*p_Nq, p_Nblock, 1);

  float nothingElapsed = 0;
  {

    // time kernel that does nothing

#if USE_GRAPH==1
    // cuda stream capture sequence for nothingKernel
    cudaGraph_t nothingGraph;
  
    cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
    
    for(int test=0;test<Ntests;++test){
      nothingKernel <<< 1, 1, 0, stream >>> ();
    }
    
    cudaStreamEndCapture(stream, &nothingGraph);

    // time graph sequence for nothing
    cudaGraphExec_t nothingInstance;
    cudaGraphInstantiate(&nothingInstance, nothingGraph, NULL, NULL, 0);
    
    cudaEventRecord(start, stream);
    
    cudaGraphLaunch(nothingInstance, stream);
    
    cudaEventRecord(end, stream);
#else
    
    cudaEventRecord(start, stream);
    
    for(int test=0;test<Ntests;++test)
    nothingKernel <<< 1, 1, 0, stream >>> ();
    
    cudaEventRecord(end, stream);
    
#endif
    
    cudaDeviceSynchronize();
    
    cudaEventElapsedTime(&nothingElapsed, start, end);
    nothingElapsed /= 1000.;
    nothingElapsed /= (double) Ntests;
  
  }
  
  int cacheSwitchNq = 8;
  // warm up call
  if(p_Nq<=cacheSwitchNq) 
    massMatrixMultiplyRegisterKernel <<< G, B, 0, stream >>>(Nelements, c_elementIds, c_cubvgeo, c_oddI, c_evenI, c_q, c_qnew);
  else
    massMatrixMultiplyConstantKernel <<< G, B, 0, stream >>>(Nelements, c_elementIds, c_cubvgeo, c_oddI, c_evenI, c_q, c_qnew);

#if USE_GRAPH==1
  // cuda stream capture
  cudaGraph_t graph;
  
  cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);

  for(int test=0;test<Ntests;++test){

    if(p_Nq<=cacheSwitchNq) 
      massMatrixMultiplyRegisterKernel <<< G, B, 0, stream >>>(Nelements, c_elementIds, c_cubvgeo, c_oddI, c_evenI, c_q, c_qnew);
    else
      massMatrixMultiplyConstantKernel <<< G, B, 0, stream >>>(Nelements, c_elementIds, c_cubvgeo, c_oddI, c_evenI, c_q, c_qnew);
    
  }

  cudaStreamEndCapture(stream, &graph);
  
  cudaGraphExec_t instance;
  cudaGraphInstantiate(&instance, graph, NULL, NULL, 0);
#endif
  
  cudaDeviceSynchronize();

  {
    cudaEventRecord(start, stream);
    
#if USE_GRAPH==0
    for(int test=0;test<Ntests;++test){

      if(p_Nq<=cacheSwitchNq) 
	massMatrixMultiplyRegisterKernel <<< G, B, 0, stream >>>(Nelements, c_elementIds, c_cubvgeo, c_oddI, c_evenI, c_q, c_qnew);
      else
	massMatrixMultiplyConstantKernel <<< G, B, 0, stream >>>(Nelements, c_elementIds, c_cubvgeo, c_oddI, c_evenI, c_q, c_qnew);
      
    }
#else
    cudaGraphLaunch(instance, stream);
#endif

    cudaEventRecord(end, stream);
    
    cudaEventSynchronize(end);
    
    float elapsed;
    cudaEventElapsedTime(&elapsed, start, end);
    elapsed /= 1000.;
    elapsed /= (double) Ntests;
    
    printf("%d %d %d %lg %lg %lg %%%% [MassMatrixMultiply: N, Nelements, Ndofs, elapsed, dofsPerSecond, nothingElapsed]\n", p_Nq-1, Nelements, p_Np*Nelements, elapsed, Nelements*(p_Np/elapsed), nothingElapsed);
  }

  // check output is correct
  massMatrixMultiplyHost(Nelements, h_elementIds, h_cubvgeo, h_I, h_q, h_qnew);

  // copy device version to host old q
  dfloat *fromDevice = (dfloat*) calloc(Nelements*p_Np, sizeof(dfloat));
  cudaMemcpy(fromDevice, c_qnew, Nelements*p_Np*sizeof(dfloat), cudaMemcpyDeviceToHost);

  dfloat maxDiff = 0;
  
  for(int e=0;e<Nelements;++e){
    for(int n=0;n<p_Np;++n){
      int id = e*p_Np + n;
      dfloat diff = fabs(h_qnew[id]-fromDevice[id]);
      maxDiff = (diff>maxDiff) ? diff:maxDiff;
    }
  }
  printf("|| Mq_{host} - Mq_{device} ||_linf = %lg\n", maxDiff);
  
  cudaEventDestroy(start);
  cudaEventDestroy(end);	
  
  return 0;

}

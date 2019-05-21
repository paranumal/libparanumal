/*

The MIT License (MIT)

Copyright (c) 2019 Tim Warburton

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

#include <hip/hip_runtime.h>

#define dfloat_t double

__forceinline__ __device__ __host__  int ijN(const int i, const int j, const int N){

  return i + j*N;

}

__forceinline__ __device__ __host__ int ijkN(const int i, const int j, const int k, const int N){

  return i + j*N + k*N*N;

}

__forceinline__ __device__ __host__ int ijklN(const int i, const int j, const int k, const int l, const int N){

  return i + j*N + k*N*N + l*N*N*N;

}

// switch:
// 1 to use HIP 10.0 stream recording
// 0 to use traditional enqueing of kernels
#define USE_GRAPH 0

#define MAX_HALF_DOFS_1D 14
#define MAX_HALF_QUAD_1D 16

#define HALF_DOFS_1D ((NUM_DOFS_1D+1)/2)
#define HALF_QUAD_1D ((NUM_QUAD_1D+1)/2)

#define p_padCubNq  ((NUM_QUAD_1D%4) ? 0:1)

#define NUM_DOFS_2D (NUM_DOFS_1D*NUM_DOFS_1D)
#define NUM_DOFS_3D (NUM_DOFS_1D*NUM_DOFS_1D*NUM_DOFS_1D)

#define NUM_QUAD_2D (NUM_QUAD_1D*NUM_QUAD_1D)
#define NUM_QUAD_3D (NUM_QUAD_1D*NUM_QUAD_1D*NUM_QUAD_1D)

#define p_Nvgeo 1
#define p_JWID 0



__constant__ dfloat_t const_oddDofToQuad[MAX_HALF_QUAD_1D*MAX_HALF_DOFS_1D];
__constant__ dfloat_t const_evenDofToQuad[MAX_HALF_QUAD_1D*MAX_HALF_DOFS_1D];

void randAlloc(int N, dfloat_t **h_a, dfloat_t **c_a){

  *h_a = (dfloat_t*) calloc(N, sizeof(dfloat_t));

  for(int n=0;n<N;++n)
    h_a[0][n] = drand48();

  hipMalloc(c_a, N*sizeof(dfloat_t));

  hipMemcpy(c_a[0], h_a[0], N*sizeof(dfloat_t), hipMemcpyHostToDevice);

}

__global__ void nothingKernel(){  }

template <int NUM_DOFS_1D, int NUM_QUAD_1D, int p_Nblock >
  __forceinline__ __device__ 
  void massMatrixMultiplyDevice(const int numElements,
				const int element,
				const dfloat_t * __restrict__ op,
				const dfloat_t * __restrict__ oddDofToQuad,
				const dfloat_t * __restrict__ evenDofToQuad,
				dfloat_t s_Ap[p_Nblock][NUM_QUAD_1D][NUM_QUAD_1D][NUM_QUAD_1D+p_padCubNq],
				dfloat_t * __restrict__ r_Ap){

  dfloat_t r_tmpOdd[HALF_QUAD_1D];
  dfloat_t r_tmpEven[HALF_QUAD_1D];

  const int t   = hipThreadIdx_x;
  const int blk = hipThreadIdx_y;
  
  // assumes barrier before s_Ap was used last
  
  // transform in 'c'
  {
    const int a = t%NUM_DOFS_1D;
    const int b = t/NUM_DOFS_1D;
    
#pragma unroll
    for(int c=0;c<HALF_DOFS_1D;++c){
      r_tmpOdd[c]  = r_Ap[c] + r_Ap[NUM_DOFS_1D-1-c];
      r_tmpEven[c] = r_Ap[c] - r_Ap[NUM_DOFS_1D-1-c];
    }
    
#pragma unroll
    for(int k=0;k<HALF_QUAD_1D;++k){
      dfloat_t resOdd = 0, resEven = 0;
      
#pragma unroll
      for(int c=0;c<HALF_DOFS_1D;++c){
	int kc = ijN(c,k,HALF_DOFS_1D);		
	resOdd  += oddDofToQuad[kc]*r_tmpOdd[c];
	resEven += evenDofToQuad[kc]*r_tmpEven[c];
      }
      
      s_Ap[blk][k][b][a]               = resOdd + resEven;
      s_Ap[blk][NUM_QUAD_1D-1-k][b][a] = resOdd - resEven;
    }
    
  }
  
  __syncthreads();

  // transform in 'b'
  {
    for(int n=t;n<NUM_DOFS_1D*NUM_QUAD_1D;n+=NUM_DOFS_2D){
      const int a = n%NUM_DOFS_1D;
      const int k = n/NUM_DOFS_1D;

#pragma unroll
      for(int b=0;b<HALF_DOFS_1D;++b){
	dfloat_t ApOdd  = s_Ap[blk][k][b][a];
	dfloat_t ApEven = s_Ap[blk][k][NUM_DOFS_1D-1-b][a];
	r_tmpOdd[b]  = ApOdd + ApEven;
	r_tmpEven[b] = ApOdd - ApEven;
      }      
      
#pragma unroll
      for(int j=0;j<HALF_QUAD_1D;++j){
	dfloat_t resOdd = 0, resEven = 0;
	
#pragma unroll
	for(int b=0;b<HALF_DOFS_1D;++b){
	  int jb = ijN(b,j,HALF_DOFS_1D);
	  resOdd  += oddDofToQuad[jb]*r_tmpOdd[b];
	  resEven += evenDofToQuad[jb]*r_tmpEven[b];
	}
	
	s_Ap[blk][k][j][a]               = resOdd+resEven;
	s_Ap[blk][k][NUM_QUAD_1D-1-j][a] = resOdd-resEven;
	
      }
    }
  }
  
  __syncthreads();

  // transform in 'a'
  {
    for(int n=t;n<NUM_QUAD_2D;n+=NUM_DOFS_2D){
      const int j = n%NUM_QUAD_1D;
      const int k = n/NUM_QUAD_1D;
      
#pragma unroll
      for(int a=0;a<HALF_DOFS_1D;++a){
	dfloat_t ApOdd  = s_Ap[blk][k][j][a];
	dfloat_t ApEven = s_Ap[blk][k][j][NUM_DOFS_1D-1-a];
	r_tmpOdd[a]  = ApOdd + ApEven;
	r_tmpEven[a] = ApOdd - ApEven;
      }
      
#pragma unroll
      for(int i=0;i<HALF_QUAD_1D;++i){
	dfloat_t resOdd = 0, resEven = 0;
	
#pragma unroll
	for(int a=0;a<HALF_DOFS_1D;++a){
	  int ia = ijN(a,i,HALF_DOFS_1D);
	  resOdd  += oddDofToQuad[ia]*r_tmpOdd[a];
	  resEven += evenDofToQuad[ia]*r_tmpEven[a];
	}

	int gid1 = ijklN(i,j,k,element, NUM_QUAD_1D);
	int gid2 = ijklN(NUM_QUAD_1D-1-i,j,k,element, NUM_QUAD_1D); 
	
	dfloat_t WJ1 = (element<numElements) ? op[gid1]: 0;
	dfloat_t WJ2 = (element<numElements) ? op[gid2]: 0;
	
	dfloat_t ApOdd  = WJ1*(resOdd + resEven);
	dfloat_t ApEven = WJ2*(resOdd - resEven);

	r_Ap[i] = ApOdd + ApEven;
	r_Ap[NUM_QUAD_1D-1-i] = ApOdd - ApEven;
      }

#pragma unroll
      for(int a=0;a<HALF_DOFS_1D;++a){
	dfloat_t resOdd = 0, resEven = 0;
	
#pragma unroll
	for(int i=0;i<HALF_QUAD_1D;++i){
	  int ia = ijN(a,i,HALF_DOFS_1D);
	  resOdd  += oddDofToQuad[ia]*r_Ap[i];
	  resEven += evenDofToQuad[ia]*r_Ap[NUM_QUAD_1D-1-i];
	}
	
	s_Ap[blk][k][j][a]               = resOdd + resEven;
	s_Ap[blk][k][j][NUM_DOFS_1D-1-a] = resOdd - resEven;
      }
    }
  }
  
  __syncthreads();

  
  // test in 'b'
  {

    for(int n=t;n<NUM_DOFS_1D*NUM_QUAD_1D;n+=NUM_DOFS_2D){
      const int a = n%NUM_DOFS_1D;
      const int k = n/NUM_DOFS_1D;

      for(int j=0;j<HALF_QUAD_1D;++j){
	dfloat_t ApOdd  = s_Ap[blk][k][j][a];
	dfloat_t ApEven = s_Ap[blk][k][NUM_QUAD_1D-1-j][a];
	r_tmpOdd[j]  = ApOdd + ApEven;
	r_tmpEven[j] = ApOdd - ApEven;
      }

#pragma unroll
      for(int b=0;b<HALF_DOFS_1D;++b){
	dfloat_t resOdd = 0, resEven = 0;
	
#pragma unroll
	for(int j=0;j<HALF_QUAD_1D;++j){
	  int jb = ijN(b,j,HALF_DOFS_1D);
	  resOdd  += oddDofToQuad[jb]*r_tmpOdd[j];
	  resEven += evenDofToQuad[jb]*r_tmpEven[j];
	}
	
	s_Ap[blk][k][b][a]               = resOdd + resEven;
	s_Ap[blk][k][NUM_DOFS_1D-1-b][a] = resOdd - resEven;
      }
    }
  }
  
  __syncthreads();

  // test in 'c'
  {
    const int a = t%NUM_DOFS_1D;
    const int b = t/NUM_DOFS_1D;

    for(int k=0;k<HALF_QUAD_1D;++k){
      dfloat_t ApOdd  = s_Ap[blk][k][b][a];
      dfloat_t ApEven = s_Ap[blk][NUM_QUAD_1D-1-k][b][a];
      r_tmpOdd[k]  = ApOdd + ApEven;
      r_tmpEven[k] = ApOdd - ApEven;
    }
    
#pragma unroll
    for(int c=0;c<HALF_DOFS_1D;++c){
      dfloat_t resOdd = 0, resEven = 0;
      
#pragma unroll
      for(int k=0;k<HALF_QUAD_1D;++k){
	int kc = ijN(c,k,HALF_DOFS_1D);
	resOdd  += oddDofToQuad[kc]*r_tmpOdd[k];
	resEven += evenDofToQuad[kc]*r_tmpEven[k];
      }
      
      r_Ap[c]               = resOdd + resEven;
      r_Ap[NUM_DOFS_1D-1-c] = resOdd - resEven;
      
    }
  }

}

template <int NUM_DOFS_1D, int NUM_QUAD_1D, int p_Nblock >
  __global__ void massMatrixMultiplyRegisterKernel(const int numElements,
						   const dfloat_t * __restrict__ op,
						   const dfloat_t * __restrict__ oddDofToQuad,
						   const dfloat_t * __restrict__ evenDofToQuad,
						   const dfloat_t * __restrict__ solIn,
						   dfloat_t * __restrict__ solOut){
  
  __shared__ dfloat_t s_tmp1[p_Nblock][NUM_QUAD_1D][NUM_QUAD_1D][NUM_QUAD_1D+p_padCubNq];
  __shared__ dfloat_t s_oddDofToQuad[HALF_DOFS_1D*HALF_QUAD_1D];
  __shared__ dfloat_t s_evenDofToQuad[HALF_QUAD_1D*HALF_DOFS_1D];

  dfloat_t r_oddDofToQuad[HALF_QUAD_1D*HALF_DOFS_1D];
  dfloat_t r_evenDofToQuad[HALF_QUAD_1D*HALF_DOFS_1D];

  dfloat_t r_Aq[NUM_QUAD_1D];

  const unsigned int t = hipThreadIdx_x;
  const int blk = hipThreadIdx_y;
  
  const int element = hipBlockIdx_x*p_Nblock + blk;
  
  const unsigned int a = t%NUM_DOFS_1D;
  const unsigned int b = t/NUM_DOFS_1D;

  if(element<numElements){
    for(int c=0;c<NUM_DOFS_1D;++c){
      
      int id = ijklN(a,b,c,element, NUM_DOFS_1D); 
      
      r_Aq[c] = solIn[id];
    }
  }
  
  for(int n=t;n<HALF_DOFS_1D*HALF_QUAD_1D;n+=NUM_DOFS_2D){
    s_oddDofToQuad[n] = oddDofToQuad[n];
    s_evenDofToQuad[n] = evenDofToQuad[n];
    n+=NUM_DOFS_2D;
  }

  __syncthreads();

  // now copy shared data to thread local register arrays
  for(int n=0;n<HALF_DOFS_1D*HALF_QUAD_1D;++n){
    r_oddDofToQuad[n] = s_oddDofToQuad[n];
    r_evenDofToQuad[n] = s_evenDofToQuad[n];
  }
  
  massMatrixMultiplyDevice <NUM_DOFS_1D, NUM_QUAD_1D, p_Nblock>
    (numElements, element, op, r_oddDofToQuad, r_evenDofToQuad, s_tmp1, r_Aq);
  
  if(element<numElements){
#pragma unroll
    for(int c=0;c<NUM_DOFS_1D;++c){
      int id = ijklN(a,b,c,element,NUM_DOFS_1D);
      solOut[id] = r_Aq[c];
    }
  }
}


template <int NUM_DOFS_1D, int NUM_QUAD_1D, int p_Nblock >
  __global__ void massMatrixMultiplySharedKernel(const int numElements,
						 const dfloat_t * __restrict__ op,
						 const dfloat_t * __restrict__ oddDofToQuad,
						 const dfloat_t * __restrict__ evenDofToQuad,
						 const dfloat_t * __restrict__ solIn,
						 dfloat_t * __restrict__ solOut){
  
  __shared__ dfloat_t s_tmp1[p_Nblock][NUM_QUAD_1D][NUM_QUAD_1D][NUM_QUAD_1D+p_padCubNq];
  __shared__ dfloat_t s_oddDofToQuad[HALF_QUAD_1D*HALF_DOFS_1D];
  __shared__ dfloat_t s_evenDofToQuad[HALF_QUAD_1D*HALF_DOFS_1D];

  dfloat_t r_Aq[NUM_QUAD_1D];

  const unsigned int t = hipThreadIdx_x;
  const int blk = hipThreadIdx_y;
  
  const int element = hipBlockIdx_x*p_Nblock + blk;
  
  const unsigned int a = t%NUM_DOFS_1D;
  const unsigned int b = t/NUM_DOFS_1D;

  if(element<numElements){
    for(int c=0;c<NUM_DOFS_1D;++c){
      int id = ijklN(a,b,c,element,NUM_DOFS_1D);
      
      r_Aq[c] = solIn[id];
    }
  }

  for(int n=t;n<HALF_DOFS_1D*HALF_QUAD_1D;n+=NUM_DOFS_1D*NUM_DOFS_1D){
    s_oddDofToQuad[n] = oddDofToQuad[n];
    s_evenDofToQuad[n] = evenDofToQuad[n];
    n+=NUM_DOFS_2D;
  }

  __syncthreads();
  
  massMatrixMultiplyDevice  <NUM_DOFS_1D, NUM_QUAD_1D, p_Nblock>
    (numElements, element, op, s_oddDofToQuad, s_evenDofToQuad, s_tmp1, r_Aq);
  
  if(element<numElements){
#pragma unroll
    for(int c=0;c<NUM_DOFS_1D;++c){
      int id = ijklN(a,b,c,element,NUM_DOFS_1D);
      solOut[id] = r_Aq[c];
    }
  }
}

template <int NUM_DOFS_1D, int NUM_QUAD_1D, int p_Nblock >
__global__ void massMatrixMultiplyConstantKernel(const int numElements,
						 const dfloat_t * __restrict__ op,
						 const dfloat_t * __restrict__ oddDofToQuad,
						 const dfloat_t * __restrict__ evenDofToQuad,
						 const dfloat_t * __restrict__ solIn,
						 dfloat_t * __restrict__ solOut){
  
  __shared__ dfloat_t s_tmp1[p_Nblock][NUM_QUAD_1D][NUM_QUAD_1D][NUM_QUAD_1D+p_padCubNq];

  dfloat_t r_Aq[NUM_QUAD_1D];

  const unsigned int t = hipThreadIdx_x;
  const int blk = hipThreadIdx_y;
  
  const int element = hipBlockIdx_x*p_Nblock + blk;
  
  const unsigned int a = t%NUM_DOFS_1D;
  const unsigned int b = t/NUM_DOFS_1D;

  if(element<numElements){
    for(int c=0;c<NUM_DOFS_1D;++c){
      
      int id = ijklN(a,b,c,element,NUM_DOFS_1D);
      
      r_Aq[c] = solIn[id];
    }
  }

  __syncthreads();
  
  massMatrixMultiplyDevice  <NUM_DOFS_1D, NUM_QUAD_1D, p_Nblock>
    (numElements, element, op, const_oddDofToQuad, const_evenDofToQuad, s_tmp1, r_Aq);
  
  if(element<numElements){
#pragma unroll
    for(int c=0;c<NUM_DOFS_1D;++c){
      int id = ijklN(a,b,c,element,NUM_DOFS_1D);
      solOut[id] = r_Aq[c];
    }
  }
}

void massMatrixMultiplyHost(int NUM_DOFS_1D, int NUM_QUAD_1D, const int numElements,
			    const dfloat_t * __restrict__ op,
			    const dfloat_t * __restrict__ cubDofToQuad,
			    const dfloat_t * __restrict__ solIn,
			    dfloat_t * __restrict__ solOut){


  dfloat_t qXXX[NUM_DOFS_1D][NUM_DOFS_1D][NUM_DOFS_1D];
  dfloat_t qIXX[NUM_QUAD_1D][NUM_DOFS_1D][NUM_DOFS_1D];
  dfloat_t qIIX[NUM_QUAD_1D][NUM_QUAD_1D][NUM_DOFS_1D];
  dfloat_t qIII[NUM_QUAD_1D][NUM_QUAD_1D][NUM_QUAD_1D];
    
  for(int e=0;e<numElements;++e){

    for(int c=0;c<NUM_DOFS_1D;++c){
      for(int b=0;b<NUM_DOFS_1D;++b){
	for(int a=0;a<NUM_DOFS_1D;++a){
	  int id = ijklN(a,b,c,e,NUM_DOFS_1D);
	  qXXX[c][b][a] = solIn[id];
	}
      }
    }
    
    for(int k=0;k<NUM_QUAD_1D;++k){
      for(int b=0;b<NUM_DOFS_1D;++b){
	for(int a=0;a<NUM_DOFS_1D;++a){
	  
	  dfloat_t res = 0;
	  
	  for(int c=0;c<NUM_DOFS_1D;++c){
	    int kc = ijN(c,k,NUM_DOFS_1D);
	    dfloat_t Ikc = cubDofToQuad[kc];
	    res += Ikc*qXXX[c][b][a];
	  }
	  
	  qIXX[k][b][a] = res;
	}
      }
    }
    
    // interpolate in b
    for(int k=0;k<NUM_QUAD_1D;++k){
      for(int j=0;j<NUM_QUAD_1D;++j){
	for(int a=0;a<NUM_DOFS_1D;++a){
	  
	  dfloat_t res = 0;
	  
	  for(int b=0;b<NUM_DOFS_1D;++b){
	    int jb = ijN(b,j,NUM_DOFS_1D);
	    dfloat_t Ijb = cubDofToQuad[jb];
	    res += Ijb*qIXX[k][b][a];
	  }
	  
	  qIIX[k][j][a] = res;
	}
      }
    }

    // interpolate in a
    for(int k=0;k<NUM_QUAD_1D;++k){
      for(int j=0;j<NUM_QUAD_1D;++j){
	for(int i=0;i<NUM_QUAD_1D;++i){

	  dfloat_t res = 0;
	  
	  for(int a=0;a<NUM_DOFS_1D;++a){
	    int ia = ijN(a,i,NUM_DOFS_1D);
	    dfloat_t Iia = cubDofToQuad[ia];
	    res += Iia*qIIX[k][j][a];
	  }
	  
	  int gid = ijklN(i,j,k,e,NUM_QUAD_1D);
	  
	  dfloat_t JW = op[gid];

	  qIII[k][j][i] = res*JW;
	}
      }
    }


    // project in a
    for(int k=0;k<NUM_QUAD_1D;++k){
      for(int j=0;j<NUM_QUAD_1D;++j){
	for(int a=0;a<NUM_DOFS_1D;++a){

	  dfloat_t res = 0;
	  
	  for(int i=0;i<NUM_QUAD_1D;++i){
	    int ia = ijN(a,i,NUM_DOFS_1D);
	    dfloat_t Iia = cubDofToQuad[ia];
	    res += Iia*qIII[k][j][i];
	  }

	  qIIX[k][j][a] = res;
	}
      }
    }


    // project in b
    for(int k=0;k<NUM_QUAD_1D;++k){
      for(int b=0;b<NUM_DOFS_1D;++b){
	for(int a=0;a<NUM_DOFS_1D;++a){

	  dfloat_t res = 0;

	  for(int j=0;j<NUM_QUAD_1D;++j){
	    int jb = ijN(b,j,NUM_DOFS_1D);
	    dfloat_t Ijb = cubDofToQuad[j*NUM_DOFS_1D+b];
	    res += Ijb*qIIX[k][j][a];
	  }
	  
	  qIXX[k][b][a] = res;

	}
      }
    }


    // project in c
    for(int c=0;c<NUM_DOFS_1D;++c){
      for(int b=0;b<NUM_DOFS_1D;++b){
	for(int a=0;a<NUM_DOFS_1D;++a){

	  dfloat_t res = 0;

	  for(int k=0;k<NUM_QUAD_1D;++k){
	    int kc = ijN(c,k,NUM_DOFS_1D);
	    dfloat_t Ikc = cubDofToQuad[kc];
	    res += Ikc*qIXX[k][b][a];
	  }

	  int id = ijklN(a,b,c,e,NUM_DOFS_1D);
	  solOut[id] = res;
	}
      }
    }
  }
  
}

void buildInterpMatrices(int NUM_DOFS_1D, int NUM_QUAD_1D, dfloat_t *h_DofToQuad,  dfloat_t **c_oddDofToQuad, dfloat_t **c_evenDofToQuad){

  dfloat_t *X = (dfloat_t*) calloc(NUM_DOFS_1D*NUM_DOFS_1D, sizeof(dfloat_t));
  dfloat_t *invX = (dfloat_t*) calloc(NUM_DOFS_1D*NUM_DOFS_1D, sizeof(dfloat_t));

  dfloat_t *cubX = (dfloat_t*) calloc(NUM_QUAD_1D*NUM_QUAD_1D, sizeof(dfloat_t));
  dfloat_t *cubInvX = (dfloat_t*) calloc(NUM_QUAD_1D*NUM_QUAD_1D, sizeof(dfloat_t));

  for(int n=0;n<NUM_QUAD_1D;++n){
    cubX[n*NUM_QUAD_1D + n] = 1;
    cubInvX[n*NUM_QUAD_1D + n] = 0.5;

    if(n<NUM_QUAD_1D/2){
      cubX[n*NUM_QUAD_1D + NUM_QUAD_1D-1-n] = -1;
      cubInvX[n*NUM_QUAD_1D + NUM_QUAD_1D-1-n] = +0.5;
    }
    
    if(n>=(NUM_QUAD_1D/2)){
      cubX[n*NUM_QUAD_1D + NUM_QUAD_1D-1-n] = +1;
      cubInvX[n*NUM_QUAD_1D + NUM_QUAD_1D-1-n] = -0.5;
    }
  }

  for(int n=0;n<NUM_DOFS_1D;++n){
    X[n*NUM_DOFS_1D + n] = 1;
    invX[n*NUM_DOFS_1D + n] = 0.5;

    if(n<NUM_DOFS_1D/2){
      X[n*NUM_DOFS_1D + NUM_DOFS_1D-1-n] = 1;
      invX[n*NUM_DOFS_1D + NUM_DOFS_1D-1-n] = -0.5;
    }
    
    if(n>=NUM_DOFS_1D/2){
      X[n*NUM_DOFS_1D + NUM_DOFS_1D-1-n] = -1;
      invX[n*NUM_DOFS_1D + NUM_DOFS_1D-1-n] = 0.5;
    }
  }

  if(NUM_DOFS_1D%2) invX[(NUM_DOFS_1D)*(NUM_DOFS_1D)/2] = 1;
  if(NUM_QUAD_1D%2) cubInvX[(NUM_QUAD_1D+1)*(NUM_QUAD_1D+1)/2] = 1;
  
  dfloat_t *IinvX = (dfloat_t*) calloc(NUM_DOFS_1D*NUM_QUAD_1D, sizeof(dfloat_t));
  dfloat_t *cubInvXIinvX = (dfloat_t*) calloc(NUM_DOFS_1D*NUM_QUAD_1D, sizeof(dfloat_t));

  // post multiply by invX
  for(int i=0;i<NUM_QUAD_1D;++i){
    for(int a=0;a<NUM_DOFS_1D;++a){
      dfloat_t res = 0;
      for(int n=0;n<NUM_DOFS_1D;++n){
	res += h_DofToQuad[i*NUM_DOFS_1D+n]*invX[n*NUM_DOFS_1D+a];
      }
      IinvX[i*NUM_DOFS_1D+a] = res;
    }
  }

  // pre multiply by invX
  for(int i=0;i<NUM_QUAD_1D;++i){
    for(int a=0;a<NUM_DOFS_1D;++a){
      dfloat_t res = 0;
      for(int n=0;n<NUM_QUAD_1D;++n){
	res += cubInvX[i*NUM_QUAD_1D+n]*IinvX[n*NUM_DOFS_1D + a];
      }
      cubInvXIinvX[i*NUM_DOFS_1D+a] = res;
    }
  }

  // now interleave the two non-zero blocks
  // [ A 0 ]  => [ A[0][0] B[0][0] A[0][1] B[0][1] .. A[0][HALF_DOFS_1D-1] B[0][HALF_DOFS_1D-1] .. 
  // [ 0 B ] 

  dfloat_t *oddDofToQuad  = (dfloat_t*) calloc(NUM_QUAD_1D*HALF_QUAD_1D, sizeof(dfloat_t));
  dfloat_t *evenDofToQuad = (dfloat_t*) calloc(NUM_QUAD_1D*HALF_QUAD_1D, sizeof(dfloat_t));
  
  for(int i=0;i<HALF_QUAD_1D;++i){
    for(int a=0;a<HALF_DOFS_1D;++a){

      oddDofToQuad[i*HALF_DOFS_1D+a] = cubInvXIinvX[i*NUM_DOFS_1D+a];
      evenDofToQuad[i*HALF_DOFS_1D+a] = cubInvXIinvX[(NUM_QUAD_1D-1-i)*NUM_DOFS_1D + NUM_DOFS_1D-1-a];
      
    }
  }
      
  int NoddDofToQuad = HALF_QUAD_1D*HALF_DOFS_1D;
  int NevenDofToQuad = HALF_QUAD_1D*HALF_DOFS_1D;
  
  hipMalloc(c_oddDofToQuad, NoddDofToQuad*sizeof(dfloat_t));
  hipMalloc(c_evenDofToQuad, NevenDofToQuad*sizeof(dfloat_t));
  
  hipMemcpy(*c_oddDofToQuad,  oddDofToQuad,  NoddDofToQuad*sizeof(dfloat_t),  hipMemcpyHostToDevice);
  hipMemcpy(*c_evenDofToQuad, evenDofToQuad, NoddDofToQuad*sizeof(dfloat_t), hipMemcpyHostToDevice);
  
  hipMemcpyToSymbol(const_oddDofToQuad,  oddDofToQuad,  NoddDofToQuad*sizeof(dfloat_t));
  hipMemcpyToSymbol(const_evenDofToQuad, evenDofToQuad, NoddDofToQuad*sizeof(dfloat_t));
}


void runMassMatrixMultiplyKernel(hipStream_t stream, int Nq, int cubNq, int numElements,
				 dfloat_t *c_op, dfloat_t *c_oddDofToQuad, dfloat_t *c_evenDofToQuad,
				 dfloat_t *c_solIn, dfloat_t *c_solOut){
  
#define massMatrixMultiplyKernel(Nq,cubNq,Nblock)			\
  {									\
    if(Nq<=10) hipLaunchKernelGGL(massMatrixMultiplyRegisterKernel<Nq,cubNq,Nblock>, G, B, 0, stream, numElements, c_op, c_oddDofToQuad, c_evenDofToQuad, c_solIn, c_solOut); \
    else       hipLaunchKernelGGL(massMatrixMultiplySharedKernel<Nq,cubNq,Nblock>, G, B, 0, stream, numElements, c_op, c_oddDofToQuad, c_evenDofToQuad, c_solIn, c_solOut); \
  }

  // Nq=6 hipLaunchKernelGGL(massMatrixMultiplyRegisterKernel<Nq,cubNq,Nblock>, G, B, 0, stream, numElements, c_op, c_oddDofToQuad, c_evenDofToQuad, c_solIn, c_solOut); \
  //   Nq==4; hipLaunchKernelGGL(massMatrixMultiplySharedKernel<Nq,cubNq,Nblock>, G, B, 0, stream, numElements, c_op, c_oddDofToQuad, c_evenDofToQuad, c_solIn, c_solOut); \
  
#define ERR printf("massMatrixMultiplyRegister with Nq=%d, cubNq=%d not available", Nq, cubNq); exit(-1)

  int Nblock = 1;
  if(Nq==2){
    Nblock = 16;
    dim3 G((numElements+Nblock-1)/Nblock, 1, 1);
    dim3 B(Nq*Nq, Nblock, 1);

    // TUNED FOR SIMD WIDTH OF 64
    
    switch(cubNq){
    case 2: massMatrixMultiplyKernel(2,2,16); break;
    case 4: massMatrixMultiplyKernel(2,4,16); break;
    case 6: massMatrixMultiplyKernel(2,6,16); break;
    default: ERR;
    }
    return;
  }

  if(Nq==4){
    Nblock = 4;
    dim3 G((numElements+Nblock-1)/Nblock, 1, 1);
    dim3 B(Nq*Nq, Nblock, 1);
    
    switch(cubNq){
    case 4: massMatrixMultiplyKernel(4,4,4); break;
    case 6: massMatrixMultiplyKernel(4,6,4); break;
    case 8: massMatrixMultiplyKernel(4,8,4); break;
    default: ERR;
    }
    return;
  }

  if(Nq==6){

    Nblock = 3;
    dim3 G((numElements+Nblock-1)/Nblock, 1, 1);
    dim3 B(Nq*Nq, Nblock, 1);
    
    switch(cubNq){
    case 6:  massMatrixMultiplyKernel(6, 6,3); break;
    case 8:  massMatrixMultiplyKernel(6, 8,3); break;
    case 10: massMatrixMultiplyKernel(6,10,3); break;
    default: ERR;
    }
    return;
  }

  if(Nq==8){

    Nblock = 2;
    dim3 G((numElements+Nblock-1)/Nblock, 1, 1);
    dim3 B(Nq*Nq, Nblock, 1)  ;
    
    switch(cubNq){
    case 8:  massMatrixMultiplyKernel(8, 8,2); break;
    case 10: massMatrixMultiplyKernel(8,10,2); break;
    case 12: massMatrixMultiplyKernel(8,12,2); break;
    default: ERR;
    }
    return;
  }

  if(Nq==10){

    Nblock = 1;
    dim3 G((numElements+Nblock-1)/Nblock, 1, 1);
    dim3 B(Nq*Nq, Nblock, 1);    
    
    switch(cubNq){
    case 10: massMatrixMultiplyKernel(10,10,1); break;
    case 12: massMatrixMultiplyKernel(10,12,1); break;
    case 14: massMatrixMultiplyKernel(10,14,1); break;
    default: ERR;
    }
    return;
  }

  
  Nblock = 1;
  dim3 G((numElements+Nblock-1)/Nblock, 1, 1);
  dim3 B(Nq*Nq, Nblock, 1);    
    
  if(Nq==12){
    switch(cubNq){
    case 12: massMatrixMultiplyKernel(12,12,1); break;
    case 14: massMatrixMultiplyKernel(12,14,1); break;
    case 16: massMatrixMultiplyKernel(12,16,1); break;
    default: ERR;
    }
    return;
  }
  
  ERR;
}


dfloat_t nothingTest(hipStream_t stream, int Ntests){

  hipEvent_t start, end;
  hipEventCreate(&start);
  hipEventCreate(&end);	

  hipDeviceSynchronize();
  
  float nothingElapsed = 0;
  {
    
    // time kernel that does nothing
    
#if USE_GRAPH==1
    // hip stream capture sequence for nothingKernel
    hipGraph_t nothingGraph;
    
    hipStreamBeginCapture(stream, hipStreamCaptureModeGlobal);
    
    for(int test=0;test<Ntests;++test){
      hipLaunchKernelGGL(nothingKernel, 1, 1, 0, stream);
    }
    
    hipStreamEndCapture(stream, &nothingGraph);
    
    // time graph sequence for nothing
    hipGraphExec_t nothingInstance;
    hipGraphInstantiate(&nothingInstance, nothingGraph, NULL, NULL, 0);
    
    hipEventRecord(start, stream);
    
    hipGraphLaunch(nothingInstance, stream);
    
    hipEventRecord(end, stream);
#else
    
    hipEventRecord(start, stream);
    
    for(int test=0;test<Ntests;++test)
      hipLaunchKernelGGL(nothingKernel, 1, 1, 0, stream);
    
    hipEventRecord(end, stream);
    
#endif
    
    hipDeviceSynchronize();
    
    hipEventElapsedTime(&nothingElapsed, start, end);
    nothingElapsed /= 1000.;
    nothingElapsed /= (double) Ntests;
    
  }

  hipEventDestroy(start);
  hipEventDestroy(end);
  
  return nothingElapsed;
}


int main(int argc, char **argv){

  hipStream_t stream;
  hipStreamCreate(&stream);
  
  if(argc!=4){
    printf("Usage: ./massMatrixMultiplyVT Nq cubNq numElements\n");
    exit(-1);
  }

  // read number of elements
  int        Nq = atoi(argv[1]);
  int     cubNq = atoi(argv[2]);
  int numElements = atoi(argv[3]);

  printf("Running: Nq=%d, cubNq=%d, numElements=%d\n", Nq, cubNq, numElements);
  
  if(Nq%2){
    printf("Nq must be even\n");
    exit(-1);
  }

  if(cubNq%2){
    printf("cubNq must be even\n");
    exit(-1);
  }
  
  int   Np = Nq*Nq*Nq;
  int   cubNp = cubNq*cubNq*cubNq;

  int halfCubNq = ((cubNq+1)/2);

  int    Ntotal = numElements*Np;
  int cubNtotal = numElements*cubNp;

  dfloat_t *h_op,      *c_op;
  dfloat_t *h_solOut,       *c_solOut;
  dfloat_t *h_solIn,        *c_solIn;
  dfloat_t *h_DofToQuad,    *c_DofToQuad;
  dfloat_t *c_oddDofToQuad, *c_evenDofToQuad;

  // float fields
  randAlloc(cubNtotal*p_Nvgeo, &h_op, &c_op);

  for(int e=0;e<numElements;++e){
    for(int n=0;n<cubNp;++n){
      h_op[e*cubNp+n] = drand48();
    }
  }
  
  hipMemcpy(c_op, h_op, p_Nvgeo*numElements*cubNp*sizeof(dfloat_t), hipMemcpyHostToDevice);
  
  randAlloc(Ntotal, &h_solIn, &c_solIn);
  randAlloc(Ntotal, &h_solOut, &c_solOut);
  
  randAlloc(Nq*cubNq, &h_DofToQuad, &c_DofToQuad);
  
  // give I the correct symmetry
  for(int i=0;i<halfCubNq;++i){
    for(int a=0;a<Nq;++a){
      h_DofToQuad[(cubNq-1-i)*Nq + Nq-1-a] = h_DofToQuad[i*Nq+a];
    }
  }

  // create Odd-even packed storage for I and transpose(I) and push to constant memory
  buildInterpMatrices (Nq,cubNq, h_DofToQuad, &c_oddDofToQuad, &c_evenDofToQuad);

  hipEvent_t start, end;
  hipEventCreate(&start);
  hipEventCreate(&end);	

  int Ntests = 100;
  // KERNEL GRID
  // do nothing kernel test
  dfloat_t nothingElapsed = nothingTest(stream, Ntests);
  nothingElapsed = nothingTest(stream, Ntests);
  
  // warm up call
  runMassMatrixMultiplyKernel (stream, Nq, cubNq, numElements, c_op, c_oddDofToQuad, c_evenDofToQuad, c_solIn, c_solOut);

#if USE_GRAPH==1
  // hip stream capture
  hipGraph_t graph;
  
  hipStreamBeginCapture(stream, hipStreamCaptureModeGlobal);

  for(int test=0;test<Ntests;++test){

    runMassMatrixMultiplyKernel (stream, Nq, cubNq, numElements, c_op, c_oddDofToQuad, c_evenDofToQuad, c_solIn, c_solOut);
    
  }

  hipStreamEndCapture(stream, &graph);
  
  hipGraphExec_t instance;
  hipGraphInstantiate(&instance, graph, NULL, NULL, 0);
#endif
  
  hipDeviceSynchronize();

  {
    hipEventRecord(start, stream);
    
#if USE_GRAPH==0
    for(int test=0;test<Ntests;++test){

      runMassMatrixMultiplyKernel (stream, Nq, cubNq, numElements, c_op, c_oddDofToQuad, c_evenDofToQuad, c_solIn, c_solOut);
      
    }
#else
    hipGraphLaunch(instance, stream);
#endif

    hipEventRecord(end, stream);
    
    hipEventSynchronize(end);
    
    float elapsed;
    hipEventElapsedTime(&elapsed, start, end);
    elapsed /= 1000.;
    elapsed /= (double) Ntests;

    int bytesMoved = (2*Np+cubNp)*sizeof(dfloat_t); // x, Mx, opa
    double bw = (bytesMoved*numElements/elapsed)/1.e9;

    double effectiveFlops =
      numElements*(2*( Nq*Nq*Nq*cubNq*2 + Nq*Nq*cubNq*cubNq*2 + Nq*cubNq*cubNq*cubNq*2)/elapsed)/1.e9;
    
    printf("%2d %8d %8d %e %e %e %e %e %%%% [MassMatrixMultiply: N, numElements, Ndofs,"
	   " elapsed, dofsPerSecond, nothingElapsed, BW in GB/s, monolithic GFLOPS]\n",
	   Nq-1, numElements, Np*numElements, elapsed, numElements*(Np/elapsed), nothingElapsed, bw, effectiveFlops);
  }

  // check output is correct
  massMatrixMultiplyHost (Nq,cubNq,numElements, h_op, h_DofToQuad, h_solIn, h_solOut);

  // copy device version to host old q
  dfloat_t *fromDevice = (dfloat_t*) calloc(numElements*Np, sizeof(dfloat_t));
  hipMemcpy(fromDevice, c_solOut, numElements*Np*sizeof(dfloat_t), hipMemcpyDeviceToHost);

  dfloat_t maxDiff = 0;
  
  for(int e=0;e<numElements;++e){
    for(int n=0;n<Np;++n){
      int id = e*Np + n;
      dfloat_t diff = fabs(h_solOut[id]-fromDevice[id]);
      maxDiff = (diff>maxDiff) ? diff:maxDiff;
    }
  }
  printf("|| Mq_{host} - Mq_{device} ||_linf = %lg\n", maxDiff);
  
  hipEventDestroy(start);
  hipEventDestroy(end);	
  
  return 0;

}

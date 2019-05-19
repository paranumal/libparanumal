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

#define hipCHECK(code,mess)\
  if(code!=hipSuccess) printf("%s generated error code %d\n", mess, code);

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

#define MAX_QUAD_1D 10
#define MAX_DOFS_1D 10
#define MAX_HALF_QUAD_1D 10
#define MAX_HALF_DOFS_1D 10


#define HALF_DOFS_1D ((NUM_DOFS_1D+1)/2)
#define HALF_QUAD_1D ((NUM_QUAD_1D+1)/2)

#define p_padCubNq 0
//((NUM_QUAD_1D%4) ? 0:1)

#define NUM_DOFS_2D (NUM_DOFS_1D*NUM_DOFS_1D)
#define NUM_DOFS_3D (NUM_DOFS_1D*NUM_DOFS_1D*NUM_DOFS_1D)

#define NUM_QUAD_2D (NUM_QUAD_1D*NUM_QUAD_1D)
#define NUM_QUAD_3D (NUM_QUAD_1D*NUM_QUAD_1D*NUM_QUAD_1D)

#define p_Nggeo 7

#define p_G00ID 0
#define p_G01ID 1
#define p_G02ID 2
#define p_G11ID 3
#define p_G12ID 4
#define p_G22ID 5
#define p_GWJID 6

__constant__ dfloat_t const_DofToQuad[MAX_DOFS_1D*MAX_QUAD_1D];
__constant__ dfloat_t const_QuadToQuadD[MAX_QUAD_1D*MAX_QUAD_1D];

__constant__ dfloat_t const_oddDofToQuad[MAX_HALF_QUAD_1D*MAX_HALF_DOFS_1D];
__constant__ dfloat_t const_evenDofToQuad[MAX_HALF_QUAD_1D*MAX_HALF_DOFS_1D];


__constant__ dfloat_t const_oddQuadToQuadD[MAX_HALF_QUAD_1D*MAX_HALF_QUAD_1D];
__constant__ dfloat_t const_evenQuadToQuadD[MAX_HALF_QUAD_1D*MAX_HALF_QUAD_1D];


void randAlloc(int N, dfloat_t **h_a, dfloat_t **c_a){

  *h_a = (dfloat_t*) calloc(N, sizeof(dfloat_t));

  for(int n=0;n<N;++n)
    h_a[0][n] = drand48();

  hipCHECK(hipMalloc(c_a, N*sizeof(dfloat_t)), "randAlloc:hipMalloc");

  hipMemcpy(c_a[0], h_a[0], N*sizeof(dfloat_t), hipMemcpyHostToDevice);

}

__global__ void nothingKernel(){  }


template <int NUM_DOFS_1D, int NUM_QUAD_1D, int p_Nblock >
  __forceinline__ __device__ 
  void stiffnessMatrixMultiplyDevice(const int numElements,
				     const int element,
				     const dfloat_t lambda,
				     const dfloat_t * __restrict__ op,
				     const dfloat_t * __restrict__ oddDofToQuad,
				     const dfloat_t * __restrict__ evenDofToQuad,
				     const dfloat_t * __restrict__ QuadToQuadD,
				     const dfloat_t * __restrict__ oddQuadToQuadD,
				     const dfloat_t * __restrict__ evenQuadToQuadD,
				     dfloat_t s_Ap[p_Nblock][NUM_QUAD_1D][NUM_QUAD_1D][NUM_QUAD_1D+p_padCubNq],
				     dfloat_t * __restrict__ r_Ap){
  
  dfloat_t r_tmpOdd[HALF_QUAD_1D];
  dfloat_t r_tmpEven[HALF_QUAD_1D];
  
  __shared__ dfloat_t s_Gpr[p_Nblock][NUM_QUAD_1D][NUM_QUAD_1D];
  __shared__ dfloat_t s_Gps[p_Nblock][NUM_QUAD_1D][NUM_QUAD_1D];
  
  dfloat_t r_p[NUM_QUAD_1D];
  
  const int t   = hipThreadIdx_x;
  const int blk = hipThreadIdx_y;
  
  // assumes barrier before s_Ap was used last
  
  // transform in 'c'
  if(t<NUM_DOFS_2D){

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
	resOdd  += const_oddDofToQuad[kc]*r_tmpOdd[c];
	resEven += const_evenDofToQuad[kc]*r_tmpEven[c];
      }
      
      s_Ap[blk][k][b][a]               = resOdd + resEven;
      s_Ap[blk][NUM_QUAD_1D-1-k][b][a] = resOdd - resEven;
    }
    
  }
  
  __syncthreads();

  // transform in 'b'
  if(t<NUM_DOFS_1D*NUM_QUAD_1D){
    const int a = t%NUM_DOFS_1D;
    const int k = t/NUM_DOFS_1D;
    
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
	resOdd  += const_oddDofToQuad[jb]*r_tmpOdd[b];
	resEven += const_evenDofToQuad[jb]*r_tmpEven[b];
      }
      
      s_Ap[blk][k][j][a]               = resOdd+resEven;
      s_Ap[blk][k][NUM_QUAD_1D-1-j][a] = resOdd-resEven;
      
    }
  }
  
  __syncthreads();

  // transform in 'a'
  {
    const int j = t%NUM_QUAD_1D;
    const int k = t/NUM_QUAD_1D;
    
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
	resOdd  += const_oddDofToQuad[ia]*r_tmpOdd[a];
	resEven += const_evenDofToQuad[ia]*r_tmpEven[a];
      }

      
      s_Ap[blk][k][j][i] = resOdd + resEven;
      s_Ap[blk][k][j][NUM_QUAD_1D-1-i] = resOdd - resEven;
    }
  }

  __syncthreads();

  {
    const int i = t%NUM_QUAD_1D;
    const int j = t/NUM_QUAD_1D;
    
    for(int k = 0; k < NUM_QUAD_1D; k++) {
      r_p[k] = s_Ap[blk][k][j][i];; // prefetch operation
      r_Ap[k] = 0.f; // zero the accumulator
    }
    
    // Layer by layer
#pragma unroll
    for(int k = 0; k < NUM_QUAD_1D; k++) {
      
      dfloat_t G00 = 0, G01 =0, G02 =0, G11 =0, G12 =0, G22 =0, GWJ =0;
      
      // prefetch geometric factors
      const int gbase = element*p_Nggeo*NUM_QUAD_3D + ijkN(i,j,k,NUM_QUAD_1D);
      
      if(element<numElements){
	G00 = op[gbase+p_G00ID*NUM_QUAD_3D];
	G01 = op[gbase+p_G01ID*NUM_QUAD_3D];
	G02 = op[gbase+p_G02ID*NUM_QUAD_3D];
	G11 = op[gbase+p_G11ID*NUM_QUAD_3D];
	G12 = op[gbase+p_G12ID*NUM_QUAD_3D];
	G22 = op[gbase+p_G22ID*NUM_QUAD_3D];
	GWJ = op[gbase+p_GWJID*NUM_QUAD_3D];
      }
      
      dfloat_t pr = 0.f, ps = 0.f, pt = 0.f;
      
      __syncthreads();
      
#pragma unroll
      for(int m = 0; m < NUM_QUAD_1D; m++) {
	int im = ijN(m,i,NUM_QUAD_1D);
	int jm = ijN(m,j,NUM_QUAD_1D);
	int km = ijN(m,k,NUM_QUAD_1D);
	pr += QuadToQuadD[im]*s_Ap[blk][k][j][m];
	ps += QuadToQuadD[jm]*s_Ap[blk][k][m][i];
	pt += QuadToQuadD[km]*r_p[m];
      }
      
      s_Gpr[blk][j][i] = (G00*pr + G01*ps + G02*pt);
      s_Gps[blk][j][i] = (G01*pr + G11*ps + G12*pt);
      
      dfloat_t Gpt = (G02*pr + G12*ps + G22*pt);
      
      dfloat_t Apk = GWJ*lambda*r_p[k];
      
      __syncthreads();
      
#pragma unroll
      for(int m = 0; m < NUM_QUAD_1D; m++){
	int mi = ijN(i,m,NUM_QUAD_1D);
	int mj = ijN(j,m,NUM_QUAD_1D);
	int km = ijN(m,k,NUM_QUAD_1D);
	Apk     += QuadToQuadD[mi]*s_Gpr[blk][j][m];
	Apk     += QuadToQuadD[mj]*s_Gps[blk][m][i];
	r_Ap[m] += QuadToQuadD[km]*Gpt; // DT(m,k)*ut(i,j,k,e)
      }
      
      r_Ap[k] += Apk;
    }

    __syncthreads();
    
    for(int k=0;k<NUM_QUAD_1D;++k){
      s_Ap[blk][k][j][i] = r_Ap[k];     
    }
  }

  __syncthreads();
  
  // test in 'a'
  {
    const int j = t%NUM_QUAD_1D;
    const int k = t/NUM_QUAD_1D;
    
    // need to load from s_Ap into r_Ap
    
#pragma unroll
    for(int i=0;i<HALF_QUAD_1D;++i){
      dfloat_t ApOdd  = s_Ap[blk][k][j][i];
      dfloat_t ApEven = s_Ap[blk][k][j][NUM_QUAD_1D-1-i];
      r_tmpOdd[i]  = ApOdd + ApEven;
      r_tmpEven[i] = ApOdd - ApEven;
    }      
    
    
#pragma unroll
    for(int a=0;a<HALF_DOFS_1D;++a){
      dfloat_t resOdd = 0, resEven = 0;
      
#pragma unroll
      for(int i=0;i<HALF_QUAD_1D;++i){
	int ia = ijN(a,i,HALF_DOFS_1D);
	resOdd  += oddDofToQuad[ia]*r_tmpOdd[i];
	resEven += evenDofToQuad[ia]*r_tmpEven[i];
      }
      
      s_Ap[blk][k][j][a]               = resOdd + resEven;
      s_Ap[blk][k][j][NUM_DOFS_1D-1-a] = resOdd - resEven;
    }
  }
  
  __syncthreads();

  
  // test in 'b'
  if(t<NUM_DOFS_1D*NUM_QUAD_1D){
    const int a = t%NUM_DOFS_1D;
    const int k = t/NUM_DOFS_1D;
    
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
	resOdd  += const_oddDofToQuad[jb]*r_tmpOdd[j];
	resEven += const_evenDofToQuad[jb]*r_tmpEven[j];
      }
      
      s_Ap[blk][k][b][a]               = resOdd + resEven;
      s_Ap[blk][k][NUM_DOFS_1D-1-b][a] = resOdd - resEven;
    }
  }
  
  __syncthreads();

  // test in 'c'
  if(t<NUM_DOFS_2D){
    
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
	resOdd  += const_oddDofToQuad[kc]*r_tmpOdd[k];
	resEven += const_evenDofToQuad[kc]*r_tmpEven[k];
      }
      
      r_Ap[c]               = resOdd + resEven;
      r_Ap[NUM_DOFS_1D-1-c] = resOdd - resEven;
      
    }
  }

}

template <int NUM_DOFS_1D, int NUM_QUAD_1D, int p_Nblock >
__global__ void stiffnessMatrixMultiplyConstantKernel(const int numElements,
						 const dfloat_t lambda,
						 const dfloat_t * __restrict__ op,
						 const dfloat_t * __restrict__ oddDofToQuad,
						 const dfloat_t * __restrict__ evenDofToQuad,
						 const dfloat_t * __restrict__ QuadToQuadD,
						 const dfloat_t * __restrict__ oddQuadToQuadD,
						 const dfloat_t * __restrict__ evenQuadToQuadD,
						 const dfloat_t * __restrict__ solIn,
						 dfloat_t * __restrict__ solOut){
  
  __shared__ dfloat_t s_tmp1[p_Nblock][NUM_QUAD_1D][NUM_QUAD_1D][NUM_QUAD_1D+p_padCubNq];
  __shared__ dfloat_t s_QuadToQuadD[NUM_QUAD_2D];
  
  dfloat_t r_Aq[NUM_QUAD_1D];

  const unsigned int t = hipThreadIdx_x;
  const int blk = hipThreadIdx_y;
  
  const int element = hipBlockIdx_x*p_Nblock + blk;
  
  const unsigned int a = t%NUM_DOFS_1D;
  const unsigned int b = t/NUM_DOFS_1D;

  s_QuadToQuadD[t] = QuadToQuadD[t];
  
  if(element < numElements){
    if(t<NUM_DOFS_2D){
      for(int c=0;c<NUM_DOFS_1D;++c){
	
	int id = ijklN(a,b,c,element,NUM_DOFS_1D);
	
	r_Aq[c] = solIn[id];
      }
    }
  }

  __syncthreads();
  
  stiffnessMatrixMultiplyDevice  <NUM_DOFS_1D, NUM_QUAD_1D, p_Nblock>
    (numElements, element, lambda, op,
     oddDofToQuad, evenDofToQuad, s_QuadToQuadD, oddQuadToQuadD, evenQuadToQuadD, s_tmp1, r_Aq);
  
  if(element<numElements){
    if(t<NUM_DOFS_2D){
#pragma unroll
      for(int c=0;c<NUM_DOFS_1D;++c){
	int id = ijklN(a,b,c,element,NUM_DOFS_1D);
	solOut[id] = r_Aq[c];
      }
    }
  }
}

void stiffnessElementalMatrixMultiplyHost(int NUM_QUAD_1D, int element, dfloat_t lambda,
					  const dfloat_t * __restrict__ op,
					  const dfloat_t * __restrict__ QuadToQuadD,
					  const dfloat_t * qIII,
					  dfloat_t *lapqIII){
  
  dfloat_t Gqr[NUM_QUAD_1D][NUM_QUAD_1D][NUM_QUAD_1D];
  dfloat_t Gqs[NUM_QUAD_1D][NUM_QUAD_1D][NUM_QUAD_1D];
  dfloat_t Gqt[NUM_QUAD_1D][NUM_QUAD_1D][NUM_QUAD_1D];

  for(int k=0;k<NUM_QUAD_1D;++k){
    for(int j=0;j<NUM_QUAD_1D;++j){
      for(int i=0;i<NUM_QUAD_1D;++i){
	
	dfloat_t qr = 0;
	dfloat_t qs = 0;
	dfloat_t qt = 0;
	
	for(int n=0;n<NUM_QUAD_1D;++n){
	  int in = ijN(n,i,NUM_QUAD_1D);
	  int jn = ijN(n,j,NUM_QUAD_1D);
	  int kn = ijN(n,k,NUM_QUAD_1D);
	  
	  int kjn = ijkN(n,j,k,NUM_QUAD_1D);
	  int kni = ijkN(i,n,k,NUM_QUAD_1D);
	  int nji = ijkN(i,j,n,NUM_QUAD_1D);
	  
	  qr += QuadToQuadD[in]*qIII[kjn];
	  qs += QuadToQuadD[jn]*qIII[kni];
	  qt += QuadToQuadD[kn]*qIII[nji];	  
	}

	const int gbase = element*p_Nggeo*NUM_QUAD_3D + ijkN(i,j,k,NUM_QUAD_1D);
	
	dfloat_t G00 = op[gbase+p_G00ID*NUM_QUAD_3D];
	dfloat_t G01 = op[gbase+p_G01ID*NUM_QUAD_3D];
	dfloat_t G02 = op[gbase+p_G02ID*NUM_QUAD_3D];
	dfloat_t G11 = op[gbase+p_G11ID*NUM_QUAD_3D];
	dfloat_t G12 = op[gbase+p_G12ID*NUM_QUAD_3D];
	dfloat_t G22 = op[gbase+p_G22ID*NUM_QUAD_3D];
	
	Gqr[k][j][i] = (G00*qr + G01*qs + G02*qt);
	Gqs[k][j][i] = (G01*qr + G11*qs + G12*qt);
	Gqt[k][j][i] = (G02*qr + G12*qs + G22*qt);
      }
    }
  }


  for(int k=0;k<NUM_QUAD_1D;++k){
    for(int j=0;j<NUM_QUAD_1D;++j){
      for(int i=0;i<NUM_QUAD_1D;++i){
  
	int kji = ijkN(i,j,k,NUM_QUAD_1D);
	
	const int gbase = element*p_Nggeo*NUM_QUAD_3D + ijkN(i,j,k,NUM_QUAD_1D);

	dfloat_t GWJ = op[gbase+p_GWJID*NUM_QUAD_3D];
	dfloat_t lapq = lambda*GWJ*qIII[kji];
	
	for(int n=0;n<NUM_QUAD_1D;++n){
	  int ni = ijN(i,n,NUM_QUAD_1D);
	  int nj = ijN(j,n,NUM_QUAD_1D);
	  int nk = ijN(k,n,NUM_QUAD_1D);

	  lapq += QuadToQuadD[ni]*Gqr[k][j][n];
	  lapq += QuadToQuadD[nj]*Gqs[k][n][i];
	  lapq += QuadToQuadD[nk]*Gqt[n][j][i];	  
	}
	
	lapqIII[kji] = lapq;
      }
    }
  }
}

void stiffnessMatrixMultiplyHost(int NUM_DOFS_1D, int NUM_QUAD_1D, const int numElements, dfloat_t lambda,
				 const dfloat_t * __restrict__ op,
				 const dfloat_t * __restrict__ DofToQuad,
				 const dfloat_t * __restrict__ QuadToQuadD,
				 const dfloat_t * __restrict__ solIn,
				 dfloat_t * __restrict__ solOut){
  

  dfloat_t qXXX[NUM_DOFS_1D][NUM_DOFS_1D][NUM_DOFS_1D];
  dfloat_t qIXX[NUM_QUAD_1D][NUM_DOFS_1D][NUM_DOFS_1D];
  dfloat_t qIIX[NUM_QUAD_1D][NUM_QUAD_1D][NUM_DOFS_1D];
  dfloat_t qIII[NUM_QUAD_3D];
  dfloat_t lapqIII[NUM_QUAD_3D];

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
	    dfloat_t Ikc = DofToQuad[kc];
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
	    dfloat_t Ijb = DofToQuad[jb];
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
	    dfloat_t Iia = DofToQuad[ia];
	    res += Iia*qIIX[k][j][a];
	  }

	  int kji = ijkN(i,j,k,NUM_QUAD_1D);
	  qIII[kji] = res;
	}
      }
    }

    stiffnessElementalMatrixMultiplyHost(NUM_QUAD_1D, e, lambda, op, QuadToQuadD, qIII, lapqIII);

#if 0
    for(int k=0;k<NUM_QUAD_1D;++k){
      for(int j=0;j<NUM_QUAD_1D;++j){
	for(int i=0;i<NUM_QUAD_1D;++i){
	  lapqIII[k][j][i] = qIII[k][j][i];
	}
      }
    }
#endif
    
    // project in a
    for(int k=0;k<NUM_QUAD_1D;++k){
      for(int j=0;j<NUM_QUAD_1D;++j){
	for(int a=0;a<NUM_DOFS_1D;++a){

	  dfloat_t res = 0;
	  
	  for(int i=0;i<NUM_QUAD_1D;++i){
	    int ia = ijN(a,i,NUM_DOFS_1D);
	    dfloat_t Iia = DofToQuad[ia];
	    int kji = ijkN(i,j,k,NUM_QUAD_1D);
	    res += Iia*lapqIII[kji];
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
	    dfloat_t Ijb = DofToQuad[jb];
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
	    dfloat_t Ikc = DofToQuad[kc];
	    res += Ikc*qIXX[k][b][a];
	  }

	  int id = ijklN(a,b,c,e,NUM_DOFS_1D);
	  solOut[id] = res;
	}
      }
    }
  }
  
}

void buildOddEvenMatrices(int NUM_COLS_OP, int NUM_ROWS_OP,
			  dfloat_t  *h_OP, dfloat_t  *h_oddOP, dfloat_t  *h_evenOP,
			  dfloat_t **c_OP, dfloat_t **c_oddOP, dfloat_t **c_evenOP){

  int HALF_COLS_OP = ((NUM_COLS_OP+1)/2);
  int HALF_ROWS_OP = ((NUM_ROWS_OP+1)/2);
  
  dfloat_t *X = (dfloat_t*) calloc(NUM_COLS_OP*NUM_COLS_OP, sizeof(dfloat_t));
  dfloat_t *invX = (dfloat_t*) calloc(NUM_COLS_OP*NUM_COLS_OP, sizeof(dfloat_t));

  dfloat_t *cubX = (dfloat_t*) calloc(NUM_ROWS_OP*NUM_ROWS_OP, sizeof(dfloat_t));
  dfloat_t *cubInvX = (dfloat_t*) calloc(NUM_ROWS_OP*NUM_ROWS_OP, sizeof(dfloat_t));

  for(int n=0;n<NUM_ROWS_OP;++n){
    cubX[n*NUM_ROWS_OP + n] = 1;
    cubInvX[n*NUM_ROWS_OP + n] = 0.5;

    if(n<NUM_ROWS_OP/2){
      cubX[n*NUM_ROWS_OP + NUM_ROWS_OP-1-n] = -1;
      cubInvX[n*NUM_ROWS_OP + NUM_ROWS_OP-1-n] = +0.5;
    }
    
    if(n>=(NUM_ROWS_OP/2)){
      cubX[n*NUM_ROWS_OP + NUM_ROWS_OP-1-n] = +1;
      cubInvX[n*NUM_ROWS_OP + NUM_ROWS_OP-1-n] = -0.5;
    }
  }

  for(int n=0;n<NUM_COLS_OP;++n){
    X[n*NUM_COLS_OP + n] = 1;
    invX[n*NUM_COLS_OP + n] = 0.5;

    if(n<NUM_COLS_OP/2){
      X[n*NUM_COLS_OP + NUM_COLS_OP-1-n] = 1;
      invX[n*NUM_COLS_OP + NUM_COLS_OP-1-n] = -0.5;
    }
    
    if(n>=NUM_COLS_OP/2){
      X[n*NUM_COLS_OP + NUM_COLS_OP-1-n] = -1;
      invX[n*NUM_COLS_OP + NUM_COLS_OP-1-n] = 0.5;
    }
  }

  if(NUM_COLS_OP%2) invX[(NUM_COLS_OP)*(NUM_COLS_OP)/2] = 1;
  if(NUM_ROWS_OP%2) cubInvX[(NUM_ROWS_OP+1)*(NUM_ROWS_OP+1)/2] = 1;
  
  dfloat_t *IinvX = (dfloat_t*) calloc(NUM_COLS_OP*NUM_ROWS_OP, sizeof(dfloat_t));
  dfloat_t *cubInvXIinvX = (dfloat_t*) calloc(NUM_COLS_OP*NUM_ROWS_OP, sizeof(dfloat_t));

  // post multiply by invX
  for(int i=0;i<NUM_ROWS_OP;++i){
    for(int a=0;a<NUM_COLS_OP;++a){
      dfloat_t resI = 0;
      for(int n=0;n<NUM_COLS_OP;++n){
	resI += h_OP [i*NUM_COLS_OP+n]*invX[n*NUM_COLS_OP+a];
      }
      IinvX[i*NUM_COLS_OP+a] = resI;
    }
  }
  
  // pre multiply by invX
  for(int i=0;i<NUM_ROWS_OP;++i){
    for(int a=0;a<NUM_COLS_OP;++a){
      dfloat_t resI = 0;
      for(int n=0;n<NUM_ROWS_OP;++n){
	resI += cubInvX[i*NUM_ROWS_OP+n]*IinvX[n*NUM_COLS_OP + a];
      }
      cubInvXIinvX[i*NUM_COLS_OP+a] = resI;
    }
  }
  
  // now interleave the two non-zero blocks
  // [ A 0 ]  => [ A[0][0] B[0][0] A[0][1] B[0][1] .. A[0][HALF_DOFS_1D-1] B[0][HALF_DOFS_1D-1] .. 
  // [ 0 B ] 

  //  dfloat_t *oddOP  = (dfloat_t*) calloc(NUM_ROWS_OP*HALF_ROWS_OP, sizeof(dfloat_t));
  //  dfloat_t *evenOP = (dfloat_t*) calloc(NUM_ROWS_OP*HALF_ROWS_OP, sizeof(dfloat_t));
  
  for(int i=0;i<HALF_ROWS_OP;++i){
    for(int a=0;a<HALF_COLS_OP;++a){

      h_oddOP[i*HALF_COLS_OP+a]  = cubInvXIinvX[i*NUM_COLS_OP+a];
      h_evenOP[i*HALF_COLS_OP+a]  = cubInvXIinvX[(NUM_ROWS_OP-1-i)*NUM_COLS_OP + NUM_COLS_OP-1-a];
    }
  }
      
  int NoddOP  = HALF_ROWS_OP*HALF_COLS_OP;
  int NevenOP = HALF_ROWS_OP*HALF_COLS_OP;
  
  hipMalloc(c_oddOP,  NoddOP*sizeof(dfloat_t));
  hipMalloc(c_evenOP, NevenOP*sizeof(dfloat_t));
  
  hipMemcpy(*c_oddOP,  h_oddOP,   NoddOP*sizeof(dfloat_t), hipMemcpyHostToDevice);
  hipMemcpy(*c_evenOP, h_evenOP, NevenOP*sizeof(dfloat_t), hipMemcpyHostToDevice);

  hipMemcpy(*c_OP, h_OP,  NUM_COLS_OP*NUM_ROWS_OP*sizeof(dfloat_t),  hipMemcpyHostToDevice);

}


void runStiffnessMatrixMultiplyKernel(hipStream_t stream, int Nq, int cubNq, int numElements, dfloat_t lambda,
				      dfloat_t *c_op,
				      dfloat_t *c_oddDofToQuad, dfloat_t *c_evenDofToQuad,
				      dfloat_t *c_QuadToQuadD, dfloat_t *c_oddQuadToQuadD, dfloat_t *c_evenQuadToQuadD,
				      dfloat_t *c_solIn, dfloat_t *c_solOut){
  
#define stiffnessMatrixMultiplyKernel(Nq,cubNq,Nblock)			\
  {									\
    dim3 G((numElements+Nblock-1)/Nblock, 1, 1);			\
    dim3 B(cubNq*cubNq, Nblock, 1);					\
    hipLaunchKernelGGL(stiffnessMatrixMultiplyConstantKernel<Nq,cubNq,Nblock>, G, B, 0, stream, numElements, lambda, c_op, c_oddDofToQuad, c_evenDofToQuad, c_QuadToQuadD, c_oddQuadToQuadD,c_evenQuadToQuadD, c_solIn, c_solOut); \
  }
  
#define ERR printf("stiffnessMatrixMultiply with Nq=%d, cubNq=%d not available", Nq, cubNq); exit(-1)

  if(Nq==2){
    switch(cubNq){
    case 2: stiffnessMatrixMultiplyKernel(2,2,1); break;
    case 4: stiffnessMatrixMultiplyKernel(2,4,1); break;
    case 6: stiffnessMatrixMultiplyKernel(2,6,1); break;
    default: ERR;
    }
    return;
  }

  if(Nq==4){
    switch(cubNq){
    case 4: stiffnessMatrixMultiplyKernel(4,4,2); break;
    case 6: stiffnessMatrixMultiplyKernel(4,6,3); break;
    case 8: stiffnessMatrixMultiplyKernel(4,8,1); break;
    default: ERR;
    }
    return;
  }

  if(Nq==6){
    switch(cubNq){
    case 6:  stiffnessMatrixMultiplyKernel(6, 6,3); break;
    case 8:  stiffnessMatrixMultiplyKernel(6, 8,1); break;
    case 10: stiffnessMatrixMultiplyKernel(6,10,1); break;
    default: ERR;
    }
    return;
  }
  
  if(Nq==8){
    switch(cubNq){
    case 8:  stiffnessMatrixMultiplyKernel(8, 8,1); break;
    case 10: stiffnessMatrixMultiplyKernel(8,10,1); break;
    case 12: stiffnessMatrixMultiplyKernel(8,12,1); break;
    default: ERR;
    }
    return;
  }

  if(Nq==10){
    switch(cubNq){
    case 10: stiffnessMatrixMultiplyKernel(10,10,1); break;
    case 12: stiffnessMatrixMultiplyKernel(10,12,1); break;
    case 14: stiffnessMatrixMultiplyKernel(10,14,1); break;
    default: ERR;
    }
    return;
  }

  if(Nq==12){
    switch(cubNq){
    case 12: stiffnessMatrixMultiplyKernel(12,12,1); break;
    case 14: stiffnessMatrixMultiplyKernel(12,14,1); break;
    case 16: stiffnessMatrixMultiplyKernel(12,16,1); break;
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
    hipEventRecord(start, stream);
    
    for(int test=0;test<Ntests;++test)
      hipLaunchKernelGGL(nothingKernel, 1, 1, 0, stream);
    
    hipEventRecord(end, stream);

    hipDeviceSynchronize();
    
    hipEventElapsedTime(&nothingElapsed, start, end);
    nothingElapsed /= 1000.;
    nothingElapsed /= (double) Ntests;
    
  }

  return nothingElapsed;
}


int main(int argc, char **argv){

  //  hipSetDevice(0);
  
  hipStream_t stream;
  hipStreamCreate(&stream);
  
  if(argc!=4){
    printf("Usage: ./stiffnessMatrixMultiplyVT Nq cubNq numElements\n");
    exit(-1);
  }

  // read number of elements
  int        Nq = atoi(argv[1]);
  int     cubNq = atoi(argv[2]);
  int numElements = atoi(argv[3]);

  dfloat_t lambda = 1;
  
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

  int halfNq = ((Nq+1)/2);
  int halfCubNq = ((cubNq+1)/2);

  int    Ntotal = numElements*Np;
  int cubNtotal = numElements*cubNp;

  dfloat_t *h_op,      *c_op;
  dfloat_t *h_solOut,       *c_solOut;
  dfloat_t *h_solIn,        *c_solIn;
  dfloat_t *h_DofToQuad,    *c_DofToQuad;

  dfloat_t *h_oddDofToQuad, *h_evenDofToQuad;
  dfloat_t *c_oddDofToQuad, *c_evenDofToQuad;

  dfloat_t *h_QuadToQuadD,    *c_QuadToQuadD;
  dfloat_t *c_oddQuadToQuadD, *c_evenQuadToQuadD;
  dfloat_t *h_oddQuadToQuadD, *h_evenQuadToQuadD;

  // float fields
  randAlloc(cubNtotal*p_Nggeo, &h_op, &c_op);
  
  randAlloc(Ntotal, &h_solIn,  &c_solIn);
  randAlloc(Ntotal, &h_solOut, &c_solOut);
  
  randAlloc(Nq*cubNq,    &h_DofToQuad,   &c_DofToQuad);
  randAlloc(cubNq*cubNq, &h_QuadToQuadD, &c_QuadToQuadD);
  
  randAlloc(halfNq*halfCubNq,    &h_oddDofToQuad,   &c_oddDofToQuad);
  randAlloc(halfCubNq*halfCubNq, &h_oddQuadToQuadD, &c_oddQuadToQuadD);
  
  randAlloc(halfNq*halfCubNq,    &h_evenDofToQuad,   &c_evenDofToQuad);
  randAlloc(halfCubNq*halfCubNq, &h_evenQuadToQuadD, &c_evenQuadToQuadD);
  
  // give I the correct symmetry
  for(int i=0;i<halfCubNq;++i){
    for(int a=0;a<Nq;++a){
      h_DofToQuad[(cubNq-1-i)*Nq + Nq-1-a] = h_DofToQuad[i*Nq+a];
    }
  }

  // give D the correct symmetry
  for(int i=0;i<halfCubNq;++i){
    for(int a=0;a<Nq;++a){
      h_QuadToQuadD[(cubNq-1-i)*Nq + Nq-1-a] = -h_QuadToQuadD[i*Nq+a];
    }
  }

  // create Odd-even packed storage for I and transpose(I) and push to constant memory
  buildOddEvenMatrices (   Nq,cubNq, h_DofToQuad,   h_oddDofToQuad, h_evenDofToQuad, &c_DofToQuad,   &c_oddDofToQuad,   &c_evenDofToQuad  );
  buildOddEvenMatrices (cubNq,cubNq, h_QuadToQuadD, h_oddQuadToQuadD, h_evenQuadToQuadD, &c_QuadToQuadD, &c_oddQuadToQuadD, &c_evenQuadToQuadD);

  printf("halfNq=%d, halfCubNq=%d\n", halfNq, halfCubNq);
  
  hipCHECK(hipMemcpyToSymbol(HIP_SYMBOL(const_oddDofToQuad),  h_oddDofToQuad, halfNq*halfCubNq*sizeof(dfloat_t), 0, hipMemcpyHostToDevice), "symbol copy oddDofToQuad");
  hipCHECK(hipMemcpyToSymbol(HIP_SYMBOL(const_DofToQuad),     h_DofToQuad,    cubNq*Nq*sizeof(dfloat_t), 0, hipMemcpyHostToDevice), "symbol copy DofToQuad");
  hipCHECK(hipMemcpyToSymbol(HIP_SYMBOL(const_evenDofToQuad), h_evenDofToQuad,halfNq*halfCubNq*sizeof(dfloat_t), 0, hipMemcpyHostToDevice), "symbol copy even DofToQuad");

#if 0
  hipMemcpyToSymbol(HIP_SYMBOL(const_QuadToQuadD),     h_QuadToQuadD,     cubNq*cubNq*sizeof(dfloat_t), 0, hipMemcpyHostToDevice);
  hipMemcpyToSymbol(HIP_SYMBOL(const_oddQuadToQuadD),  h_oddQuadToQuadD,  halfCubNq*halfCubNq*sizeof(dfloat_t), 0, hipMemcpyHostToDevice);
  hipMemcpyToSymbol(HIP_SYMBOL(const_evenQuadToQuadD), h_evenQuadToQuadD, halfCubNq*halfCubNq*sizeof(dfloat_t), 0, hipMemcpyHostToDevice);
#endif
  
  hipEvent_t start, end;
  hipEventCreate(&start);
  hipEventCreate(&end);	

  hipDeviceSynchronize();
  
  int Ntests = 10;
  // KERNEL GRID
  // do nothing kernel test
  dfloat_t nothingElapsed = nothingTest(stream, Ntests);
  nothingElapsed = nothingTest(stream, Ntests);


  hipDeviceSynchronize();
  //  dfloat_t nothingElapsed =0;
  
  // warm up call
  runStiffnessMatrixMultiplyKernel (stream, Nq, cubNq, numElements, lambda,
				    c_op,
				    c_oddDofToQuad, c_evenDofToQuad,
				    c_QuadToQuadD, c_oddQuadToQuadD, c_evenQuadToQuadD,
				    c_solIn, c_solOut);
  
  hipDeviceSynchronize();

  {
    hipEventRecord(start, stream);
    
    for(int test=0;test<Ntests;++test){

      runStiffnessMatrixMultiplyKernel (stream, Nq, cubNq, numElements, lambda,
					c_op,
					c_oddDofToQuad, c_evenDofToQuad,
					c_QuadToQuadD, c_oddQuadToQuadD, c_evenQuadToQuadD,
					c_solIn, c_solOut);
      
    }

    hipEventRecord(end, stream);
    
    hipEventSynchronize(end);
    
    float elapsed;
    hipEventElapsedTime(&elapsed, start, end);
    elapsed /= 1000.;
    elapsed /= (double) Ntests;
    
    printf("%2d %8d %8d %e %e %e %%%% [MassMatrixMultiply: N, numElements, Ndofs,"
	   " elapsed, dofsPerSecond, nothingElapsed]\n",
	   Nq-1, numElements, Np*numElements, elapsed, numElements*(Np/elapsed), nothingElapsed);
  }

  // check output is correct
  stiffnessMatrixMultiplyHost (Nq, cubNq, numElements, lambda, h_op, h_DofToQuad, h_QuadToQuadD, h_solIn, h_solOut);

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

#if 0
  for(int e=0;e<numElements;++e){
    for(int n=0;n<cubNp;++n){
      int gid = e*cubNp*p_Nggeo + n;
      printf("vgeo[%d] = %e %e %e %e %e %e %e\n",
	     gid,
	     h_op[gid+p_G00ID*cubNp],
	     h_op[gid+p_G01ID*cubNp],
	     h_op[gid+p_G02ID*cubNp],
	     h_op[gid+p_G11ID*cubNp],
	     h_op[gid+p_G12ID*cubNp],
	     h_op[gid+p_G22ID*cubNp],
	     h_op[gid+p_GWJID*cubNp]);
    }
  }
#endif

#if 0
  for(int j=0;j<cubNq;++j){
    for(int i=0;i<cubNq;++i){
      printf(" %e ", h_QuadToQuadD[ijN(i,j,cubNq)]);
    }
    printf("\n");
  }
#endif
  
  hipEventDestroy(start);
  hipEventDestroy(end);	
  
  return 0;

}

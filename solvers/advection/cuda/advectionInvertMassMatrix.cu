#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

static const int p_Nq = 8;
static const int p_cubNq = 8;

#define p_Nq2 (p_Nq*p_Nq)
#define p_Np  (p_Nq*p_Nq*p_Nq)

#define p_cubNq2 (p_cubNq*p_cubNq)
#define p_cubNp  (p_cubNq*p_cubNq*p_cubNq)

#define p_Nvgeo 1
#define p_JWID 0

#define p_Nwarps ((p_Nq2+32-1)/32)

#define p_MAX_ITERATIONS 4

#define dlong int
#define hlong dlong
#define dfloat double


__constant__ dfloat const_I[p_cubNq][p_Nq];

//__device__ void advectionMassMatrixMultiply(const dlong element, const dfloat * __restrict__  cubvgeo, dfloat s_Ap[p_cubNq][p_cubNq][p_cubNq]){
__device__ void advectionMassMatrixMultiply_v00(const dlong element, const dfloat * __restrict__  r_WJ, dfloat s_Ap[p_cubNq][p_cubNq][p_cubNq]){

  dfloat r_tmp[p_Nq];
  
  const int t = threadIdx.x;

  __syncthreads();

  // transform in 'a'
  {
    const int b = t%p_cubNq;
    const int c = t/p_cubNq;
    
    if(b<p_Nq && c<p_Nq){
      for(int a=0;a<p_Nq;++a){
	r_tmp[a] = s_Ap[c][b][a];
      }

#pragma unroll p_cubNq
      for(int i=0;i<p_cubNq;++i){
	dfloat res = 0;

#pragma unroll p_Nq
	for(int a=0;a<p_Nq;++a){
	  res += const_I[i][a]*r_tmp[a];
	}
	s_Ap[c][b][i] = res;
      }
    }

  }
  
  __syncthreads();

  // transform in 'b'
  {
    const int i = t%p_cubNq;
    const int c = t/p_cubNq;
    
    if(c<p_Nq){
      for(int b=0;b<p_Nq;++b){
	r_tmp[b] = s_Ap[c][b][i];
      }

#pragma unroll p_cubNq
      for(int j=0;j<p_cubNq;++j){
	dfloat res = 0;

#pragma unroll p_Nq
	for(int b=0;b<p_Nq;++b){
	  res += const_I[j][b]*r_tmp[b];
	}
	s_Ap[c][j][i] = res;
      }
    }
  }


  __syncthreads();

  // transform in 'c'
  {
    const int i = t%p_cubNq;
    const int j = t/p_cubNq;

    for(int c=0;c<p_Nq;++c){
      r_tmp[c] = s_Ap[c][j][i];
    }

#pragma unroll p_cubNq
    for(int k=0;k<p_cubNq;++k){
      dfloat res = 0;
#pragma unroll p_Nq
      for(int c=0;c<p_Nq;++c){
	res += const_I[k][c]*r_tmp[c];
      }

      //      const dlong gid = element*p_cubNp*p_Nvgeo + t + k*p_cubNq*p_cubNq + p_JWID*p_cubNp;
      //      const dfloat WJ = cubvgeo[gid];
      const dfloat WJ = r_WJ[k];
      
      s_Ap[k][j][i] = WJ*res;
    }
  }

  __syncthreads();

  // test in 'c'
  {
    const int i = t%p_cubNq;
    const int j = t/p_cubNq;
    
    for(int k=0;k<p_cubNq;++k){
      r_tmp[k] = s_Ap[k][j][i];
    }

#pragma unroll p_Nq
    for(int c=0;c<p_Nq;++c){
      dfloat res = 0;
#pragma unroll p_cubNq
      for(int k=0;k<p_cubNq;++k){
	res += const_I[k][c]*r_tmp[k];
      }
      
      s_Ap[c][j][i] = res;
    }
  }

  __syncthreads();

  // test in 'b'
  {
    const int i = t%p_cubNq;
    const int c = t/p_cubNq;
    
    if(c<p_Nq){
      for(int j=0;j<p_cubNq;++j){
	r_tmp[j] = s_Ap[c][j][i];
      }

#pragma unroll p_Nq
      for(int b=0;b<p_Nq;++b){
	dfloat res = 0;
#pragma unroll  9
	for(int j=0;j<p_cubNq;++j){
	  res += const_I[j][b]*r_tmp[j];
	}
	s_Ap[c][b][i] = res;
      }
    }
  }
  
  __syncthreads();

  // test in 'a'
  {
    const int b = t%p_cubNq;
    const int c = t/p_cubNq;
    
    if(b<p_Nq && c<p_Nq){
      for(int i=0;i<p_cubNq;++i){
	r_tmp[i] = s_Ap[c][b][i];
      }

#pragma unroll p_Nq
      for(int a=0;a<p_Nq;++a){
	dfloat res = 0;
#pragma unroll p_cubNq
	for(int i=0;i<p_cubNq;++i){
	  res += const_I[i][a]*r_tmp[i];
	}
	s_Ap[c][b][a] = res;
      }
    }
  }
}

// use L1 for I matrix
__device__ void advectionMassMatrixMultiply_v01(const dlong element, const dfloat * __restrict__ c_I, const dfloat * __restrict__  r_WJ, dfloat s_Ap[p_cubNq][p_cubNq][p_cubNq]){

  dfloat r_tmp[p_Nq];
  
  const int t = threadIdx.x;

  __syncthreads();

  // transform in 'a'
  {
    const int b = t%p_cubNq;
    const int c = t/p_cubNq;
    
    if(b<p_Nq && c<p_Nq){
#pragma unroll p_Nq
      for(int a=0;a<p_Nq;++a){
	r_tmp[a] = s_Ap[c][b][a];
      }

#pragma unroll p_cubNq
      for(int i=0;i<p_cubNq;++i){
	dfloat res = 0;

#pragma unroll p_Nq
	for(int a=0;a<p_Nq;++a){
	  res += c_I[i*p_Nq+a]*r_tmp[a];
	}
	s_Ap[c][b][i] = res;
      }
    }

  }
  
  __syncthreads();

  // transform in 'b'
  {
    const int i = t%p_cubNq;
    const int c = t/p_cubNq;
    
    if(c<p_Nq){
      for(int b=0;b<p_Nq;++b){
	r_tmp[b] = s_Ap[c][b][i];
      }

#pragma unroll p_cubNq
      for(int j=0;j<p_cubNq;++j){
	dfloat res = 0;

#pragma unroll p_Nq
	for(int b=0;b<p_Nq;++b){
	  res += c_I[j*p_Nq+b]*r_tmp[b];
	}
	s_Ap[c][j][i] = res;
      }
    }
  }


  __syncthreads();

  // transform in 'c'
  {
    const int i = t%p_cubNq;
    const int j = t/p_cubNq;
    
    for(int c=0;c<p_Nq;++c){
      r_tmp[c] = s_Ap[c][j][i];
    }

#pragma unroll p_cubNq
    for(int k=0;k<p_cubNq;++k){
      dfloat res = 0;
#pragma unroll p_Nq
      for(int c=0;c<p_Nq;++c){
	res += c_I[k*p_Nq+c]*r_tmp[c];
      }

      //      const dlong gid = element*p_cubNp*p_Nvgeo + t + k*p_cubNq*p_cubNq + p_JWID*p_cubNp;
      //      const dfloat WJ = cubvgeo[gid];
      const dfloat WJ = r_WJ[k];
      
      s_Ap[k][j][i] = WJ*res;
    }
  }

  __syncthreads();

  // test in 'c'
  {
    const int i = t%p_cubNq;
    const int j = t/p_cubNq;
    
    for(int k=0;k<p_cubNq;++k){
      r_tmp[k] = s_Ap[k][j][i];
    }

#pragma unroll p_Nq
    for(int c=0;c<p_Nq;++c){
      dfloat res = 0;
#pragma unroll p_cubNq
      for(int k=0;k<p_cubNq;++k){
	res += c_I[k*p_Nq+c]*r_tmp[k];
      }
      
      s_Ap[c][j][i] = res;
    }
  }

  __syncthreads();

  // test in 'b'
  {
    const int i = t%p_cubNq;
    const int c = t/p_cubNq;
    
    if(c<p_Nq){
      for(int j=0;j<p_cubNq;++j){
	r_tmp[j] = s_Ap[c][j][i];
      }

#pragma unroll p_Nq
      for(int b=0;b<p_Nq;++b){
	dfloat res = 0;
#pragma unroll p_cubNq
	for(int j=0;j<p_cubNq;++j){
	  res += c_I[j*p_Nq+b]*r_tmp[j];
	}
	s_Ap[c][b][i] = res;
      }
    }
  }
  
  __syncthreads();

  // test in 'a'
  {
    const int b = t%p_cubNq;
    const int c = t/p_cubNq;
    
    if(b<p_Nq && c<p_Nq){
      for(int i=0;i<p_cubNq;++i){
	r_tmp[i] = s_Ap[c][b][i];
      }

#pragma unroll p_Nq
      for(int a=0;a<p_Nq;++a){
	dfloat res = 0;
#pragma unroll p_cubNq
	for(int i=0;i<p_cubNq;++i){
	  res += c_I[i*p_Nq+a]*r_tmp[i];
	}
	s_Ap[c][b][a] = res;
      }
    }
  }
}


// do C op first
__device__ void advectionMassMatrixMultiply(const dlong element, const dfloat * __restrict__  r_WJ, dfloat * __restrict__ r_p, dfloat s_Ap[p_cubNq][p_cubNq][p_cubNq]){

  dfloat r_tmp[p_Nq];
  
  const int t = threadIdx.x;

  //  __syncthreads();

  // transform in 'c'
  {
    const int a = t%p_Nq;
    const int b = t/p_Nq;
    
    if(b<p_Nq){
#pragma unroll p_cubNq
      for(int k=0;k<p_cubNq;++k){
	dfloat res = 0;

#pragma unroll p_Nq
	for(int c=0;c<p_Nq;++c){
	  res += const_I[k][c]*r_p[c];
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
      
#pragma unroll p_cubNq
      for(int j=0;j<p_cubNq;++j){
	dfloat res = 0;

#pragma unroll p_Nq
	for(int b=0;b<p_Nq;++b){
	  res += const_I[j][b]*r_tmp[b];
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

#pragma unroll p_cubNq
    for(int i=0;i<p_cubNq;++i){
      dfloat res = 0;
#pragma unroll p_Nq
      for(int a=0;a<p_Nq;++a){
	res += const_I[i][a]*r_tmp[a];
      }

      //      const dlong gid = element*p_cubNp*p_Nvgeo + t + k*p_cubNq*p_cubNq + p_JWID*p_cubNp;
      //      const dfloat WJ = cubvgeo[gid];
      const dfloat WJ = r_WJ[i];
      
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

#pragma unroll p_Nq
    for(int a=0;a<p_Nq;++a){
      dfloat res = 0;
#pragma unroll p_cubNq
      for(int i=0;i<p_cubNq;++i){
	res += const_I[i][a]*r_tmp[i];
      }
      
      s_Ap[k][j][a] = res;
    }
  }

  __syncthreads();

  // test in 'b'
  {
    const int a = t%p_Nq;
    const int k = t/p_Nq;
    
    if(k<p_cubNq){
      for(int j=0;j<p_cubNq;++j){
	r_tmp[j] = s_Ap[k][j][a];
      }

#pragma unroll p_Nq
      for(int b=0;b<p_Nq;++b){
	dfloat res = 0;
#pragma unroll  9
	for(int j=0;j<p_cubNq;++j){
	  res += const_I[j][b]*r_tmp[j];
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

#pragma unroll p_Nq
      for(int c=0;c<p_Nq;++c){
	dfloat res = 0;
#pragma unroll p_cubNq
	for(int k=0;k<p_cubNq;++k){
	  res += const_I[k][c]*r_tmp[k];
	}
	s_Ap[c][b][a] = res;
      }
    }
  }
}






__device__ dfloat dotProduct(dfloat a, volatile dfloat *s_a, volatile dfloat *s_warpa){
  
  const int t = threadIdx.x;
  int w = t/32;
  int n = t%32;
  
  __syncthreads();

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
__global__ void advectionInvertMassMatrix(const dlong Nelements,
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

  __shared__ dfloat s_tmp1[p_cubNq][p_cubNq][p_cubNq];
  volatile __shared__ dfloat s_tmp2[p_Nq2];
  volatile __shared__ dfloat s_tmpWarp[p_Nwarps];

  dfloat r_r[p_Nq], r_z[p_Nq], r_p[p_Nq], r_x[p_Nq], r_WJ[p_cubNq], r_precon[p_Nq];
  
  const dlong e = blockIdx.x;

  const dlong element = elementIds[e];

  const int t = threadIdx.x;
  
  const int a = t%p_Nq;
  const int b = t/p_Nq;

#if 0
  if(a==0 && b==0 && e==0){
    for(int n=0;n<p_cubNq;++n){
      for(int m=0;m<p_Nq;++m){
	const dfloat cI = const_I[n][m];
	printf("%f ", cI);
      }
      printf("\n");
    }
  }
#endif
  
  dfloat rdotz_ab = 0;
  dfloat rdotr_ab = 0;

  for(int k=0;k<p_cubNq;++k){
    const dlong gid = element*p_cubNp*p_Nvgeo + t + k*p_cubNq*p_cubNq + p_JWID*p_cubNp;
    r_WJ[k] = cubvgeo[gid];
  }
  
  if(t<p_Nq2){    
    for(int c=0;c<p_Nq;++c){

      dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;

      r_precon[c] = precon[id];
      
      r_r[c] = rhsq[id];
      //      r_z[c] = precon[id]*r_r[c];
      r_z[c] = r_precon[c]*r_r[c];
      r_p[c] = r_z[c];

      rdotz_ab += r_r[c]*r_z[c];
    }
  }

  dfloat rdotz = dotProduct(rdotz_ab, s_tmp2, s_tmpWarp);

  for(int it=0;it<p_MAX_ITERATIONS;++it){

#if 0
    // Ap
    if(t<p_Nq2){
#pragma unroll p_Nq
      for(int c=0;c<p_Nq;++c){
	s_tmp1[c][b][a] = r_p[c];
      }
    }
#endif
    
    //    advectionMassMatrixMultiply(element, cubvgeo, s_tmp1);
    advectionMassMatrixMultiply(element, r_p, r_WJ, s_tmp1);
    //    advectionMassMatrixMultiply(element, cubI, r_WJ, s_tmp1);

    dfloat pAp_ab = 0;
    if(t<p_Nq2){
#pragma unroll p_Nq
      for(int c=0;c<p_Nq;++c)
	pAp_ab += s_tmp1[c][b][a]*r_p[c];
    }
    dfloat pAp = dotProduct(pAp_ab, s_tmp2, s_tmpWarp);
  
    dfloat alpha = rdotz/pAp;

    rdotz_ab = 0;
    rdotr_ab = 0;
    
    if(t<p_Nq2){
#pragma unroll p_Nq
      for(int c=0;c<p_Nq;++c){
	r_x[c] += alpha*r_p[c];
	r_r[c] -= alpha*s_tmp1[c][b][a];
	
	//const dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;
	//	r_z[c] = precon[id]*r_r[c];
	r_z[c] = r_precon[c]*r_r[c];
	
	rdotz_ab += r_r[c]*r_z[c];
	rdotr_ab += r_r[c]*r_r[c];
      }
    }
    
    // r.z
    const dfloat rdotz_new = dotProduct(rdotz_ab, s_tmp2, s_tmpWarp);
    const dfloat rdotr     = dotProduct(rdotr_ab, s_tmp2, s_tmpWarp);

    if(rdotr<tol) break;
    
    const dfloat beta = rdotz_new/rdotz;

    rdotz = rdotz_new;

    if(t<p_Nq2){
#pragma unroll p_Nq
      for(int c=0;c<p_Nq;++c){
	r_p[c] = r_z[c] + beta*r_p[c];
      }
    }
  } // end iterations
  
  if(t<p_Nq2){
#pragma unroll p_Nq
    for(int c=0;c<p_Nq;++c){
      dlong id = a + b*p_Nq + c*p_Nq2 + element*p_Np;

      dfloat r_qcba = q[id];
      dfloat r_resq = resq[id];
      r_resq = rka*r_resq + dt*r_x[c];
      r_qcba += rkb*r_resq;
      
      resq[id] = r_resq;
      qnew[id] = r_qcba;
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
  cudaMemcpyToSymbol(const_I, h_I, p_Nq*p_cubNq*sizeof(dfloat));

  // LSERK fake constants
  dfloat rka = 1.3, rkb = 3.2, dt = 0.01;
  
  // call matrix inverse
  dim3 G(Nelements, 1, 1);
  dim3 B(p_cubNq*p_cubNq, 1, 1);

  int Ntests = 10;
  for(int test=0;test<Ntests;++test)
    advectionInvertMassMatrix <<< G, B >>>
      (Nelements, c_elementIds, dt, rka, rkb, tol, c_cubvgeo, c_I, c_precon, c_resq, c_rhsq, c_q, c_qnew);
      //      (Nelements, c_elementIds, c_cubvgeo, c_precon, c_rhsq, c_q);

  cudaDeviceSynchronize();

  return 0;

}

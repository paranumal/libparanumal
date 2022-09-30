#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

#include <cublas_v2.h>
#include <curand.h>


#if 0
#define dfloat double
#define dfloatString "double"
#else
#define dfloat float
#define dfloatString "float"
#endif

#define p_dim 3
#define p_Nfields 4

// scraped from recent 


#define p_Nvgeo 9

void generateRandArray(int sz, dfloat * a){
  //  a = (dfloat*) calloc(sz, sizeof(dfloat)); 
  for (int n=0; n<sz;++n){
    a[n] = drand48()-0.5;
    //printf("a[%d] = %f \n ", n, a[n]);
  }
}

//coalested reads
__global__ void geofactorsKernelv1(const int Nelements, 
				   const dfloat * __restrict__ vgeo,
				   const dfloat * __restrict__ Dq,
				   dfloat * __restrict__ Aq  ){

  int p_Np = blockDim.x; // define this
  int e  = blockIdx.x; 
  const int t  = threadIdx.x;

  const dfloat rx = vgeo[e*p_Nvgeo + 0];
  const dfloat sx = vgeo[e*p_Nvgeo + 1];
  const dfloat tx = vgeo[e*p_Nvgeo + 2];
  const dfloat ry = vgeo[e*p_Nvgeo + 3];
  const dfloat sy = vgeo[e*p_Nvgeo + 4];
  const dfloat ty = vgeo[e*p_Nvgeo + 5];
  const dfloat rz = vgeo[e*p_Nvgeo + 6];
  const dfloat sz = vgeo[e*p_Nvgeo + 7];
  const dfloat tz = vgeo[e*p_Nvgeo + 8];

  const int Dbase = e*p_dim*p_Np*p_Nfields + t;
  
  dfloat ur = Dq[Dbase + p_Np*p_dim*0 + p_Np*0];
  dfloat us = Dq[Dbase + p_Np*p_dim*0 + p_Np*1];
  dfloat ut = Dq[Dbase + p_Np*p_dim*0 + p_Np*2];
  dfloat vr = Dq[Dbase + p_Np*p_dim*1 + p_Np*0];
  dfloat vs = Dq[Dbase + p_Np*p_dim*1 + p_Np*1];
  dfloat vt = Dq[Dbase + p_Np*p_dim*1 + p_Np*2];
  dfloat wr = Dq[Dbase + p_Np*p_dim*2 + p_Np*0];
  dfloat ws = Dq[Dbase + p_Np*p_dim*2 + p_Np*1];
  dfloat wt = Dq[Dbase + p_Np*p_dim*2 + p_Np*2];
  dfloat pr = Dq[Dbase + p_Np*p_dim*3 + p_Np*0];
  dfloat ps = Dq[Dbase + p_Np*p_dim*3 + p_Np*1];
  dfloat pt = Dq[Dbase + p_Np*p_dim*3 + p_Np*2];

  const int base = e*p_Np*p_Nfields + t;

  Aq[base + p_Np*0] = rx*pr + sx*ps + tx*pt;
  Aq[base + p_Np*1] = ry*pr + sy*ps + ty*pt;
  Aq[base + p_Np*2] = rz*pr + sz*ps + tz*pt;
  

  dfloat divU = rx*ur + sx*us + tx*ut;
  divU += ry*vr + sy*vs + ty*vt;
  divU += rz*wr + sz*ws + tz*wt;

  Aq[base + p_Np*3] = divU;
}

void gpuFillRand(int N, dfloat **h_v, dfloat **c_v){

  *h_v = (dfloat*) calloc(N, sizeof(dfloat));

  for(int n=0;n<N;++n) h_v[0][n] = drand48();

  cudaMalloc(c_v, N*sizeof(dfloat));

  if(c_v==NULL){
    printf("gpuFillRand(%d,..,..) failed\n", N);
  }
  
  cudaMemcpy(*c_v, *h_v, N*sizeof(dfloat), cudaMemcpyHostToDevice);


}

void gpuBlasGemm(cublasHandle_t &handle, const dfloat *A, const dfloat *B, 
		 dfloat *C, const int m, const int n, const int k, const int option) {
 
  if (option == 1){
    //tranpose the matrices
    //7* p_Np, E,  p_Np
    //m, n, k

    //m: number of rows of matrix op(A) -- first matrix --  and C.
    // op(B) = trans(q) = trans(p_Np*E) - > E
    //n: number of columns in op(b) -- second matrix - and C 7*p_Np
    // op(A) = trans(D) = trans(7*p_Np x p_Np) -> [p_Np x 7*p_Np] -> 7p_Np
    // k the remaining number which is p_Np
    // result is C^T so size [7*Np x E] 
    // int lda=m,ldb=k,ldc=m;
    //C [p_Np*7 x E]
    //NW int lda=m,ldb=n,ldc=k;
    int lda = m, ldb =k, ldc=n; 
    //NW int lda = k, ldb = m, ldc =n;
    //NW int lda = k, ldb =n, ldc = m;
    //NW int lda = n, ldb =k, ldc =m;
    //NW int lda =n, ldb =m, ldc = k;



    const dfloat alf = 1;
    const dfloat bet = 0;
    const dfloat *alpha = &alf;
    const dfloat *beta = &bet;


    // Do the actual multiplication
    if(sizeof(dfloat)==8)
      cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_T, n, m, k, 
		  (double*)alpha, (double*)B, lda, 
		  (double*)A, lda, 
		  (double*)beta, 
		  (double*)C, ldc);
    else
      cublasSgemm(handle, CUBLAS_OP_T, CUBLAS_OP_T, n, m, k, (float*)alpha, (float*)B, ldb, (float*)A, lda, (float*)beta, (float*)C, ldc);
  }
  else{
    //default

    int lda=m,ldb=k,ldc=m;
    const dfloat alf = 1;
    const dfloat bet = 0;
    const dfloat *alpha = &alf;
    const dfloat *beta = &bet;


    // Do the actual multiplication
    if(sizeof(dfloat)==8)
      cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, 
		  (double*)alpha, (double*)A, lda, 
		  (double*)B, ldb, 
		  (double*)beta, 
		  (double*)C, ldc);
    else
      cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, (float*)alpha, (float*)A, lda, (float*)B, ldb, (float*)beta, (float*)C, ldc);
  }

}



int main(int argc, char **argv){

  size_t aux;

  int NX = 20, NY = 20, NZ = 20, p_N = 4;

  if(argc==5){
    NX = atoi(argv[1]);
    NY = atoi(argv[2]);
    NZ = atoi(argv[3]);
    p_N = atoi(argv[4]);
  }

  int E = NX*NY*NZ*6;
  int p_Np = ((p_N+1)*(p_N+2)*(p_N+3))/6;
  
  printf("E= %d p_N = %d p_Np = %d \n", E, p_N, p_Np);

  // number of geometric factors
  
  int Niter = 10;
  //(2*p_Np-1)*7*p_Np + 14*p_Np;
  unsigned long long int gflops= p_dim*p_Nfields*p_Np*p_Np*2 + p_Np*2*15;

  dfloat  *h_D, *h_q, *h_Dq, *h_Aq, *h_vgeo;
  dfloat  *d_D, *d_q, *d_Dq, *d_Aq, *d_vgeo;

  srand48(12345);  

  // allocate D
  gpuFillRand( (p_dim*p_Np)*p_Np, &h_D, &d_D); 

  //allocate q
  gpuFillRand( p_Nfields*p_Np*E,   &h_q,    &d_q); 

  //allocate Dq
  gpuFillRand(p_Nfields*p_dim*p_Np*E,   &h_Dq,    &d_Dq);
  
  //allocate Aq
  gpuFillRand(p_Nfields*p_Np*E,   &h_Aq,    &d_Aq); 

  //allocate geofactors
  gpuFillRand(p_Nvgeo*E,   &h_vgeo,    &d_vgeo);


  size_t free, total;

  // Create a handle for CUBLA
  cublasHandle_t handle;
  cublasCreate(&handle);

  // warm up
  gpuBlasGemm(handle, d_D, d_q, d_Dq, p_dim*p_Np, p_Nfields*E,  p_Np, 0);

  // create events
  cudaEvent_t start, stop;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);


  cudaEventRecord(start);

  for(int it=0;it<Niter;++it){
    //  [Dr;Ds;Dt] x [r,u,v,w,p]
    gpuBlasGemm(handle, d_D, d_q, d_Dq, p_dim*p_Np, p_Nfields*E,  p_Np, 0);

    //printf("done! %d \n", it);
    geofactorsKernelv1<<< E, p_Np >>> (E, d_vgeo, d_Dq, d_Aq);
  }

  cudaEventRecord(stop);
  cudaEventSynchronize(stop);

  float  elapsedCublas =0.0f;

  cudaEventElapsedTime(&elapsedCublas, start, stop);
  elapsedCublas /= (1000.*Niter);

  printf("[ DDD %d  %17.15f %17.15f ]\n", p_Np*E,  gflops*E/(elapsedCublas*1.e9), p_Np*E/(elapsedCublas*1.e9)); 

  //full version
  printf("TIME %5.4e flops = %llu  GFLOPS: %17.17f, NFLOP: %lld \n",elapsedCublas,gflops, (((dfloat)gflops*(dfloat)E/1.e9)/elapsedCublas), gflops*E);

  cudaMemcpy(h_Aq, d_Aq, p_Nfields*E*p_Np*sizeof(dfloat), cudaMemcpyDeviceToHost);

  // Destroy the handle
  cublasDestroy(handle);

  exit(0);
  return 0;

}


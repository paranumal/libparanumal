// to test for orders 1 to 10:
// for N in `seq 1 10` ; do nvcc -Dp_N=$N -arch=sm_60 --use_fast_math -o dgemm dgemm.cu -lcublas -lm; ./dgemm ; done

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

#include <cublas_v2.h>
#include <curand.h>

#define dfloat double

#ifndef p_N
#define p_N 4
#endif

#define p_Np ((int)((p_N+1)*(p_N+2))/2) 

// scraped from recent 
#if p_N==1
#define p_cubNp 6
#endif

#if p_N==2
#define p_cubNp 12
#endif

#if p_N==3
#define p_cubNp 19
#endif

#if p_N==4
#define p_cubNp 36
#endif

#if p_N==5
#define p_cubNp 54
#endif

#if p_N==6
#define p_cubNp 73
#endif

#if p_N==7
#define p_cubNp 93
#endif

#if p_N==8
#define p_cubNp 118
#endif

#if p_N==9
#define p_cubNp 145
#endif

#if p_N==10
#define p_cubNp 256
#endif

#define p_Nvgeo 7
#define p_RXID 0
#define p_RYID 1
#define p_SXID 2
#define p_SYID 3

__global__ void volumeFlux(const int Nelements, 
			   const dfloat * __restrict__ vgeo,
			   const dfloat * __restrict__ q,
			   dfloat * __restrict__ rhsq  ){

   const int e  = blockIdx.x; 
   const int t  = threadIdx.x;
   const int id = t + e*p_cubNp*4;
  
   const dfloat rx = vgeo[e*p_Nvgeo + p_RXID];
   const dfloat ry = vgeo[e*p_Nvgeo + p_RYID];
   const dfloat sx = vgeo[e*p_Nvgeo + p_SXID];
   const dfloat sy = vgeo[e*p_Nvgeo + p_SYID];
   
   const dfloat un  = q[id + 0*p_cubNp];
   const dfloat vn  = q[id + 1*p_cubNp];
   const dfloat udn = q[id + 2*p_cubNp];
   const dfloat vdn = q[id + 3*p_cubNp];

    const dfloat f11 = un*udn;
    const dfloat f12 = vn*udn;

    const dfloat f21 = un*vdn;
    const dfloat f22 = vn*vdn;
    
    rhsq[id + 0*p_cubNp] = rx*f11 + ry*f12;
    rhsq[id + 1*p_cubNp] = sx*f11 + sy*f12;
    rhsq[id + 2*p_cubNp] = rx*f21 + ry*f22;
    rhsq[id + 3*p_cubNp] = sx*f21 + sy*f22;
}

void gpuFillRand(int N, dfloat **h_v, dfloat **c_v){
  
  *h_v = (dfloat*) calloc(N, sizeof(dfloat));

  for(int n=0;n<N;++n) h_v[0][n] = drand48();
  
  cudaMalloc(c_v, N*sizeof(dfloat));
  
  cudaMemcpy(*c_v, *h_v, N*sizeof(dfloat), cudaMemcpyHostToDevice);

  
}

void gpuBlasGemm(cublasHandle_t &handle, const dfloat *A, const dfloat *B, dfloat *C, const int m, const int k, const int n) {
  int lda=m,ldb=k,ldc=m;
  const dfloat alf = 1;
  const dfloat bet = 0;
  const dfloat *alpha = &alf;
  const dfloat *beta = &bet;


  // Do the actual multiplication
  if(sizeof(dfloat)==8)
    cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, (double*)alpha, (double*)A, lda, (double*)B, ldb, (double*)beta, (double*)C, ldc);
  else
    cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, (float*)alpha, (float*)A, lda, (float*)B, ldb, (float*)beta, (float*)C, ldc);

}



int main(int argc, char **argv){

  int Nelements = (argc==1) ? 10000:atoi(argv[1]); // Write exact element number
  int Np        = p_Np;
  int Ncub      = p_cubNp;

  // fields q = (u,v,ud,vd)
  dfloat *h_q, *h_cq; 
  dfloat *d_q, *d_cq; 

  // geofacs
  dfloat *h_vgeo;
  dfloat *d_vgeo;
  
  // matrices
  dfloat *h_cI, *h_Div;
  dfloat *d_cI, *d_Div;

  // results
  dfloat *h_flux, *h_rhs;
  dfloat *d_flux, *d_rhs;

  // allocate geofacs
  gpuFillRand(p_Nvgeo*Nelements, &h_vgeo, &d_vgeo);

  // allocate arrays for matrices
  gpuFillRand(Ncub*Np,          &h_cI,   &d_cI); 
  gpuFillRand(2*Ncub*Np,        &h_Div,  &d_Div); 

  // allocate arrays for data
  gpuFillRand(4*Np*Nelements,   &h_q,    &d_q);
  gpuFillRand(4*Ncub*Nelements, &h_cq,   &d_cq);
  gpuFillRand(4*Ncub*Nelements, &h_flux, &d_flux); 
  gpuFillRand(2*Np*Nelements,   &h_rhs,  &d_rhs); 

  // Create a handle for CUBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);

  // create events
  cudaEvent_t start, stop;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);

  int Niterations = 100;

  for(int it=0;it<Niterations;++it){
    
    // interpolate from nodes to cubature
    gpuBlasGemm(handle, d_cI, d_q, d_cq, Ncub, Np, 4*Nelements);
    
    // compute volume fluxes
    dim3 G(Nelements,1,1);
    dim3 B(p_cubNp,1,1);
    
    volumeFlux<<< G, B >>> (Nelements, d_vgeo, d_cq, d_flux);
    
    // compute divergence
    gpuBlasGemm(handle, d_Div, d_flux, d_rhs, Np,  2*Ncub, 2*Nelements);
  }

  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  
  float elapsed; 
  
  cudaEventElapsedTime(&elapsed, start, stop);
  elapsed /= (Niterations*1000.);

  // minimal amount of data that could have moved (excluding matrices)
  long long int minData  = (4*Np + 2*Np + p_Nvgeo )*sizeof(dfloat);
  long long int actData  = (4*Np + 4*Ncub + 4*Ncub + 4*Ncub + 4*Ncub + 2*Np)*sizeof(dfloat);
  //long long int minFlops = (2*Np*Ncub*4 + Ncub*16 + 2*Np*Ncub*4);
  long long int minFlops = (2*Np*Ncub*4 + Ncub*4 + 2*Np*Ncub*4 + 8*Np);
  
  double GIG = 1024*1024*1024;
  double minBW  = Nelements*(minData/elapsed)/GIG;
  double actBW  = Nelements*(actData/elapsed)/GIG;
  double gflops = Nelements*(minFlops/elapsed)/GIG;
  
  printf("N=%d, K=%d, elapsed = %5.7E, minBW = %5.7E, actBW (est) = %5.7E, estGF = %5.7E\n", p_N, Nelements, elapsed, minBW, actBW, gflops);

  // Destroy the handle
  cublasDestroy(handle);


  exit(0);
  return 0;
  
}
  

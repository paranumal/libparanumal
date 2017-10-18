#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

#include <cublas_v2.h>
#include <curand.h>
#define dfloat double

#define p_N 5
#define p_Np ((int)((p_N+1)*(p_N+2))/2) 
#define p_cubNp ((int)(3*p_Np))
#define p_Nvgeo 7
#define p_RXID 0
#define p_RYID 1
#define p_SXID 2
#define p_SYID 3

__global__ void volumeFlux(const int Nelements, 
           const dfloat * __restrict__ vgeo,
           const dfloat * __restrict__ u,
           const dfloat * __restrict__ v,
           const dfloat * __restrict__ ud,
           const dfloat * __restrict__ vd,
           dfloat * __restrict__ rhsu  ){

   const int e  = blockIdx.x; 
   const int t  = threadIdx.x;
   const int id = t + e*blockDim.x; 
  
   const dfloat rx = vgeo[e*p_Nvgeo + p_RXID];
   const dfloat ry = vgeo[e*p_Nvgeo + p_RYID];
   const dfloat sx = vgeo[e*p_Nvgeo + p_SXID];
   const dfloat sy = vgeo[e*p_Nvgeo + p_SYID];

   // now have u,v,ur,us,vr,vs at cubature node c
    const dfloat un  = u [id + 0*p_cubNp];
    const dfloat vn  = v [id + 1*p_cubNp];
    const dfloat udn = ud[id + 2*p_cubNp];
    const dfloat vdn = vd[id + 3*p_cubNp];
    

    const dfloat f11 = un*udn;
    const dfloat f12 = vn*udn;

    const dfloat f21 = un*vdn;
    const dfloat f22 = vn*vdn;


    rhsu[id + 0*p_cubNp] = rx*f11 + ry*f21;
    rhsu[id + 1*p_cubNp] = sx*f11 + sy*f21;
    rhsu[id + 2*p_cubNp] = rx*f12 + ry*f12;
    rhsu[id + 3*p_cubNp] = rx*f22 + ry*f22;
  }

void gpuFillRand(int N, dfloat **h_v, dfloat **c_v){
  
  *h_v = (dfloat*) calloc(N, sizeof(dfloat));

  for(int n=0;n<N;++n) h_v[0][n] = drand48();
  
  cudaMalloc(c_v, N*sizeof(dfloat));
  
  cudaMemcpy(*c_v, *h_v, N*sizeof(dfloat), cudaMemcpyHostToDevice);

  
}

void gpuBlasGemm(const dfloat *A, const dfloat *B, dfloat *C, const int m, const int k, const int n) {
  int lda=m,ldb=k,ldc=m;
  const dfloat alf = 1;
  const dfloat bet = 0;
  const dfloat *alpha = &alf;
  const dfloat *beta = &bet;

  // Create a handle for CUBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);

  // Do the actual multiplication
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

  // Destroy the handle
  cublasDestroy(handle);
}



int main(int argc, char **argv){

  int N         = atoi(argv[1]);
  int Nelements = 10000; // Write exact element number
  int Np        = (N+1)*(N+2)/2;
  int Ncub      = 3*Np; // (Use exact values later)
  int nrI, ncI, ncB; 
  // Velocity
  dfloat *h_u, *h_v, *h_ud, *h_vd; 
  dfloat *d_u, *d_v, *d_ud, *d_vd; 
  // Velocity Vector;
  dfloat *h_I, *h_P, *h_B, *h_C, *h_D, *h_vgeo; 
  dfloat *d_I, *d_P, *d_B, *d_C, *d_D, *d_vgeo; 
 
  
   // C = I*B // Interpolate
   nrI = Ncub; ncI = Np; ncB = 4*Nelements; 
   gpuFillRand(nrI*ncI, &h_I, &d_I); 
   gpuFillRand(ncI*ncB, &h_B, &d_B);
   gpuFillRand(nrI*ncB, &h_C, &d_C);

    #if 1
   //
   gpuBlasGemm(d_I, d_B, d_C, nrI, ncI, ncB);

   #endif
   
   


   gpuFillRand(p_cubNp*Nelements, &h_u, &d_u); 
   gpuFillRand(p_cubNp*Nelements, &h_v, &d_v); 
   gpuFillRand(p_cubNp*Nelements, &h_ud, &d_ud); 
   gpuFillRand(p_cubNp*Nelements, &h_vd, &d_vd); 
   gpuFillRand(p_Nvgeo*Nelements, &h_vgeo, &d_vgeo); 
   gpuFillRand(p_cubNp*4*Nelements, &h_D, &d_D); 


   dim3 G(Nelements,1,1);
   dim3 B(p_cubNp,1,1);

   volumeFlux<<< G, B >>> (Nelements, d_vgeo, d_u, d_v, d_ud, d_vd, d_D);





   // K = P*D // Project
    int nrP, ncP, ncC; 
    nrP = Np;   ncP = 2*Ncub; ncC = 2*Nelements;
    gpuFillRand(nrP*ncP, &h_P, &d_P); 
    gpuFillRand(nrP*ncC, &h_D, &d_D); 
  
    gpuBlasGemm(d_P, d_C, d_D, nrP, ncP, ncC);
  
   

  // printf("Nbytes = %llu\n", Nbytes);
  
  // dim3 G(Nelements,1,1);
  // dim3 B(p_NSIMD,1,1);

  // printf("p_Np = %d\n", p_Np);
  // printf("p_cubNp = %d\n", p_cubNp);
  // printf("p_NSIMD = %d\n", p_NSIMD);
  // printf("p_CSIMD = %d\n", p_CSIMD);
  // printf("p_BSIMD = %d\n", p_BSIMD);
  // printf("G.x = %d, B.x = %d\n", G.x, B.x);
  
  
  // experimentalVolumeKernel <<< G, B >>> (Nelements, c_vgeo, c_cI, c_cDr, c_cDs, c_cProj, c_u, c_v, c_Nu, c_Nv);

  exit(0);
  return 0;
  
}
  
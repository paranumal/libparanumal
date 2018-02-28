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


void generateRandArray(int sz, dfloat * a){
  a = (dfloat*) calloc(sz, sizeof(dfloat)); 
  dfloat sum = 0;
  for (int n=0; n<sz;++n){
    a[n] = drand48()-0.5;
    sum += pow((*pt)[n],2);
  }
}

__global__ void geofactorsKernel(const int Nelements, 
    const dfloat * __restrict__ ggeo,
    const dfloat * __restrict__ q,
    dfloat * __restrict__ Aq  ){

  const int e  = blockIdx.x; 
  const int t  = threadIdx.x;
  const int id = t + e*p_cubNp*4;

  const int gid = e*p_Nggeo;
  const datafloat Grr = ggeo[gid + p_G00ID];
  const datafloat Grs = ggeo[gid + p_G01ID];
  const datafloat Grt = ggeo[gid + p_G02ID];
  const datafloat Gss = ggeo[gid + p_G11ID];
  const datafloat Gst = ggeo[gid + p_G12ID];
  const datafloat Gtt = ggeo[gid + p_G22ID];
  const datafloat J   = ggeo[gid + p_GWJID];

  const int id = n + e*p_Np;

  Aq[id] = 
    Grr*q[e*p_Np + 0*7 + t]+
    Grs*q[e*p_Np + 1*7 + t]+
    Grt*q[e*p_Np + 2*7 + t]+
    Gss*q[e*p_Np + 3*7 + t]+
    Gst*q[e*p_Np + 4*7 + t]+
    Gtt*q[e*p_Np + 5*7 + t]+ 
    +J*lambda*q[e*p_Np + 6*7 + t];


}


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


  int E = (argc>=2) ? atoi(argv[1]):512;
  int p_N = (argc>=3) ? atoi(argv[2]):5;
  int option = (argc>=4) ? atoi(argv[3]):1;  
  int p_Ne = (argc>=5) ? atoi(argv[4]):1;
  int p_Nb = (argc>=6) ? atoi(argv[5]):1;

  int p_Np = ((p_N+1)*(p_N+2)*(p_N+3))/6;
  int p_Nfp = ((p_N+1)*(p_N+2))/2;
  int p_Nfaces = 4;
  int p_NfacesNfp = p_Nfaces*p_Nfp;

  int BSIZE  = p_Np;


  // number of geometric factors
  int Nvgeo =  10;
  int Nsgeo = 6;

  int Niter = 10, it;

  int  gflops = p_Np*20*(1+p_Np);

  gflops *= Niter;

  dfloat  *d_Dcombined, *d_q, *d_Dq, *d_Aq;

  dfloat *h_Dcombined, *h_q, *h_Dq, *h_Aq;

  //other - needed for setup
  dfloat *SrrT,  *SrsT,  *SrtT,  *SssT,  *SsrT, *SstT,  *StsT, *StrT,  *SttT,  *MM; 
  dfloat *g_ggeo;


  //order of operations 
  //1) random calloc SrrT, etc
  generateRandArray(p_Np*p_Np, SrrT); 
  generateRandArray(p_Np*p_Np, SrsT); 
  generateRandArray(p_Np*p_Np, SrtT); 
  generateRandArray(p_Np*p_Np, SssT); 
  generateRandArray(p_Np*p_Np, SstT); 
  generateRandArray(p_Np*p_Np, SttT); 
  generateRandArray(p_Np*p_Np, MM); 

  //allocate q
  gpuFillRand(p_Np*E,   &h_q,    &d_q); 

  //allocate Dq

  gpuFillRand(7*p_Np*E,   &h_Dq,    &d_Dq); 
//allocate Aq
  gpuFillRand(p_Np*E,   &h_Aq,    &d_Aq); 


  //2) put Dcombined together
  h_Dcombined = (dfloat*) calloc(p_Np*p_Np*7, sizeof(dfloat)); 

  for (int n=0; n<p_Np; ++n){

    for (int k=0; k<7; ++k){
      if (k ==0){
        //copy Np entries from SrrT
        for (int m=0; m<p_Np; ++m){
          h_Dcombined[p_Np*n*7 + k*Np +m] = SrrT[p_Np*n+m]; 
        }
      }

      if (k ==1){
        //copy Np entries from SrsT
        for (int m=0; m<p_Np; ++m){
          h_Dcombined[p_Np*n*7 + k*Np +m] = SrsT[p_Np*n+m]; 
        }
      }

      if (k ==2){
        //copy Np entries from SrtT
        for (int m=0; m<p_Np; ++m){
          h_Dcombined[p_Np*n*7 + k*Np +m] = SrtT[p_Np*n+m]; 
        }
      }

      if (k ==3){
        //copy Np entries from SssT
        for (int m=0; m<p_Np; ++m){
          h_Dcombined[p_Np*n*7 + k*Np +m] = SssT[p_Np*n+m]; 
        }
      }

      if (k ==4){
        //copy Np entries from SstT
        for (int m=0; m<p_Np; ++m){
          h_Dcombined[p_Np*n*7 + k*Np +m] = SstT[p_Np*n+m]; 
        }
      }

      if (k == 5){
        //copy Np entries from SttT
        for (int m=0; m<p_Np; ++m){
          h_Dcombined[p_Np*n*7 + k*Np +m] = SttT[p_Np*n+m]; 
        }
      }

      if (k ==6){
        //copy Np entries from MM
        for (int m=0; m<p_Np; ++m){
          h_Dcombined[p_Np*n*7 + k*Np +m] = MM[p_Np*n+m]; 
        }
      }


    }
  }


  //copy Dcombined to device
  cudaMalloc(d_Dcombined, p_Np*p_Np*7*sizeof(dfloat)); 

  cudaMemcpy(d_Dcombined, h_Dcombined, p_Np*p_Np*7*sizeof(dfloat), cudaMemcpyHostToDevice);


  //3) use cublasDgemm -> outputs 7NpxNel array

  // Create a handle for CUBLA
  cublasHandle_t handle;
  cublasCreate(&handle);

  // create events
  cudaEvent_t start, stop;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /*
     m
     number of rows of matrix op(A) and rows of matrix C; m must be at least zero.
     n
     number of columns of matrix op(B) and number of columns of C; n must be at least zero.
     k
     number of columns of matrix op(A) and number of rows of op(B);k must be at least zero.
   */
  cudaEventRecord(start);


  for(int it=0;it<Niter;++it){

    // D*q
    gpuBlasGemm(handle, d_Dcombined, d_q, d_Dq, 7*p_Np, E, p_Np);
    //geo-factors

    dim3 G(Nelements,1,1);
    dim3 B(BSIZE,1,1);

    geofactorsKernel<<< G, B >>> (Nelements, d_ggeo, d_Dq, d_Aq);
}

  cudaEventRecord(stop);
  cudaEventSynchronize(stop);

  float elapseCublas

  cudaEventElapsedTime(&elapsedCublas, start, stop);
  elapsedCublas /= (Niterations*1000.);

  //6) use the original occa kernel to check for correctness
  //6.1) allocate the data for occa
  //6.2) run the kernel
  //6.3) copy the data back to CPU
  //6.4) compute a norm of the difference

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

  // Create a handle for CUBLA
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
  printf("%d %d %5.7E %5.7E %5.7E %5.7E\n", p_N, Nelements, elapsed, minBW, actBW, gflops);

  // Destroy the handle
  cublasDestroy(handle);


  exit(0);
  return 0;

}


/// for N in `seq 1 10` ; do nvcc -Dp_N=$N -arch=sm_60 --use_fast_math -o dgemm dgemm.cu -lcublas -lm; ./dgemm ; done

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



#ifndef p_N
#define p_N 2
#endif
 

// scraped from recent 


#define p_Nggeo 7

#define p_G00ID 0
#define p_G01ID 1
#define p_G02ID 2
#define p_G11ID 3
#define p_G12ID 4
#define p_G22ID 5
#define p_GWJID 6



void generateRandArray(int sz, dfloat * a){
//  a = (dfloat*) calloc(sz, sizeof(dfloat)); 
  for (int n=0; n<sz;++n){
    a[n] = drand48()-0.5;
//printf("a[%d] = %f \n ", n, a[n]);
  }
}





__global__ void geofactorsKernel(const int Nelements, 
    const dfloat * __restrict__ ggeo,
    const dfloat * __restrict__ q,
const dfloat lambda,
    dfloat * __restrict__ Aq  ){

int p_Np = blockDim.x;
  const int e  = blockIdx.x; 
  const int t  = threadIdx.x;
  const int id = t + e*p_Np;

  const int gid = e*p_Nggeo;
  const dfloat Grr = ggeo[gid + p_G00ID];
  const dfloat Grs = ggeo[gid + p_G01ID];
  const dfloat Grt = ggeo[gid + p_G02ID];
  const dfloat Gss = ggeo[gid + p_G11ID];
  const dfloat Gst = ggeo[gid + p_G12ID];
  const dfloat Gtt = ggeo[gid + p_G22ID];
  const dfloat J   = ggeo[gid + p_GWJID];

 

  Aq[id] = 
    Grr*q[e*7*p_Np + 0*7 + t]+
    Grs*q[e*7*p_Np + 1*7 + t]+
    Grt*q[e*7*p_Np + 2*7 + t]+
    Gss*q[e*7*p_Np + 3*7 + t]+
    Gst*q[e*7*p_Np + 4*7 + t]+
    Gtt*q[e*7*p_Np + 5*7 + t]+ 
    +J*lambda*q[e*7*p_Np + 6*7 + t];


}



void gpuFillRand(int N, dfloat **h_v, dfloat **c_v){

  *h_v = (dfloat*) calloc(N, sizeof(dfloat));

  for(int n=0;n<N;++n) h_v[0][n] = drand48();

  cudaMalloc(c_v, N*sizeof(dfloat));

  cudaMemcpy(*c_v, *h_v, N*sizeof(dfloat), cudaMemcpyHostToDevice);


}

void gpuBlasGemm(cublasHandle_t &handle, const dfloat *A, const dfloat *B, 
dfloat *C, const int m, const int n, const int k) {
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



int main(int argc, char **argv){


  int E = (argc>=2) ? atoi(argv[1]):512;
  //int p_N = (argc>=3) ? atoi(argv[2]):5;
  //int option = (argc>=4) ? atoi(argv[3]):1;  
  //int p_Ne = (argc>=5) ? atoi(argv[4]):1;
  //int p_Nb = (argc>=6) ? atoi(argv[5]):1;

  int p_Np = ((p_N+1)*(p_N+2)*(p_N+3))/6;
  int p_Nfp = ((p_N+1)*(p_N+2))/2;
  int p_Nfaces = 4;
  int p_NfacesNfp = p_Nfaces*p_Nfp;

  int BSIZE  = p_Np;
printf("E= %d p_N = %d p_Np = %d \n", E, p_N, p_Np);

  // number of geometric factors
  
  int Niter = 10, it;

  unsigned long long int  gflopsOLD = p_Np*20*(1+p_Np);

unsigned long long int gflops=(2*p_Np-1)*7*p_Np + 14*p_Np;
// gflops *= Niter;
//gflopsOLD *= Niter;
  dfloat  *d_Dcombined, *d_q, *d_Dq, *d_Aq, *d_ggeo;

  dfloat *h_Dcombined, *h_q, *h_Dq, *h_Aq, *h_ggeo;

  //other - needed for setup
dfloat *SrrT,  *SrsT,  *SrtT,  *SssT,  *SsrT, *SstT,  *StsT, *StrT,  *SttT,  *MM; 


SrrT = (dfloat*) calloc( p_Np* p_Np, sizeof(dfloat));
SrsT = (dfloat*) calloc( p_Np* p_Np, sizeof(dfloat));
SrtT = (dfloat*) calloc( p_Np* p_Np, sizeof(dfloat));
SssT = (dfloat*) calloc( p_Np* p_Np, sizeof(dfloat));
SstT = (dfloat*) calloc( p_Np* p_Np, sizeof(dfloat));
SttT = (dfloat*) calloc( p_Np* p_Np, sizeof(dfloat));
MM = (dfloat*) calloc( p_Np* p_Np, sizeof(dfloat));


 srand48(12345);  

 

//printf("check 1, about to allocate \n");
  //order of operations 
  //1) random calloc SrrT, etc
//printf("SrrT \n");
  generateRandArray( p_Np* p_Np, SrrT); 
//printf("\n\n");
 for (int n=0; n<p_Np*p_Np;++n){
    
//printf("SrrT[%d] = %f \n ", n, SrrT[n]);
  }


//printf("SrsT \n");
  generateRandArray( p_Np* p_Np, SrsT); 
//printf("SrtT \n");
  generateRandArray( p_Np* p_Np, SrtT); 
  
//printf("SssT \n");
generateRandArray( p_Np* p_Np, SssT); 

//printf("SstT \n");
  generateRandArray( p_Np* p_Np, SstT); 
//printf("SttT \n");
  generateRandArray( p_Np* p_Np, SttT); 
//printf("MM \n");
  generateRandArray( p_Np* p_Np, MM); 
printf("check 2, allocateD \n");
  //allocate q
  gpuFillRand( p_Np*E,   &h_q,    &d_q); 


  //allocate Dq

  gpuFillRand(7* p_Np*E,   &h_Dq,    &d_Dq); 
//allocate Aq
  gpuFillRand( p_Np*E,   &h_Aq,    &d_Aq); 
//allocate geofactors
  gpuFillRand(p_Nggeo*E,   &h_ggeo,    &d_ggeo);
printf("check 3, allocateD some more \n");

  //2) put Dcombined together
  h_Dcombined = (dfloat*) calloc( p_Np* p_Np*7, sizeof(dfloat)); 
//printf("check 4, Dcombined allocateD \n");
  for (int n=0; n< p_Np; ++n){

    for (int k=0; k<7; ++k){

      if (k ==0){
//printf("check 4.1, copying the contents of SrrT \n");
        //copy k*p_Np entries from SrrT
        for (int m=0; m< p_Np; ++m){

          h_Dcombined[ p_Np*n*7 + k*k*p_Np +m] = SrrT[ p_Np*n+m]; 
//printf("SrrT[%d] = %f \n", p_Np*n+m,SrrT[ n*p_Np +m] );
        }
      }

      if (k ==1){
        //copy k*p_Np entries from SrsT
        for (int m=0; m< p_Np; ++m){
          h_Dcombined[ p_Np*n*7 + k* p_Np +m] = SrsT[ p_Np*n+m]; 
        }
      }

      if (k ==2){
        //copy k*p_Np entries from SrtT
        for (int m=0; m< p_Np; ++m){
          h_Dcombined[ p_Np*n*7 + k*k*p_Np +m] = SrtT[ p_Np*n+m]; 
        }
      }

      if (k ==3){
        //copy k*p_Np entries from SssT
        for (int m=0; m< p_Np; ++m){
          h_Dcombined[ p_Np*n*7 + k*k*p_Np +m] = SssT[ p_Np*n+m]; 
        }
      }

      if (k ==4){
        //copy k*p_Np entries from SstT
        for (int m=0; m< p_Np; ++m){
          h_Dcombined[ p_Np*n*7 + k*k*p_Np +m] = SstT[ p_Np*n+m]; 
        }
      }

      if (k == 5){
        //copy k*p_Np entries from SttT
        for (int m=0; m< p_Np; ++m){
          h_Dcombined[ p_Np*n*7 + k*k*p_Np +m] = SttT[ p_Np*n+m]; 
        }
      }

      if (k ==6){
        //copy k*p_Np entries from MM
        for (int m=0; m< p_Np; ++m){
          h_Dcombined[ p_Np*n*7 + k*k*p_Np +m] = MM[ p_Np*n+m]; 
        }
      }


    }
  }
//printf("check 5, Dcombined ready\n");

  //copy Dcombined to device
  cudaMalloc(&d_Dcombined,  p_Np* p_Np*7*sizeof(dfloat)); 


  cudaMemcpy(d_Dcombined, h_Dcombined,  p_Np* p_Np*7*sizeof(dfloat), cudaMemcpyHostToDevice);


  //3) use cublasDgemm -> outputs 7k*p_NpxNel array

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
    gpuBlasGemm(handle, d_Dcombined, d_q, d_Dq, 7* p_Np, E,  p_Np);
    //geo-factors

    dim3 G(E,1,1);
    dim3 B(BSIZE,1,1);
dfloat lambda = 1.0;
    geofactorsKernel<<< G, B >>> (E, d_ggeo, d_Dq,lambda, d_Aq);
}

  cudaEventRecord(stop);
  cudaEventSynchronize(stop);

  float elapsedCublas;

  cudaEventElapsedTime(&elapsedCublas, start, stop);
  elapsedCublas /= (1000.*Niter);
printf("TIME %17.16f flops = %llu  GFLOPS: %17.17f \n",elapsedCublas,gflops, 
(((dfloat)gflops*(dfloat)E)/elapsedCublas)/1000000000.0);

cudaMemcpy(h_Aq, d_Aq, E*p_Np*sizeof(dfloat), cudaMemcpyDeviceToHost);


      dfloat normAq = 0;

      for(int n=0;n<E*BSIZE;++n)
        normAq += (h_Aq[n]*h_Aq[n]);
      normAq = sqrt(normAq);

      printf("CUDA-CUBLAS: error Aq = %17.15lf\n", normAq);
#if 0
  //6) use the original occa kernel to check for correctness
occa::device device;
 occa::kernel Tet3Dkernel;  

occa::kernel Tet3Dkernel[8];
  occa::kernel correctRes;
  occa::kernelInfo kernelInfo;
int p_Nb = 1;
int p_Ne = 1;

 kernelInfo.addDefine("dfloat", dfloatString);
  kernelInfo.addDefine("p_N", p_N);
  kernelInfo.addDefine("p_Np", p_Np);
  kernelInfo.addDefine("p_Nfp", p_Nfp);
  kernelInfo.addDefine("p_Nfaces", p_Nfaces);
  kernelInfo.addDefine("p_NfacesNfp", p_Nfaces*p_Nfp);
  kernelInfo.addDefine("BSIZE", (int)BSIZE);
  kernelInfo.addDefine("p_Nvgeo", 11);
  kernelInfo.addDefine("p_Nsgeo", 7);
  kernelInfo.addDefine("p_Nggeo", 7);
  kernelInfo.addDefine("p_Ne", p_Ne);
  kernelInfo.addDefine("p_Nb", p_Nb);

    kernelInfo.addDefine("p_G00ID", 0);
    kernelInfo.addDefine("p_G01ID", 1);
    kernelInfo.addDefine("p_G02ID", 2);
    kernelInfo.addDefine("p_G11ID", 3);
    kernelInfo.addDefine("p_G12ID", 4);
    kernelInfo.addDefine("p_G22ID", 5);
    kernelInfo.addDefine("p_GWJID", 6);


    int  *elementList;
dfloat *Aq;
    occa::memory o_elementList, o_ggeo, o_SrrT, o_SrsT, 
o_SrtT, o_SsrT,o_SssT,o_SstT, o_StrT, o_StsT, o_SttT, o_MM, o_q, o_Aq; 

  //6.1) allocate the data for occa 


o_ggeo = device.malloc(E*7*sizeof(dfloat), h_ggeo);
o_SrrT = device.malloc(p_Np*p_Np*sizeof(dfloat), SrrT);
o_SrsT = device.malloc(p_Np*p_Np*sizeof(dfloat), SrsT); 
o_SrtT = device.malloc(p_Np*p_Np*sizeof(dfloat), SrtT); 
o_SssT = device.malloc(p_Np*p_Np*sizeof(dfloat), SssT); 
o_SstT = device.malloc(p_Np*p_Np*sizeof(dfloat), SstT); 
o_SttT = device.malloc(p_Np*p_Np*sizeof(dfloat), SttT); 
o_MM = device.malloc(p_Np*p_Np*sizeof(dfloat), MM); 

  o_q = device.malloc(p_Np*E*sizeof(dfloat), h_q);
o_Aq =  device.malloc(p_Np*E*sizeof(dfloat), h_q);

Aq = (dfloat *) calloc(E*p_Np,sizeof(dfloat))
 elementList = (int *) calloc(E, sizeof(int));


    for(int j=0;j<E;++j){
      elementList[j] = j;
    }
    o_elementList = device.malloc(E*sizeof(int), elementList);
    char buf[200];
      sprintf(buf, "ellipticPartialAxTet3D_Ref3"); 
      Tet3Dkernel = device.buildKernelFromSource("ellipticAxTet3D.okl", buf, kernelInfo);
    
    occa::initTimer(device);



//6.2) run the kernel

 Tet3Dkernel(E, o_elementList,
            o_ggeo,
            o_SrrT,
            o_SrsT,
            o_SrtT,
            o_SsrT,
            o_SssT,
            o_SstT,
            o_StrT,
            o_StsT,
            o_SttT,
            o_MM,
            lambda,
            o_q,       
            o_Aq);


  //6.3) copy the data back to CPU
  o_Aq.copyTo(Aq);
cudaMemcpy(*h_Aq, *d_Aq, E*p_Np*sizeof(dfloat), cudaMemcpyDeviceToHost);


      datafloat normAq = 0;

      for(int n=0;n<E*BSIZE;++n)
        normAq += (Aq[n]*Aq[n]-h_Aq[n]);
      normAq = sqrt(normAq);

      printf("OCCA: error Aq = %17.15lf\n", normAq);


  //6.4) compute a norm of the difference
#endif
  // Destroy the handle
  cublasDestroy(handle);

  exit(0);
  return 0;

}


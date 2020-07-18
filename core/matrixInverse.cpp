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

#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

extern "C" {
  void dgesv_ ( int     *N, int     *NRHS, double  *A,
                int     *LDA,
                int     *IPIV,
                double  *B,
                int     *LDB,
                int     *INFO );

  void sgesv_(int *N, int *NRHS,float  *A, int *LDA, int *IPIV, float  *B, int *LDB,int *INFO);

  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

  void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);
  void sgetri_(int* N, float* A, int* lda, int* IPIV, float* WORK, int* lwork, int* INFO);

  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
              double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );

  double dlange_(char *NORM, int *M, int *N, double *A, int *LDA, double *WORK);
  void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM,
                double *RCOND, double *WORK, int *IWORK, int *INFO );

  void dgesvd_( char *JOBU, char *JOBVT, int *M, int *N, double *A,
                int *LDA, double *S, double *U, int *LDU, double *VT, int* LDVT,
                double *WORK, int *LWORK, int *INFO );

  void sgesvd_( char *JOBU, char *JOBVT, int *M, int *N, float *A,
                int *LDA, float *S, float *U, int *LDU, float *VT, int *LDVT,
                float *WORK, int *LWORK, int *INFO);
}

void matrixInverse(int N, double *A){
  int lwork = N*N;
  int info;

  // compute inverse mass matrix
  double *tmpInvA = (double*) calloc(N*N, sizeof(double));

  int *ipiv = (int*) calloc(N, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));

  for(int n=0;n<N*N;++n){
    tmpInvA[n] = A[n];
  }

  dgetrf_ (&N, &N, tmpInvA, &N, ipiv, &info);
  dgetri_ (&N, tmpInvA, &N, ipiv, work, &lwork, &info);

  if(info)
    printf("inv: dgetrf/dgetri reports info = %d when inverting matrix\n", info);

  for(int n=0;n<N*N;++n)
    A[n] = tmpInvA[n];

  free(work);
  free(ipiv);
  free(tmpInvA);
}

void matrixInverse(int N, float *A){
  int lwork = N*N;
  int info;

  // compute inverse mass matrix
  float *tmpInvA = (float*) calloc(N*N, sizeof(float));

  int *ipiv = (int*) calloc(N, sizeof(int));
  float *work = (float*) calloc(lwork, sizeof(float));

  for(int n=0;n<N*N;++n){
    tmpInvA[n] = A[n];
  }

  //NC: are we missing these?
  // sgetrf_ (&N, &N, tmpInvA, &N, ipiv, &info);
  // sgetri_ (&N, tmpInvA, &N, ipiv, work, &lwork, &info);
  info =1; //NC: throw an error for now

  if(info)
    printf("inv: sgetrf/sgetri reports info = %d when inverting matrix\n", info);

  for(int n=0;n<N*N;++n)
    A[n] = tmpInvA[n];

  free(work);
  free(ipiv);
  free(tmpInvA);
}

// Assumes row major ordering
void matrixPInverse(int M, int N, double *A){
#if 1
  double *tmpA = (double *)calloc(M*N, sizeof(double));

  for(int m=0; m<M*N; m++)
    tmpA[m] = A[m]; 

double *V = (double *) malloc(M*N*sizeof(double));
// First compute  inv(trans(A)*A)
for(int n=0; n<N; n++){
  for(int m=0; m<N; m++){
    double tmp = 0; 
    for(int i = 0; i<M; i++){
     tmp += tmpA[i*N +n]*tmpA[i*N +m];
    }
    V[n*N + m] = tmp; 
  }
}
matrixInverse(N, V); 

// (P^T*P)^-1 * P^T
for(int n=0; n<N; n++){
  for(int m=0; m<M; m++){
    double tmp = 0; 
    for(int i=0; i<N; i++){
    tmp += V[n*N + i]*tmpA[m*N + i];
    }
    A[n*M + m] = tmp; 
  }
}

free(V); free(tmpA); 

#else  

// transpose
double *tmpA = (double *)calloc(M*N, sizeof(double));
  for(int m=0; m<M; m++){
    for(int n=0; n<N; n++){
    tmpA[n*M + m] = A[m*N + n]; 
  }
}
  int INFO;  
  char opt = 'A'; 
  int LWORK  = M*N;
  int minN = M>N ? N : M; 

  double *S    = (double *)calloc(minN, sizeof(double)); 
  double *UT   = (double *)calloc(M*M, sizeof(double)); 
  double *V    = (double *)calloc(N*N, sizeof(double)); 
  double *WORK = (double *)calloc(LWORK, sizeof(double)); 

  dgesvd_(&opt,&opt, &M, &N, tmpA, &M, S, UT, &M, V, &N, WORK, &LWORK, &INFO);
  
  if(INFO)
    printf("Singular Value Decomp: dgesvd reports info = %d when decomposing matrix\n", INFO);

  double *SM = (double *)calloc(M*N, sizeof(double));
  double *VS = (double *)calloc(M*N, sizeof(double));
  
  // SM = transpose(1.0/S(S>0))
  for(int n=0; n<minN; n++){
    SM[n*M +n]= 1.0/S[n];
  }

  // Pinv = V*transpose(S^-1)*transpose(U)
   for(int n=0; n<N; n++){
    for(int m=0; m<M; m++){
      double tmp = 0; 
      for(int i=0; i<N; i++){
          tmp += V[n*N + i]*SM[i*M + m]; 
      }
      VS[n*M + m] = tmp; 
    }
  }

  for(int n=0; n<N; n++){
    for(int m=0; m<M; m++){
      double tmp = 0; 
      for(int i=0; i<M; i++){
          tmp += VS[n*M + i]*UT[i*M + m]; 
      }
      A[n*M + m] = tmp; 
    }
  }

free(S); free(WORK); free(V); free(UT); 
free(tmpA); free(SM); free(VS); 
#endif
  
}



// Not tested yet
void matrixPInverse(int M, int N, float *A){
// // transpose
// float *tmpA = (float *)calloc(M*N, sizeof(float));
//   for(int m=0; m<M; m++){
//     for(int n=0; n<N; n++){
//     tmpA[n*M + m] = A[m*N + n]; 
//   }
// }
//   int INFO;  
//   char opt = 'A'; 
//   int LWORK  = M*N;
//   int minN = M>N ? N : M; 

//   float *S    = (float *)calloc(minN, sizeof(float)); 
//   float *UT   = (float *)calloc(M*M, sizeof(float)); 
//   float *V    = (float *)calloc(N*N, sizeof(float)); 
//   float *WORK = (float *)calloc(LWORK, sizeof(float)); 

//   // sgesvd_(&opt,&opt, &M, &N, tmpA, &M, S, UT, &M, V, &N, WORK, &LWORK, &INFO); !!!!
  
//   if(INFO)
//     printf("Singular Value Decomp: dgesvd reports info = %d when decomposing matrix\n", INFO);

//   float *SM = (float *)calloc(M*N, sizeof(float));
//   float *VS = (float *)calloc(M*N, sizeof(float));
  
//   // SM = transpose(1.0/S(S>0))
//   for(int n=0; n<minN; n++){
//     SM[n*M +n]= 1.0/S[n];
//   }

//   // Pinv = V*transpose(S^-1)*transpose(U)
//    for(int n=0; n<N; n++){
//     for(int m=0; m<M; m++){
//       float tmp = 0; 
//       for(int i=0; i<N; i++){
//           tmp += V[n*N + i]*SM[i*M + m]; 
//       }
//       VS[n*M + m] = tmp; 
//     }
//   }

//   for(int n=0; n<N; n++){
//     for(int m=0; m<M; m++){
//       float tmp = 0; 
//       for(int i=0; i<M; i++){
//           tmp += VS[n*M + i]*UT[i*M + m]; 
//       }
//       A[n*M + m] = tmp; 
//     }
//   }

// free(S); free(WORK); free(V); free(UT); 
// free(tmpA); free(SM); free(VS); 
  
}
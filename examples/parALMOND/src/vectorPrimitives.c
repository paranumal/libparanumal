#include "agmg.h"

dfloat norm(iint n, dfloat *a){
  dfloat result = 0.;
  #pragma omp parallel for reduction(+:result)
  for(iint i=0; i<n; i++){
    result += a[i]*a[i];
  }
  return sqrt(result);
}

dfloat innerProd(iint n, dfloat *a, dfloat *b){
  dfloat result = 0.;
  #pragma omp parallel for reduction(+:result)
  for(iint i=0; i<n; i++)
    result += a[i]*b[i];
  return result;
}

void doubleInnerProd(iint n, dfloat *aDotbc, dfloat *a, dfloat *b, dfloat *c) {
  dfloat aDotb = 0.;
  dfloat aDotc = 0.;
  #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc)
  for(iint i=0; i<n; i++) {
    aDotb += a[i]*b[i];
    aDotc += a[i]*c[i];
  }
  aDotbc[0] = aDotb;
  aDotbc[1] = aDotc;
}

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(iint n, dfloat *aDotbc, dfloat *a, dfloat *b, dfloat *c) {
  dfloat aDotb = 0.;
  dfloat aDotc = 0.;
  dfloat bDotb = 0.;
  #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:bDotb)
  for(iint i=0; i<n; i++) {
    aDotb += a[i]*b[i];
    aDotc += a[i]*c[i];
    bDotb += b[i]*b[i];
  }
  aDotbc[0] = aDotb;
  aDotbc[1] = aDotc;
  aDotbc[2] = bDotb;
}

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(iint n, dfloat *aDotbcd, dfloat *a, dfloat *b, dfloat *c, dfloat* d) {
  dfloat aDotb = 0.;
  dfloat aDotc = 0.;
  dfloat aDotd = 0.;
  #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:aDotd)
  for(iint i=0; i<n; i++) {
    aDotb += a[i]*b[i];
    aDotc += a[i]*c[i];
    aDotd += a[i]*d[i];
  }
  aDotbcd[0] = aDotb;
  aDotbcd[1] = aDotc;
  aDotbcd[2] = aDotd;
}

// y = beta*y + alpha*x
void vectorAdd(iint n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y){
  #pragma omp parallel for
  for(iint i=0; i<n; i++)
    y[i] = beta*y[i] + alpha*x[i];
}

// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(iint n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y){
  dfloat result = 0.;
  #pragma omp parallel for reduction(+:result)
  for(iint i=0; i<n; i++) {
    y[i] = beta*y[i] + alpha*x[i];
    result += y[i]*y[i];
  }
  return result;
}

void dotStar(iint m, dfloat *a, dfloat *b){
  #pragma omp parallel for
  for(iint i=0; i<m; i++)
    b[i] *= a[i];
}

void scaleVector(iint m, dfloat *a, dfloat alpha){
  #pragma omp parallel for
  for(iint i=0; i<m; i++)
    a[i] *= alpha;
}

void setVector(iint m, dfloat *a, dfloat alpha){
  #pragma omp parallel for
  for(iint i=0; i<m; i++)
    a[i] = alpha;
}

void addScalar(iint m, dfloat alpha, dfloat *a){
  #pragma omp parallel for
  for(iint i=0; i<m; i++)
    a[i] += alpha;
}

dfloat sumVector(iint m, dfloat *a){
  dfloat alpha = 0.;

  #pragma omp parallel for reduction(+:alpha)
  for (iint i=0; i<m; i++) {
    alpha += a[i];
  }
  return alpha;
}
    
void randomize(iint m, dfloat *a){
  for(iint i=0; i<m; i++)
    a[i] = (dfloat) drand48();
}

dfloat maxEntry(iint n, dfloat *a){
  if(n == 0)
    return 0;

  dfloat maxVal = 0.;
  //  #pragma omp parallel for reduction(max:maxVal)
  for(iint i=0; i<n; i++){
    dfloat a2 = (a[i] < 0) ? -a[i] : a[i];
    if(maxVal < a2){
      maxVal = a2;
    }
  }
  return maxVal;
}

#define RDIMX 32
#define RDIMY 8
#define RLOAD 1

void scaleVector(parAlmond_t *parAlmond, iint N, occa::memory o_a, dfloat alpha){
  if (N) parAlmond->scaleVectorKernel(N, alpha, o_a);
}

void setVector(parAlmond_t *parAlmond, iint N, occa::memory o_a, dfloat alpha){
  if (N) parAlmond->setVectorKernel(N, alpha, o_a);
}

dfloat sumVector(parAlmond_t *parAlmond, iint N, occa::memory o_a){
  iint numBlocks = ((N+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  if(!numBlocks) numBlocks = 1;

  dfloat alpha =0., zero = 0.;

  parAlmond->setVectorKernel(1, zero, parAlmond->o_rho);
  if (N) parAlmond->sumVectorKernel(numBlocks,N,o_a,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(&alpha,1*sizeof(dfloat),0);

  return alpha;
}

void addScalar(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_a){
  if (N) parAlmond->addScalarKernel(N, alpha, o_a);
}

void dotStar(parAlmond_t *parAlmond, iint N, occa::memory o_a, occa::memory o_b){
  if (N) parAlmond->simpleDotStarKernel(N, o_a, o_b);
}

void dotStar(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_a,
	           occa::memory o_b, dfloat beta, occa::memory o_c){
  if (N) parAlmond->dotStarKernel(N, alpha, beta, o_a, o_b, o_c);
}
dfloat innerProd(parAlmond_t *parAlmond, iint N,
                  occa::memory o_x, occa::memory o_y){

  iint numBlocks = ((N+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  if(!numBlocks) numBlocks = 1;

#if 0
  dfloat result =0.;

  parAlmond->o_rho.copyFrom(&result,1*sizeof(dfloat));
  parAlmond->innerProdKernel(numBlocks,N,o_x,o_y,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(&result,1*sizeof(dfloat));

#else
  dfloat zero = 0, result;
  parAlmond->setVectorKernel(1, zero, parAlmond->o_rho);
  parAlmond->innerProdKernel(numBlocks,N,o_x,o_y,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(&result,1*sizeof(dfloat),0);
#endif
  return result;
}

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(parAlmond_t *parAlmond, iint N, dfloat *aDotbc, occa::memory o_a,
                                        occa::memory o_b, occa::memory o_c) {
  iint numBlocks = ((N+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  if(!numBlocks) numBlocks = 1;
#if 0
  aDotbc[0] = 0.;
  aDotbc[1] = 0.;
  aDotbc[2] = 0.;

  parAlmond->o_rho.copyFrom(aDotbc);
  parAlmond->kcycleCombinedOp1Kernel(numBlocks,N,o_a,o_b,o_c,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(aDotbc);
#else
  dfloat zero = 0;
  parAlmond->setVectorKernel(3, zero, parAlmond->o_rho);
  parAlmond->kcycleCombinedOp1Kernel(numBlocks,N,o_a,o_b,o_c,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(aDotbc);
#endif

}

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(parAlmond_t *parAlmond, iint N, dfloat *aDotbcd, occa::memory o_a,
                                              occa::memory o_b, occa::memory o_c, occa::memory o_d) {

  iint numBlocks = ((N+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  if(!numBlocks) numBlocks = 1;

#if 0
  aDotbcd[0] = 0.;
  aDotbcd[1] = 0.;
  aDotbcd[2] = 0.;

  parAlmond->o_rho.copyFrom(aDotbcd);
  parAlmond->kcycleCombinedOp2Kernel(numBlocks,N,o_a,o_b,o_c,o_d,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(aDotbcd);
#else
  dfloat zero = 0;
  parAlmond->setVectorKernel(3, zero, parAlmond->o_rho);
  parAlmond->kcycleCombinedOp2Kernel(numBlocks,N,o_a,o_b,o_c,o_d,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(aDotbcd);
#endif
}

// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_x,
                                                          dfloat beta, occa::memory o_y){
  iint numBlocks = ((N+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  if(!numBlocks) numBlocks = 1;

  dfloat zero = 0, result;
#if 0
  parAlmond->o_rho.copyFrom(&result,1*sizeof(dfloat));
  parAlmond->vectorAddInnerProdKernel(numBlocks,N,alpha,beta,o_x,o_y,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(&result,1*sizeof(dfloat));
#else
  parAlmond->setVectorKernel(1, zero, parAlmond->o_rho);
  parAlmond->vectorAddInnerProdKernel(numBlocks,N,alpha,beta,o_x,o_y,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(&result,1*sizeof(dfloat),0);
#endif
  return result;
}


void vectorAdd(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y){
  parAlmond->vectorAddKernel(N, alpha, beta, o_x, o_y);
}

void vectorAdd(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_x,
	 dfloat beta, occa::memory o_y, occa::memory o_z){
  parAlmond->vectorAddKernel2(N, alpha, beta, o_x, o_y, o_z);
}

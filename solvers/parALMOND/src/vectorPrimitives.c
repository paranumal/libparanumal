/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "agmg.h"

dfloat norm(dlong n, dfloat *a){
  dfloat result = 0.;
  #pragma omp parallel for reduction(+:result)
  for(dlong i=0; i<n; i++){
    result += a[i]*a[i];
  }
  return sqrt(result);
}

dfloat innerProd(dlong n, dfloat *a, dfloat *b){
  dfloat result = 0.;
  #pragma omp parallel for reduction(+:result)
  for(dlong i=0; i<n; i++)
    result += a[i]*b[i];
  return result;
}

void doubleInnerProd(dlong n, dfloat *aDotbc, dfloat *a, dfloat *b, dfloat *c) {
  dfloat aDotb = 0.;
  dfloat aDotc = 0.;
  #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc)
  for(dlong i=0; i<n; i++) {
    aDotb += a[i]*b[i];
    aDotc += a[i]*c[i];
  }
  aDotbc[0] = aDotb;
  aDotbc[1] = aDotc;
}

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(dlong n, dfloat *aDotbc, dfloat *a, 
                      dfloat *b, dfloat *c, dfloat* w, bool weighted) {
  dfloat aDotb = 0.;
  dfloat aDotc = 0.;
  dfloat bDotb = 0.;
  if (weighted) {
    #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:bDotb)
    for(dlong i=0; i<n; i++) {
      aDotb += w[i]*a[i]*b[i];
      aDotc += w[i]*a[i]*c[i];
      bDotb += w[i]*b[i]*b[i];
    }
  } else {
    #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:bDotb)
    for(dlong i=0; i<n; i++) {
      aDotb += a[i]*b[i];
      aDotc += a[i]*c[i];
      bDotb += b[i]*b[i];
    }
  }
  aDotbc[0] = aDotb;
  aDotbc[1] = aDotc;
  aDotbc[2] = bDotb;
}

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(dlong n, dfloat *aDotbcd, dfloat *a, dfloat *b, 
                        dfloat *c, dfloat* d, dfloat *w, bool weighted) {
  dfloat aDotb = 0.;
  dfloat aDotc = 0.;
  dfloat aDotd = 0.;
  if (weighted) {
    #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:aDotd)
    for(dlong i=0; i<n; i++) {
      aDotb += w[i]*a[i]*b[i];
      aDotc += w[i]*a[i]*c[i];
      aDotd += w[i]*a[i]*d[i];
    }  
  } else {
    #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:aDotd)
    for(dlong i=0; i<n; i++) {
      aDotb += a[i]*b[i];
      aDotc += a[i]*c[i];
      aDotd += a[i]*d[i];
    }  
  }
  
  aDotbcd[0] = aDotb;
  aDotbcd[1] = aDotc;
  aDotbcd[2] = aDotd;
}

// y = beta*y + alpha*x
void vectorAdd(dlong n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y){
  #pragma omp parallel for
  for(dlong i=0; i<n; i++)
    y[i] = beta*y[i] + alpha*x[i];
}

// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(dlong n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y,
                          dfloat *w, bool weighted){
  dfloat result = 0.;
  if (weighted) {
    #pragma omp parallel for reduction(+:result)
    for(dlong i=0; i<n; i++) {
      y[i] = beta*y[i] + alpha*x[i];
      result += w[i]*y[i]*y[i];
    }  
  } else {
    #pragma omp parallel for reduction(+:result)
    for(dlong i=0; i<n; i++) {
      y[i] = beta*y[i] + alpha*x[i];
      result += y[i]*y[i];
    }  
  }
  return result;
}

void dotStar(dlong m, dfloat *a, dfloat *b){
  #pragma omp parallel for
  for(dlong i=0; i<m; i++)
    b[i] *= a[i];
}

void scaleVector(dlong m, dfloat *a, dfloat alpha){
  #pragma omp parallel for
  for(dlong i=0; i<m; i++)
    a[i] *= alpha;
}

void setVector(dlong m, dfloat *a, dfloat alpha){
  #pragma omp parallel for
  for(dlong i=0; i<m; i++)
    a[i] = alpha;
}

void addScalar(dlong m, dfloat alpha, dfloat *a){
  #pragma omp parallel for
  for(dlong i=0; i<m; i++)
    a[i] += alpha;
}

dfloat sumVector(dlong m, dfloat *a){
  dfloat alpha = 0.;

  #pragma omp parallel for reduction(+:alpha)
  for (dlong i=0; i<m; i++) {
    alpha += a[i];
  }
  return alpha;
}
    
void randomize(dlong m, dfloat *a){
  for(dlong i=0; i<m; i++)
    a[i] = (dfloat) drand48();
}

dfloat maxEntry(dlong n, dfloat *a){
  if(n == 0)
    return 0;

  dfloat maxVal = 0.;
  //  #pragma omp parallel for reduction(max:maxVal)
  for(dlong i=0; i<n; i++){
    dfloat a2 = (a[i] < 0) ? -a[i] : a[i];
    if(maxVal < a2){
      maxVal = a2;
    }
  }
  return maxVal;
}




void scaleVector(parAlmond_t *parAlmond, dlong N, occa::memory o_a, dfloat alpha){
  if (N) parAlmond->scaleVectorKernel(N, alpha, o_a);
}

void setVector(parAlmond_t *parAlmond, dlong N, occa::memory o_a, dfloat alpha){
  if (N) parAlmond->setVectorKernel(N, alpha, o_a);
}

dfloat sumVector(parAlmond_t *parAlmond, dlong N, occa::memory o_a){
  dlong numBlocks = ((N+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  if(!numBlocks) numBlocks = 1;

  if (N) parAlmond->sumVectorKernel(numBlocks,N,o_a,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(parAlmond->rho,numBlocks*sizeof(dfloat),0);
  
  dfloat alpha =0.;
  #pragma omp parallel for reduction(+:alpha)
  for (dlong i=0; i<numBlocks; i++) {
    alpha += parAlmond->rho[i];
  }

  return alpha;
}

void addScalar(parAlmond_t *parAlmond, dlong N, dfloat alpha, occa::memory o_a){
  if (N) parAlmond->addScalarKernel(N, alpha, o_a);
}

void dotStar(parAlmond_t *parAlmond, dlong N, occa::memory o_a, occa::memory o_b){
  if (N) parAlmond->simpleDotStarKernel(N, o_a, o_b);
}

void dotStar(parAlmond_t *parAlmond, dlong N, dfloat alpha, occa::memory o_a,
	           occa::memory o_b, dfloat beta, occa::memory o_c){
  if (N) parAlmond->dotStarKernel(N, alpha, beta, o_a, o_b, o_c);
}

dfloat innerProd(parAlmond_t *parAlmond, dlong N,
                  occa::memory o_x, occa::memory o_y){

  dlong numBlocks = ((N+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  if(!numBlocks) numBlocks = 1;

  parAlmond->innerProdKernel(numBlocks,N,o_x,o_y,parAlmond->o_rho);
  parAlmond->o_rho.copyTo(parAlmond->rho,numBlocks*sizeof(dfloat),0);
  
  dfloat result =0.;
  #pragma omp parallel for reduction(+:result)
  for (dlong i=0; i<numBlocks; i++) {
    result += parAlmond->rho[i];
  }

  return result;
}

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(parAlmond_t *parAlmond, dlong N, dfloat *aDotbc, 
                        occa::memory o_a, occa::memory o_b, 
                        occa::memory o_c, occa::memory o_w, bool weighted) {
  dlong numBlocks = ((N+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  if(!numBlocks) numBlocks = 1;

  if (weighted) {
    parAlmond->kcycleWeightedCombinedOp1Kernel(numBlocks,N,o_a,o_b,o_c,o_w,parAlmond->o_rho);
  } else {
    parAlmond->kcycleCombinedOp1Kernel(numBlocks,N,o_a,o_b,o_c,parAlmond->o_rho);
  }
  parAlmond->o_rho.copyTo(parAlmond->rho,3*numBlocks*sizeof(dfloat),0);
  
  dfloat aDotb = 0., aDotc = 0., bDotb = 0.;
  #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:bDotb)
  for(dlong i=0; i<numBlocks; i++) {
    aDotb += parAlmond->rho[3*i+0];
    aDotc += parAlmond->rho[3*i+1];
    bDotb += parAlmond->rho[3*i+2];
  }  

  aDotbc[0] = aDotb;
  aDotbc[1] = aDotc;
  aDotbc[2] = bDotb;
}

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(parAlmond_t *parAlmond, dlong N, dfloat *aDotbcd, 
                        occa::memory o_a, occa::memory o_b, 
                        occa::memory o_c, occa::memory o_d,
                        occa::memory o_w, bool weighted) {

  dlong numBlocks = ((N+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  if(!numBlocks) numBlocks = 1;

  if (weighted) {
    parAlmond->kcycleWeightedCombinedOp2Kernel(numBlocks,N,o_a,o_b,o_c,o_d,o_w,parAlmond->o_rho);
  } else {
    parAlmond->kcycleCombinedOp2Kernel(numBlocks,N,o_a,o_b,o_c,o_d,parAlmond->o_rho);
  }
  parAlmond->o_rho.copyTo(parAlmond->rho,3*numBlocks*sizeof(dfloat),0);
  
  dfloat aDotb = 0., aDotc = 0., aDotd = 0.;
  #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:aDotd)
  for(dlong i=0; i<numBlocks; i++) {
    aDotb += parAlmond->rho[3*i+0];
    aDotc += parAlmond->rho[3*i+1];
    aDotd += parAlmond->rho[3*i+2];
  }  

  aDotbcd[0] = aDotb;
  aDotbcd[1] = aDotc;
  aDotbcd[2] = aDotd;
}

// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(parAlmond_t *parAlmond, dlong N, dfloat alpha, occa::memory o_x,
                                                          dfloat beta, occa::memory o_y,
                                                          occa::memory o_w, bool weighted){
  dlong numBlocks = ((N+RDIMX*RDIMY-1)/(RDIMX*RDIMY))/RLOAD;
  if(!numBlocks) numBlocks = 1;

  if (weighted) {
    parAlmond->vectorAddWeightedInnerProdKernel(numBlocks,N,alpha,beta,o_x,o_y,o_w,parAlmond->o_rho);
  } else {
    parAlmond->vectorAddInnerProdKernel(numBlocks,N,alpha,beta,o_x,o_y,parAlmond->o_rho);
  }
  parAlmond->o_rho.copyTo(parAlmond->rho,numBlocks*sizeof(dfloat),0);
  
  dfloat result =0.;
  #pragma omp parallel for reduction(+:result)
  for (dlong i=0; i<numBlocks; i++) {
    result += parAlmond->rho[i];
  }

  return result;
}


void vectorAdd(parAlmond_t *parAlmond, dlong N, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y){
  parAlmond->vectorAddKernel(N, alpha, beta, o_x, o_y);
}

void vectorAdd(parAlmond_t *parAlmond, dlong N, dfloat alpha, occa::memory o_x,
	 dfloat beta, occa::memory o_y, occa::memory o_z){
  parAlmond->vectorAddKernel2(N, alpha, beta, o_x, o_y, o_z);
}

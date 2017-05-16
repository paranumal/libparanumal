#include "parAlmond.h"


dfloat norm(iint n, dfloat *a){

  occa::tic("norm");
  dfloat result = 0.;
  for(iint i=0; i<n; i++){
    result += a[i]*a[i];
  }

  result = sqrt(result);
  occa::toc("norm");
  return result;
}


dfloat innerProd(iint n, dfloat *a, dfloat *b){
  
  occa::tic("innerProd");
  dfloat result = 0.;
  for(iint i=0; i<n; i++)
    result += a[i]*b[i];
  occa::toc("innerProd");

  return result;
}


// y = beta*y + alpha*x
void vectorAdd(iint n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y){
  occa::tic("vectorAdd");
  for(iint i=0; i<n; i++)
    y[i] = beta*y[i] + alpha*x[i];
  occa::toc("vectorAdd");
}

void dotStar(iint m, dfloat *a, dfloat *b){
  occa::tic("dotStar");
  for(iint i=0; i<m; i++)
    b[i] *= a[i];
  occa::toc("dotStar");
}

void scaleVector(iint m, dfloat *a, dfloat alpha){
  occa::tic("scaleVector");
  for(iint i=0; i<m; i++)
    a[i] *= alpha;
  occa::toc("scaleVector");
}

void randomize(iint m, dfloat *a){
  for(iint i=0; i<m; i++)
    a[i] = (dfloat) drand48();
}


dfloat maxEntry(iint n, dfloat *a){
  if(n == 0)
    return 0;

  dfloat maxVal = 0.;

  for(iint i=0; i<n; i++){
    dfloat a2 = (a[i] < 0) ? -a[i] : a[i];

    if(maxVal < a2){
      maxVal = a2;
    }
  }

  return maxVal;
}

void copyVector(parAlmond_t *parAlmond, iint N, occa::memory o_a, occa::memory o_b){

  occaTimerTic(parAlmond->device,"copyKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;

  parAlmond->copyKernel(numBlocks, AGMGBDIM, N, o_a, o_b);
  occaTimerToc(parAlmond->device,"copyKernel");
}


void scaleVector(parAlmond_t *parAlmond, iint N, occa::memory o_a, dfloat alpha){

  occaTimerTic(parAlmond->device,"scaleKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;

  parAlmond->scaleVectorKernel(numBlocks, AGMGBDIM, N, alpha, o_a);
  occaTimerToc(parAlmond->device,"scaleKernel");
}

void dotStar(parAlmond_t *parAlmond, iint N, occa::memory o_a, occa::memory o_b){

  occaTimerTic(parAlmond->device,"dotStarKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;

  parAlmond->simpleDotStarKernel(numBlocks, AGMGBDIM, N, o_a, o_b);
  occaTimerToc(parAlmond->device,"dotStarKernel");
}

void dotStar(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_a,
	           occa::memory o_b, dfloat beta, occa::memory o_c){

  occaTimerTic(parAlmond->device,"dotStarKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;

  parAlmond->dotStarKernel(numBlocks, AGMGBDIM, N, alpha, beta, o_a, o_b, o_c);
  occaTimerToc(parAlmond->device,"dotStarKernel");
}

dfloat innerProd(parAlmond_t *parAlmond, iint N, occa::memory o_a, occa::memory o_b){

  occaTimerTic(parAlmond->device,"innerProdKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;

  dfloat *hostRed = (dfloat *) calloc(numBlocks, sizeof(dfloat));

  occa::memory o_dRed =
    parAlmond->device.malloc(numBlocks*sizeof(dfloat), hostRed);

  parAlmond->partialInnerProdKernel(numBlocks, AGMGBDIM, N, o_a, o_b, o_dRed);

  o_dRed.copyTo(hostRed);

  dfloat result = 0.;
  for(iint i=0; i<numBlocks; i++)
    result += hostRed[i];

  free(hostRed);
  o_dRed.free();
  occaTimerToc(parAlmond->device,"innerProdKernel");

  return result;
}

void vectorAdd(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y){

  occaTimerTic(parAlmond->device,"vectorAddKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;
  
  parAlmond->vectorAddKernel(numBlocks, AGMGBDIM, N, alpha, beta, o_x, o_y);
  occaTimerToc(parAlmond->device,"vectorAddKernel");
}

void vectorAdd(parAlmond_t *parAlmond, iint N, dfloat alpha, occa::memory o_x,
	 dfloat beta, occa::memory o_y, occa::memory o_z){

  occaTimerTic(parAlmond->device,"vectorAddKernel2");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;
  
  parAlmond->vectorAddKernel2(numBlocks, AGMGBDIM, N, alpha, beta, o_x, o_y, o_z);
  occaTimerToc(parAlmond->device,"vectorAddKernel2");
}

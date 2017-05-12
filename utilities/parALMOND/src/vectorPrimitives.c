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

void copyVector(almond_t *almond, iint N, occa::memory o_a, occa::memory o_b){

  occaTimerTic(almond->device,"copyKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;

  almond->copyKernel(numBlocks, AGMGBDIM, N, o_a, o_b);
  occaTimerToc(almond->device,"copyKernel");
}


void scaleVector(almond_t *almond, iint N, occa::memory o_a, dfloat alpha){

  occaTimerTic(almond->device,"scaleKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;

  almond->scaleVectorKernel(numBlocks, AGMGBDIM, N, alpha, o_a);
  occaTimerToc(almond->device,"scaleKernel");
}

void dotStar(almond_t *almond, iint N, occa::memory o_a, occa::memory o_b){

  occaTimerTic(almond->device,"dotStarKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;

  almond->simpleDotStarKernel(numBlocks, AGMGBDIM, N, o_a, o_b);
  occaTimerToc(almond->device,"dotStarKernel");
}

void dotStar(almond_t *almond, iint N, dfloat alpha, occa::memory o_a,
	           occa::memory o_b, dfloat beta, occa::memory o_c){

  occaTimerTic(almond->device,"dotStarKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;

  almond->dotStarKernel(numBlocks, AGMGBDIM, N, alpha, beta, o_a, o_b, o_c);
  occaTimerToc(almond->device,"dotStarKernel");
}

dfloat innerProd(almond_t *almond, iint N, occa::memory o_a, occa::memory o_b){

  occaTimerTic(almond->device,"innerProdKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;

  dfloat *hostRed = (dfloat *) calloc(numBlocks, sizeof(dfloat));

  occa::memory o_dRed =
    almond->device.malloc(numBlocks*sizeof(dfloat), hostRed);

  almond->partialInnerProdKernel(numBlocks, AGMGBDIM, N, o_a, o_b, o_dRed);

  o_dRed.copyTo(hostRed);

  dfloat result = 0.;
  for(iint i=0; i<numBlocks; i++)
    result += hostRed[i];

  free(hostRed);
  o_dRed.free();
  occaTimerToc(almond->device,"innerProdKernel");

  return result;
}

void vectorAdd(almond_t *almond, iint N, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y){

  occaTimerTic(almond->device,"vectorAddKernel");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;
  
  almond->vectorAddKernel(numBlocks, AGMGBDIM, N, alpha, beta, o_x, o_y);
  occaTimerToc(almond->device,"vectorAddKernel");
}

void vectorAdd(almond_t *almond, iint N, dfloat alpha, occa::memory o_x,
	 dfloat beta, occa::memory o_y, occa::memory o_z){

  occaTimerTic(almond->device,"vectorAddKernel2");
  const iint numBlocks = (N+AGMGBDIM-1)/AGMGBDIM;
  
  almond->vectorAddKernel2(numBlocks, AGMGBDIM, N, alpha, beta, o_x, o_y, o_z);
  occaTimerToc(almond->device,"vectorAddKernel2");
}

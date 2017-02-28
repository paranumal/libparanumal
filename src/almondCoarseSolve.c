
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include <map>
#include <vector>
#include <almondHeaders.hpp>

#pragma message("WARNING : HARD CODED TO FLOAT/INT\n")

#define dfloat double
#define iint int

typedef struct {
  almond::csr<dfloat> *A;
  almond::agmg<dfloat> M;
  std::vector<dfloat>  nullA;
  std::vector<dfloat>  rhs;
  std::vector<dfloat>  sol;
  uint numLocalRows;
  uint nnz;

  char* iintType;
  char* dfloatType;

} almond_t;

void * almondSetup(uint  numLocalRows, 
		   void* rowIds,
		   uint  nnz, 
		   void* Ai,
		   void* Aj,
		   void* Avals,
		   int   nullSpace,
		   const char* iintType, 
		   const char* dfloatType) {

  int n;
  almond_t *almond = (almond_t*) calloc(1, sizeof(almond_t));

  std::vector<int>    vAj(nnz);
  std::vector<dfloat> vAvals(nnz);
  std::vector<int>    vRowStarts(numLocalRows+1);

  int *iAi = (int*) Ai;
  int *iAj = (int*) Aj;
  dfloat *dAvals = (dfloat*) Avals;

  // assumes presorted
  int cnt = 0;
  for(n=0;n<nnz;++n){
    if(iAi[n]>=numLocalRows || iAj[n]>=numLocalRows)
      printf("errant nonzero %d,%d,%g\n", iAi[n], iAj[n], dAvals[n]);
    if(n==0 || (iAi[n]!=iAi[n-1])){
      //      printf("*\n");
      vRowStarts[cnt] = n;
      ++cnt;
    }else{
      //      printf("\n");
    }
    vAj[n] = iAj[n];
    vAvals[n] = dAvals[n];
  }
  vRowStarts[cnt] = n;
  printf("cnt=%d, numLocalRows=%d\n", cnt, numLocalRows);
  
  almond->A = new almond::csr<dfloat>(vRowStarts, vAj, vAvals);

  almond->rhs.resize(numLocalRows);
  almond->sol.resize(numLocalRows);
  almond->nullA.resize(numLocalRows);
  for (int i=0;i<numLocalRows;i++)almond->nullA[i] = 1;
  
  almond->M.setup(*(almond->A), almond->nullA);
  almond->M.report();
  almond->M.ktype = almond::PCG;

  almond->numLocalRows = numLocalRows;
  
  almond->iintType = strdup(iintType);
  almond->dfloatType = strdup(dfloatType);

  return (void *) almond;
}

int almondSolve(void* x,
		void* A,
		void* rhs) {

  almond_t *almond = (almond_t*) A;
  dfloat *dx = (dfloat*)x;
  dfloat *drhs = (dfloat*)rhs;
  int N = almond->numLocalRows;
  for(iint i=0;i<N;++i) almond->rhs[i] = drhs[i];

  almond->M.solve(almond->rhs, almond->sol);

  for(iint i=0;i<N;++i) dx[i] = almond->sol[i];
  
  return 0;
}

int almondFree(void* A) {
#if 0
  almond_t *almond = (almond_t *) A;

  crs_free(almond->A);

  if (!strcmp(almond->dfloatType,"float")) { 
    free(almond->Avals);
    free(almond->x);  
    free(almond->rhs);
  }

  if (!strcmp(almond->iintType,"int")) { 
    free(almond->rowIds);
  }
#endif
  return 0;
}

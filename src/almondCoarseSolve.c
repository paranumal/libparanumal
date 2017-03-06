
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include <map>
#include <vector>

#include "occa.hpp"
#include "almondHeaders.hpp"
#include "mesh.h"
#include "mpi.h"

#pragma message("WARNING : HARD CODED TO FLOAT/INT\n")

#define amgFloat double

typedef struct {
  almond::csr<amgFloat> *A;
  almond::agmg<amgFloat> M;
  std::vector<amgFloat>  rhs;
  std::vector<amgFloat>  x;
  std::vector<amgFloat>  nullA;

  uint numLocalRows;
  uint Nnum;
  uint recvNnum;
  uint nnz;

  int *globalSortId, *compressId;

  amgFloat *xUnassembled;

  amgFloat *xSort;
  amgFloat *rhsSort;

  char* iintType;
  char* dfloatType;

} almond_t;

void * almondSetup(uint  Nnum,
		   int* rowStarts, 
		   void* rowIds,
		   uint  nnz, 
		   void* Ai,
		   void* Aj,
		   void* Avals,
		   int  *globalSortId, 
		   int  *compressId,
		   int   nullSpace,
		   const char* iintType, 
		   const char* dfloatType) {

  int n;
  almond_t *almond = (almond_t*) calloc(1, sizeof(almond_t));

  int *iAi = (int*) Ai;
  int *iAj = (int*) Aj;
  dfloat *dAvals = (dfloat*) Avals;

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );

  almond->numLocalRows = rowStarts[myid+1]-rowStarts[myid];
  int globalOffset = rowStarts[myid];
  printf("rank %d, numLocalRows %d, offset %d, \n", myid, almond->numLocalRows, globalOffset);


  std::vector<int>    vAj(nnz);
  std::vector<amgFloat> vAvals(nnz);
  std::vector<int>    vRowStarts(almond->numLocalRows+1);

  // assumes presorted
  int cnt = 0;
  for(n=0;n<nnz;++n){
    if(  (iAi[n] >= (almond->numLocalRows + globalOffset)) || (iAj[n] >= (almond->numLocalRows + globalOffset))
	 || (iAi[n] <  globalOffset)                  || (iAj[n] < globalOffset) ) 
      printf("errant nonzero %d,%d,%g, rank %d \n", iAi[n], iAj[n], dAvals[n], myid);
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
  printf("cnt=%d, numLocalRows=%d\n", cnt, almond->numLocalRows);

  almond->Nnum = Nnum;
  //almond->numLocalRows = numLocalRows;
  almond->recvNnum = compressId[almond->numLocalRows];

  almond->globalSortId = (int*) calloc(almond->recvNnum,sizeof(int));
  almond->compressId  = (int*) calloc(almond->numLocalRows+1,sizeof(int));

  for (n=0;n<Nnum;n++) almond->globalSortId[n] = globalSortId[n];
  for (n=0;n<almond->numLocalRows+1;n++) almond->compressId[n] = compressId[n];
  
  almond->xUnassembled = (amgFloat*) calloc(Nnum,sizeof(amgFloat));
  almond->xSort = (amgFloat*) calloc(Nnum,sizeof(amgFloat));
  almond->rhsSort = (amgFloat*) calloc(Nnum,sizeof(amgFloat));

  almond->rhs.resize(almond->numLocalRows);
  almond->x.resize(almond->numLocalRows);
  
  almond->A = new almond::csr<amgFloat>(vRowStarts, vAj, vAvals);

  almond->nullA.resize(almond->numLocalRows);
  for (int i=0;i<almond->numLocalRows;i++)almond->nullA[i] = 1;
  
  almond->M.setup(*(almond->A), almond->nullA);
  almond->M.report();
  almond->M.ktype = almond::PCG;
  
  almond->iintType = strdup(iintType);
  almond->dfloatType = strdup(dfloatType);

  return (void *) almond;
}

int almondSolve(void* x,
		void* A,
		void* rhs) {
  
  almond_t *almond = (almond_t*) A;

  dfloat *dx = (dfloat*) x;
  dfloat *drhs = (dfloat*) rhs;
  
  //sort by globalid 
  for (iint n=0;n<almond->Nnum;n++) 
    almond->rhsSort[n] = drhs[almond->globalSortId[n]];

  //gather
  for (iint n=0;n<almond->numLocalRows;++n){
    almond->rhs[n] = 0.;
    almond->x[n] = 0.;
    for (iint id=almond->compressId[n];id<almond->compressId[n+1];id++) 
      almond->rhs[n] += almond->rhsSort[id];
  }

  if(1){
    almond->M.solve(almond->rhs, almond->x);
  }
  else{
    int maxIt = 40;
    amgFloat tol = 1e-1;
    almond::pcg<amgFloat>(almond->A[0],
			almond->rhs,
			almond->x,
			almond->M,
			maxIt,
			tol);
  }

  //scatter
  for (iint n=0;n<almond->numLocalRows;++n){
    for (iint id = almond->compressId[n];id<almond->compressId[n+1];id++) 
      almond->xSort[id] = almond->x[n]; 
  }

  //sort by to original numbering
  for (iint n=0;n<almond->Nnum;n++) 
    almond->xUnassembled[almond->globalSortId[n]] = almond->xSort[n];

  for(iint i=0;i<almond->Nnum;++i) dx[i] = almond->xUnassembled[i];
  
  return 0;
}

int almondFree(void* A) {
  return 0;
}

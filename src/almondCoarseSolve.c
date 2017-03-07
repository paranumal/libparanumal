
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

  int *sendSortId, *globalSortId, *compressId;
  int *sendCounts, *sendOffsets,  *recvCounts, *recvOffsets;

  amgFloat *xUnassembled;
  amgFloat *rhsUnassembled;

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
		   int    *sendSortId, 
		   int    *globalSortId, 
		   int    *compressId,
		   int    *sendCounts, 
		   int    *sendOffsets, 
		   int    *recvCounts, 
		   int    *recvOffsets,
		   int   nullSpace,
		   const char* iintType, 
		   const char* dfloatType) {

  int n;
  almond_t *almond = (almond_t*) calloc(1, sizeof(almond_t));

  int *iAi = (int*) Ai;
  int *iAj = (int*) Aj;
  dfloat *dAvals = (dfloat*) Avals;

  int num_procs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs );
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );

  almond->numLocalRows = rowStarts[myid+1]-rowStarts[myid];
  int globalOffset = rowStarts[myid];
  printf("rank %d, numLocalRows %d, offset %d, \n", myid, almond->numLocalRows, globalOffset);


  //std::vector<int>    vAj(nnz);
  //std::vector<amgFloat> vAvals(nnz);
  std::vector<int>    vRowStarts(almond->numLocalRows+1);

  // assumes presorted
  int cnt = 0; //nnz counter
  int cnt2 =0; //row start counter
  for(n=0;n<nnz;++n){
    if(  (iAi[n] >= (almond->numLocalRows + globalOffset)) || (iAi[n] < globalOffset))
      printf("errant nonzero %d,%d,%g, rank %d \n", iAi[n], iAj[n], dAvals[n], myid);

    if(cnt2==0 || (iAi[n]!=iAi[n-1])) vRowStarts[cnt2++] = cnt;
      
    if (iAj[n] >= globalOffset && iAj[n] < almond->numLocalRows + globalOffset) cnt++;
  }

  vRowStarts[cnt2] = cnt;
  printf("cnt2=%d, numLocalRows=%d\n", cnt2, almond->numLocalRows);

  std::vector<int>    vAj(cnt);
  std::vector<amgFloat> vAvals(cnt);

  cnt = 0;
  for(n=0;n<nnz;++n){
    if (iAj[n] >= globalOffset && iAj[n] < almond->numLocalRows + globalOffset) {
      vAj[cnt] = iAj[n]-globalOffset;
      vAvals[cnt++] = dAvals[n];  
    }
  }

  almond->Nnum = Nnum;
  //almond->numLocalRows = numLocalRows;
  almond->recvNnum = compressId[almond->numLocalRows];

  almond->sendSortId = (int*) calloc(Nnum,sizeof(int));
  for (n=0;n<Nnum;n++)  almond->sendSortId[n] = sendSortId[n];

  almond->sendCounts = (int*) calloc(num_procs,sizeof(int));
  almond->recvCounts = (int*) calloc(num_procs,sizeof(int));
  for (n=0;n<num_procs;n++){
    almond->sendCounts[n] = sendCounts[n]*sizeof(amgFloat);
    almond->recvCounts[n] = recvCounts[n]*sizeof(amgFloat);
  }

  almond->sendOffsets = (int*) calloc(num_procs+1,sizeof(int));
  almond->recvOffsets = (int*) calloc(num_procs+1,sizeof(int));
  for (n=0;n<num_procs+1;n++){
    almond->sendOffsets[n] = sendOffsets[n]*sizeof(amgFloat);
    almond->recvOffsets[n] = recvOffsets[n]*sizeof(amgFloat);
  }

  almond->globalSortId = (int*) calloc(almond->recvNnum,sizeof(int));
  for (n=0;n<almond->recvNnum;n++) almond->globalSortId[n] = globalSortId[n];

  almond->compressId  = (int*) calloc(almond->numLocalRows+1,sizeof(int));
  for (n=0;n<almond->numLocalRows+1;n++) almond->compressId[n] = compressId[n];

  int Nmax =  (Nnum >almond->recvNnum)? Nnum : almond->recvNnum;
  almond->xUnassembled   = (amgFloat *) malloc(Nmax*sizeof(amgFloat));
  almond->rhsUnassembled = (amgFloat *) malloc(Nmax*sizeof(amgFloat));
  almond->xSort   = (amgFloat *) malloc(Nmax*sizeof(amgFloat));
  almond->rhsSort = (amgFloat *) malloc(Nmax*sizeof(amgFloat));

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

  for (iint n=0;n<almond->Nnum;n++)
    almond->rhsUnassembled[n] = ((dfloat*) rhs)[n];
  
  //sort by owner
  for (iint n=0;n<almond->Nnum;n++) 
    almond->rhsSort[n] = almond->rhsUnassembled[almond->sendSortId[n]];

  //Scatter nodes to their owners
  MPI_Alltoallv(almond->rhsSort, almond->sendCounts, almond->sendOffsets, MPI_CHAR,
                almond->rhsUnassembled, almond->recvCounts, almond->recvOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  //sort by globalid 
  for (iint n=0;n<almond->recvNnum;n++) 
    almond->rhsSort[n] = almond->rhsUnassembled[almond->globalSortId[n]];

  //gather
  for(iint n=0;n<almond->numLocalRows;++n){
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
  for(iint n=0;n<almond->numLocalRows;++n){
    for (iint id = almond->compressId[n];id<almond->compressId[n+1];id++) 
      almond->xSort[id] = almond->x[n]; 
  }

  //sort by original rank
  for (iint n=0;n<almond->recvNnum;n++) 
    almond->xUnassembled[almond->globalSortId[n]] = almond->xSort[n];

  //Scatter nodes back to their original rank
  MPI_Alltoallv(almond->xUnassembled, almond->recvCounts, almond->recvOffsets, MPI_CHAR,
                almond->xSort, almond->sendCounts, almond->sendOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  //sort back to original ordering
  for (iint n=0;n<almond->Nnum;n++) 
    almond->xUnassembled[almond->sendSortId[n]] = almond->xSort[n];


  for(iint i=0;i<almond->Nnum;++i) ((dfloat *) x)[i] = almond->xUnassembled[i];
  
  return 0;
}

int almondFree(void* A) {
  return 0;
}

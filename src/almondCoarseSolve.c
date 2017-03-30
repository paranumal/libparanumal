
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

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  almond->numLocalRows = rowStarts[rank+1]-rowStarts[rank];
  int globalOffset = rowStarts[rank];

;
  std::vector<int>    vRowStarts(almond->numLocalRows+1);

  // assumes presorted
  int cnt = 0; //nnz counter
  int cnt2 =0; //row start counter
  for(n=0;n<nnz;++n){
    if(  (iAi[n] >= (almond->numLocalRows + globalOffset)) || (iAi[n] < globalOffset))
      printf("errant nonzero %d,%d,%g, rank %d \n", iAi[n], iAj[n], dAvals[n], rank);

    if(cnt2==0 || (iAi[n]!=iAi[n-1])) vRowStarts[cnt2++] = cnt;
      
    if (iAj[n] >= globalOffset && iAj[n] < almond->numLocalRows + globalOffset) cnt++;
  }

  vRowStarts[cnt2] = cnt;

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

  almond->sendCounts = (int*) calloc(size,sizeof(int));
  almond->recvCounts = (int*) calloc(size,sizeof(int));
  for (n=0;n<size;n++){
    almond->sendCounts[n] = sendCounts[n]*sizeof(amgFloat);
    almond->recvCounts[n] = recvCounts[n]*sizeof(amgFloat);
  }

  almond->sendOffsets = (int*) calloc(size+1,sizeof(int));
  almond->recvOffsets = (int*) calloc(size+1,sizeof(int));
  for (n=0;n<size+1;n++){
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
  
  almond->M.setup(*(almond->A), almond->nullA, NULL);
  for (iint r=0;r<size;r++) {
    if (r==rank) {
      printf("----------Rank %d ------------------------\n", rank);
      almond->M.report();
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);  
  }
  almond->M.ktype = almond::PCG;
  
  almond->iintType = strdup(iintType);
  almond->dfloatType = strdup(dfloatType);

  return (void *) almond;
}


void * almondGlobalSetup(uint  Nnum,
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

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  almond->numLocalRows = rowStarts[rank+1]-rowStarts[rank];
  int globalOffset = rowStarts[rank];

  int numGlobalRows = rowStarts[size];

  std::vector<int>    vRowStarts(numGlobalRows+1);  
  std::vector<int>    vAj(nnz);
  std::vector<amgFloat> vAvals(nnz);

  // assumes presorted
  int cnt2 =0; //row start counter
  for(n=0;n<nnz;++n) {
    if(cnt2==0 || (iAi[n]!=iAi[n-1])) vRowStarts[cnt2++] = n;      
    vAj[n] = iAj[n];
    vAvals[n] = dAvals[n];
  }
  vRowStarts[cnt2] = nnz;
  
  almond->A = new almond::csr<amgFloat>(vRowStarts, vAj, vAvals);

  almond->nullA.resize(numGlobalRows);
  for (int i=0;i<numGlobalRows;i++)almond->nullA[i] = 1;
  
  almond->M.setup(*(almond->A), almond->nullA, rowStarts);
  for (iint r=0;r<size;r++) {
    if (r==rank) {
      printf("----------Rank %d ------------------------\n", rank);
      almond->M.report();
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);  
  }
  almond->M.ktype = almond::PCG;


  almond->Nnum = Nnum;
  almond->recvNnum = compressId[almond->numLocalRows];

  almond->sendSortId = (int*) calloc(Nnum,sizeof(int));
  for (n=0;n<Nnum;n++)  almond->sendSortId[n] = sendSortId[n];

  almond->sendCounts = (int*) calloc(size,sizeof(int));
  almond->recvCounts = (int*) calloc(size,sizeof(int));
  for (n=0;n<size;n++){
    almond->sendCounts[n] = sendCounts[n]*sizeof(amgFloat);
    almond->recvCounts[n] = recvCounts[n]*sizeof(amgFloat);
  }

  almond->sendOffsets = (int*) calloc(size+1,sizeof(int));
  almond->recvOffsets = (int*) calloc(size+1,sizeof(int));
  for (n=0;n<size+1;n++){
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
  
  almond->iintType = strdup(iintType);
  almond->dfloatType = strdup(dfloatType);

  return (void *) almond;
}


int almondSolve(void* x,
		void* A,
		void* rhs,
    void (*coarseSolve)(void *x, void *A, void *rhs),
    void *coarseA,
    int coarseTotal,
    int coarseOffset) {
  

  almond_t *almond = (almond_t*) A;

  iint rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

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
    almond->M.solve(almond->rhs, almond->x,coarseSolve,coarseA,coarseTotal,coarseOffset);
  }
  else{
    int maxIt = 40;
    amgFloat tol = 1e-2;
    almond::pcg<amgFloat>(almond->A[0],
			  almond->rhs,
			  almond->x,
			  almond->M,
			  maxIt,
			  tol,
        coarseSolve,coarseA,
        coarseTotal,coarseOffset);
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


void almondProlongateCoarseProblem(void *ALMOND, int *coarseNp, int *coarseOffsets, void **B) { 

  almond_t *almond = (almond_t*) ALMOND;

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  int numLevels = almond->M.levels.size();
  
  //number of coarse basis vectors
  int localCoarseNp = almond->M.levels[numLevels-1]->A.nrows();

  MPI_Allgather(&localCoarseNp, 1, MPI_INT, coarseNp, 1, MPI_INT, MPI_COMM_WORLD);

  for (iint n=0;n<size;n++) coarseOffsets[n+1] = coarseOffsets[n] + coarseNp[n];

  int coarseTotal = coarseOffsets[size];

  *B = (amgFloat *) calloc(localCoarseNp*almond->Nnum, sizeof(amgFloat));

  //allocate storage for prolanged basis vector at all levels
  std::vector<amgFloat> *b = (std::vector<amgFloat> *) calloc(numLevels,sizeof(std::vector<amgFloat>));
  for (int n=0; n<numLevels;n++)
    b[n].resize(almond->M.levels[n]->A.nrows());

  for (int k=0;k<localCoarseNp;k++) {
    //fill coarsest vector
    for (int n=0;n<localCoarseNp;n++) b[numLevels-1][n] = 0.0;
    b[numLevels-1][k] = 1.0;

    //prolongate
    for (int m=numLevels-1;m>0;m--) 
      almond->M.levels[m-1]->interpolate(b[m], b[m-1]);
    
    //scatter
    for(int n=0;n<almond->numLocalRows;++n){
      for (iint id = almond->compressId[n];id<almond->compressId[n+1];id++) 
        almond->xSort[id] = b[0][n]; 
    }

    //sort by original rank
    for (int n=0;n<almond->recvNnum;n++) 
      almond->xUnassembled[almond->globalSortId[n]] = almond->xSort[n];

    //Fake scatter nodes back to their original rank
    int recvCount = almond->recvCounts[rank]/sizeof(amgFloat);
    int recvOffset = almond->recvOffsets[rank]/sizeof(amgFloat);
    int sendOffset = almond->sendOffsets[rank]/sizeof(amgFloat);
    for (int n=0; n<almond->Nnum; n++) almond->xSort[n] = 0.;
    for (int n=0; n<recvCount;n++) almond->xSort[n+sendOffset] = almond->xUnassembled[n+recvOffset];

    //sort back to original ordering
    for (int n=0;n<almond->Nnum;n++) 
      almond->xUnassembled[almond->sendSortId[n]] = almond->xSort[n];
    
    //save
    for (int n=0;n<almond->Nnum;n++) 
      ((amgFloat *)*B)[n+k*almond->Nnum] = almond->xUnassembled[n];
  }
}

void almondGlobalCoarseSetup(void *ALMOND, int *coarseNp, int *coarseOffsets, int **globalNumbering,
                    int *nnz, int **rows, int **cols, void **vals) {

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  almond_t *almond = (almond_t*) ALMOND;

  int numLevels = almond->M.levels.size();

  int Np = almond->M.levels[numLevels-1]->A.nrows();

  MPI_Allgather(&Np, 1, MPI_INT, coarseNp, 1, MPI_INT, MPI_COMM_WORLD);

  coarseOffsets[0] = 0;
  for (int r=0;r<size;r++)
    coarseOffsets[r+1] = coarseOffsets[r] + coarseNp[r];

  *globalNumbering = (int *) calloc(coarseOffsets[size],sizeof(int));
  for (int n=0;n<coarseOffsets[size];n++)
    (*globalNumbering)[n] = n;

  *nnz = almond->M.levels[numLevels-1]->A.nnz();
  amgFloat *dvals;
  if (*nnz) {
    *rows = (int *) calloc(*nnz,sizeof(int));
    *cols = (int *) calloc(*nnz,sizeof(int));
    dvals = (amgFloat *) calloc(*nnz,sizeof(amgFloat));
  }

  for (int n=0;n<Np;n++) {
    for (int m=almond->M.levels[numLevels-1]->A.rowStarts[n];
             m<almond->M.levels[numLevels-1]->A.rowStarts[n+1];m++) {
      (*rows)[m] = n + coarseOffsets[rank];
      int col = almond->M.levels[numLevels-1]->A.cols[m];
      (*cols)[m] = almond->M.levels[numLevels-1]->A.colMap[col];
      dvals[m] = almond->M.levels[numLevels-1]->A.coefs[m];
    }
  }

  *vals = dvals;
}

